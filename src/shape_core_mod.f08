!> Shared invalidation engine and status codes for shape parameterization
!! libraries (the shape_core deliverable of the shape parameterization
!! contract, atlas contracts/shape-parameterization-api.md).
!!
!! The engine is bookkeeping-only: it never calls library code. A library
!! embeds one shape_engine_t per thread inside its cache type, calls
!! shape_engine_begin_s with each incoming parameter vector, queries
!! shape_engine_needs_f per cached intermediate, recomputes stale ones
!! itself (producers before consumers), and reports each recompute via
!! shape_engine_note_computed_s.
!!
!! Parameters are diffed as ikl bit patterns via transfer, never by FP
!! comparison: deterministic under -ffast-math; -0.0 vs +0.0 counts as a
!! change; identical NaN bit patterns count as unchanged.
module shape_core_mod

    use precision_utilities_mod, only: rk, ik, ikl

    implicit none
    private

    public :: shape_engine_t
    public :: shape_engine_init_s
    public :: shape_engine_begin_s
    public :: shape_engine_needs_f
    public :: shape_engine_note_computed_s
    public :: shape_engine_invalidate_all_s
    public :: shape_engine_recompute_count_f
    public :: SHAPE_CACHE_MAX_PARAMS
    public :: SHAPE_STANDALONE_MAX_PARAMS
    public :: SHAPE_MAX_INTERMEDIATES
    public :: SHAPE_VALID
    public :: SHAPE_ERROR_TOO_MANY_PARAMS
    public :: SHAPE_ERROR_CACHE_NOT_INITIALIZED
    public :: SHAPE_ERROR_INVALID_GRID
    public :: SHAPE_ERROR_WRONG_PARAM_COUNT
    public :: SHAPE_ERROR_INVALID_INIT
    public :: SHAPE_ERROR_TABLES_NOT_INITIALIZED

    !> Maximum parameter count for the cached tier (tiers 2/3). Binds the
    !! public cache_init boundary of every adopting library.
    integer(kind = ik), parameter :: SHAPE_CACHE_MAX_PARAMS = 8_ik
    !> Maximum parameter count for one-shot (tier-1) calls; effective cap is
    !! min(64, N_max) per library.
    integer(kind = ik), parameter :: SHAPE_STANDALONE_MAX_PARAMS = 64_ik
    !> Maximum number of cached intermediates one engine can track. Engine
    !! sizing, not a contract constant: ~2x the largest planned dependency
    !! map; raising it later is a backward-compatible minor bump.
    integer(kind = ik), parameter :: SHAPE_MAX_INTERMEDIATES = 16_ik

    ! Shared status codes (contract range 0-99, append-only, never
    ! renumbered; 7-99 reserved). INVALID_GRID and TABLES_NOT_INITIALIZED
    ! are raised by libraries, not the engine - they live here so every
    ! library reports the same code for the same failure.
    !> Success.
    integer(kind = ik), parameter :: SHAPE_VALID = 0_ik
    !> n_params above SHAPE_CACHE_MAX_PARAMS at engine/cache init.
    integer(kind = ik), parameter :: SHAPE_ERROR_TOO_MANY_PARAMS = 1_ik
    !> Compute call before init.
    integer(kind = ik), parameter :: SHAPE_ERROR_CACHE_NOT_INITIALIZED = 2_ik
    !> Grid resolution below the library's documented floor (library-raised).
    integer(kind = ik), parameter :: SHAPE_ERROR_INVALID_GRID = 3_ik
    !> params length differs from the n_params fixed at init.
    integer(kind = ik), parameter :: SHAPE_ERROR_WRONG_PARAM_COUNT = 4_ik
    !> Invalid init arguments: n_params < 1, intermediates outside [1, 16],
    !! or a mask using bits >= n_params.
    integer(kind = ik), parameter :: SHAPE_ERROR_INVALID_INIT = 5_ik
    !> cache_init_shared with an uninitialized tables_t (library-raised).
    integer(kind = ik), parameter :: SHAPE_ERROR_TABLES_NOT_INITIALIZED = 6_ik

    !> Compile-time assert: the bit-pattern diff transfers rk reals into ikl
    !! integers, so their widths must match exactly. If rk ever grows past
    !! ikl (e.g. real128), the initializer below divides by zero and the
    !! module fails to compile - a truncated diff would silently leave stale
    !! intermediates marked valid. Value is 1 when the widths match; it is
    !! folded into param_bits' extent so the constant is always referenced
    !! (-Wunused-parameter runs under -Werror in this toolchain).
    integer(kind = ik), parameter :: SHAPE_RK_IKL_WIDTH_ASSERT = &
            1_ik / merge(1_ik, 0_ik, storage_size(1.0_rk) == bit_size(0_ikl))

    !> Per-thread invalidation engine. Embed one inside a library cache
    !! type; all components are private - interact only through the
    !! shape_engine_* procedures. Default initialization is the
    !! uninitialized state (every procedure is safe on it).
    type :: shape_engine_t
        private
        integer(kind = ik)  :: n_params        = 0_ik      !! Fixed at init
        integer(kind = ik)  :: n_intermediates = 0_ik      !! Fixed at init (= size(masks))
        logical             :: initialized     = .false.   !! Set only by a successful init
        logical             :: has_params      = .false.   !! .false. => next begin_s is cold
        !> Bit patterns of the last accepted vector. The extent folds in the
        !! width assert (value 1), keeping the assert referenced and tying it
        !! to the exact array whose width it guarantees.
        integer(kind = ikl) :: param_bits(SHAPE_CACHE_MAX_PARAMS * SHAPE_RK_IKL_WIDTH_ASSERT) = 0_ikl
        integer(kind = ik)  :: masks(SHAPE_MAX_INTERMEDIATES)           = 0_ik   !! Bit i-1 set => depends on parameter i
        logical             :: valid(SHAPE_MAX_INTERMEDIATES)           = .false. !! Per-intermediate validity stamps
        integer(kind = ikl) :: recompute_count(SHAPE_MAX_INTERMEDIATES) = 0_ikl  !! Test-facing recompute counters
    end type shape_engine_t

contains

    !> Initialize the engine: fix the parameter count and the per-intermediate
    !! dependency masks; reset all validity stamps and counters.
    !!
    !! intent(out) on engine means entry resets every component to its
    !! default initializer, so a failed init always leaves an engine with
    !! initialized = .false. - no half-configured state exists.
    !!
    !! Validation order is normative: n_params, then size(masks), then mask
    !! bits - the first failure determines the status.
    pure subroutine shape_engine_init_s(engine, n_params, masks, status)
        type(shape_engine_t), intent(out) :: engine    !! Engine to configure (reset on entry)
        integer(kind = ik), intent(in) :: n_params     !! Parameter count, in [1, SHAPE_CACHE_MAX_PARAMS]
        integer(kind = ik), intent(in) :: masks(:)     !! Dependency mask per intermediate; bit i-1 = parameter i
        integer(kind = ik), intent(out) :: status      !! SHAPE_VALID, SHAPE_ERROR_TOO_MANY_PARAMS, or SHAPE_ERROR_INVALID_INIT

        integer(kind = ik) :: n_intermediates
        integer(kind = ik) :: allowed_bits
        integer(kind = ik) :: j

        if (n_params > SHAPE_CACHE_MAX_PARAMS) then
            status = SHAPE_ERROR_TOO_MANY_PARAMS
            return
        end if
        if (n_params < 1_ik) then
            status = SHAPE_ERROR_INVALID_INIT
            return
        end if

        n_intermediates = size(masks, kind = ik)
        if (n_intermediates < 1_ik .or. n_intermediates > SHAPE_MAX_INTERMEDIATES) then
            status = SHAPE_ERROR_INVALID_INIT
            return
        end if

        ! Bits 0 .. n_params-1 set; any mask bit outside this range is a
        ! dependency on a parameter that does not exist
        allowed_bits = maskr(n_params, ik)
        do j = 1_ik, n_intermediates
            if (iand(masks(j), not(allowed_bits)) /= 0_ik) then
                status = SHAPE_ERROR_INVALID_INIT
                return
            end if
        end do

        engine%n_params = n_params
        engine%n_intermediates = n_intermediates
        engine%masks(1:n_intermediates) = masks
        engine%initialized = .true.
        status = SHAPE_VALID
    end subroutine shape_engine_init_s

    !> Present a parameter vector for the next compute call. Diffs it
    !! against the previous vector by ikl bit pattern and clears the
    !! validity stamp of every intermediate whose mask intersects the
    !! changed set; a cold engine (no previous vector) clears all stamps.
    !!
    !! On SHAPE_ERROR_WRONG_PARAM_COUNT the engine self-invalidates (all
    !! stamps cleared, next begin_s cold) - contract failure semantics: any
    !! nonzero status means the next call runs cold. The library's half is
    !! zero-filling its outputs and returning.
    pure subroutine shape_engine_begin_s(engine, params, status)
        type(shape_engine_t), intent(inout) :: engine  !! Engine tracking this cache
        real(kind = rk), intent(in) :: params(:)       !! Parameter vector, length must equal init n_params
        integer(kind = ik), intent(out) :: status      !! SHAPE_VALID, SHAPE_ERROR_CACHE_NOT_INITIALIZED, or SHAPE_ERROR_WRONG_PARAM_COUNT

        integer(kind = ikl) :: bits(SHAPE_CACHE_MAX_PARAMS)
        integer(kind = ik) :: changed
        integer(kind = ik) :: i
        integer(kind = ik) :: j

        if (.not. engine%initialized) then
            status = SHAPE_ERROR_CACHE_NOT_INITIALIZED
            return
        end if
        if (size(params, kind = ik) /= engine%n_params) then
            engine%valid = .false.
            engine%has_params = .false.
            status = SHAPE_ERROR_WRONG_PARAM_COUNT
            return
        end if

        do i = 1_ik, engine%n_params
            bits(i) = transfer(params(i), 0_ikl)
        end do

        if (engine%has_params) then
            changed = 0_ik
            do i = 1_ik, engine%n_params
                if (bits(i) /= engine%param_bits(i)) then
                    changed = ibset(changed, i - 1_ik)
                end if
            end do
            do j = 1_ik, engine%n_intermediates
                if (iand(engine%masks(j), changed) /= 0_ik) then
                    engine%valid(j) = .false.
                end if
            end do
        else
            engine%valid = .false.
        end if

        engine%param_bits(1:engine%n_params) = bits(1:engine%n_params)
        engine%has_params = .true.
        status = SHAPE_VALID
    end subroutine shape_engine_begin_s

    !> Does this intermediate need recomputation? Fail-safe: pre-init and
    !! out-of-range queries return .true. - the worst case is a redundant
    !! recompute, never a stale result.
    pure function shape_engine_needs_f(engine, intermediate) result(needs)
        type(shape_engine_t), intent(in) :: engine     !! Engine to query
        integer(kind = ik), intent(in) :: intermediate !! 1-based intermediate index
        logical :: needs                               !! .true. when recomputation is required

        if (.not. engine%initialized) then
            needs = .true.
        else if (intermediate < 1_ik .or. intermediate > engine%n_intermediates) then
            needs = .true.
        else
            needs = .not. engine%valid(intermediate)
        end if
    end function shape_engine_needs_f

    !> Report that an intermediate was just recomputed for the current
    !! parameter vector: sets its validity stamp and increments its
    !! recompute counter. Infallible: pre-init or out-of-range is a no-op
    !! (the contract's minimality tests make misuse visible).
    pure subroutine shape_engine_note_computed_s(engine, intermediate)
        type(shape_engine_t), intent(inout) :: engine  !! Engine to update
        integer(kind = ik), intent(in) :: intermediate !! 1-based intermediate index

        if (.not. engine%initialized) return
        if (intermediate < 1_ik .or. intermediate > engine%n_intermediates) return
        engine%valid(intermediate) = .true.
        engine%recompute_count(intermediate) = engine%recompute_count(intermediate) + 1_ikl
    end subroutine shape_engine_note_computed_s

    !> Clear every validity stamp and forget the stored parameter vector:
    !! the next begin_s runs cold. Infallible; counters are NOT reset (they
    !! accumulate until the next init). Libraries MUST call this on any
    !! failure inside a cached compute - contract failure semantics
    !! (zero-filled outputs are the library's half).
    pure subroutine shape_engine_invalidate_all_s(engine)
        type(shape_engine_t), intent(inout) :: engine  !! Engine to reset to cold

        engine%valid = .false.
        engine%has_params = .false.
    end subroutine shape_engine_invalidate_all_s

    !> Recompute counter accessor for the contract's mandated minimality
    !! tests. Returns 0 for pre-init or out-of-range queries.
    pure function shape_engine_recompute_count_f(engine, intermediate) result(count)
        type(shape_engine_t), intent(in) :: engine     !! Engine to query
        integer(kind = ik), intent(in) :: intermediate !! 1-based intermediate index
        integer(kind = ikl) :: count                   !! Times this intermediate was recomputed since init

        if (.not. engine%initialized) then
            count = 0_ikl
        else if (intermediate < 1_ik .or. intermediate > engine%n_intermediates) then
            count = 0_ikl
        else
            count = engine%recompute_count(intermediate)
        end if
    end function shape_engine_recompute_count_f

end module shape_core_mod
