!> Test program for shape_core_mod: invalidation engine, shared constants,
!! shared status codes.
!!
!! Implements the coverage list of
!! docs/superpowers/specs/2026-08-02-shape-core-design.md (tests 1-9).
program test_shape_core_p

    use precision_utilities_mod, only: rk, ik, ikl
    use shape_core_mod, only: shape_engine_t, shape_engine_init_s, &
            shape_engine_needs_f, shape_engine_note_computed_s, &
            shape_engine_recompute_count_f, shape_engine_begin_s, &
            shape_engine_invalidate_all_s, &
            SHAPE_CACHE_MAX_PARAMS, SHAPE_STANDALONE_MAX_PARAMS, &
            SHAPE_MAX_INTERMEDIATES, &
            SHAPE_VALID, SHAPE_ERROR_TOO_MANY_PARAMS, &
            SHAPE_ERROR_CACHE_NOT_INITIALIZED, SHAPE_ERROR_INVALID_GRID, &
            SHAPE_ERROR_WRONG_PARAM_COUNT, SHAPE_ERROR_INVALID_INIT, &
            SHAPE_ERROR_TABLES_NOT_INITIALIZED

    implicit none

    integer(kind = ik) :: failure_count
    type(shape_engine_t) :: engine
    type(shape_engine_t) :: virgin_engine
    integer(kind = ik) :: status
    integer(kind = ik) :: no_masks(0)
    integer(kind = ik) :: one_mask(1)
    integer(kind = ik) :: three_masks(3)
    integer(kind = ik) :: sixteen_masks(16)
    integer(kind = ik) :: seventeen_masks(17)
    real(kind = rk) :: params3(3)
    integer(kind = ikl), parameter :: NEG_ZERO_BITS = int(z'8000000000000000', kind = ikl)  !! IEEE754 -0.0
    integer(kind = ikl), parameter :: QNAN_BITS = int(z'7FF8000000000000', kind = ikl)      !! IEEE754 quiet NaN
    !> volatile: under -fno-signed-zeros (Release -ffast-math) the compiler
    !! may elide a store of -0.0 over a location it knows holds +0.0 (equal
    !! values under that flag), so the bit pattern would never reach the
    !! engine; volatile forces every store to be materialized.
    real(kind = rk), volatile :: special_value(1)
    real(kind = rk) :: params2(2)

    failure_count = 0_ik

    ! --- Constants and status codes are the contract's shared values ---
    call check_condition("SHAPE_CACHE_MAX_PARAMS is 8", SHAPE_CACHE_MAX_PARAMS == 8_ik)
    call check_condition("SHAPE_STANDALONE_MAX_PARAMS is 64", SHAPE_STANDALONE_MAX_PARAMS == 64_ik)
    call check_condition("SHAPE_MAX_INTERMEDIATES is 16", SHAPE_MAX_INTERMEDIATES == 16_ik)
    call check_status("SHAPE_VALID", SHAPE_VALID, 0_ik)
    call check_status("SHAPE_ERROR_TOO_MANY_PARAMS", SHAPE_ERROR_TOO_MANY_PARAMS, 1_ik)
    call check_status("SHAPE_ERROR_CACHE_NOT_INITIALIZED", SHAPE_ERROR_CACHE_NOT_INITIALIZED, 2_ik)
    call check_status("SHAPE_ERROR_INVALID_GRID", SHAPE_ERROR_INVALID_GRID, 3_ik)
    call check_status("SHAPE_ERROR_WRONG_PARAM_COUNT", SHAPE_ERROR_WRONG_PARAM_COUNT, 4_ik)
    call check_status("SHAPE_ERROR_INVALID_INIT", SHAPE_ERROR_INVALID_INIT, 5_ik)
    call check_status("SHAPE_ERROR_TABLES_NOT_INITIALIZED", SHAPE_ERROR_TABLES_NOT_INITIALIZED, 6_ik)

    ! --- Test 1: init validation, every error path + legal boundaries ---
    one_mask = [1_ik]

    ! n_params bounds
    call shape_engine_init_s(engine, 9_ik, one_mask, status)
    call check_status("init n_params=9 rejected", status, SHAPE_ERROR_TOO_MANY_PARAMS)
    call shape_engine_init_s(engine, 0_ik, one_mask, status)
    call check_status("init n_params=0 rejected", status, SHAPE_ERROR_INVALID_INIT)

    ! Validation precedence (spec: n_params, then size(masks), then mask bits):
    ! n_params=9 with an EMPTY mask array must report TOO_MANY_PARAMS, not INVALID_INIT
    call shape_engine_init_s(engine, 9_ik, no_masks, status)
    call check_status("init precedence: n_params before masks", status, SHAPE_ERROR_TOO_MANY_PARAMS)

    ! intermediates bounds
    call shape_engine_init_s(engine, 1_ik, no_masks, status)
    call check_status("init 0 intermediates rejected", status, SHAPE_ERROR_INVALID_INIT)
    seventeen_masks = 0_ik
    call shape_engine_init_s(engine, 1_ik, seventeen_masks, status)
    call check_status("init 17 intermediates rejected", status, SHAPE_ERROR_INVALID_INIT)

    ! mask uses bits >= n_params: bit 2 with n_params = 2
    call shape_engine_init_s(engine, 2_ik, [4_ik], status)
    call check_status("init mask bit >= n_params rejected", status, SHAPE_ERROR_INVALID_INIT)

    ! Legal boundaries: 1 param / 1 intermediate; 8 params / 16 intermediates
    ! (all 8 bits used); mask 0 legal
    call shape_engine_init_s(engine, 1_ik, one_mask, status)
    call check_status("init 1 param, 1 intermediate", status, SHAPE_VALID)
    sixteen_masks = 255_ik    ! all 8 parameter bits, every intermediate
    call shape_engine_init_s(engine, 8_ik, sixteen_masks, status)
    call check_status("init 8 params, 16 intermediates", status, SHAPE_VALID)
    call shape_engine_init_s(engine, 3_ik, [0_ik], status)
    call check_status("init mask 0 legal", status, SHAPE_VALID)

    ! --- Test 9 (edges): defined pre-init and out-of-range behavior ---
    ! Pre-init: needs is fail-safe .true., counter reads 0
    call check_condition("needs_f pre-init is true", shape_engine_needs_f(virgin_engine, 1_ik))
    call check_count("count pre-init is 0", shape_engine_recompute_count_f(virgin_engine, 1_ik), 0_ikl)

    ! Fresh init: every declared intermediate stale, counters 0
    three_masks = [1_ik, 6_ik, 7_ik]     ! deps: {1}, {2,3}, {1,2,3}
    call shape_engine_init_s(engine, 3_ik, three_masks, status)
    call check_status("init 3 params, 3 intermediates", status, SHAPE_VALID)
    call check_condition("post-init all stale", &
            shape_engine_needs_f(engine, 1_ik) .and. &
            shape_engine_needs_f(engine, 2_ik) .and. &
            shape_engine_needs_f(engine, 3_ik))

    ! note_computed_s sets the stamp and counts
    call shape_engine_note_computed_s(engine, 1_ik)
    call check_condition("note clears needs", .not. shape_engine_needs_f(engine, 1_ik))
    call check_count("note increments counter", shape_engine_recompute_count_f(engine, 1_ik), 1_ikl)

    ! Out-of-range queries: needs true, count 0; out-of-range note is a no-op
    call check_condition("needs_f(0) is true", shape_engine_needs_f(engine, 0_ik))
    call check_condition("needs_f(4) is true (beyond n_intermediates)", shape_engine_needs_f(engine, 4_ik))
    call check_condition("needs_f(17) is true", shape_engine_needs_f(engine, 17_ik))
    call check_count("count out-of-range is 0", shape_engine_recompute_count_f(engine, 0_ik), 0_ikl)
    call shape_engine_note_computed_s(engine, 0_ik)
    call shape_engine_note_computed_s(engine, 4_ik)
    call shape_engine_note_computed_s(engine, 17_ik)
    call check_condition("no-op notes leave stamps unchanged", &
            (.not. shape_engine_needs_f(engine, 1_ik)) .and. &
            shape_engine_needs_f(engine, 2_ik) .and. &
            shape_engine_needs_f(engine, 3_ik))
    ! All counters, not just I1: full observable state must be untouched
    ! (components are private, so "bitwise unchanged" is verified through
    ! the complete public surface; Task 4's selective-invalidation tests
    ! extend this to mask/bit integrity once begin_s exists)
    call check_count("no-op notes leave counter I1", shape_engine_recompute_count_f(engine, 1_ik), 1_ikl)
    call check_count("no-op notes leave counter I2", shape_engine_recompute_count_f(engine, 2_ik), 0_ikl)
    call check_count("no-op notes leave counter I3", shape_engine_recompute_count_f(engine, 3_ik), 0_ikl)

    ! note_computed_s pre-init is a no-op too (nothing observable changes)
    call shape_engine_note_computed_s(virgin_engine, 1_ik)
    call check_condition("pre-init note is a no-op", shape_engine_needs_f(virgin_engine, 1_ik))

    ! --- Tests 2/3/4: cold start, selective invalidation, union, no-op ---
    three_masks = [1_ik, 6_ik, 7_ik]     ! deps: {1}, {2,3}, {1,2,3}
    call shape_engine_init_s(engine, 3_ik, three_masks, status)
    call check_status("re-init for begin tests", status, SHAPE_VALID)

    ! Test 2: first begin_s is cold - every intermediate stale
    params3 = [1.0_rk, 2.0_rk, 3.0_rk]
    call shape_engine_begin_s(engine, params3, status)
    call check_status("cold begin status", status, SHAPE_VALID)
    call check_condition("cold begin marks all stale", &
            shape_engine_needs_f(engine, 1_ik) .and. &
            shape_engine_needs_f(engine, 2_ik) .and. &
            shape_engine_needs_f(engine, 3_ik))
    call shape_engine_note_computed_s(engine, 1_ik)
    call shape_engine_note_computed_s(engine, 2_ik)
    call shape_engine_note_computed_s(engine, 3_ik)

    ! Test 4: identical params - nothing stale, counters unchanged
    call shape_engine_begin_s(engine, params3, status)
    call check_status("no-op begin status", status, SHAPE_VALID)
    call check_condition("no-op begin keeps all valid", &
            .not. (shape_engine_needs_f(engine, 1_ik) .or. &
                   shape_engine_needs_f(engine, 2_ik) .or. &
                   shape_engine_needs_f(engine, 3_ik)))
    ! Spec test 4 requires ALL counters unchanged, not a sample
    call check_count("no-op keeps counter I1", shape_engine_recompute_count_f(engine, 1_ik), 1_ikl)
    call check_count("no-op keeps counter I2", shape_engine_recompute_count_f(engine, 2_ik), 1_ikl)
    call check_count("no-op keeps counter I3", shape_engine_recompute_count_f(engine, 3_ik), 1_ikl)

    ! Test 3: change parameter 1 only - exactly I1 (mask 1) and I3 (mask 7) stale
    params3 = [1.5_rk, 2.0_rk, 3.0_rk]
    call shape_engine_begin_s(engine, params3, status)
    call check_condition("param-1 change: I1 stale", shape_engine_needs_f(engine, 1_ik))
    call check_condition("param-1 change: I2 still valid", .not. shape_engine_needs_f(engine, 2_ik))
    call check_condition("param-1 change: I3 stale", shape_engine_needs_f(engine, 3_ik))
    call shape_engine_note_computed_s(engine, 1_ik)
    call shape_engine_note_computed_s(engine, 3_ik)
    ! Minimality via counters: I1 and I3 recomputed twice, I2 once
    call check_count("minimality counter I1", shape_engine_recompute_count_f(engine, 1_ik), 2_ikl)
    call check_count("minimality counter I2", shape_engine_recompute_count_f(engine, 2_ik), 1_ikl)
    call check_count("minimality counter I3", shape_engine_recompute_count_f(engine, 3_ik), 2_ikl)

    ! Test 3 (union): change params 1 AND 2 in one begin_s - dirty set is the
    ! union of both masks (catches overwrite-instead-of-union bugs: a wrong
    ! implementation that keeps only the last changed param's bit would
    ! leave I1 valid)
    params3 = [2.5_rk, 4.0_rk, 3.0_rk]
    call shape_engine_begin_s(engine, params3, status)
    call check_condition("two-param change: I1 stale", shape_engine_needs_f(engine, 1_ik))
    call check_condition("two-param change: I2 stale", shape_engine_needs_f(engine, 2_ik))
    call check_condition("two-param change: I3 stale", shape_engine_needs_f(engine, 3_ik))
    call shape_engine_note_computed_s(engine, 1_ik)
    call shape_engine_note_computed_s(engine, 2_ik)
    call shape_engine_note_computed_s(engine, 3_ik)

    ! --- Test 7: errors ---
    ! begin_s pre-init
    call shape_engine_begin_s(virgin_engine, params3, status)
    call check_status("begin pre-init", status, SHAPE_ERROR_CACHE_NOT_INITIALIZED)

    ! Wrong param count self-invalidates (contract: any nonzero status =>
    ! next call runs cold), even when the next vector equals the previous one
    call shape_engine_begin_s(engine, params3(1:2), status)
    call check_status("begin wrong count", status, SHAPE_ERROR_WRONG_PARAM_COUNT)
    call shape_engine_begin_s(engine, params3, status)
    call check_status("begin after wrong count", status, SHAPE_VALID)
    call check_condition("wrong count forced cold restart", &
            shape_engine_needs_f(engine, 1_ik) .and. &
            shape_engine_needs_f(engine, 2_ik) .and. &
            shape_engine_needs_f(engine, 3_ik))

    ! --- Test 5: bit semantics under fast-math ---
    call shape_engine_init_s(engine, 1_ik, [1_ik], status)
    call check_status("init for bit-semantics tests", status, SHAPE_VALID)

    ! +0.0 then -0.0: different bit patterns => invalidates
    special_value(1) = 0.0_rk
    call shape_engine_begin_s(engine, special_value, status)
    call shape_engine_note_computed_s(engine, 1_ik)
    special_value(1) = transfer(NEG_ZERO_BITS, 1.0_rk)
    call shape_engine_begin_s(engine, special_value, status)
    call check_condition("-0.0 vs +0.0 invalidates", shape_engine_needs_f(engine, 1_ik))
    call shape_engine_note_computed_s(engine, 1_ik)

    ! Same -0.0 again: identical bits => no invalidation
    special_value(1) = transfer(NEG_ZERO_BITS, 1.0_rk)
    call shape_engine_begin_s(engine, special_value, status)
    call check_condition("identical -0.0 does not invalidate", .not. shape_engine_needs_f(engine, 1_ik))

    ! NaN with identical bit pattern: counts as unchanged (bit diff, not FP
    ! comparison - NaN /= NaN would invalidate every call)
    special_value(1) = transfer(QNAN_BITS, 1.0_rk)
    call shape_engine_begin_s(engine, special_value, status)
    call check_condition("first NaN invalidates (bits changed)", shape_engine_needs_f(engine, 1_ik))
    call shape_engine_note_computed_s(engine, 1_ik)
    special_value(1) = transfer(QNAN_BITS, 1.0_rk)
    call shape_engine_begin_s(engine, special_value, status)
    call check_condition("identical NaN bits do not invalidate", .not. shape_engine_needs_f(engine, 1_ik))

    ! --- Test 6: interleaving - needs/note on subsets ---
    ! Two intermediates, both depending on the single param: output A uses
    ! I1, output B uses I2. Compute A's intermediate under vector p, then
    ! ask for B under the SAME p: the diff is empty but I2 was never
    ! computed - it must still read stale.
    call shape_engine_init_s(engine, 1_ik, [1_ik, 1_ik], status)
    call check_status("init for interleaving tests", status, SHAPE_VALID)
    params2(1) = 4.0_rk
    call shape_engine_begin_s(engine, params2(1:1), status)
    call shape_engine_note_computed_s(engine, 1_ik)          ! output A only
    call shape_engine_begin_s(engine, params2(1:1), status)  ! same vector: empty diff
    call check_condition("interleaving: A valid under empty diff", .not. shape_engine_needs_f(engine, 1_ik))
    call check_condition("interleaving: B stale despite empty diff", shape_engine_needs_f(engine, 2_ik))
    call shape_engine_note_computed_s(engine, 2_ik)

    ! --- Test 8: recovery via invalidate_all_s ---
    call shape_engine_invalidate_all_s(engine)
    call shape_engine_begin_s(engine, params2(1:1), status)  ! same vector again
    call check_condition("invalidate_all: next begin is cold", &
            shape_engine_needs_f(engine, 1_ik) .and. shape_engine_needs_f(engine, 2_ik))
    call shape_engine_note_computed_s(engine, 1_ik)
    call shape_engine_note_computed_s(engine, 2_ik)
    ! Counters accumulate across invalidate_all (reset only by init):
    ! I1 computed at cold start + after invalidate_all = 2; I2 likewise
    call check_count("counters survive invalidate_all I1", shape_engine_recompute_count_f(engine, 1_ik), 2_ikl)
    call check_count("counters survive invalidate_all I2", shape_engine_recompute_count_f(engine, 2_ik), 2_ikl)

    ! --- Summary ---
    if (failure_count > 0_ik) then
        print "(a, i0, a)", "FAILED: ", failure_count, " test(s)"
        error stop 1
    else
        print "(a)", "All tests passed"
    end if

contains

    !> Check that a logical condition holds.
    subroutine check_condition(test_name, condition)
        character(len = *), intent(in) :: test_name  !! Test label for the report
        logical, intent(in) :: condition             !! Must be .true. to pass

        if (condition) then
            print "(a, a)", "PASS: ", test_name
        else
            print "(a, a)", "FAIL: ", test_name
            failure_count = failure_count + 1_ik
        end if
    end subroutine check_condition

    !> Check that a status code has the expected value.
    subroutine check_status(test_name, actual, expected)
        character(len = *), intent(in) :: test_name  !! Test label for the report
        integer(kind = ik), intent(in) :: actual     !! Returned status
        integer(kind = ik), intent(in) :: expected   !! Required status

        if (actual == expected) then
            print "(a, a)", "PASS: ", test_name
        else
            print "(a, a, a, i0, a, i0)", "FAIL: ", test_name, &
                    " status = ", actual, " expected = ", expected
            failure_count = failure_count + 1_ik
        end if
    end subroutine check_status

    !> Check that a recompute counter has the expected value.
    subroutine check_count(test_name, actual, expected)
        character(len = *), intent(in) :: test_name  !! Test label for the report
        integer(kind = ikl), intent(in) :: actual    !! Returned counter
        integer(kind = ikl), intent(in) :: expected  !! Required counter

        if (actual == expected) then
            print "(a, a)", "PASS: ", test_name
        else
            print "(a, a, a, i0, a, i0)", "FAIL: ", test_name, &
                    " count = ", actual, " expected = ", expected
            failure_count = failure_count + 1_ik
        end if
    end subroutine check_count

end program test_shape_core_p
