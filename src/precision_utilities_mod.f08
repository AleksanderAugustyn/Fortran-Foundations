!> Module for precision and tolerance utilities, used to compare floating point numbers
!!
!! Defines kinds for integers and reals, as well as tolerance constants and helper functions for comparing floating point numbers.
!! If you need to change precision (for example real64 -> real128), change only rk parameter, the change will propagate to all reals
!!
!! @note
!! Uses iso_fortran_env for guaranteed bit-width kinds independent of C ABI.
!! For C-interoperable type aliases, see c_bindings_mod.
!! @endnote
module precision_utilities_mod

    use, intrinsic :: iso_fortran_env, only: int32, int64, real32, real64, real128
#ifdef DEBUG
    use, intrinsic :: ieee_arithmetic, only: ieee_is_nan, ieee_is_finite
#endif

    implicit none

    private

    ! -- Integer Kinds --
    integer, parameter, public :: ik = int32    !! Default integer kind (32-bit)
    integer, parameter, public :: ikl = int64   !! 64-bit integers, use only when needed

    ! -- Real Kind --
    integer, parameter, public :: rk_single = real32    !! Single precision real kind (32-bit), use only when needed
    integer, parameter, public :: rk = real64           !! Default real kind (64-bit double precision)
    integer, parameter, public :: rk_high = real128     !! High precision real kind (128-bit), use only when needed

    ! -- Tolerance Constants and Helper Functions for Reals --
    real(kind = rk), parameter, public :: REAL_TOL = sqrt(epsilon(1.0_rk))  !! Tolerance for comparing two reals, can be adjusted if needed

    ! Public functions
    public :: are_close_f, is_zero_f

contains

    !> A helper function for comparing two reals, handling NaNs and Infinities correctly
    !!
    !! Compares two real numbers to determine if they are "close enough" to be considered equal, taking into account special cases like NaN and Infinity.
    !! NaN is not considered close to anything, including itself.
    !! Two infinities are considered close only if they are of the same sign.
    !! Exactly equal finite numbers (in the IEEE 754 sense) are considered close. This correctly treats +0.0 and -0.0 as equal.
    !! For all other finite numbers, a relative/absolute tolerance check is performed.
    !!
    !! @note
    !! In Release builds without DEBUG defined, NaN/Inf checks are disabled for performance.
    !! @endnote
    pure logical function are_close_f(val1, val2)
        real(kind = rk), intent(in) :: val1 !! First real number to compare
        real(kind = rk), intent(in) :: val2 !! Second real number to compare

#ifdef DEBUG
        ! First, handle NaNs. According to the standard, NaN is not equal to anything,
        ! including itself. Therefore, if either value is NaN, they are not close.
        if (ieee_is_nan(val1) .or. ieee_is_nan(val2)) then
            are_close_f = .false.
            return
        end if

        ! Second, handle the case where both numbers are infinite.
        ! "Close" only makes sense if they are the *same* infinity.
        ! This check is true if both are +Inf or both are -Inf.
        if (.not. ieee_is_finite(val1) .and. .not. ieee_is_finite(val2)) then
            are_close_f = (val1 > 0.0_rk .and. val2 > 0.0_rk) .or. (val1 < 0.0_rk .and. val2 < 0.0_rk)
            return
        end if
#endif

        ! Third, handle the most common case of identical finite numbers.
        ! This also correctly handles +0.0 vs -0.0_rk
        if (abs(val1 - val2) <= 0.0_rk) then
            are_close_f = .true.
            return
        end if

        ! Finally, for all remaining cases (which must be non-identical, finite numbers),
        ! perform the standard relative/absolute tolerance check.
        are_close_f = (abs(val1 - val2) <= max(REAL_TOL * max(abs(val1), abs(val2)), REAL_TOL))

        return
    end function are_close_f

    !> A helper function for checking "close to zero"
    !!
    !! Checks if the absolute value of a real number is less than a predefined tolerance. Is used to avoid direct equality comparisons between floating point numbers and zero.
    pure logical function is_zero_f(val)
        real(kind = rk), intent(in) :: val  !! Input real number to check
        is_zero_f = (abs(val) < REAL_TOL)

        return
    end function is_zero_f

end module precision_utilities_mod