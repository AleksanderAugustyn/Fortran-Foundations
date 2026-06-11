!> Test program for Fortran-Foundations.
!!
!! Covers the complete elliptic integral routines (both entry points, the
!! degenerate k² = 1 handling with status codes, and the near-degenerate
!! regime that regressed in the 2026-06-11 PES bug) and is_zero_f.
!!
!! Reference values computed with mpmath at 25 significant digits
!! (ellipk/ellipe take m = k²).
program test_fortran_foundations_p

    use precision_utilities_mod, only: rk, ik, is_zero_f
    use mathematical_utilities_mod, only: complete_elliptic_integral_generalized_s, &
            complete_elliptic_integral_complementary_s

    implicit none

    ! Reference values (mpmath at mp.dps = 25, rounded to 16 significant digits
    ! so the literals stay within double precision)
    real(kind = rk), parameter :: K_AT_HALF = 1.685750354812596_rk     !! K(k = 0.5)
    real(kind = rk), parameter :: E_AT_HALF = 1.467462209339427_rk     !! E(k = 0.5)
    real(kind = rk), parameter :: K_AT_06 = 1.750753802915753_rk       !! K(k = 0.6)
    real(kind = rk), parameter :: E_AT_06 = 1.418083394448724_rk       !! E(k = 0.6)
    real(kind = rk), parameter :: PI_OVER_TWO = 1.570796326794897_rk   !! K(0) = E(0)
    real(kind = rk), parameter :: K_AT_C_1E11 = 14.05051237261977_rk   !! K at 1 - k² = 1e-11
    real(kind = rk), parameter :: K_AT_C_2E11 = 13.70393878237069_rk   !! K at 1 - k² = 2e-11

    integer(kind = ik) :: failure_count
    real(kind = rk) :: result_value
    integer(kind = ik) :: status

    failure_count = 0_ik

    ! --- Moderate-modulus values via the primary (modulus) entry point ---
    call complete_elliptic_integral_generalized_s(result_value, 0.5_rk, 1.0_rk, 1.0_rk)
    call check_close("K(0.5) primary entry", result_value, K_AT_HALF, 1.0e-13_rk)

    call complete_elliptic_integral_generalized_s(result_value, 0.5_rk, 1.0_rk, 0.75_rk)
    call check_close("E(0.5) primary entry", result_value, E_AT_HALF, 1.0e-13_rk)

    call complete_elliptic_integral_generalized_s(result_value, 0.0_rk, 1.0_rk, 1.0_rk)
    call check_close("K(0) = pi/2", result_value, PI_OVER_TWO, 1.0e-14_rk)

    ! Success path must report status = 0 when the argument is supplied
    call complete_elliptic_integral_generalized_s(result_value, 0.5_rk, 1.0_rk, 1.0_rk, status)
    call check_status("K(0.5) success status", status, 0_ik)

    ! --- Complementary entry point agrees with the primary one ---
    call complete_elliptic_integral_complementary_s(result_value, 0.64_rk, 1.0_rk, 1.0_rk)
    call check_close("K(0.6) complementary entry", result_value, K_AT_06, 1.0e-13_rk)

    call complete_elliptic_integral_complementary_s(result_value, 0.64_rk, 1.0_rk, 0.64_rk, status)
    call check_close("E(0.6) complementary entry", result_value, E_AT_06, 1.0e-13_rk)
    call check_status("E(0.6) success status", status, 0_ik)

    ! --- Regression for the 2026-06-11 bug: near-degenerate 1 - k² far below ---
    ! --- the old sqrt(epsilon) window must follow the log singularity, not  ---
    ! --- jump to a sentinel                                                 ---
    call complete_elliptic_integral_complementary_s(result_value, 1.0e-11_rk, 1.0_rk, 1.0_rk, status)
    call check_close("K at 1 - k**2 = 1e-11", result_value, K_AT_C_1E11, 1.0e-12_rk)
    call check_status("K at 1 - k**2 = 1e-11 status", status, 0_ik)

    call complete_elliptic_integral_complementary_s(result_value, 2.0e-11_rk, 1.0_rk, 1.0_rk, status)
    call check_close("K at 1 - k**2 = 2e-11", result_value, K_AT_C_2E11, 1.0e-12_rk)
    call check_status("K at 1 - k**2 = 2e-11 status", status, 0_ik)

    ! --- Exact degeneracy: divergent K-type integral reports status = 1 ---
    call complete_elliptic_integral_complementary_s(result_value, 0.0_rk, 1.0_rk, 1.0_rk, status)
    call check_condition("degenerate K returns +huge", result_value >= huge(1.0_rk))
    call check_status("degenerate K status", status, 1_ik)

    call complete_elliptic_integral_complementary_s(result_value, 0.0_rk, 1.0_rk, -1.0_rk, status)
    call check_condition("degenerate K (b0 < 0) returns -huge", result_value <= -huge(1.0_rk))
    call check_status("degenerate K (b0 < 0) status", status, 1_ik)

    ! --- Exact degeneracy: E-type integral stays finite, E(1) = coeff_a0 ---
    call complete_elliptic_integral_complementary_s(result_value, 0.0_rk, 1.0_rk, 0.0_rk, status)
    call check_close("degenerate E returns coeff_a0", result_value, 1.0_rk, 1.0e-15_rk)
    call check_status("degenerate E status", status, 0_ik)

    ! Primary entry with k = 1 exactly: 1 - 1*1 == 0 exactly, same behavior
    call complete_elliptic_integral_generalized_s(result_value, 1.0_rk, 1.0_rk, 0.0_rk, status)
    call check_close("E(k=1) primary entry", result_value, 1.0_rk, 1.0e-15_rk)
    call check_status("E(k=1) primary status", status, 0_ik)

    call complete_elliptic_integral_generalized_s(result_value, 1.0_rk, 1.0_rk, 1.0_rk, status)
    call check_status("K(k=1) primary status", status, 1_ik)

    ! --- is_zero_f: default and explicit tolerance ---
    call check_condition("is_zero_f(1e-9) default", is_zero_f(1.0e-9_rk))
    call check_condition("is_zero_f(1e-7) default is false", .not. is_zero_f(1.0e-7_rk))
    call check_condition("is_zero_f(1e-7, tol=1e-6)", is_zero_f(1.0e-7_rk, 1.0e-6_rk))
    call check_condition("is_zero_f(1e-13, tol=1e-14) is false", .not. is_zero_f(1.0e-13_rk, 1.0e-14_rk))

    ! --- Summary ---
    if (failure_count > 0_ik) then
        print "(a, i0, a)", "FAILED: ", failure_count, " test(s)"
        error stop 1
    else
        print "(a)", "All tests passed"
    end if

contains

    !> Check that a value matches a reference to a relative tolerance.
    subroutine check_close(test_name, actual, expected, relative_tolerance)
        character(len = *), intent(in) :: test_name           !! Test label for the report
        real(kind = rk), intent(in) :: actual                 !! Computed value
        real(kind = rk), intent(in) :: expected               !! Reference value
        real(kind = rk), intent(in) :: relative_tolerance     !! Allowed |actual - expected| / |expected|

        if (abs(actual - expected) <= relative_tolerance * abs(expected)) then
            print "(a, a)", "PASS: ", test_name
        else
            print "(a, a, a, es24.16, a, es24.16)", "FAIL: ", test_name, &
                    " actual = ", actual, " expected = ", expected
            failure_count = failure_count + 1_ik
        end if
    end subroutine check_close

    !> Check that a status code has the expected value.
    subroutine check_status(test_name, actual, expected)
        character(len = *), intent(in) :: test_name      !! Test label for the report
        integer(kind = ik), intent(in) :: actual         !! Returned status
        integer(kind = ik), intent(in) :: expected       !! Expected status

        if (actual == expected) then
            print "(a, a)", "PASS: ", test_name
        else
            print "(a, a, a, i0, a, i0)", "FAIL: ", test_name, &
                    " status = ", actual, " expected = ", expected
            failure_count = failure_count + 1_ik
        end if
    end subroutine check_status

    !> Check that a logical condition holds.
    subroutine check_condition(test_name, condition)
        character(len = *), intent(in) :: test_name  !! Test label for the report
        logical, intent(in) :: condition             !! Condition that must be true

        if (condition) then
            print "(a, a)", "PASS: ", test_name
        else
            print "(a, a)", "FAIL: ", test_name
            failure_count = failure_count + 1_ik
        end if
    end subroutine check_condition

end program test_fortran_foundations_p
