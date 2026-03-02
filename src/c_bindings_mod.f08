!> Module for C-interoperability type aliases and conversion utilities
!!
!! Provides the centralized boundary layer between pure Fortran internals and
!! C-ABI consumers (Python/ctypes, Rust FFI, C++ extern "C"). Downstream API
!! wrapper modules should use this module instead of importing iso_c_binding directly.
!!
!! @note
!! Internal Fortran mathematics should use precision_utilities_mod (iso_fortran_env).
!! This module is only for code that crosses the language boundary.
!! @endnote
module c_bindings_mod

    use, intrinsic :: iso_c_binding, only: c_int, c_int64_t, c_float, c_double, c_bool, &
            c_ptr, c_char, c_null_ptr, c_null_char, c_f_pointer, c_loc, c_associated

    implicit none

    private

    ! -- C-equivalent type aliases --
    integer, parameter, public :: ik_c = c_int             !! C int
    integer, parameter, public :: ikl_c = c_int64_t        !! C int64_t
    integer, parameter, public :: rk_single_c = c_float    !! C float
    integer, parameter, public :: rk_c = c_double          !! C double
    integer, parameter, public :: bool_c = c_bool          !! C _Bool

    ! -- Re-exported iso_c_binding types for downstream convenience --
    public :: c_ptr, c_char, c_null_ptr, c_null_char

    ! -- String conversion utilities --
    public :: c_string_to_fortran, fortran_string_to_c

    ! -- Safety validation --
    public :: verify_c_bindings

contains

    !> Convert a null-terminated C string pointer to a Fortran allocatable string
    !!
    !! Walks the C memory starting at c_string_ptr until c_null_char is found,
    !! then copies the characters into a Fortran allocatable string.
    !!
    !! @warning
    !! The caller must ensure c_string_ptr points to valid, null-terminated memory.
    !! Passing c_null_ptr or a non-null-terminated buffer is undefined behaviour.
    !! @endwarning
    subroutine c_string_to_fortran(c_string_ptr, f_string)

        ! Dummy arguments
        type(c_ptr), intent(in) :: c_string_ptr                       !! Pointer to a null-terminated C string
        character(len = :), allocatable, intent(out) :: f_string      !! Resulting Fortran string

        ! Local variables
        integer, parameter :: MAX_SCAN_LENGTH = 1048576  !! 1 MiB safety limit
        character(kind = c_char), pointer :: char_array(:)
        integer :: str_length, i

        ! Guard against null pointer
        if (.not. c_associated(c_string_ptr)) then
            allocate(character(len = 0) :: f_string)
            return
        end if

        ! First pass: find the string length by scanning for null terminator
        call c_f_pointer(c_string_ptr, char_array, [MAX_SCAN_LENGTH])
        str_length = 0
        i = 1
        do while (i <= MAX_SCAN_LENGTH)
            if (char_array(i) == c_null_char) exit
            str_length = str_length + 1
            i = i + 1
        end do

        ! Second pass: copy characters into Fortran string
        allocate(character(len = str_length) :: f_string)
        do i = 1, str_length
            f_string(i:i) = char_array(i)
        end do

    end subroutine c_string_to_fortran

    !> Convert a Fortran string to a null-terminated C character array
    !!
    !! The output array has length len(f_string) + 1 to accommodate the null terminator.
    !! Suitable for passing string data back to C/C++/Rust callers.
    pure subroutine fortran_string_to_c(f_string, c_string)

        ! Dummy arguments
        character(len = *), intent(in) :: f_string                    !! Input Fortran string
        character(kind = c_char), intent(out) :: c_string(len(f_string) + 1)  !! Output null-terminated C array

        ! Local variables
        integer :: i, n

        n = len(f_string)
        do i = 1, n
            c_string(i) = f_string(i:i)
        end do
        c_string(n + 1) = c_null_char

    end subroutine fortran_string_to_c

    !> Verify that internal Fortran kinds and C-boundary kinds have matching storage sizes
    !!
    !! Asserts at runtime that rk == rk_c and ik == ik_c in byte width.
    !! Call this once during shared-library initialization to catch ABI mismatches
    !! caused by cross-compilation or non-standard compiler configurations.
    subroutine verify_c_bindings()

        use precision_utilities_mod, only: rk, ik

        if (storage_size(1.0_rk) /= storage_size(1.0_rk_c)) then
            error stop "FATAL: Fortran rk and C rk_c storage size mismatch"
        end if

        if (storage_size(1_ik) /= storage_size(1_ik_c)) then
            error stop "FATAL: Fortran ik and C ik_c storage size mismatch"
        end if

    end subroutine verify_c_bindings

end module c_bindings_mod
