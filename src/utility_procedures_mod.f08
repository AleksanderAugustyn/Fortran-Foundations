!> Module containing utility procedures for general operations
!!
!! This module provides utility procedures for sorting arrays and related operations.
!! This module has no dependencies on higher-level calculation or data modules.
module utility_procedures_mod

    use precision_utilities_mod, only: ik, rk

    implicit none

    private

    ! Public subroutines
    public :: sort_array_ascending_s, sort_ascending_with_indices_s

contains

    pure subroutine sort_array_ascending_s(array)
        implicit none
        real(kind = rk), intent(inout) :: array(:)

        integer(kind = ik) :: n, i
        real(kind = rk) :: temp

        n = size(array, kind = ik)
        if (n <= 1_ik) return

        ! Build max-heap (heapify)
        do i = n / 2_ik, 1_ik, -1_ik
            call sift_down(array, i, n)
        end do

        ! Extract elements from heap
        do i = n, 2_ik, -1_ik
            temp = array(1_ik)
            array(1_ik) = array(i)
            array(i) = temp
            call sift_down(array, 1_ik, i - 1_ik)
        end do

    contains

        pure subroutine sift_down(arr, start, end_idx)
            real(kind = rk), intent(inout) :: arr(:)
            integer(kind = ik), intent(in) :: start, end_idx

            integer(kind = ik) :: root, child, swap
            real(kind = rk) :: temp

            root = start
            do while (2_ik * root <= end_idx)
                child = 2_ik * root
                swap = root

                if (arr(swap) < arr(child)) swap = child
                if (child + 1_ik <= end_idx) then
                    if (arr(swap) < arr(child + 1_ik)) swap = child + 1_ik
                end if

                if (swap == root) return

                temp = arr(root)
                arr(root) = arr(swap)
                arr(swap) = temp
                root = swap
            end do
        end subroutine sift_down

    end subroutine sort_array_ascending_s

    !> Sorts a real array in ascending order with two associated integer arrays
    !!
    !! Performs an in-place ascending sort of array `a` using selection sort,
    !! simultaneously reordering two associated integer arrays (`k` and `it`) to
    !! maintain correspondence. This is commonly used for sorting eigenvalues while
    !! preserving quantum number labels.
    !!
    !! ## Algorithm
    !!
    !! Uses selection sort: for each position i, finds the minimum element in the
    !! unsorted portion [i+1, n] and swaps it with position i if necessary. All three
    !! arrays are reordered identically to maintain element-to-element correspondence.
    !!
    !! ## Complexity
    !!
    !! - Time: O(n²) comparisons, O(n²) worst-case swaps
    !! - Space: O(1) auxiliary storage
    !! - Suitable for small to moderate arrays (n ≲ 1000)
    !!
    !! @note
    !! For large arrays, consider using quicksort or merge sort for O(n log n) performance.
    !! However, selection sort has minimal memory overhead and is stable for small n.
    !! @endnote
    pure subroutine sort_ascending_with_indices_s(n, primary_values, associated_indices_1, associated_indices_2)

        implicit none

        ! Dummy arguments
        integer(kind = ik), intent(in), value :: n !! Number of elements in the arrays
        real(kind = rk), dimension(n), intent(inout) :: primary_values !! Primary array to sort (input: unsorted, output: ascending) [dimensionless]
        integer(kind = ik), dimension(n), intent(inout) :: associated_indices_1 !! First associated integer array (reordered to match sorted values)
        integer(kind = ik), dimension(n), intent(inout) :: associated_indices_2 !! Second associated integer array (reordered to match sorted values)

        ! Local variables
        integer(kind = ik) :: i !! Outer loop index (current position being filled)
        integer(kind = ik) :: j !! Inner loop index (search for minimum in unsorted portion)
        integer(kind = ik) :: temp_int_1 !! Temporary storage for swapping elements in `associated_indices_1`
        integer(kind = ik) :: temp_int_2 !! Temporary storage for swapping elements in `associated_indices_2`
        real(kind = rk) :: temp_real !! Temporary storage for swapping elements in `primary_values`

        ! Selection sort: for each position i = 1 to n-1
        ! Find minimum in [i, n] and swap to position i
        do i = 1, n - 1
            do j = i + 1, n
                ! If current element is greater than j-th element, swap them
                if (primary_values(i) > primary_values(j)) then
                    ! Swap primary values
                    temp_real = primary_values(i)
                    primary_values(i) = primary_values(j)
                    primary_values(j) = temp_real

                    ! Swap first associated indices
                    temp_int_1 = associated_indices_1(i)
                    associated_indices_1(i) = associated_indices_1(j)
                    associated_indices_1(j) = temp_int_1

                    ! Swap second associated indices
                    temp_int_2 = associated_indices_2(i)
                    associated_indices_2(i) = associated_indices_2(j)
                    associated_indices_2(j) = temp_int_2
                end if
            end do
        end do

    end subroutine sort_ascending_with_indices_s

end module utility_procedures_mod