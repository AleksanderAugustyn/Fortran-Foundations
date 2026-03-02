module mathematical_utilities_mod

    use precision_utilities_mod, only: ik, rk

    implicit none

    private

    ! Public functions
    public :: integrate_simpson_grid_f

    ! Public subroutines
    public :: complete_elliptic_integral_generalized_s
    public :: compute_gauss_legendre_quadrature_s
    public :: compute_spherical_harmonics_normalization_constants_s
    public :: compute_strutinsky_curvature_coefficients_s
    public :: compute_hermite_polynomials_s
    public :: compute_laguerre_polynomials_s
    public :: compute_gauss_hermite_quadrature_s
    public :: compute_gauss_laguerre_quadrature_s
    public :: compute_legendre_polynomials_s

contains

    !> Integrate a function sampled on a uniform grid using Simpson's 1/3 rule.
    !!
    !! Computes the definite integral of f(x) given discrete samples y_i = f(x_i)
    !! on a uniform grid with spacing dx, using composite Simpson's rule:
    !!
    !!   ∫f(x)dx ≈ (dx/3) × [y₁ + 4y₂ + 2y₃ + 4y₄ + ... + 4y_{n-1} + y_n]
    !!
    !! @warning
    !! Requires an odd number of grid points (even number of intervals).
    !! Validation checks are included in Debug builds.
    !! @endwarning
    pure function integrate_simpson_grid_f(y, dx) result(integral)

        ! Dummy arguments
        real(kind = rk), intent(in) :: y(:)    !! Function values at uniformly-spaced grid points. Array length must be odd and >= 3.
        real(kind = rk), intent(in) :: dx      !! Grid spacing (must be positive)

        ! Result variable
        real(kind = rk) :: integral    !! Approximation to the definite integral

        ! Local variables
        integer(kind = ik) :: n !! Number of grid points

        ! Get number of points
        n = size(y)

#ifdef DEBUG
    ! Validate inputs in Debug builds
    if (n < 3) error stop "integrate_simpson_grid_f: need at least 3 points"
    if (mod(n, 2) == 0) error stop "integrate_simpson_grid_f: n must be odd"
    if (dx <= 0.0_rk) error stop "integrate_simpson_grid_f: dx must be positive"
#endif

        ! Using array slices with stride 2:
        ! y(2:n-1:2) extracts y2, y4, y6...
        ! y(3:n-2:2) extracts y3, y5, y7...

        ! Compute integral using Simpson's 1/3 rule
        integral = (y(1) + y(n) + 4.0_rk * sum(y(2:n - 1:2)) + 2.0_rk * sum(y(3:n - 2:2))) * (dx / 3.0_rk)

    end function integrate_simpson_grid_f

    !> Computes a generalized complete elliptic integral using the AGM algorithm.
    !!
    !! Evaluates integrals of the form:
    !!
    !! \[ \int_0^{\pi/2} \frac{a_0 \cos^2\phi + b_0 \sin^2\phi}
    !!    {\sqrt{\cos^2\phi + (1-k^2)\sin^2\phi}} d\phi \]
    !!
    !! using the arithmetic-geometric mean iteration which converges quadratically.
    !!
    !! ## Special Cases
    !!
    !! - For \( a_0 = b_0 = 1 \): Returns \( K(k) \), the complete elliptic integral of the first kind
    !! - For \( a_0 = 1, b_0 = 1 - k^2 \): Returns \( E(k) \), the complete elliptic integral of the second kind
    !!
    !! ## Algorithm
    !!
    !! Uses the AGM (Arithmetic-Geometric Mean) iteration:
    !! 1. Initialize: \( a_0 = 1 \), \( g_0 = \sqrt{1-k^2} \)
    !! 2. Iterate: \( a_{n+1} = (a_n + g_n)/2 \), \( g_{n+1} = \sqrt{a_n g_n} \)
    !! 3. Converges when \( |a_n - g_n| < \epsilon \)
    !!
    !! @note
    !! Original name: cel2
    !! @endnote
    pure subroutine complete_elliptic_integral_generalized_s(result_integral, elliptic_modulus, coeff_a0, coeff_b0)

        use mathematical_and_physical_constants_mod, only: PI_OVER_FOUR_C
        use precision_utilities_mod, only: is_zero_f

        implicit none

        ! Arguments
        real(kind = rk), intent(out) :: result_integral !! The computed elliptic integral value
        real(kind = rk), intent(in), value :: elliptic_modulus !! The elliptic modulus k (must satisfy 0 ≤ k² ≤ 1)
        real(kind = rk), intent(in), value :: coeff_a0         !! First coefficient in the integral combination
        real(kind = rk), intent(in), value :: coeff_b0         !! Second coefficient in the integral combination

        ! Named constants
        real(kind = rk), parameter :: CONVERGENCE_TOLERANCE = 1.0e-8_rk !! Relative tolerance for AGM convergence
        real(kind = rk), parameter :: LARGE_VALUE = 1.0e75_rk           !! Representation of effectively infinite result

        ! Local variables for AGM iteration
        real(kind = rk) :: complementary_modulus_sq !! 1 - k², the complementary modulus squared (geo in original)
        real(kind = rk) :: agm_geometric            !! Geometric mean sequence value
        real(kind = rk) :: agm_arithmetic           !! Arithmetic mean sequence value (ari in original)
        real(kind = rk) :: agm_arithmetic_prev      !! Previous arithmetic mean (aari in original)
        real(kind = rk) :: coeff_current            !! Current a-coefficient in iteration (aao in original)
        real(kind = rk) :: coeff_combined           !! Combined coefficient (ano in original)
        real(kind = rk) :: weight_accumulator       !! Weighted sum accumulator (w in original)

        complementary_modulus_sq = 1.0_rk - elliptic_modulus * elliptic_modulus

        ! Validate: complementary modulus squared must be non-negative (i.e., |k| ≤ 1)
        if (complementary_modulus_sq < 0.0_rk) then
            error stop "complete_elliptic_integral_generalized_s: invalid elliptic modulus (|k| > 1)"
        else if (is_zero_f(complementary_modulus_sq)) then
            ! Degenerate case: k² = 1, integral has singularity
            if (coeff_b0 < 0.0_rk) then
                result_integral = -LARGE_VALUE
            else if (is_zero_f(coeff_b0)) then
                result_integral = coeff_a0
            else
                result_integral = LARGE_VALUE
            end if
            return
        end if

        ! Initialize AGM iteration
        agm_geometric = sqrt(complementary_modulus_sq)
        agm_arithmetic = 1.0_rk
        coeff_current = coeff_a0
        coeff_combined = coeff_a0 + coeff_b0
        weight_accumulator = coeff_b0

        ! AGM iteration loop - converges quadratically
        do
            weight_accumulator = weight_accumulator + coeff_current * agm_geometric
            weight_accumulator = weight_accumulator + weight_accumulator  ! Double the accumulator
            coeff_current = coeff_combined
            agm_arithmetic_prev = agm_arithmetic
            agm_arithmetic = agm_geometric + agm_arithmetic
            coeff_combined = weight_accumulator / agm_arithmetic + coeff_combined

            ! Check for convergence: relative change in arithmetic mean
            if (agm_arithmetic_prev - agm_geometric <= CONVERGENCE_TOLERANCE * agm_arithmetic_prev) exit

            ! Update geometric mean for next iteration
            agm_geometric = sqrt(agm_geometric * agm_arithmetic_prev)
            agm_geometric = agm_geometric + agm_geometric  ! Double for the AGM step
        end do

        ! Final result: π/4 * (combined coefficient) / (arithmetic mean)
        ! Note: The original uses π/4 directly, giving the specific integral normalization
        result_integral = PI_OVER_FOUR_C * coeff_combined / agm_arithmetic

    end subroutine complete_elliptic_integral_generalized_s

    !> Computes normalization constants for axially-symmetric spherical harmonics
    !!
    !! Calculates the normalization coefficients C_λ for spherical harmonics
    !! Y_{λ,0}(θ,φ) used in nuclear shape parameterization:
    !!
    !! C_λ = sqrt[(2λ + 1)/(4π)]
    !!
    !! These constants ensure proper normalization according to the convention:
    !!
    !! ∫₀²π dφ ∫₋₁¹ d(cosθ) |Y_{λ,μ}(θ,φ)|² = 1
    !!
    !! # Mathematical Background
    !!
    !! In the nuclear deformation expansion:
    !!
    !! R(θ,φ) = R₀[1 + Σ_{λ,μ} β_{λ,μ} Y_{λ,μ}(θ,φ)]
    !!
    !! the normalization constants ensure dimensional consistency and proper
    !! scaling of deformation parameters.
    !!
    !! @note
    !! Array index k corresponds to multipolarity λ = k
    !! @endnote
    pure subroutine compute_spherical_harmonics_normalization_constants_s(&
            spherical_harmonics_normalization_constants, &
            max_multipole_order)

        use mathematical_and_physical_constants_mod, only: PI_C
        implicit none

        ! Dummy arguments
        ! Input arguments
        integer(kind = ik), intent(in), value :: max_multipole_order                                                !! Maximum multipole order λ_max
        ! Output arguments
        real(kind = rk), dimension(max_multipole_order), intent(out) :: spherical_harmonics_normalization_constants !! Normalization constants C_λ for λ = 1, ..., max_multipole

        ! Local variables
        integer(kind = ik) :: lambda

        ! Input validation
        if (max_multipole_order < 1_ik .or. max_multipole_order > size(spherical_harmonics_normalization_constants)) then
            error stop "compute_spherical_harmonics_normalization_constants_s: invalid max_multipole_order"
        end if

        ! Compute C_λ = √[(2λ+1)/(4π)] for λ = 1, 2, ..., λ_max
        do lambda = 1, max_multipole_order
            spherical_harmonics_normalization_constants(lambda) = sqrt(real(2_ik * lambda + 1_ik, rk)) * (1.0_rk / sqrt(4.0_rk * PI_C))
        end do

    end subroutine compute_spherical_harmonics_normalization_constants_s

    !> Computes Gauss-Legendre quadrature nodes and weights for arbitrary order
    !!
    !! Uses Newton-Raphson iteration to find roots of the Legendre polynomial P_n(x),
    !! then computes the corresponding quadrature weights.
    !!
    !! ## Algorithm
    !!
    !! 1. Initial guess for root i: x_i ≈ cos(π(i - 0.25)/(n + 0.5))
    !! 2. Newton-Raphson: x_{k+1} = x_k - P_n(x_k) / P'_n(x_k)
    !! 3. Legendre polynomials computed via three-term recurrence
    !! 4. Weight: w_i = 2 / [(1 - x_i²)(P'_n(x_i))²]
    !!
    !! ## Symmetry
    !!
    !! Legendre polynomial roots are symmetric about x = 0, so only n/2 roots
    !! are computed; the rest are obtained by symmetry.
    !!
    !! @note
    !! Convergence tolerance is 3×10⁻¹⁴, suitable for double precision.
    !! @endnote
    pure subroutine compute_gauss_legendre_quadrature_s(n, nodes, weights)

        use mathematical_and_physical_constants_mod, only: PI_C
        implicit none

        ! Dummy arguments
        integer(kind = ik), intent(in), value :: n          !! Number of quadrature points
        real(kind = rk), intent(out) :: nodes(n)            !! Quadrature nodes on [-1, 1]
        real(kind = rk), intent(out) :: weights(n)          !! Quadrature weights

        ! Local constants
        real(kind = rk), parameter :: CONVERGENCE_TOL = 3.0e-14_rk  !! Newton-Raphson tolerance
        integer(kind = ik), parameter :: MAX_ITERATIONS = 100_ik   !! Max Newton iterations

        ! Local variables
        integer(kind = ik) :: i, j, m, iter
        real(kind = rk) :: z, z_old, p1, p2, p3, pp

        ! Number of roots to compute (exploit symmetry)
        m = (n + 1_ik) / 2_ik

        ! Find roots of P_n(x) using Newton-Raphson
        do i = 1_ik, m
            ! Initial approximation for root i
            z = cos(PI_C * (real(i, rk) - 0.25_rk) / (real(n, rk) + 0.5_rk))

            ! Newton-Raphson iteration
            do iter = 1_ik, MAX_ITERATIONS
                ! Compute P_n(z) and P'_n(z) using recurrence
                p1 = 1.0_rk    ! P_0(z)
                p2 = 0.0_rk    ! P_{-1}(z) (not used, but needed for recurrence start)

                do j = 1_ik, n
                    p3 = p2
                    p2 = p1
                    ! Three-term recurrence: (j+1)P_{j+1} = (2j+1)xP_j - jP_{j-1}
                    p1 = ((2.0_rk * real(j, rk) - 1.0_rk) * z * p2 - (real(j, rk) - 1.0_rk) * p3) / real(j, rk)
                end do

                ! p1 is now P_n(z)
                ! Derivative: P'_n(z) = n(zP_n - P_{n-1}) / (z² - 1)
                pp = real(n, rk) * (z * p1 - p2) / (z * z - 1.0_rk)

                z_old = z
                z = z_old - p1 / pp

                if (abs(z - z_old) < CONVERGENCE_TOL) exit
            end do

            ! Store symmetric roots and weights
            nodes(i) = z
            nodes(n + 1_ik - i) = -z

            ! Weight formula: w_i = 2 / [(1 - x_i²)(P'_n(x_i))²]
            weights(i) = 2.0_rk / ((1.0_rk - z * z) * pp * pp)
            weights(n + 1_ik - i) = weights(i)
        end do

    end subroutine compute_gauss_legendre_quadrature_s

    !> Compute Strutinsky curvature coefficients for shell correction averaging
    !!
    !! Calculates coefficients cm(k) = (-1)^(k/2) / (2^k * (k/2)!) for even k = 2, 4, 6, ...
    !! These coefficients are used in Hermite polynomial corrections for
    !! Strutinsky energy smoothing in nuclear structure calculations.
    !!
    !! ## Mathematical Background
    !!
    !! The Strutinsky shell correction method smooths single-particle level densities
    !! using Gaussian averaging. The curvature coefficients correct for the finite
    !! smoothing width via Hermite polynomial expansions:
    !!
    !!   n_i(λ) = 0.5*(1 + erf(u)) - exp(-u²)*P(u)/√π
    !!
    !! where P(u) = Σ cm(k)*H_k(u) is the curvature correction polynomial.
    !!
    !! ## Array Indexing
    !!
    !! Coefficients are stored at even indices: cm(2), cm(4), cm(6), ...
    !! Odd indices remain zero.
    pure subroutine compute_strutinsky_curvature_coefficients_s(max_order, coefficients)

        implicit none

        ! Dummy arguments
        integer(kind = ik), intent(in), value :: max_order            !! Maximum order (coefficients computed for even k up to max_order)
        real(kind = rk), intent(out) :: coefficients(max_order)       !! Curvature coefficients cm(k) for k = 1, ..., max_order

        ! Local variables
        real(kind = rk) :: factorial_product    !! Running product for (k/2)! computation
        integer(kind = ik) :: k                 !! Loop index over even orders

        ! Initialize all coefficients to zero (odd indices stay zero)
        coefficients = 0.0_rk

        ! Compute cm(k) = (-1)^(k/2) / (2^k * (k/2)!) for even k
        ! Use running factorial product: (k/2)! = 1 * 2 * 3 * ... * (k/2)
        ! which is equivalent to: k/2 * (k-2)/2 * (k-4)/2 * ... = k!! / 2^(k/2)
        factorial_product = 1.0_rk

        do k = 2, max_order, 2
            ! Update factorial: multiply by (k/2)
            ! This builds (k/2)! incrementally: 1! -> 2! -> 3! -> ...
            factorial_product = factorial_product * real(k, rk) / 2.0_rk

            ! cm(k) = (-1)^(k/2) / (2^k * (k/2)!)
            coefficients(k) = ((-1.0_rk)**(k / 2)) / (2.0_rk**k * factorial_product)
        end do

    end subroutine compute_strutinsky_curvature_coefficients_s

    !> Compute Hermite polynomials H_k(x) and optionally their derivatives up to order n
    !!
    !! Uses the recurrence relation for physicist's Hermite polynomials:
    !!   H_{k+1}(x) = 2x H_k(x) - 2k H_{k-1}(x)
    !! and the derivative relation:
    !!   H'_k(x) = 2k H_{k-1}(x)
    !!
    !! @note
    !! Original name: dhep
    !! @endnote
    pure subroutine compute_hermite_polynomials_s(x_eval, max_order, hermite_vals, hermite_derivs)
        real(kind = rk), intent(in), value :: x_eval                    !! Evaluation point
        integer(kind = ik), intent(in), value :: max_order              !! Maximum polynomial order (n >= 0)
        real(kind = rk), intent(out) :: hermite_vals(:)                 !! H_k(x) for k = 0, 1, ..., max_order (size >= max_order + 1)
        real(kind = rk), intent(out), optional :: hermite_derivs(:)     !! Optional H'_k(x) for k = 0, 1, ..., max_order (size >= max_order + 1)

        integer(kind = ik) :: i_order
        real(kind = rk) :: recurrence_term, derivative_term
        real(kind = rk) :: two_x
        logical :: compute_derivatives

        ! Check if derivatives are requested
        compute_derivatives = present(hermite_derivs)

        ! H_0(x) = 1, H'_0(x) = 0
        hermite_vals(1) = 1.0_rk
        if (compute_derivatives) hermite_derivs(1) = 0.0_rk
        if (max_order <= 0_ik) return

        ! H_1(x) = 2x, H'_1(x) = 2
        two_x = x_eval + x_eval
        hermite_vals(2) = two_x
        if (compute_derivatives) hermite_derivs(2) = 2.0_rk
        if (max_order <= 1_ik) return

        ! Recurrence: H_{k+1}(x) = 2x H_k(x) - 2k H_{k-1}(x)
        ! Derivative:  H'_{k+1}(x) = 2(k+1) H_k(x)
        do i_order = 2_ik, max_order
            recurrence_term = x_eval * hermite_vals(i_order) - real(i_order - 1_ik, rk) * hermite_vals(i_order - 1)
            hermite_vals(i_order + 1) = recurrence_term + recurrence_term

            if (compute_derivatives) then
                derivative_term = real(i_order, rk) * hermite_vals(i_order)
                hermite_derivs(i_order + 1) = derivative_term + derivative_term
            end if
        end do

    end subroutine compute_hermite_polynomials_s

    !> Compute generalized Laguerre polynomials L^(alpha)_k(x) and their weighted derivatives
    !!
    !! Uses the recurrence relation:
    !!   L^(alpha)_{k+1}(x) = [(2k + alpha + 1 - x) L^(alpha)_k(x) - (k + alpha) L^(alpha)_{k-1}(x)] / (k + 1)
    !! The derivative output follows a modified form:
    !!   dxlag_{k+1} = (k+1) L^(alpha)_{k+1}(x) - (k + 1 + alpha) L^(alpha)_k(x)
    !!
    !! @note
    !! Original name: dlap
    !! @endnote
    pure subroutine compute_laguerre_polynomials_s(x_eval, laguerre_vals, laguerre_derivs, max_order, alpha_param)
        real(kind = rk), intent(in), value :: x_eval           !! Evaluation point
        real(kind = rk), intent(out) :: laguerre_vals(:)       !! L^(alpha)_k(x) for k = 0, ..., max_order (size >= max_order + 1)
        real(kind = rk), intent(out) :: laguerre_derivs(:)     !! Weighted derivative values (size >= max_order + 1)
        integer(kind = ik), intent(in), value :: max_order     !! Maximum polynomial order (n >= 0)
        integer(kind = ik), intent(in), value :: alpha_param   !! Laguerre parameter alpha (integer)

        integer(kind = ik) :: i_order
        real(kind = rk) :: recurrence_term
        real(kind = rk) :: alpha_real
        real(kind = rk) :: coeff_linear, coeff_prev

        alpha_real = real(alpha_param, rk)

        ! L^(alpha)_0(x) = 1
        laguerre_vals(1) = 1.0_rk
        laguerre_derivs(1) = 0.0_rk
        if (max_order <= 0_ik) return

        ! L^(alpha)_1(x) = 1 + alpha - x
        laguerre_vals(2) = 1.0_rk + alpha_real - x_eval
        laguerre_derivs(2) = -x_eval
        if (max_order <= 1_ik) return

        ! Recurrence relation for k = 2, ..., max_order
        do i_order = 2_ik, max_order
            ! Coefficient: (2k + alpha - 1 - x) where k = i_order
            coeff_linear = real(i_order + i_order + alpha_param - 1_ik, rk) - x_eval
            ! Coefficient: (k - 1 + alpha) where k = i_order
            coeff_prev = real(i_order - 1_ik + alpha_param, rk)

            recurrence_term = coeff_linear * laguerre_vals(i_order) - coeff_prev * laguerre_vals(i_order - 1)
            laguerre_vals(i_order + 1) = recurrence_term / real(i_order, rk)

            ! Weighted derivative formula
            laguerre_derivs(i_order + 1) = real(i_order, rk) * laguerre_vals(i_order + 1) &
                    - real(i_order + alpha_param, rk) * laguerre_vals(i_order)
        end do

    end subroutine compute_laguerre_polynomials_s

    !> Compute nodes and weights for Gauss-Hermite quadrature using Golub-Welsch algorithm
    !!
    !! Computes nodes and weights for Gauss-Hermite quadrature
    !! with weight function w(x) = e^{-x²} on (-∞, ∞).
    !!
    !! The Golub-Welsch algorithm constructs a symmetric tridiagonal Jacobi matrix
    !! from the three-term recurrence relation of Hermite polynomials:
    !!
    !!   H_{k+1}(x) = 2x H_k(x) - 2k H_{k-1}(x)
    !!
    !! For monic (probabilist's) normalization, the recurrence becomes:
    !!   P_{k+1}(x) = x P_k(x) - k P_{k-1}(x)
    !!
    !! The Jacobi matrix has:
    !!   - Diagonal: all zeros (Hermite polynomials are symmetric)
    !!   - Off-diagonal(k): sqrt(k)
    !!
    !! The quadrature approximates: ∫_{-∞}^{∞} f(x) e^{-x²} dx ≈ Σ w_i f(x_i)
    !!
    !! Complexity: O(n²) for eigenvalue computation
    !! Memory: O(n²) for eigenvector storage during computation
    !!
    !! @note
    !! Replaces the original hermit/hroot/hrecur implementation.
    !! @endnote
    pure subroutine compute_gauss_hermite_quadrature_s(n_points, nodes, weights, info)

        use mathematical_and_physical_constants_mod, only: SQRT_PI_C

        implicit none

        integer(kind = ik), intent(in) :: n_points          !! Number of quadrature points (n ≥ 1)
        real(kind = rk), intent(out) :: nodes(n_points)     !! Quadrature nodes, sorted ascending
        real(kind = rk), intent(out) :: weights(n_points)   !! Corresponding quadrature weights
        integer(kind = ik), intent(out) :: info             !! Status flag: 0 = success, >0 = QR iteration failed

        integer(kind = ik), parameter :: MAX_QR_ITER = 30
        real(kind = rk), parameter :: EPS = epsilon(1.0_rk)

        real(kind = rk) :: diagonal(n_points)
        real(kind = rk) :: off_diag(n_points)
        real(kind = rk) :: Q(n_points, n_points)
        real(kind = rk) :: mu_0, c, s, r, p, g, f, b, dd, shift, temp
        integer(kind = ik) :: i, j, k, l, m, iter

        if (n_points < 1) then
            info = -1
            return
        end if

        info = 0

        !---------------------------------------------------------------------------
        ! Build symmetric tridiagonal Jacobi matrix for Hermite polynomials
        !---------------------------------------------------------------------------
        ! For the "physicist's" Hermite polynomials with weight e^{-x²}:
        !   H_{k+1}(x) = 2x H_k(x) - 2k H_{k-1}(x)
        !
        ! The monic orthogonal polynomials satisfy:
        !   P_k(x) = x P_{k-1}(x) - b_k P_{k-2}(x)
        ! where b_k = k/2 for physicist's Hermite
        !
        ! Jacobi matrix:
        !   diagonal: 0 (all zeros due to symmetry of Hermite polynomials)
        !   off_diagonal(k) = sqrt(b_{k+1}) = sqrt((k+1)/2) = sqrt(k/2) for k=1..n-1

        diagonal = 0.0_rk

        off_diag(n_points) = 0.0_rk
        do k = 1, n_points - 1
            off_diag(k) = sqrt(real(k, rk) / 2.0_rk)
        end do

        !---------------------------------------------------------------------------
        ! Implicit QR algorithm with Wilkinson shift for symmetric tridiagonal matrix
        !---------------------------------------------------------------------------
        ! Initialize Q to identity matrix
        Q = 0.0_rk
        do i = 1, n_points
            Q(i, i) = 1.0_rk
        end do

        ! Process each eigenvalue
        do l = 1, n_points
            iter = 0

            qr_iteration: do
                ! Find smallest m such that off_diag(m) is negligible
                do m = l, n_points - 1
                    dd = abs(diagonal(m)) + abs(diagonal(m + 1))
                    ! Handle case where both diagonal elements are zero
                    if (dd < EPS) dd = 1.0_rk
                    if (abs(off_diag(m)) <= EPS * dd) exit
                end do

                ! Converged for eigenvalue l
                if (m == l) exit qr_iteration

                iter = iter + 1
                if (iter > MAX_QR_ITER) then
                    info = l  ! Failed to converge
                    return
                end if

                ! Wilkinson shift
                g = (diagonal(l + 1) - diagonal(l)) / (2.0_rk * off_diag(l))
                r = sqrt(g * g + 1.0_rk)
                if (g >= 0.0_rk) then
                    shift = diagonal(m) - diagonal(l) + off_diag(l) / (g + r)
                else
                    shift = diagonal(m) - diagonal(l) + off_diag(l) / (g - r)
                end if

                c = 1.0_rk
                s = 1.0_rk
                p = 0.0_rk

                ! Chase the bulge from bottom to top
                do i = m - 1, l, -1
                    f = s * off_diag(i)
                    b = c * off_diag(i)

                    ! Compute Givens rotation
                    if (abs(f) >= abs(shift)) then
                        c = shift / f
                        r = sqrt(c * c + 1.0_rk)
                        off_diag(i + 1) = f * r
                        s = 1.0_rk / r
                        c = c * s
                    else
                        s = f / shift
                        r = sqrt(s * s + 1.0_rk)
                        off_diag(i + 1) = shift * r
                        c = 1.0_rk / r
                        s = s * c
                    end if

                    ! Update diagonal elements
                    shift = diagonal(i + 1) - p
                    r = (diagonal(i) - shift) * s + 2.0_rk * c * b
                    p = s * r
                    diagonal(i + 1) = shift + p
                    shift = c * r - b

                    ! Accumulate rotation into Q
                    do j = 1, n_points
                        f = Q(j, i + 1)
                        Q(j, i + 1) = s * Q(j, i) + c * f
                        Q(j, i) = c * Q(j, i) - s * f
                    end do
                end do

                diagonal(l) = diagonal(l) - p
                off_diag(l) = shift
                off_diag(m) = 0.0_rk

            end do qr_iteration
        end do

        !---------------------------------------------------------------------------
        ! Sort eigenvalues and eigenvectors by ascending eigenvalue
        !---------------------------------------------------------------------------
        do i = 1, n_points - 1
            k = i
            p = diagonal(i)

            ! Find minimum in remaining elements
            do j = i + 1, n_points
                if (diagonal(j) < p) then
                    k = j
                    p = diagonal(j)
                end if
            end do

            ! Swap if necessary
            if (k /= i) then
                diagonal(k) = diagonal(i)
                diagonal(i) = p

                ! Swap corresponding eigenvector columns
                do j = 1, n_points
                    temp = Q(j, i)
                    Q(j, i) = Q(j, k)
                    Q(j, k) = temp
                end do
            end if
        end do

        !---------------------------------------------------------------------------
        ! Compute output
        !---------------------------------------------------------------------------
        nodes = diagonal

        ! Weights = μ_0 × (first component of eigenvector)²
        ! where μ_0 = ∫_{-∞}^{∞} e^{-x²} dx = √π
        mu_0 = SQRT_PI_C
        do k = 1, n_points
            weights(k) = mu_0 * Q(1, k)**2
        end do

    end subroutine compute_gauss_hermite_quadrature_s

    !> Compute nodes and weights for Gauss-Laguerre quadrature
    !!
    !! Uses Golub-Welsch algorithm: eigenvalues of the Jacobi matrix give the nodes,
    !! and weights are computed from the first components of the eigenvectors.
    !!
    !! The three-term recurrence for monic generalized Laguerre polynomials is:
    !!   P_k(x) = (x - a_k) P_{k-1}(x) - b_k P_{k-2}(x)
    !! where:
    !!   a_k = 2k + α - 1  (k = 1, 2, ..., n)
    !!   b_k = (k-1)(k-1+α) (k = 2, 3, ..., n)
    !!
    !! The symmetric Jacobi matrix has:
    !!   diagonal(k)     = a_k
    !!   off_diagonal(k) = sqrt(b_{k+1})
    !!
    !! Complexity: O(n²) for eigenvalue computation
    !! Memory: O(n²) for eigenvector storage during computation
    !!
    !! @note
    !! Replaces the original laguer/lgroot/lgrecr implementation.
    !! This approach eliminates the empirically-tuned "magic constants" required
    !! by Newton-Raphson based methods for initial root guesses.
    !! @endnote
    !!
    !! @note
    !! For very large n (> 1000), consider asymptotic methods for better performance.
    !! @endnote
    pure subroutine compute_gauss_laguerre_quadrature_s(n_points, alpha_param, nodes, weights, info)

        implicit none

        ! Arguments
        integer(kind = ik), intent(in) :: n_points          !! Number of quadrature points (n ≥ 1)
        real(kind = rk), intent(in) :: alpha_param          !! Laguerre parameter α (≥ 0)
        real(kind = rk), intent(out) :: nodes(n_points)     !! Quadrature nodes, sorted ascending
        real(kind = rk), intent(out) :: weights(n_points)   !! Corresponding quadrature weights
        integer(kind = ik), intent(out) :: info             !! Status flag: 0 = success, >0 = QR iteration failed

        integer(kind = ik), parameter :: MAX_QR_ITER = 30
        real(kind = rk), parameter :: EPS = epsilon(1.0_rk)

        real(kind = rk) :: diagonal(n_points)
        real(kind = rk) :: off_diag(n_points)
        real(kind = rk) :: Q(n_points, n_points)
        real(kind = rk) :: mu_0, c, s, r, p, g, f, b, dd, shift, temp
        integer(kind = ik) :: i, j, k, l, m, iter

        if (n_points < 1) then
            info = -1
            return
        end if

        info = 0

        !---------------------------------------------------------------------------
        ! Build symmetric tridiagonal Jacobi matrix
        !---------------------------------------------------------------------------
        ! Diagonal: a_k = 2k + α - 1
        do j = 1, n_points
            diagonal(j) = real(2 * j - 1, rk) + alpha_param
        end do

        ! Off-diagonal: sqrt(b_{k+1}) = sqrt(k(k+α))
        off_diag(n_points) = 0.0_rk
        do j = 1, n_points - 1
            off_diag(j) = sqrt(real(j, rk) * (real(j, rk) + alpha_param))
        end do

        !---------------------------------------------------------------------------
        ! Implicit QR algorithm with Wilkinson shift for symmetric tridiagonal matrix
        !---------------------------------------------------------------------------
        ! Initialize Q to identity matrix
        Q = 0.0_rk
        do i = 1, n_points
            Q(i, i) = 1.0_rk
        end do

        ! Process each eigenvalue
        do l = 1, n_points
            iter = 0

            qr_iteration: do
                ! Find smallest m such that off_diag(m) is negligible
                do m = l, n_points - 1
                    dd = abs(diagonal(m)) + abs(diagonal(m + 1))
                    if (abs(off_diag(m)) <= EPS * dd) exit
                end do

                ! Converged for eigenvalue l
                if (m == l) exit qr_iteration

                iter = iter + 1
                if (iter > MAX_QR_ITER) then
                    info = l  ! Failed to converge
                    return
                end if

                ! Wilkinson shift: eigenvalue of trailing 2×2 closer to diagonal(m)
                g = (diagonal(l + 1) - diagonal(l)) / (2.0_rk * off_diag(l))
                r = sqrt(g * g + 1.0_rk)
                if (g >= 0.0_rk) then
                    shift = diagonal(m) - diagonal(l) + off_diag(l) / (g + r)
                else
                    shift = diagonal(m) - diagonal(l) + off_diag(l) / (g - r)
                end if

                c = 1.0_rk
                s = 1.0_rk
                p = 0.0_rk

                ! Chase the bulge from bottom to top
                do i = m - 1, l, -1
                    f = s * off_diag(i)
                    b = c * off_diag(i)

                    ! Compute Givens rotation
                    if (abs(f) >= abs(shift)) then
                        c = shift / f
                        r = sqrt(c * c + 1.0_rk)
                        off_diag(i + 1) = f * r
                        s = 1.0_rk / r
                        c = c * s
                    else
                        s = f / shift
                        r = sqrt(s * s + 1.0_rk)
                        off_diag(i + 1) = shift * r
                        c = 1.0_rk / r
                        s = s * c
                    end if

                    ! Update diagonal elements
                    shift = diagonal(i + 1) - p
                    r = (diagonal(i) - shift) * s + 2.0_rk * c * b
                    p = s * r
                    diagonal(i + 1) = shift + p
                    shift = c * r - b

                    ! Accumulate rotation into Q
                    do j = 1, n_points
                        f = Q(j, i + 1)
                        Q(j, i + 1) = s * Q(j, i) + c * f
                        Q(j, i) = c * Q(j, i) - s * f
                    end do
                end do

                diagonal(l) = diagonal(l) - p
                off_diag(l) = shift
                off_diag(m) = 0.0_rk

            end do qr_iteration
        end do

        !---------------------------------------------------------------------------
        ! Sort eigenvalues and eigenvectors by ascending eigenvalue
        !---------------------------------------------------------------------------
        do i = 1, n_points - 1
            k = i
            p = diagonal(i)

            ! Find minimum in remaining elements
            do j = i + 1, n_points
                if (diagonal(j) < p) then
                    k = j
                    p = diagonal(j)
                end if
            end do

            ! Swap if necessary
            if (k /= i) then
                diagonal(k) = diagonal(i)
                diagonal(i) = p

                ! Swap corresponding eigenvector columns
                do j = 1, n_points
                    temp = Q(j, i)
                    Q(j, i) = Q(j, k)
                    Q(j, k) = temp
                end do
            end if
        end do

        !---------------------------------------------------------------------------
        ! Compute output
        !---------------------------------------------------------------------------
        nodes = diagonal

        ! Weights = Γ(α+1) × (first component of eigenvector)²
        mu_0 = gamma(alpha_param + 1.0_rk)
        do k = 1, n_points
            weights(k) = mu_0 * Q(1, k)**2
        end do

    end subroutine compute_gauss_laguerre_quadrature_s

    !> Computes Legendre polynomials and optionally their derivatives using recurrence
    !!
    !! Uses three-term recurrence relations for efficient computation:
    !!
    !! \[ P_n(x) = \frac{(2n-1)xP_{n-1}(x) - (n-1)P_{n-2}(x)}{n} \]
    !!
    !! \[ P'_n(x) = xP'_{n-1}(x) + nP_{n-1}(x) \]
    !!
    !! with initial conditions \( P_0(x) = 1 \), \( P_1(x) = x \).
    !!
    !! ## Array Indexing Convention
    !! Arrays use 1-based Fortran indexing where `p(i)` stores \( P_{i-1}(x) \):
    !!
    !! - `p(1)` contains \( P_0(x) \)
    !! - `p(2)` contains \( P_1(x) \)
    !! - `p(n)` contains \( P_{n-1}(x) \)
    !!
    !! ## Numerical Properties
    !! - Complexity: O(n) time, O(n) space
    !! - Numerically stable for \( |x| \leq 1 \)
    !! - No overflow protection for large n
    !!
    !! @warning
    !! For \( |x| > 1 \), polynomials may grow unbounded
    !! @endwarning
    !!
    !! @note
    !! This routine does not normalize the polynomials
    !! @endnote
    pure subroutine compute_legendre_polynomials_s(n, x, p, dp)

        implicit none

        ! Dummy arguments
        ! Input arguments
        integer(kind = ik), intent(in), value :: n          !! Order index: computes polynomials up to \( P_{n-1}(x) \) (must be n >= 1)
        real(kind = rk), intent(in), value :: x             !! Evaluation point, typically \( -1 \leq x \leq 1 \)
        ! Output arguments
        real(kind = rk), intent(out) :: p(:)                !! Array of length n storing \( P_0(x) \) through \( P_{n-1}(x) \)
        real(kind = rk), intent(out), optional :: dp(:)     !! Optional array of length n storing \( P'_0(x) \) through \( P'_{n-1}(x) \)

        ! Local variables
        integer(kind = ik) :: l                     !! Loop index
        real(kind = rk) :: xl_minus_1, xl_minus_2   !! Cache polynomial values
        real(kind = rk) :: dxl_minus_1              !! Cache derivative value
        real(kind = rk) :: inv_l_minus_1            !! Precompute division
        logical :: compute_derivatives              !! Flag for derivative computation

        ! Check if derivatives are requested
        compute_derivatives = present(dp)

        ! Base case: P_0(x) = 1, P'_0(x) = 0
        p(1) = 1.0_rk
        if (compute_derivatives) dp(1) = 0.0_rk
        if (n == 1_ik) return

        ! Base case: P_1(x) = x, P'_1(x) = 1
        p(2) = x
        if (compute_derivatives) dp(2) = 1.0_rk
        if (n == 2_ik) return

        ! Initialize cached values
        xl_minus_2 = 1.0_rk     ! P_0
        xl_minus_1 = x          ! P_1
        dxl_minus_1 = 1.0_rk    ! P'_1

        ! Recursive computation
        do l = 3, n
            ! Precompute division
            inv_l_minus_1 = 1.0_rk / real(l - 1_ik, rk)

            ! P_l = ((2l-3)*x*P_{l-1} - (l-2)*P_{l-2}) / (l-1)
            p(l) = (real(2_ik * l - 3_ik, rk) * x * xl_minus_1 - real(l - 2_ik, rk) * xl_minus_2) * inv_l_minus_1

            if (compute_derivatives) then
                ! P'_l = x*P'_{l-1} + (l-1)*P_{l-1}
                dp(l) = x * dxl_minus_1 + real(l - 1_ik, rk) * xl_minus_1
                dxl_minus_1 = dp(l)
            end if

            ! Update cached values for next iteration
            xl_minus_2 = xl_minus_1
            xl_minus_1 = p(l)
        end do

    end subroutine compute_legendre_polynomials_s

end module mathematical_utilities_mod