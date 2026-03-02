# Fortran-Foundations

A foundational Fortran library providing precision management, physical constants, and numerical algorithms for scientific computing. Designed as a self-contained base layer for nuclear structure and quantum mechanics codes.

## Modules

### `precision_utilities_mod`

Kind parameters and floating-point comparison utilities.

- **Kind parameters:** `ik` (int), `ikl` (int64), `rk` (double), `rk_single` (float), `rk_high` (long double) — all C-interoperable via `iso_c_binding`
- **Comparison functions:** `are_close_f` (robust relative/absolute tolerance with NaN/Inf handling in Debug builds), `is_zero_f` (tolerance-based zero check)

### `mathematical_and_physical_constants_mod`

Mathematical constants and 2022 CODATA physical constants.

- **Mathematical:** &pi;, 2&pi;, 4&pi;, &pi;/4, &pi;&sup2;, &radic;&pi;, 1/&radic;2, *e*
- **Physical:** &#x210F;c, fine-structure constant &alpha;, particle masses (proton, neutron, electron), atomic mass unit, nuclear radius constant r&#x2080;, Coulomb energy unit, mass excesses

### `mathematical_utilities_mod`

Numerical integration, special functions, and quadrature rules.

| Procedure | Description | Algorithm |
|---|---|---|
| `integrate_simpson_grid_f` | Definite integral on uniform grid | Composite Simpson's 1/3 rule |
| `complete_elliptic_integral_generalized_s` | Complete elliptic integrals K(k), E(k) | Arithmetic-geometric mean |
| `compute_hermite_polynomials_s` | Physicist's Hermite polynomials + derivatives | 3-term recurrence |
| `compute_laguerre_polynomials_s` | Generalized Laguerre polynomials + derivatives | 3-term recurrence |
| `compute_legendre_polynomials_s` | Legendre polynomials + derivatives | 3-term recurrence |
| `compute_gauss_legendre_quadrature_s` | Gauss-Legendre nodes and weights | Newton-Raphson root finding |
| `compute_gauss_hermite_quadrature_s` | Gauss-Hermite nodes and weights | Golub-Welsch + implicit QR |
| `compute_gauss_laguerre_quadrature_s` | Gauss-Laguerre nodes and weights | Golub-Welsch + implicit QR |
| `compute_spherical_harmonics_normalization_constants_s` | Axially-symmetric spherical harmonics normalization | Direct evaluation |
| `compute_strutinsky_curvature_coefficients_s` | Strutinsky shell-correction coefficients | Hermite polynomial expansion |

### `utility_procedures_mod`

Array sorting algorithms.

- `sort_array_ascending_s` — in-place heap sort (O(n log n))
- `sort_ascending_with_indices_s` — selection sort preserving two associated index arrays (for eigenvalue/quantum-number correspondence)

## Building

**Requirements:** CMake 3.20+, a Fortran compiler (gfortran recommended)

```bash
cmake -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build
```

Build types: `Debug`, `Release`, `RelWithDebInfo`. Debug builds enable IEEE 754 NaN/Inf validation checks.

The library builds as a static library (`fortran_foundations`). To use it in a downstream CMake project:

```cmake
add_subdirectory(Fortran-Foundations)
target_link_libraries(your_target PRIVATE FortranFoundations::fortran_foundations)
```

## Dependencies

- [GCC-Compiler-Options](https://github.com/AleksanderAugustyn/GCC-Compiler-Options) — automatically fetched by CMake at configure time

No external numerical libraries (BLAS, LAPACK, etc.) are required.

## Design

- Fortran 2018 standard (`.f08` extension)
- `implicit none` in every scope
- All procedures are `pure`
- Explicit `intent` on all arguments
- Module exports use `only:` imports
- C-interoperable kind parameters for potential mixed-language use
