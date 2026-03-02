!> Module for defining key mathematical and physical constants.
!!
!! This module centralizes the definitions of important constants used throughout
!! the calculations, ensuring consistency and ease of maintenance.
!!
!! # References
!!
!! - Eite Tiesinga, Peter J. Mohr, David B. Newell, and Barry N. Taylor (2024), "The 2022 CODATA Recommended Values of the Fundamental Physical Constants" (Web Version 9.0). Database developed by J. Baker, M. Douma, and S. Kotochigova.
!! Available at https://physics.nist.gov/constants, National Institute of Standards and Technology, Gaithersburg, MD 20899.
module mathematical_and_physical_constants_mod

    use precision_utilities_mod, only: rk

    implicit none

    private

    ! Mathematical constants
    real(kind = rk), parameter, public :: PI_C = 3.141592653589793_rk               !! Value of \( \pi \)
    real(kind = rk), parameter, public :: TWO_PI_C = 2.0_rk * PI_C                  !! Value of \( 2\pi \)
    real(kind = rk), parameter, public :: FOUR_PI_C = 4.0_rk * PI_C                 !! Value of \( 4\pi \)
    real(kind = rk), parameter, public :: PI_OVER_FOUR_C = PI_C / 4.0_rk            !! Value of \( \pi/4 \)
    real(kind = rk), parameter, public :: PI_SQUARED_C = PI_C**2                    !! Value of \( \pi^2 \)
    real(kind = rk), parameter, public :: SQRT_PI_C = sqrt(PI_C)                    !! Value of \( \sqrt{\pi} \)
    real(kind = rk), parameter, public :: SQRT_TWO_INV_C = 1.0_rk / sqrt(2.0_rk)    !! Value of \( 1/\sqrt{2} \)
    real(kind = rk), parameter, public :: EULER_E_C = 2.718281828459045_rk          !! Euler's number \( e \)

    ! Physical constants
    real(kind = rk), parameter, public :: HBARC_C = 197.3269804_rk                                                  !! \( \hbar c \) in MeV·fm
    real(kind = rk), parameter, public :: ALPHA_FINE_C = 7.2973525643e-3_rk                                         !! Fine-structure constant \( \alpha \)
    real(kind = rk), parameter, public :: ALPHA_FINE_INV_C = 1.0_rk / ALPHA_FINE_C                                  !! Inverse fine-structure constant \( 1/\alpha \)
    real(kind = rk), parameter, public :: AMU_MEV_C = 931.49410372_rk                                               !! Atomic mass unit in MeV/c^2
    real(kind = rk), parameter, public :: MASS_PROTON_MEV_C = 938.27208943_rk                                       !! Proton mass in MeV/c^2
    real(kind = rk), parameter, public :: MASS_NEUTRON_MEV_C = 939.56542194_rk                                      !! Neutron mass in MeV/c^2
    real(kind = rk), parameter, public :: MASS_ELECTRON_MEV_C = 0.51099895069_rk                                    !! Electron mass in MeV/c^2
    real(kind = rk), parameter, public :: MASS_NUCLEON_MEV_C = (MASS_PROTON_MEV_C + MASS_NEUTRON_MEV_C) / 2.0_rk    !! Average nucleon mass in MeV/c^2
    real(kind = rk), parameter, public :: MASS_HYDROGEN_MEV_C = MASS_PROTON_MEV_C + MASS_ELECTRON_MEV_C             !! Hydrogen atom mass in MeV/c^2
    real(kind = rk), parameter, public :: R0_CONSTANT_C = 1.16_rk                                                   !! Nuclear radius constant in fm
    real(kind = rk), parameter, public :: COULOMB_ENERGY_UNIT_C = ALPHA_FINE_C * HBARC_C                            !! \( e^2/(4\pi\epsilon_0) \) in MeV·fm
    real(kind = rk), parameter, public :: NUCLEON_SPIN_C = 0.5_rk                                                   !! Intrinsic spin quantum number for nucleons (s = 1/2).
    real(kind = rk), parameter, public :: NEUTRON_MASS_EXCESS_MEV_C = MASS_NEUTRON_MEV_C - AMU_MEV_C                !! Neutron mass excess in MeV
    real(kind = rk), parameter, public :: PROTON_MASS_EXCESS_MEV_C = MASS_PROTON_MEV_C - AMU_MEV_C                  !! Proton mass excess in MeV
    real(kind = rk), parameter, public :: HYDROGEN_MASS_EXCESS_MEV_C = MASS_HYDROGEN_MEV_C - AMU_MEV_C              !! Hydrogen mass excess in MeV

end module mathematical_and_physical_constants_mod