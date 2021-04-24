module trop_mam4_analytical

!-------------------------------------------------------
!  analytical solution of the following system
!  k0  [jh2o2]       H2O2 + hv ->
!  k1  [usr_HO2_HO2] HO2 + HO2 -> H2O2
!  k2                H2O2 + OH -> H2O + HO2                                           ; 2.9e-12, -160
!  k3  [usr_SO2_OH]  SO2 + OH -> H2SO4
!  k4                DMS + OH -> SO2                                                  ; 9.6e-12, -234.
!  k5  [usr_DMS_OH]  DMS + OH -> .5 * SO2 + .5 * HO2
!  k6                DMS + NO3 -> SO2 + HNO3                                          ; 1.9e-13,  520.
!
!  Fixed (constant) values for HO2, OH, NO3
!-------------------------------------------------------

   use chemical_state,            only : chemical_state_t
   use environment,               only : environment_t
   use photolysis_rate_constants, only : photolysis_rate_constants_t
   use shr_kind_mod,              only: r8 => shr_kind_r8

   implicit none

   private
   public :: advance_chemical_system
   real(r8), parameter :: k_b = 1.380649e-23 ! J/K, should be specified elsewhere

contains

   subroutine advance_chemical_system(time_step_seconds, environment, &
                                      photo_rate_constants, &
                                      initial_state, final_state)
      implicit none

!-------------------------------------------------------
!       ... dummy arguments
!-------------------------------------------------------
      real(r8),                           intent(in)  :: time_step_seconds
      class(environment_t),               intent(in)  :: environment
      class(photolysis_rate_constants_t), intent(in)  :: photo_rate_constants
      class(chemical_state_t),            intent(in)  :: initial_state
      class(chemical_state_t),            intent(out) :: final_state

!-------------------------------------------------------
!       ... local variables
!-------------------------------------------------------

      real(r8) :: k(0:6)  ! rate constants for each reaction
      real(r8) :: alpha, beta, gamma, eta, zeta, kappa, phi  ! matrix entries for this linear system
      real(r8) :: ko, kinf, fc   ! temporary parameters for rate constants

      ! rate constants
      k(0) = photo_rate_constants%H2O2__s_

      ko = 3.0e-13*exp(460./environment%temperature__K_)
      kinf = 32.1e-33*initial_state%M__molec_cm3_*exp(920./environment%temperature__K_)
      fc = 1 + 1.4e-12*initial_state%H2O__molec_cm3_*exp(2200./environment%temperature__K_)
      k(1) = (ko + kinf)*fc

      k(2) = 1.8e-12  ! different from what is in the CAM pp_trop_mam = 2.9e-12_r8 * exp( -160._r8 * itemp)

      fc = 3.0e-31*(300./environment%temperature__K_)**3.3
      ko = fc*initial_state%M__molec_cm3_/(1 + fc*initial_state%M__molec_cm3_/1.5e-12)
      k(3) = ko*.6**(1 + (log10(fc*initial_state%M__molec_cm3_/1.5e-12))**2)**(-1)

      k(4) = 9.6e-12*exp(-3.23071866E-21/(k_b*environment%temperature__K_)) ! 9.6e-12_r8 * exp( -234._r8 * itemp)

      ko = 1._r8 + 5.5e-31_r8*exp(7460/environment%temperature__K_)*initial_state%M__molec_cm3_*0.21
      k(5) = 1.7e-42*exp(7810./environment%temperature__K_)*initial_state%M__molec_cm3_*0.21/ko

      k(6) = 1.9e-13_r8*exp(7.1793748E-21/(k_b*environment%temperature__K_)) ! 1.9e-13_r8 * exp( 520._r8 * itemp)

!-------------------------------------------------------
!     Analytical Solution
!-------------------------------------------------------

      ! unsolved concentrations remain unchanged
      final_state%M__molec_cm3_   = initial_state%M__molec_cm3_
      final_state%H2O__molec_cm3_ = initial_state%H2O__molec_cm3_
      final_state%HO2__molec_cm3_ = initial_state%HO2__molec_cm3_
      final_state%OH__molec_cm3_  = initial_state%OH__molec_cm3_
      final_state%NO3__molec_cm3_ = initial_state%NO3__molec_cm3_

      ! DMS
      ! d(DMS)/dt = gamma*DMS
      gamma = -k(4)*initial_state%OH__molec_cm3_ - k(5)*initial_state%OH__molec_cm3_ - k(6)*initial_state%NO3__molec_cm3_
      final_state%DMS__molec_cm3_ = initial_state%DMS__molec_cm3_*exp(gamma*time_step_seconds)

      ! SO2
      alpha = -k(3)*initial_state%OH__molec_cm3_
      beta = k(4)*initial_state%OH__molec_cm3_ + 0.5*k(5)*initial_state%OH__molec_cm3_ + k(6)*initial_state%NO3__molec_cm3_
      !  if you want a sulfur-conserving mechanism
      !  beta = k(4)*initial_state%OH__molec_cm3_ + 1.0*k(5)*initial_state%OH__molec_cm3_ + k(6)*initial_state%NO3__molec_cm3_
      final_state%SO2__molec_cm3_ = initial_state%SO2__molec_cm3_*exp(alpha*time_step_seconds) &
              + initial_state%DMS__molec_cm3_*beta*(exp(gamma*time_step_seconds) - &
                                                    exp(alpha*time_step_seconds)) &
                                                   /(gamma - alpha)

      ! H2SO4
      kappa = k(3)*initial_state%OH__molec_cm3_
      final_state%H2SO4__molec_cm3_ = initial_state%H2SO4__molec_cm3_ &
                + kappa*initial_state%SO2__molec_cm3_*( (exp(alpha * time_step_seconds) - 1)/ alpha) &
                + kappa*initial_state%DMS__molec_cm3_*(beta/(gamma - alpha)) &
                *( &
                (exp(gamma*time_step_seconds) - 1)/gamma &
                - (exp(alpha*time_step_seconds) - 1)/alpha &
                )

      ! HNO3
      eta = k(6)*initial_state%NO3__molec_cm3_
      final_state%HNO3__molec_cm3_ = initial_state%HNO3__molec_cm3_ &
               + initial_state%DMS__molec_cm3_*eta*((exp(gamma*time_step_seconds) - 1)/gamma)

      ! H2O2
      zeta = -k(0) - k(2)*initial_state%OH__molec_cm3_
      phi = k(1)*initial_state%HO2__molec_cm3_*initial_state%HO2__molec_cm3_
      final_state%H2O2__molec_cm3_ = (initial_state%H2O2__molec_cm3_ + phi/zeta) &
               *exp(zeta*time_step_seconds) &
               - phi/zeta

   end subroutine advance_chemical_system

end module trop_mam4_analytical
