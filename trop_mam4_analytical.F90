!> \file
!> Analytic solution for trop-mam4

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!> \mainpage trop-mam4 Analytic Solution
!!
!! The trop-mam4 gas-phase mechanism consists of 7 reactions:
!!
!! \f[
!!   \ce{H2O2 + hv  ->[k_0]}
!! \f]\f[
!!   \ce{HO2 + HO2 ->[k_1] H2O2}
!! \f]\f[
!!   \ce{H2O2 + OH ->[k_2] H2O + HO2}
!! \f]\f[
!!   \ce{SO2 + OH  ->[k_3] H2SO4}
!! \f]\f[
!!   \ce{DMS + OH  ->[k_4] SO2}
!! \f]\f[
!!   \ce{DMS + OH  ->[k_5]} 0.5 \, \ce{SO2 +} 0.5 \, \ce{HO2}
!! \f]\f[
!!   \ce{DMS + NO3 ->[k_6] SO2 + HNO3}
!! \f]
!! with the concentrations of the following chemicals are specified elsewhere and unprognosed by the gas-phase chemistry step.
!! \f[ \ce{HO2, OH, NO3} \f]
!! In this solution, we ignore the resulting tendency for water vapor, \f$ \ce{H2O} \f$ as well.
!!
!! Rate Constants
!!---------------
!! The rate constant \f$ k_0 \f$ is a time-evolving value specified by the model and considered to be a constant during the chemistry time step.
!!
!! The rate constant \f$ k_1 \f$ corresponds to a sum of Arrhenius reactions
!! \f{eqnarray*}{
!! k_{\alpha} &=& (3\cdot 10^{-13}) \cdot e^{460/T} \\
!! k_{\beta} &=& (2 \cdot 10^{-33}) \cdot  \ce{M} e^{920/T} \\
!! F &=& 1.4 + (2 \cdot 10^{21}) \cdot \ce{H2O} e^{2200/T} \\
!! k_1 &=& (k_{\alpha} + k_{\beta}) \cdot F
!! \f}
!!
!! The rate constant \f$ k_2 \f$ corresponds to an Arrhenius reactions
!! \f[
!! k_2 = 1.8 \cdot 10^{-12}
!! \f]
!!
!! The rate constant \f$ k_3 \f$ corresponds to a Troe reaction
!! \f{eqnarray*}{
!! f_c &=& (3\cdot 10^{-31}) \left({\frac{300}{T}}\right)^{3.3} \\
!! f_1 &=& \frac{f_c  \cdot  \ce{M} }{ (1 + f_c \ce{M} / (1.5 \cdot 10^{-12}))}\\
!! e_{fac} &=& \left({\mbox{log}_{10}(f_c \cdot \ce{M}/{(1.5 \cdot 10^{-12})}}\right)^{2} \\
!! k_3 &=& f_1 \cdot (0.6)^{ \left(1 + {e_{fac}}^{2}   \right)^{-1}  }
!! \f}
!! 
!! The rate constant \f$ k_4 \f$ corresponds to an Arrhenius reaction
!! \f[
!! k_4 = (9.6\cdot 10^{-12}) \cdot e^{(3.23071866\cdot 10^{-21})/(k_b T)}
!! \f]
!!
!! The rate constant \f$ k_5 \f$ is specfied as follows
!! \f{eqnarray*}{
!! k_{high} &=& 1 + (5.5 \cdot 10^{-31}) \cdot e^{7460/T} \cdot \ce{M} (0.21) \\
!! k_5 &=& (1.7\cdot 10^{-42}) \cdot e^{(7810/T)} \cdot \ce{M} (0.21) / k_{high}
!! \f}
!!
!! The rate constant \f$ k_6 \f$ is a typical Arrhenius reaction
!!\f[
!! k_6 = (1.9\cdot 10^{-13}) \cdot e^{(-7.1793748 \cdot 10^{-21} / (k_b T) )}
!!\f]
!!
!! Solution
!!---------
!! The solution of the resulting differential equations can be solved analytically.  In each of the following expressions, the greek variables are constants.
!!
!! DMS
!!----
!! The time evolution of \f$ \ce{DMS}\f$ is governed by the differential equation
!! \f{eqnarray*}{
!! \frac{d \, \ce{DMS}} {dt} &=& \gamma \ce{DMS} \\
!! \gamma &=& -k_4 \ce{OH} -k_5 \ce{OH} -k_6 \ce{NO3} \\
!! \rightarrow \ce{DMS} (t) &=& \ce{DMS}(0) \,\, e^{\gamma t}
!! \f}
!!
!! HNO3
!!----
!! The time evolution of \f$ \ce{HNO3}\f$ is governed by the differential equation
!! \f{eqnarray*}{
!! \frac{d \, \ce{HNO3}} {dt} &=& \eta \,\ce{DMS} \\
!! \eta &=& k_6 \ce{NO3} \\
!! \rightarrow \ce{HNO3} (t) &=& \ce{HNO3}(0) + \frac{\eta \, \ce{DMS}(0) }{\gamma} \left( e^{\gamma t} - 1 \right) \\
!! \f}
!!
!! SO2
!!----
!! The time evolution of \f$ \ce{SO2}\f$ is governed by the differential equation
!! \f{eqnarray*}{
!! \frac{d \, \ce{SO2}} {dt} &=& \alpha \, \ce{SO2} + \beta \, \ce{DMS} \\
!! \alpha &=& -k_3 \ce{OH} \\
!! \beta &=& k_4 \ce{OH} + 0.5 \, k_5 \ce{OH} + k_6 \ce {NO3} \\
!! \rightarrow \ce{SO2} (t) &=& \ce{SO2}(0)\, e^{\alpha t} + \frac{\beta} {\gamma-\alpha} \ce{DMS}(0) \left( e^{\gamma t} - e^{\alpha t} \right) \\
!! \f}
!!
!! H2SO4
!!----
!! The time evolution of \f$ \ce{H2SO4}\f$ is governed by the differential equation
!! \f{eqnarray*}{
!! \frac{d \, \ce{H2SO4}} {dt} &=& \kappa \, \ce{SO2} \\
!! &=& \kappa \, \left( \ce{SO2}(0)\, e^{\alpha t} + \frac{\beta} {\gamma-\alpha} \ce{DMS}(0) \left( e^{\gamma t} - e^{\alpha t} \right)   \right) \\
!! \kappa &=& k_3 \ce{OH} \\
!! \rightarrow \ce{H2SO4} (t) &=& \ce{H2SO4}(0)  \\
!!     && + \kappa \, \ce{SO2}(0) \left( \frac{e^{\alpha t}-1}{\alpha}  \right) \\
!!     && + \kappa \frac{\beta}{\gamma - \alpha} \ce{DMS}(0) \left( \frac{e^{\gamma t} -1}{\gamma} - \frac{e^{\alpha t}-1}{\alpha}    \right) 
!! \f}
!!
!! H2O2
!!----
!! The time evolution of \f$ \ce{H2O2}\f$ is governed by the differential equation
!! \f{eqnarray*}{
!! \frac{d \, \ce{H2O2}} {dt} &=& \zeta \, \ce{H2O2} + \phi \\
!! \zeta &=& -k_0 -k_2 \ce{OH} \\
!! \phi &=& k_1 \ce{HO2} \ce{HO2} \\
!! \rightarrow \ce{H2O2} (t) &=& \left( \ce{H2O2}(0) + \frac{\phi} {\zeta}\right) e^{\zeta t} - \frac{\phi}{\zeta}
!! \f}



module trop_mam4_analytical

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

      ko = 3.0e-13_r8*exp(460._r8/environment%temperature__K_)
      kinf = 2.1e-33_r8*initial_state%M__molec_cm3_*exp(920._r8/environment%temperature__K_)
      fc = 1.0_r8 + 1.4e-21_r8*initial_state%H2O__molec_cm3_*exp(2200._r8/environment%temperature__K_)
      k(1) = (ko + kinf)*fc

      k(2) = 1.8e-12  ! different from what is in the CAM pp_trop_mam = 2.9e-12_r8 * exp( -160._r8 * itemp)

      fc = 3.0e-31*(300./environment%temperature__K_)**3.3
      ko = fc*initial_state%M__molec_cm3_/(1 + fc*initial_state%M__molec_cm3_/1.5e-12)
      k(3) = ko*.6**(1 + (log10(fc*initial_state%M__molec_cm3_/1.5e-12))**2)**(-1)

      k(4) = 9.6e-12*exp(3.23071866E-21/(k_b*environment%temperature__K_)) ! 9.6e-12_r8 * exp( -234._r8 * itemp)

      ko = 1._r8 + 5.5e-31_r8*exp(7460/environment%temperature__K_)*initial_state%M__molec_cm3_*0.21
      k(5) = 1.7e-42*exp(7810./environment%temperature__K_)*initial_state%M__molec_cm3_*0.21/ko

      k(6) = 1.9e-13_r8*exp(-7.1793748E-21/(k_b*environment%temperature__K_)) ! 1.9e-13_r8 * exp( 520._r8 * itemp)

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
