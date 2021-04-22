      module trop_mam4_analytical
!-------------------------------------------------------
!  k0  [jh2o2]       H2O2 + hv ->
!  k1  [usr_HO2_HO2] HO2 + HO2 -> H2O2
!  k2                H2O2 + OH -> H2O + HO2                                           ; 2.9e-12, -160
!  k3  [usr_SO2_OH]  SO2 + OH -> H2SO4
!  k4                DMS + OH -> SO2                                                  ; 9.6e-12, -234.
!  k5  [usr_DMS_OH]  DMS + OH -> .5 * SO2 + .5 * HO2
!  k6                DMS + NO3 -> SO2 + HNO3                                          ; 1.9e-13,  520.
!
!  gas_array = 'H2O2', HO2, SO2, H2SO4, DMS, NO3, HNO3
!-------------------------------------------------------

      use shr_kind_mod, only : r8 => shr_kind_r8

      private
      public :: solution
      public :: set_rates

      contains

      subroutine set_rates(temperature, pressure, rate_constant )

      implicit none

!-------------------------------------------------------
!       ... dummy arguments
!-------------------------------------------------------
      real(r8), intent(in)    :: temperature
      real(r8), intent(in)    :: pressure
      real(r8), intent(out) :: rate_constant(0:6)

!-------------------------------------------------------
!       ... local variables
!-------------------------------------------------------
      itemp = 1/temperature
      (0,1,3,5) are so far undefined, but need to be added here?
      rate_constant(2) = 2.9e-12_r8 * exp( -160._r8 * itemp(:,:) )
      rate_constant(4) = 9.6e-12_r8 * exp( -234._r8 * itemp(:,:) )
      rate_constant(6) = 1.9e-13_r8 * exp( 520._r8 * itemp(:,:) )

      end subroutine setrxt


      subroutine time_advance( rate_constant, gas_array )

      implicit none

!-------------------------------------------------------
!       ... dummy arguments
!-------------------------------------------------------
      real(r8), intent(in)    :: rate_constant
      real(r8), intent(inout) :: gas_array
      real(r8), intent(in)    :: time_step_seconds

!-------------------------------------------------------
!       ... local variables
!-------------------------------------------------------

!  gas_array = 'H2O2', HO2, SO2, H2SO4, DMS, NO3, HNO3
      ! DMS(t) given DMS_0 = DMS(t=0)
      c = -rate_constant(4)* OH - rate_constant(5)*OH-rate_constant(6)
      DMS = DMS_0 * exp( c * time_step_seconds)

      ! SO2
      alpha = -rate_constant(3) * OH
      beta = rate_constant(4) + 0.5*rate_constant(5)*OH + rate_constant(6)*NO3
      ratio = beta/(c-alpha)
      coeff_1 = (SO2_0 - ratio) exp(alpha*time_step_seconds) + ratio*exp(c * time_step_seconds)
      

      end subroutine time_advance

      end module trop_mam4_analytical
