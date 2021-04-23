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

      use shr_kind_mod, only : r8 => shr_kind_r8

      private
      public :: solution
      public :: set_rates

      contains

      subroutine time_advance( time_step_seconds, temperature, pressure, &
        j_h2o2, &
        HO2, OH, NO3, &
        H2O2_t, H2O2_0, &
        SO2_t, SO2_0, &
        H2SO4_t, H2SO4_0, &
        DMS_t, DMS_0, &
        HNO3_t, HNO3_0 &
        )
      implicit none

!-------------------------------------------------------
!       ... dummy arguments
!-------------------------------------------------------
      real(r8), intent(in)    :: time_step_seconds
      real(r8), intent(in)    :: j_h2o2         ! photodecomp rate of H2O2 -> 2*OH
      real(r8), intent(in)    :: HO2, OH, NO3   ! fixed values
      real(r8), intent(in)    :: H2O2_0, SO2_0, H2SO4_0, DMS_0, HNO3_0 ! initial value
      real(r8), intent(out)   :: H2O2_t, SO2_t, H2SO4_t, DMS_t, HNO3_t ! final value
      real(r8) :: rate_constant(0:6)


!-------------------------------------------------------
!       ... local variables
!-------------------------------------------------------

      itemp = 1/temperature
      (1,3,5) are so far undefined, but need to be added here?
      rate_constant(0) = j_h2o2
      rate_constant(1) = 
      rate_constant(2) = 2.9e-12_r8 * exp( -160._r8 * itemp)
      rate_constant(3) = 
      rate_constant(4) = 9.6e-12_r8 * exp( -234._r8 * itemp)
      rate_constant(5) = 
      rate_constant(6) = 1.9e-13_r8 * exp( 520._r8 * itemp)


      ! DMS(t) given DMS_0 = DMS(t=0)
      C = -rate_constant(4)* OH - rate_constant(5)*OH-rate_constant(6)
      DMS = DMS_0 * exp( C * time_step_seconds)

      ! SO2
      A = -rate_constant(3) * OH
      B = rate_constant(4) + 0.5*rate_constant(5)*OH + rate_constant(6)*NO3
      beta = B * DMS_0
      alpha = A
      gamma = C
      SO2_t = SO2_0 &
        + beta * ( exp(gamma*time_step_seconds) - exp( alpha * time_step_seconds ) &
          / (gamma - alpha )

      ! H2SO4
      kappa = rate_constant(3) * OH
      H2SO4_t = H2SO4_0 &
        + kappa * SO2_0 * time_step_seconds &
        + (beta / (gamma - alpha) ) * ( &
             ( exp( gamma * time_step_seconds ) - 1) / gamma &
            -(exp( alpha * time_step_seconds ) - 1) / alpha &
          )

      ! HNO3
      eta = rate_constant(6) * OH
      HNO3_t = HNO3_0 &
        + eta * DMS_0 * ( (exp(gamma*time_step_seconds) - 1) / gamma )
   
      ! H2O2
      zeta = - rate_constant(0) - rate_constant(2) * OH
      phi = rate_constant(1) * HO2 * HO2
      H2O2_t = (H2O2_0 + phi/zeta ) *exp(zeta * time_step_seconds) &
        - phi/zeta


      end subroutine time_advance

end module trop_mam4_analytical
