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
      real (r8), parameter :: k_b = ??

      contains

      subroutine time_advance( time_step_seconds, temperature, pressure, &
        j_h2o2, &
        m, h2ovmr&
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

      real(r8), intent(in)    :: j_h2o2                 ! photodecomp rate of H2O2 -> 2*OH
      real(r8), intent(in)    :: m                      ! number density of air
      real(r8), intent(in)    :: HO2, OH, NO3, h2ovmr   ! fixed values

      real(r8), intent(in)    :: H2O2_0, SO2_0, H2SO4_0, DMS_0, HNO3_0 ! initial value
      real(r8), intent(out)   :: H2O2_t, SO2_t, H2SO4_t, DMS_t, HNO3_t ! final value

!-------------------------------------------------------
!       ... local variables
!-------------------------------------------------------

      real(r8) :: k(0:6)  ! rate constants for each reaction
      real(r8) :: alpha, beta, gamma, eta, zeta, kappa, phi

      ! rate constants
      k(0) = j_h2o2
      
      ko   = 3.0e-13 * exp ( 460. / temperature)
      kinf = 32.1e-33 * m * exp( 920. / temperature)
      fc = 1 + 1.4e-12 * m * h2ovmr * exp(2200. / temperature)
      k(1) = (k0+kinf) * fc

      k(2) = 1.8e-12  ! different from what is in the CAM pp_trop_mam = 2.9e-12_r8 * exp( -160._r8 * itemp)

      fc = 3.0e-31 *(300. / temperature)**3.3
      ko = fc*m / (1 + fc*m/1.5e-12)
      k(3) = ko * .6**(1 + (log10(fc*m/1.5e-12))**2)**(-1)

      k(4) = 9.6e-12 * exp (-3.23071866E-21/k_b / temperature) ! 9.6e-12_r8 * exp( -234._r8 * itemp)

      ko = 1._r8 + 5.5e-31_r8 * exp(7460/temperature) * m * 0.21
      k(5) = 1.7e-42 * exp(7810./temperature) * m(:,k) * 0.21 / ko

      k(6) = 1.9e-13_r8 * exp( 7.1793748E-21/k_b / temperature) ! 1.9e-13_r8 * exp( 520._r8 * itemp)


!-------------------------------------------------------
!       Solutions!!!
!-------------------------------------------------------

      ! DMS(t) given DMS_0 = DMS(t=0)
      gamma = -k(4)* OH - k(5)*OH-k(6)
      DMS = DMS_0 * exp( gamma * time_step_seconds)

      ! SO2
      alpha = -k(3) * OH
      beta = k(4) + 0.5*k(5)*OH + k(6)*NO3 
      SO2_t = SO2_0 &
        + DMS_0 *beta * ( exp(gamma*time_step_seconds) - exp( alpha * time_step_seconds ) &
                      / (gamma - alpha )

      ! H2SO4
      kappa = k(3) * OH
      H2SO4_t = H2SO4_0 &
        + kappa * SO2_0 * time_step_seconds &
        + DMS_0 * (beta / (gamma - alpha) ) &
          * ( &
               ( exp( gamma * time_step_seconds ) - 1) / gamma &
              -( exp( alpha * time_step_seconds ) - 1) / alpha &
            )

      ! HNO3
      eta = k(6) * OH
      HNO3_t = HNO3_0 &
        + DMS_0 * eta * ( (exp(gamma*time_step_seconds) - 1) / gamma )
   
      ! H2O2
      zeta = - k(0) - k(2) * OH
      phi = k(1) * HO2 * HO2
      H2O2_t = (H2O2_0 + phi/zeta ) &
               * exp(zeta * time_step_seconds) &
               - phi/zeta


      end subroutine time_advance

end module trop_mam4_analytical
