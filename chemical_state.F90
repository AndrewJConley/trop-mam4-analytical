module chemical_state

   use shr_kind_mod,                   only: r8 => shr_kind_r8

   implicit none

   private
   public :: chemical_state_t
   real(r8), parameter :: k_b = 1.380649e-23 ! J/K, should be specified elsewhere

   ! chemical state
   type :: chemical_state_t
     real(r8) :: M__molec_cm3_     = 0.0_r8
     real(r8) :: H2O__molec_cm3_   = 0.0_r8
     real(r8) :: HO2__molec_cm3_   = 0.0_r8
     real(r8) :: OH__molec_cm3_    = 0.0_r8
     real(r8) :: NO3__molec_cm3_   = 0.0_r8
     real(r8) :: H2O2__molec_cm3_  = 0.0_r8
     real(r8) :: SO2__molec_cm3_   = 0.0_r8
     real(r8) :: H2SO4__molec_cm3_ = 0.0_r8
     real(r8) :: DMS__molec_cm3_   = 0.0_r8
     real(r8) :: HNO3__molec_cm3_  = 0.0_r8
   contains
     procedure :: set_ideal_M
   end type chemical_state_t

contains

   subroutine set_ideal_M( this, environment )
     use environment,                  only : environment_t
     class(chemical_state_t), intent(inout) :: this
     class(environment_t),    intent(in)    :: environment
     this%M__molec_cm3_ = environment%pressure__Pa_ &
                          / ( k_b * environment%temperature__K_ ) * 1e-6
   end subroutine set_ideal_M

end module chemical_state
