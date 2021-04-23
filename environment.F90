module environment

  use shr_kind_mod,                    only: r8 => shr_kind_r8

  implicit none

  private
  public :: environment_t

  ! environmental conditions
  type :: environment_t
    real(r8) :: temperature__K_ = 298.15_r8
    real(r8) :: pressure__Pa_   = 101325.0_r8
  end type environment_t

end module environment
