module photolysis_rate_constants

  use shr_kind_mod,                    only: r8 => shr_kind_r8

  implicit none

  private
  public :: photolysis_rate_constants_t

  ! photolysis rate constants
  type :: photolysis_rate_constants_t
    real(r8) :: H2O2__s_          = 0.0_r8
  end type photolysis_rate_constants_t

end module photolysis_rate_constants
