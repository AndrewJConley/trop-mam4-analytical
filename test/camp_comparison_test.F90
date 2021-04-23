!> \file
!> Comparison of CAMP solver output with the analytical solution for trop-mam4

!> CAMP comparison module
program camp_comparison_test

  use chemical_state,                  only : chemical_state_t
  use environment,                     only : environment_t
  use photolysis_rate_constants,       only : photolysis_rate_constants_t
  use shr_kind_mod,                    only : r8 => shr_kind_r8
  use trop_mam4_analytical,            only : advance_chemical_system

  implicit none

  real(r8),         parameter       :: time_step__s = 60.0_r8
  integer,          parameter       :: number_of_time_steps = 60 * 24
  integer,          parameter       :: file_unit = 10
  character(len=*), parameter       :: file_name = "out/analytic_solution.csv"
  type(photolysis_rate_constants_t) :: photo_rate_constants
  type(environment_t)               :: environment
  type(chemical_state_t)            :: initial_state, final_state
  integer                           :: i_time_step

  environment%temperature__K_     = 298.15_r8
  environment%pressure__Pa_       = 101325.0_r8

  photo_rate_constants%H2O2__s_   = 2.0e-5_r8

  call initial_state%set_ideal_M( environment )
  initial_state%H2O__molec_cm3_   = 1.2e16_r8
  initial_state%HO2__molec_cm3_   = 1.0e8_r8
  initial_state%OH__molec_cm3_    = 1.0e5_r8
  initial_state%NO3__molec_cm3_   = 3.61e8_r8
  initial_state%H2O2__molec_cm3_  = 5.42e14_r8
  initial_state%SO2__molec_cm3_   = 1.0e5_r8
  initial_state%H2SO4__molec_cm3_ = 1.0e5_r8
  initial_state%DMS__molec_cm3_   = 3.6e15_r8
  initial_state%HNO3__molec_cm3_  = 1.0e5_r8

  open( unit = file_unit, file = file_name, status = "replace",               &
        action = "write" )
  write( file_unit, * ) "time,H2O2,SO2,H2SO4,DMS,HNO3"
  call output_state( 0.0_r8, initial_state )
  do i_time_step = 1, number_of_time_steps
    call advance_chemical_system( time_step__s, environment,                  &
                                  photo_rate_constants,                       &
                                  initial_state, final_state )
    call output_state( i_time_step * time_step__s, final_state )
    initial_state = final_state
  end do
  close( file_unit )
  write(*,*) "Passed!"

contains

  subroutine output_state( time__s, state )
    use pmc_util,                      only : to_string
    real(r8),                intent(in) :: time__s
    class(chemical_state_t), intent(in) :: state
    write( file_unit, * ) trim( to_string( time__s          ) )//","//        &
                          trim( to_string( state%H2O2__molec_cm3_  ) )//","// &
                          trim( to_string( state%SO2__molec_cm3_   ) )//","// &
                          trim( to_string( state%H2SO4__molec_cm3_ ) )//","// &
                          trim( to_string( state%DMS__molec_cm3_   ) )//","// &
                          trim( to_string( state%HNO3__molec_cm3_  ) )
  end subroutine output_state

end program camp_comparison_test
