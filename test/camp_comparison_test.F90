!> \file
!> Comparison of CAMP solver output with the analytical solution for trop-mam4

!> CAMP comparison module
program camp_comparison_test

  use camp,                            only : camp_t
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
  integer, parameter                :: M     = 1
  integer, parameter                :: H2O   = 2
  integer, parameter                :: HO2   = 3
  integer, parameter                :: OH    = 4
  integer, parameter                :: NO3   = 5
  integer, parameter                :: H2O2  = 6
  integer, parameter                :: SO2   = 7
  integer, parameter                :: H2SO4 = 8
  integer, parameter                :: DMS   = 9
  integer, parameter                :: HNO3  = 10
  type(photolysis_rate_constants_t) :: photo_rate_constants
  type(environment_t)               :: environment
  type(chemical_state_t)            :: initial_state, final_state
  type(camp_t), pointer             :: camp
  integer                           :: camp_ids( 10 )
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
  initial_state%DMS__molec_cm3_   = 6.0e11_r8
  initial_state%HNO3__molec_cm3_  = 1.0e5_r8

  camp => camp_t( 'config.json' )
  camp_ids( M     ) = camp%get_state_index( 'M'     )
  camp_ids( H2O   ) = camp%get_state_index( 'H2O'   )
  camp_ids( HO2   ) = camp%get_state_index( 'HO2'   )
  camp_ids( OH    ) = camp%get_state_index( 'OH'    )
  camp_ids( NO3   ) = camp%get_state_index( 'NO3'   )
  camp_ids( H2O2  ) = camp%get_state_index( 'H2O2'  )
  camp_ids( SO2   ) = camp%get_state_index( 'SO2'   )
  camp_ids( H2SO4 ) = camp%get_state_index( 'H2SO4' )
  camp_ids( DMS   ) = camp%get_state_index( 'DMS'   )
  camp_ids( HNO3  ) = camp%get_state_index( 'HNO3'  )

  call camp%set_environment( environment )
  call camp%set_state( camp_ids( M     ), initial_state%M__molec_cm3_     )
  call camp%set_state( camp_ids( H2O   ), initial_state%H2O__molec_cm3_   )
  call camp%set_state( camp_ids( HO2   ), initial_state%HO2__molec_cm3_   )
  call camp%set_state( camp_ids( OH    ), initial_state%OH__molec_cm3_    )
  call camp%set_state( camp_ids( NO3   ), initial_state%NO3__molec_cm3_   )
  call camp%set_state( camp_ids( H2O2  ), initial_state%H2O2__molec_cm3_  )
  call camp%set_state( camp_ids( SO2   ), initial_state%SO2__molec_cm3_   )
  call camp%set_state( camp_ids( H2SO4 ), initial_state%H2SO4__molec_cm3_ )
  call camp%set_state( camp_ids( DMS   ), initial_state%DMS__molec_cm3_   )
  call camp%set_state( camp_ids( HNO3  ), initial_state%HNO3__molec_cm3_  )

  open( unit = file_unit, file = file_name, status = "replace",               &
        action = "write" )
  write( file_unit, * ) "time,H2O2,H2O2_camp,SO2,SO2_camp,H2SO4,H2SO4_camp,"//&
                        "DMS,DMS_camp,HNO3,HNO3_camp"
  call output_state( 0.0_r8, initial_state, camp, camp_ids )
  do i_time_step = 1, number_of_time_steps
    call advance_chemical_system( time_step__s, environment,                  &
                                  photo_rate_constants,                       &
                                  initial_state, final_state )
    call camp%solve( time_step__s )
    call output_state( i_time_step * time_step__s, final_state, camp,         &
                       camp_ids )
    initial_state = final_state
  end do
  close( file_unit )
  write(*,*) "Passed!"

contains

  subroutine output_state( time__s, state, camp, camp_ids )
    use pmc_util,                      only : to_string
    real(r8),                intent(in) :: time__s
    class(chemical_state_t), intent(in) :: state
    class(camp_t),           intent(in) :: camp
    integer,                 intent(in) :: camp_ids(:)
    write( file_unit, * ) trim( to_string( time__s           ) )//","//       &
        trim( to_string( state%H2O2__molec_cm3_              ) )//","//       &
        trim( to_string( camp%get_state( camp_ids( H2O2  ) ) ) )//","//       &
        trim( to_string( state%SO2__molec_cm3_               ) )//","//       &
        trim( to_string( camp%get_state( camp_ids( SO2   ) ) ) )//","//       &
        trim( to_string( state%H2SO4__molec_cm3_             ) )//","//       &
        trim( to_string( camp%get_state( camp_ids( H2SO4 ) ) ) )//","//       &
        trim( to_string( state%DMS__molec_cm3_               ) )//","//       &
        trim( to_string( camp%get_state( camp_ids( DMS   ) ) ) )//","//       &
        trim( to_string( state%HNO3__molec_cm3_              ) )//","//       &
        trim( to_string( camp%get_state( camp_ids( HNO3  ) ) ) )
  end subroutine output_state

end program camp_comparison_test
