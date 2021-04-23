!> CAMP wrapper
module camp

  use environment,                     only : environment_t
  use pmc_camp_core,                   only : camp_core_t
  use pmc_camp_state,                  only : camp_state_t
  use pmc_chem_spec_data,              only : chem_spec_data_t
  use pmc_constants,                   only : dp

  implicit none

  private
  public :: camp_t

  type :: camp_t
    private
    type(camp_core_t),      pointer :: core_
    type(camp_state_t),     pointer :: state_
    type(chem_spec_data_t), pointer :: chem_data_
    real(kind=dp)                   :: molec_cm3_to_ppm_ = 0.0
  contains
    procedure :: get_state_index
    procedure :: solve
    procedure :: set_environment
    procedure :: get_state
    procedure :: set_state
  end type camp_t

  interface camp_t
    module procedure :: constructor
  end interface

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  function constructor( config_file_path ) result( new_obj )
    use pmc_util,                      only : assert
    type(camp_t),     pointer    :: new_obj
    character(len=*), intent(in) :: config_file_path
    allocate( new_obj )
    new_obj%core_       => camp_core_t( config_file_path )
    call new_obj%core_%initialize( )
    call new_obj%core_%solver_initialize( )
    new_obj%state_ => new_obj%core_%new_state( )
    call assert( 125584813,                                                   &
                 new_obj%core_%get_chem_spec_data( new_obj%chem_data_ ) )
  end function constructor

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  integer function get_state_index( this, species_name )
    class(camp_t),    intent(inout) :: this
    character(len=*), intent(in)    :: species_name
    get_state_index = this%chem_data_%gas_state_id( species_name )
  end function get_state_index

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine solve( this, time_step )
    class(camp_t), intent(inout) :: this
    real(kind=dp), intent(in)    :: time_step
    call this%core_%solve( this%state_, time_step )
  end subroutine solve

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine set_environment( this, environment )
    class(camp_t),        intent(inout) :: this
    class(environment_t), intent(in)    :: environment
    real(kind=dp), parameter :: k_b = 1.380649e-23 ! J/K
    call this%state_%env_states(1)%set_temperature_K(                         &
           environment%temperature__K_ )
    call this%state_%env_states(1)%set_pressure_Pa(                           &
           environment%pressure__Pa_   )
    this%molec_cm3_to_ppm_ = k_b * environment%temperature__K_                &
                             / environment%pressure__Pa_ * 1e6
  end subroutine set_environment

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  real(kind=dp) function get_state( this, species_index )
    class(camp_t), intent(in) :: this
    integer,       intent(in) :: species_index
    get_state = this%state_%state_var( species_index )                        &
                      / this%molec_cm3_to_ppm_
  end function get_state

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine set_state( this, species_index, value )
    class(camp_t), intent(inout) :: this
    integer,       intent(in)    :: species_index
    real(kind=dp), intent(in)    :: value
    this%state_%state_var( species_index ) = value * this%molec_cm3_to_ppm_
  end subroutine set_state

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module camp
