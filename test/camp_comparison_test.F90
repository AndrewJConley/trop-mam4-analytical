!> \file
!> Comparison of CAMP solver output with the analytical solution for trop-mam4

!> CAMP comparison module
program camp_comparison_test

  use shr_kind_mod,                    only : r8 => shr_kind_r8
  use trop_mam4_analytical,            only : advance_chemical_system

  implicit none

  ! Environmental conditions
  real(r8), parameter :: boltzmann_constant__J_K = 1.380649e-23_r8
  real(r8), parameter :: temperature__K          = 298.15_r8
  real(r8), parameter :: pressure__Pa            = 101325.0_r8
  real(r8), parameter :: air_density__molec_cm3  = pressure__Pa &
                                                   / boltzmann_constant__J_K &
                                                   / temperature__K

  ! Photolysis rate constants
  real(r8), parameter :: H2O2_photo_rate_constant__s = 2.0e-5_r8

  ! Fixed chemical species concentrations
  ! (\todo come up with more realistic initial concentration values)
  real(r8), parameter :: H2O__molec_cm3   = 1.2e16_r8
  real(r8), parameter :: HO2__molec_cm3   = 1.0e8_r8
  real(r8), parameter :: OH__molec_cm3    = 1.0e5_r8
  real(r8), parameter :: NO3__molec_cm3   = 3.61e8_r8

  ! Variable chemical species concentrations
  real(8) :: H2O2__molec_cm3  = 5.42e14_r8
  real(8) :: SO2__molec_cm3   = 1.0e5_r8
  real(8) :: H2SO4__molec_cm3 = 1.0e5_r8
  real(8) :: DMS__molec_cm3   = 3.6e15_r8
  real(8) :: HNO3__molec_cm3  = 1.0e5_r8

  ! simulation conditions
  real(r8), parameter :: time_step__s = 60.0_r8
  integer,  parameter :: number_of_time_steps = 60 * 24

  ! output file properties
  integer, parameter :: file_unit = 10
  character(len=*), parameter :: file_name = "out/analytic_solution.csv"

  integer :: i_time_step

  open( unit = file_unit, file = file_name, status = "replace",               &
        action = "write" )
  write( file_unit, * ) "time,H2O2,SO2,H2SO4,DMS,HNO3"
  write( file_unit, * ) 0.0_r8,",",H2O2__molec_cm3, ",", SO2__molec_cm3, ",", &
                        H2SO4__molec_cm3, ",", DMS__molec_cm3, ",",           &
                        HNO3__molec_cm3
  do i_time_step = 1, number_of_time_steps
    call advance_chemical_system( time_step__s, temperature__K,               &
        pressure__Pa, H2O2_photo_rate_constant__s,                            &
        air_density__molec_cm3,                                               &
        H2O__molec_cm3 / air_density__molec_cm3,                              &
        HO2__molec_cm3,                                                       &
        OH__molec_cm3,                                                        &
        NO3__molec_cm3,                                                       &
        H2O2__molec_cm3,  H2O2__molec_cm3,                                    &
        SO2__molec_cm3,   SO2__molec_cm3,                                     &
        H2SO4__molec_cm3, H2SO4__molec_cm3,                                   &
        DMS__molec_cm3,   DMS__molec_cm3,                                     &
        HNO3__molec_cm3,  HNO3__molec_cm3 )
    write( file_unit, * ) i_time_step * time_step__s, ",",                    &
                          H2O2__molec_cm3, ",", SO2__molec_cm3, ",",          &
                          H2SO4__molec_cm3, ",", DMS__molec_cm3, ",",         &
                          HNO3__molec_cm3
  end do
  close( file_unit )
  write(*,*) "Passed!"

end program camp_comparison_test
