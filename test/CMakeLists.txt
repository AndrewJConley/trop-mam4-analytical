################################################################################
# Copy test data

add_custom_target(copy_test_data ALL ${CMAKE_COMMAND} -E copy_directory
  ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_BINARY_DIR}/test)

################################################################################
# Tests

add_executable(camp_comparison_test ${PROJECT_SOURCE_DIR}/chemical_state.F90
                                    ${PROJECT_SOURCE_DIR}/environment.F90
                                    ${PROJECT_SOURCE_DIR}/photolysis_rate_constants.F90
                                    ${PROJECT_SOURCE_DIR}/shr_kind_mod.F90
                                    ${PROJECT_SOURCE_DIR}/trop_mam4_analytical.F90
                                    camp_comparison_test.F90)
target_link_libraries(camp_comparison_test ${PARTMC_LIB})
add_test(NAME camp_comparison_test COMMAND camp_comparison_test.sh)

################################################################################
