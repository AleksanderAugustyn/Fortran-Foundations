add_executable(fortran_foundations_tests test_fortran_foundations_p.f08)

target_link_libraries(fortran_foundations_tests PRIVATE fortran_foundations fortran_foundations_flags)

set_target_properties(fortran_foundations_tests PROPERTIES
        Fortran_MODULE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/mod_files
)

add_test(NAME fortran_foundations_tests COMMAND fortran_foundations_tests)
