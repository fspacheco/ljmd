include(FetchContent)
FetchContent_Declare(gtest  QUIET
  URL https://github.com/google/googletest/archive/release-1.12.0.tar.gz
)
# configure build of googletest
set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)
set(BUILD_GMOCK OFF CACHE BOOL "" FORCE)
FetchContent_MakeAvailable(gtest)

add_executable(unit_tests kinetic_energy.cpp first_velverlet.cpp utilities.cpp forces.cpp input_output.cpp)
target_link_options(unit_tests PRIVATE -fopenmp)
target_link_libraries(unit_tests PRIVATE mdlib-util mdlib-mpi gtest_main)

# automatic discovery of unit tests
include(GoogleTest)
gtest_discover_tests(unit_tests PROPERTIES
    LABELS "unit"
  DISCOVERY_TIMEOUT  # how long to wait (in seconds) before crashing
    240)

add_test(NAME CheckCodeWithRef
        COMMAND bash -c "cd ..;make check")

add_custom_target(run_check_ref COMMAND ${CMAKE_CTEST_COMMAND} --verbose -R CheckCodeWithRef)

