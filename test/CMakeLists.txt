include(FetchContent)
FetchContent_Declare(gtest  QUIET
  URL https://github.com/google/googletest/archive/release-1.12.0.tar.gz
)
# configure build of googletest
set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)
set(BUILD_GMOCK OFF CACHE BOOL "" FORCE)
FetchContent_MakeAvailable(gtest)

add_executable(unit_tests kinetic_energy.cpp first_velverlet.cpp utilities.cpp forces.cpp input_output.cpp)

target_link_libraries(unit_tests PRIVATE mdlib gtest_main)

# automatic discovery of unit tests
include(GoogleTest)
gtest_discover_tests(unit_tests PROPERTIES
    LABELS "unit"
  DISCOVERY_TIMEOUT  # how long to wait (in seconds) before crashing
    240)


add_test(NAME CheckCodeWithOpenMP
	COMMAND bash -c "cd ../examples; ../ljmd-openmp.x < argon_108.inp; \
	head -10 argon_108.dat | LC_ALL=C awk '{printf(\"%d %.6f %.6f %.6f\",$1,$2,$3,$4);}' > test.dat; \
	head -10 ../../reference/argon_108.dat | LC_ALL=C awk '{printf(\"%d %.6f %.6f %.6f\",$1,$2,$3,$4);}' > ref.dat; \
	cmp test.dat ref.dat || exit 1; rm -f test.dat ref.dat"
)

add_custom_target(run_openmp COMMAND ${CMAKE_CTEST_COMMAND} --verbose -R CheckCodeWithOpenMP)

