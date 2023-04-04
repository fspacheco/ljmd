# CMake build for P2.1 code
cmake_minimum_required(VERSION 3.10)
project(ljmd VERSION 1.0 LANGUAGES C CXX) 
set(CMAKE_BUILD_TYPE RelWithDebInfo Release)

# enable OpenMP support by default?
#set(BUILD_OMP_DEFAULT OFF)
find_package(OpenMP REQUIRED)#COMPONENTS CXX QUIET)
if (OpenMP_FOUND)
	message(STATUS "OpenMP found")
	set(BUILD_OMP_DEFAULT ON)
endif()

#if(OpenMP_CXX_FOUND)
#	#check_omp_h_include()
#	#if(HAVE_OMP_H_INCLUDE)
#    set(BUILD_OMP_DEFAULT ON)
#	#endif()
#endif()
#option(BUILD_OMP "Build with OpenMP support" ${BUILD_OMP_DEFAULT})

# Add executable 
add_executable(ljmd-serial.x src/ljmd.c)
# compile options
target_compile_options(ljmd-serial.x PRIVATE -Wall -g -O0 -DLJMD_VERSION=2.0)#3 -ffast-math -fomit-frame-pointer)
# include directory for header files
include_directories(include)

# add library for force computation
add_library(mdlib src/utilities.c src/verlet_ti.c src/forces.c src/energies.c src/input.c src/output.c src/reading.c)
target_compile_options(mdlib PRIVATE -Wall -g -O0 -DLJMD_VERSION=2.0)#3 -ffast-math -fomit-frame-pointer)
target_link_libraries(mdlib m)

# link libraries
target_link_libraries(ljmd-serial.x m)
target_link_libraries(ljmd-serial.x mdlib)

# Add executable 
add_executable(ljmd-openmp.x src/ljmd.c)
# compile options
#target_compile_options(ljmd-openmp.x PRIVATE -fopenmp -Wall -g -O0 -DLJMD_VERSION=2.0)
target_compile_options(ljmd-openmp.x PRIVATE -fopenmp -O3 -DLJMD_VERSION=2.0)
# include directory for header files
include_directories(include)

# add library for force computation
add_library(mdlib-openmp src/utilities.c src/verlet_ti.c src/forces.c src/energies.c src/input.c src/output.c src/reading.c)
#target_compile_options(mdlib-openmp PRIVATE -fopenmp -Wall -g -O0 -DLJMD_VERSION=2.0)
target_compile_options(mdlib-openmp PRIVATE -fopenmp -O3 -DLJMD_VERSION=2.0)
target_link_libraries(mdlib-openmp m)

# link libraries
target_link_options(ljmd-openmp.x PRIVATE -fopenmp)
target_link_libraries(ljmd-openmp.x m)
target_link_libraries(ljmd-openmp.x mdlib-openmp)


# add subdirectories
add_subdirectory(examples)


# set target for serial version
add_custom_target(serial
    COMMAND make
)

# set target for cleaning
add_custom_target(clean_data_files
    COMMAND ${CMAKE_MAKE_PROGRAM} -C examples clean_example
)

# set target for checking
add_custom_target(check
    DEPENDS serial
    COMMAND ${CMAKE_MAKE_PROGRAM} -C examples check_example
)

enable_testing()

add_subdirectory(test)

