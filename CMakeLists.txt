############################
# MESOHH CMAKE FILE
# Supports (thusfar) GNU, Intel and IBM compilers
# To set the compiler to a specific one, set the FC environment variable
# NETCDF can be set in the config/default.cmake, or through command line options
# NETCDF_INCLUDE, NETCDF_LIB, FFTW_INCLUDE and FFTW_LIB, respectively.
# Command line options are fed to CMake using the -D switch
# Options include:
# * MPI (TRUE/FALSE) to toggle the use of MPI. True by default.
# * PROFILER (SCALASCA/MARMOT) to enable a specific profiler. Default is none.
# * CMAKE_BUILD_TYPE (RELEASE/DEBUG) to build in optimized or debug mode. Default is RELEASE
# All options are persistent in the sense that once set, they will last until changed by the user.
# For daily use "cmake .. && make" should suffice to build the code, but a command line could look like:
# export FC="gfortran" && cmake -D SYST=tornado -D MPI=TRUE -D PROFILER=SCALASCA -D CMAKE_BUILD_TYPE=DEBUG .. && make
############################
set (CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_CURRENT_SOURCE_DIR}/config)

INCLUDE(CMakeForceCompiler) #Necessary to change between MPI/parallel/profiler compilers, without having to do a make clean
cmake_minimum_required (VERSION 2.8.1)
set (CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR})


#  make sure that the default is a RELEASE
if (NOT CMAKE_BUILD_TYPE)
  set (CMAKE_BUILD_TYPE RELEASE CACHE STRING
      "Choose the type of build, options are: None Debug Release." FORCE )
else()
  string(TOUPPER ${CMAKE_BUILD_TYPE} TEMP)
  set(CMAKE_BUILD_TYPE ${TEMP} CACHE STRING
      "Choose the type of build, options are: None Debug Release."  FORCE)
endif ()

#Load system specific settings
if (NOT SYST)
  set (SYST default CACHE STRING
      "Choose the location: mpipc thunder fedora archlinux supermuc jugene" FORCE )
endif ()
include (${SYST} OPTIONAL)

#Start the project only after all the variables are set
project (uclales Fortran)

# Set MPI on or off. By default it is on.
set (MPI TRUE CACHE BOOL "Choose to enable MPI or not.")
if(MPI)
 MESSAGE(STATUS "MPI is Enabled.")
 FILE(GLOB mpifile "${CMAKE_CURRENT_SOURCE_DIR}/src/mpi/mpi_interface.f90")
else()
  MESSAGE(STATUS "MPI is Disabled.")
  FILE(GLOB mpifile "${CMAKE_CURRENT_SOURCE_DIR}/src/seq/seq_interface.f90")
endif()
execute_process(COMMAND ${CMAKE_COMMAND} -E copy_if_different ${mpifile} ${CMAKE_CURRENT_SOURCE_DIR}/src/mpi_interface.f90)

get_filename_component (Fortran_COMPILER_NAME ${CMAKE_Fortran_COMPILER_ID} NAME)

if (NOT ${Fortran_COMPILER_WRAPPER}  STREQUAL "")
  CMAKE_FORCE_Fortran_COMPILER(${Fortran_COMPILER_WRAPPER} ${Fortran_COMPILER_NAME})
endif ()
MESSAGE(STATUS "Fortran Compiler " ${CMAKE_Fortran_COMPILER})

set(CMAKE_Fortran_FLAGS ${USER_Fortran_FLAGS} )
set(CMAKE_Fortran_FLAGS_RELEASE ${USER_Fortran_FLAGS_RELEASE} )
set(CMAKE_Fortran_FLAGS_DEBUG ${USER_Fortran_FLAGS_DEBUG})
MESSAGE(STATUS "Build Type " ${CMAKE_BUILD_TYPE})
if (CMAKE_BUILD_TYPE STREQUAL "RELEASE")
  MESSAGE(STATUS "Compiler Flags " ${CMAKE_Fortran_FLAGS} ${CMAKE_Fortran_FLAGS_RELEASE})
else()
  MESSAGE(STATUS "Compiler Flags " ${CMAKE_Fortran_FLAGS} ${CMAKE_Fortran_FLAGS_DEBUG})
endif()


set(INCLUDE_DIRS ${NETCDF_INCLUDE_DIR})

option(WITH_TENSTREAM "Include the tenstream solver for radiation" OFF)
if(WITH_TENSTREAM)
  if((NOT TENSTREAM_DIR) OR (NOT TENSTREAM_BUILD_DIR))
    message( FATAL_ERROR "You asked to include the TENSTREAM solver, yet you did not specify a path to the library -- please set: TENSTREAM_DIR(root dir of tenstream solver) and TENSTREAM_BUILD_DIR(directory where you issued cmake&&make)" )
  endif((NOT TENSTREAM_DIR) OR (NOT TENSTREAM_BUILD_DIR))

  set (CMAKE_C_FLAGS "-cpp -W -Wall -Wuninitialized --std=c99" )

  if(NOT HDF5_ROOT_DIR_HINT)
    set( HDF5_ROOT_DIR_HINT "${TENSTREAM_DIR}/cmake/hdf5")
  endif(NOT HDF5_ROOT_DIR_HINT)
  FIND_PACKAGE (HDF5 NAMES hdf5 REQUIRED)

  # Petsc
  enable_language(Fortran C)
  list (APPEND CMAKE_MODULE_PATH "${TENSTREAM_DIR}/cmake/petsc_cmake_modules/")
  find_package(PETSc REQUIRED)
  add_definitions (${PETSC_DEFINITIONS})

  include_directories (${PETSC_INCLUDES} ${ITAPS_INCLUDES} ${HDF5_INCLUDE_DIRS})


  FIND_LIBRARY(TENSTREAM_LIBRARY tenstream "${TENSTREAM_BUILD_DIR}/lib/")
  #  LIST(APPEND LIBS "${TENSTREAM_LIBRARY}")
  set(INCLUDE_DIRS ${INCLUDE_DIRS} "${TENSTREAM_BUILD_DIR}/src/")

  list (APPEND LIBS "${HDF5_LIBRARIES}" "${PETSC_LIBRARIES}" "${TENSTREAM_LIBRARY}")

  add_definitions(-DHAVE_TENSTREAM)
endif(WITH_TENSTREAM)

add_subdirectory(src)


#####################
# DOCUMENTATION
#####################
add_custom_target(todo ALL)
ADD_CUSTOM_COMMAND(TARGET todo POST_BUILD
                  COMMAND echo "UCLALES TODO LIST" > TODO
                  COMMAND date  >> TODO
                  COMMAND grep -Rin \\todo  src | sed 's/!.*TODO//I' >>  TODO
                  WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
                  )

