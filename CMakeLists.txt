############################
# MESOHH CMAKE FILE
# Supports (thusfar) GNU, Intel and IBM compilers
# To set the compiler to a specific one, set the FC environment variable
# NETCDF and FFTW libraries can be set using the following environment variables:
# NETCDF_INCLUDE, NETCDF_LIB, FFTW_INCLUDE and FFTW_LIB, respectively.
# Command line options are fed to CMake using the -D switch
# Options include:
# * MPI (TRUE/FALSE) to toggle the use of MPI. True by default.
# * PROFILER (SCALASCA/MARMOT) to enable a specific profiler. Default is none.
# * CMAKE_BUILD_TYPE (RELEASE/DEBUG) to build in optimized or debug mode. Default is RELEASE
# All options are persistent in the sense that once set, they will last until changed by the user.
# For daily use "cmake .. && make" should suffice to build the code, but a command line could look like:
# export FC="gfortran" && cmake -DMPI=TRUE -PROFILER=SCALASCA -CMAKE_BUILD_TYPE=DEBUG .. && make
############################
INCLUDE(CMakeForceCompiler) #Necessary to change between MPI/parallel/profiler compilers, without having to do a make clean
cmake_minimum_required (VERSION 2.8.1)
project (uclales Fortran)

# Set the name of the executable, and build it into the build directory
set (CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR})
set(target_name uclales)

# Set the default paths for libraries
set (CMAKE_PREFIX_PATH $ENV{HOME} /usr
    CACHE STRING "Root dirs of the include and lib directories")
set (CMAKE_LIBRARY_PATH /usr/lib64 $ENV{SARA_NETCDF_LIB} $ENV{NETCDF_LIB}
    CACHE STRING "Specific library directories")
set (CMAKE_INCLUDE_PATH $ENV{SARA_NETCDF_INCLUDE} $ENV{NETCDF_INCLUDE}
    CACHE STRING "Specific include directories")
set (FIND_LIBRARY_USE_LIB64_PATHS TRUE)

#  make sure that the default is a RELEASE
if (NOT CMAKE_BUILD_TYPE)
  set (CMAKE_BUILD_TYPE RELEASE CACHE STRING
      "Choose the type of build, options are: None Debug Release." FORCE )
endif ()

# Set MPI on or off. By default it is on.
set (MPI TRUE CACHE BOOL "Choose to enable MPI or not.")
if(MPI)
 MESSAGE(STATUS "MPI is Enabled.")
 FILE(GLOB mpifile "${CMAKE_CURRENT_SOURCE_DIR}/src/mpi/mpi_interface.f90")
else()
  MESSAGE(STATUS "MPI is Disabled.")
  FILE(GLOB mpifile "${CMAKE_CURRENT_SOURCE_DIR}/mpi/seq_interface.f90")
endif()
execute_process(COMMAND ${CMAKE_COMMAND} -E copy_if_different ${mpifile} ${CMAKE_CURRENT_SOURCE_DIR}/src/mpi_interface.f90)

######################
# Retrieve the type of Fortran compiler. This is necassary to set the compiler options.
######################
get_filename_component (Fortran_COMPILER_NAME ${CMAKE_Fortran_COMPILER_ID} NAME)
if (Fortran_COMPILER_NAME STREQUAL "GNU")
  if(MPI)
    set (COMPILER_WRAPPER mpif90)
  endif()
  set (CMAKE_Fortran_FLAGS "-pedantic -fopenmp -finit-real=nan -W -Wall -fdefault-real-8 -ffree-line-length-none")
  set (CMAKE_Fortran_FLAGS_RELEASE "-funroll-all-loops -fno-f2c -O3")
  set (CMAKE_Fortran_FLAGS_DEBUG   " -Wuninitialized -O -fbounds-check -fbacktrace -fno-f2c -O0 -g -ffpe-trap=invalid,zero,overflow")
elseif (Fortran_COMPILER_NAME STREQUAL "Intel")
  if(MPI)
    set (COMPILER_WRAPPER mpif90)
  endif()
    set(CMAKE_Fortran_FLAGS "-O3 -r8 -diag-disable remark")
    set (CMAKE_Fortran_FLAGS_RELEASE "-O3")
    set (CMAKE_Fortran_FLAGS_DEBUG   "-O0 -g -ffpe-trap=invalid,zero,overflow")
elseif (Fortran_COMPILER_NAME STREQUAL "XL")
  if(MPI)
    set (COMPILER_WRAPPER mpfort)
  endif()
    set(CMAKE_Fortran_FLAGS "-qsmp=omp:noauto -qfree=F90 -qrealsize=8  -qwarn64 -qnosave -qinitauto=FFF00000 -qflttrap=en:ov:zero:inv:imp -qflag=w:e")
    set(CMAKE_Fortran_FLAGS_RELEASE "-O4 -qnoipa -qstrict=none:exceptions  ")
    set(CMAKE_Fortran_FLAGS_DEBUG "-O2 -g -qfullpath -C -qflttrp=enable:nanq:overflow:zerodivide -qsigtrap -qinitauto=ff")
endif ()


# Set Profiler
set (PROFILER CACHE STRING "Enable a profiler. Options: SCALASCA or MARMOT")
if(PROFILER STREQUAL "SCALASCA")
  set (COMPILER_WRAPPER scalf90)
elseif(PROFILER STREQUAL "MARMOT")
  set (COMPILER_WRAPPER marmotf90)
endif()

# Set the compiler name to the wrapper, if necessary
if (NOT ${COMPILER_WRAPPER} STREQUAL "")
  CMAKE_FORCE_Fortran_COMPILER(${COMPILER_WRAPPER} ${Fortran_COMPILER_NAME})
endif ()
MESSAGE(STATUS "Compiler " ${CMAKE_Fortran_COMPILER})

# FIND NETCDF
find_path(NETCDF_INCLUDE_DIR netcdf.mod NETCDF.mod
  DOC "NetCDF include directory (must contain netcdf.mod)"
)
if(NETCDF_INCLUDE_DIR)
  MESSAGE(STATUS "NETCDF module found in " ${NETCDF_INCLUDE_DIR})
else()
  MESSAGE(STATUS "WARNING: No NETCDF bindings are found.")
endif()

find_library(NETCDF_C_LIB netcdf)
if(NETCDF_C_LIB)
  set(NETCDF_LIBS ${NETCDF_C_LIB})
  MESSAGE(STATUS "NETCDF C library: " ${NETCDF_C_LIB})
else()
  MESSAGE(STATUS "WARNING: No NETCDF bindings are found.")
endif()

find_library(NETCDF_FORTRAN_LIB netcdff
  DOC "NetCDF Fortran library"
)
if(NETCDF_FORTRAN_LIB)
  set(NETCDF_LIBS ${NETCDF_LIBS} ${NETCDF_FORTRAN_LIB})
  MESSAGE(STATUS "NETCDF Fortran library: " ${NETCDF_FORTRAN_LIB})
else()
  MESSAGE(STATUS "WARNING: No Fortran NETCDF bindings are found. This is not always fatal.")
endif()



#####################
# CONCATENATE INCLUDE DIRECTORIES
#####################
include_directories(${NETCDF_INCLUDE_DIR})

#####################
# COMPILE SRC DIR
#####################
ADD_SUBDIRECTORY( src)

