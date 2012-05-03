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
INCLUDE(CMakeForceCompiler) #Necessary to change between MPI/parallel/profiler compilers, without having to do a make clean
cmake_minimum_required (VERSION 2.8.1)
project (uclales Fortran)

# Set the name of the executable, and build it into the build directory
set (CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR})
set(target_name uclales)
set (CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_CURRENT_SOURCE_DIR}/config)

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
  FILE(GLOB mpifile "${CMAKE_CURRENT_SOURCE_DIR}/src/seq/seq_interface.f90")
endif()
execute_process(COMMAND ${CMAKE_COMMAND} -E copy_if_different ${mpifile} ${CMAKE_CURRENT_SOURCE_DIR}/src/mpi_interface.f90)

#READ in local settings (e.g. paths)
set (FIND_LIBRARY_USE_LIB64_PATHS TRUE)

set (SYST default CACHE STRING "Choose the platform we are on.")

include (${SYST} OPTIONAL)
######################
# Retrieve the type of Fortran compiler. This is necassary to set the compiler options.
######################
get_filename_component (Fortran_COMPILER_NAME ${CMAKE_Fortran_COMPILER_ID} NAME)
if (Fortran_COMPILER_NAME STREQUAL "GNU")
  if(MPI)
    set (COMPILER_WRAPPER mpif90)
  endif()
  set (CMAKE_Fortran_FLAGS "-fbacktrace -finit-real=nan -fdefault-real-8  -fno-f2c -ffree-line-length-none")
  set (CMAKE_Fortran_FLAGS_RELEASE "-funroll-all-loops -O3")
  set (CMAKE_Fortran_FLAGS_DEBUG   "-W -Wall -Wuninitialized -fcheck=all -fbacktrace -O0 -g -ffpe-trap=invalid,zero,overflow")
  set (PROFILER_FLAG "-pg")
elseif (Fortran_COMPILER_NAME STREQUAL "NAG")
  if(MPI)
    set (COMPILER_WRAPPER mpif90)
  endif()
  set (CMAKE_PREFIX_PATH "")
  set (CMAKE_Fortran_FLAGS "-kind=byte -colour -C -f2003 -nan -r8 -free -w=x95 -wmismatch=mpi_allreduce ")
  set (CMAKE_Fortran_FLAGS_RELEASE "-O4")
  set (CMAKE_Fortran_FLAGS_DEBUG   "-O0 -gline -g")
  set (PROFILER_FLAG "-pg")
elseif (Fortran_COMPILER_NAME STREQUAL "Intel")
  if(MPI)
    set (COMPILER_WRAPPER mpif90)
  endif()
  set(CMAKE_Fortran_FLAGS "-r8 -ftz -extend_source")
  set (CMAKE_Fortran_FLAGS_RELEASE "-O3 -no-prec-div -xHost")
  set (CMAKE_Fortran_FLAGS_DEBUG   "-traceback -fpe0 -O0 -g -check all -check noarg_temp_created -r8 -debug extended -ftrapuv -fp-stack-check -check pointers -check bounds")
  set (PROFILER_FLAG "-p")
elseif (Fortran_COMPILER_NAME STREQUAL "XL")
  if(MPI)
    set (COMPILER_WRAPPER mpfort)
  endif()
  set(CMAKE_Fortran_FLAGS "-qthreaded -qsmp -qfree=F90 -qrealsize=8  -qwarn64 -qnosave -qinitauto=FFF00000 -qflttrap=en:ov:zero:inv:imp -qflag=w:e")
  set(CMAKE_Fortran_FLAGS_RELEASE "-O4 -qnoipa -qstrict=none:exceptions -qinitauto=ff -qsigtrap")
  set(CMAKE_Fortran_FLAGS_DEBUG "-O2 -qfullpath -C -qflttrp=enable:inexact:invalid:nanq:overflow:zerodivide -qsigtrap -qinitauto")
  set (PROFILER_FLAG "-pg")
endif ()


# Set Profiler
set (PROFILER CACHE STRING "Enable a profiler. Options: NONE, GPROF, SCALASCA or MARMOT")
if(PROFILER STREQUAL "GPROF")
  set (CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} ${PROFILER_FLAG}")
elseif(PROFILER STREQUAL "SCALASCA")
  set (COMPILER_WRAPPER scalf90)
elseif(PROFILER STREQUAL "MARMOT")
  set (COMPILER_WRAPPER marmotf90)
endif()


find_package(NETCDF REQUIRED COMPONENTS F90)
set (UCLALES_INCLUDE_DIRS ${NETCDF_INCLUDE_DIRS})
set (UCLALES_LIBRARIES ${NETCDF_LIBRARIES} ${HDF5_LIBRARIES})


# Set the compiler name to the wrapper, if necessary
if (NOT ${COMPILER_WRAPPER}  STREQUAL "")
  CMAKE_FORCE_Fortran_COMPILER(${COMPILER_WRAPPER} ${Fortran_COMPILER_NAME})
endif ()
MESSAGE(STATUS "Compiler " ${CMAKE_Fortran_COMPILER})

#####################
# COMPILE SRC DIR
#####################
ADD_SUBDIRECTORY( src)
