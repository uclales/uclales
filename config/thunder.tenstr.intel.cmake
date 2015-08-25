# For IntelMPI
#set(CMAKE_C_COMPILER "mpiicc")
#set(CMAKE_CXX_COMPILER "mpiicpc")
#set(CMAKE_Fortran_COMPILER "mpiifort")

# For OpenMPI
set(CMAKE_C_COMPILER "mpicc")
set(CMAKE_CXX_COMPILER "mpicxx")
set(CMAKE_Fortran_COMPILER "mpif90")

set(USER_C_FLAGS "-std=c99 ")
set(USER_Fortran_FLAGS " -fpp -traceback -r8 -ftz -extend-source -g ")
set(USER_Fortran_FLAGS_RELEASE " -O3 -no-prec-div -xAVX ")
set(USER_Fortran_FLAGS_DEBUG "-fpe0 -O0 -g -check all -check nopointers -check noarg_temp_created")

set(NETCDF_INCLUDE_DIR "$ENV{NETCDFFROOT}/include")
set(NETCDF_LIB_1       "$ENV{NETCDFFROOT}/lib/libnetcdff.a")
set(NETCDF_LIB_2       "$ENV{NETCDFCROOT}/lib/libnetcdf.a")
set(HDF5_LIB_1         "$ENV{HDF5ROOT}/lib/libhdf5_hl.a")
set(HDF5_LIB_2         "$ENV{HDF5ROOT}/lib/libhdf5.a")
set(SZIP_LIB           "$ENV{SZIPROOT}/lib/libsz.a")
set(LIBS ${NETCDF_LIB_1} ${NETCDF_LIB_2} ${HDF5_LIB_1} ${HDF5_LIB_2} ${SZIP_LIB} m z curl)

# TENSTREAM LIB for 3D Radiation
set(TENSTREAM_DIR       "$ENV{HOME}/tenstream/")
set(TENSTREAM_BUILD_DIR "$ENV{HOME}/tenstream/build/")
