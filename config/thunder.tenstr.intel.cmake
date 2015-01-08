# Default GCC
#set(CMAKE_Fortran_COMPILER "ifort")
#set(Fortran_COMPILER_WRAPPER mpif90)

set(CMAKE_C_COMPILER "mpicc")
set(CMAKE_CXX_COMPILER "mpicxx")
set(CMAKE_Fortran_COMPILER "mpif90")

set(USER_C_FLAGS "-std=c99 ")
set(USER_Fortran_FLAGS " -fpp -traceback -r8 -ftz -extend_source -g ")
set(USER_Fortran_FLAGS_RELEASE " -O3 -no-prec-div -xAVX -fp-model source ")
set(USER_Fortran_FLAGS_DEBUG "-fpe0 -O0 -g -check all -check nopointers -check noarg_temp_created")

#set(NETCDF_INCLUDE_DIR "/scratch/mpi/mpiaes/m300362/libs/netcdf-fortran-intel15/include")
#set(NETCDF_LIB_1       "/scratch/mpi/mpiaes/m300362/libs/netcdf-fortran-intel15/lib/libnetcdff.a")
#set(NETCDF_LIB_2       "/scratch/mpi/mpiaes/m300362/libs/netcdf-intel15/lib/libnetcdf.a")
#
#list(APPEND HDF5_LIBRARIES "/scratch/mpi/mpiaes/m300362/libs/hdf5-intel15/lib/libhdf5_hl.a")
#list(APPEND HDF5_LIBRARIES "/scratch/mpi/mpiaes/m300362/libs/hdf5-intel15/lib/libhdf5hl_fortran.a")
#list(APPEND HDF5_LIBRARIES "/scratch/mpi/mpiaes/m300362/libs/hdf5-intel15/lib/libhdf5.a")
#list(APPEND HDF5_LIBRARIES "/scratch/mpi/mpiaes/m300362/libs/hdf5-intel15/lib/libhdf5_fortran.a")
#
#set(SZIP_LIB           "/sw/squeeze-x64/szip-2.1-static/lib/libsz.a")
#set(LIBS ${NETCDF_LIB_1} ${NETCDF_LIB_2} ${HDF5_LIBRARIES} ${SZIP_LIB} m z curl dl)

set(NETCDF_INCLUDE_DIR "/sw/squeeze-x64/netcdf_fortran-4.2-static-intel15/include")
set(NETCDF_LIB_1       "/sw/squeeze-x64/netcdf_fortran-4.2-static-intel15/lib/libnetcdff.a")
set(NETCDF_LIB_2       "/sw/squeeze-x64/netcdf-4.2-static/lib/libnetcdf.a")
set(HDF5_LIB_1         "/sw/squeeze-x64/hdf5-1.8.8-static/lib/libhdf5_hl.a")
set(HDF5_LIB_2         "/sw/squeeze-x64/hdf5-1.8.8-static/lib/libhdf5.a")
set(SZIP_LIB           "/sw/squeeze-x64/szip-2.1-static/lib/libsz.a")
set(LIBS ${NETCDF_LIB_1} ${NETCDF_LIB_2} ${HDF5_LIB_1} ${HDF5_LIB_2} ${SZIP_LIB} m z curl)

# TENSTREAM LIB for 3D Radiation
set(TENSTREAM_DIR       "/home/zmaw/m300362/tenstream/")
set(TENSTREAM_BUILD_DIR "/home/zmaw/m300362/tenstream/build-intel15/")
