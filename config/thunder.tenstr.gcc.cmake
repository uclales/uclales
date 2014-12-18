# Default GCC
set(CMAKE_Fortran_COMPILER "gfortran")
set(Fortran_COMPILER_WRAPPER mpif90)

set(USER_Fortran_FLAGS "-fbacktrace -finit-real=nan -fdefault-real-8 -fdefault-double-8 -fno-f2c -ffree-line-length-none")
set(USER_Fortran_FLAGS_RELEASE "-funroll-all-loops -O3 -march=native -mtune=native")
set(USER_Fortran_FLAGS_DEBUG "-W -Wall -Wuninitialized -fcheck=all -fbacktrace -O0 -g -ffpe-trap=invalid,zero,overflow")

set(NETCDF_INCLUDE_DIR "/scratch/mpi/mpiaes/m300362/libs/netcdf-fortran-gcc48/include")
set(NETCDF_LIB_1       "/scratch/mpi/mpiaes/m300362/libs/netcdf-fortran-gcc48/lib/libnetcdff.a")
set(NETCDF_LIB_2       "/scratch/mpi/mpiaes/m300362/libs/netcdf-gcc48/lib/libnetcdf.a")

list(APPEND HDF5_LIBRARIES "/scratch/mpi/mpiaes/m300362/libs/hdf5/lib/libhdf5_hl.a")
list(APPEND HDF5_LIBRARIES "/scratch/mpi/mpiaes/m300362/libs/hdf5/lib/libhdf5hl_fortran.a")
list(APPEND HDF5_LIBRARIES "/scratch/mpi/mpiaes/m300362/libs/hdf5/lib/libhdf5.a")
list(APPEND HDF5_LIBRARIES "/scratch/mpi/mpiaes/m300362/libs/hdf5/lib/libhdf5_fortran.a")

set(SZIP_LIB           "/sw/squeeze-x64/szip-2.1-static/lib/libsz.a")
set(LIBS ${NETCDF_LIB_1} ${NETCDF_LIB_2} ${HDF5_LIBRARIES} ${SZIP_LIB} m z curl dl)

# TENSTREAM LIB for 3D Radiation
set(TENSTREAM_DIR       "/home/zmaw/m300362/tenstream/")
set(TENSTREAM_CMAKE     "/home/zmaw/m300362/tenstream/config/default.cmake")
set(TENSTREAM_BUILD_DIR "/home/zmaw/m300362/tenstream/build/")
