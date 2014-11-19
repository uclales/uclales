# Default GCC
set(CMAKE_Fortran_COMPILER "/home/opt/cosmo_tica_lib/ompi1.8.1/openmpi-1.8.1/install/bin/mpif90")
set(Fortran_COMPILER_WRAPPER "/home/opt/cosmo_tica_lib/ompi1.8.1/openmpi-1.8.1/install/bin/mpif90")

set(USER_Fortran_FLAGS "-fdefault-real-8 -ffree-line-length-none") #
set(USER_Fortran_FLAGS_RELEASE "-O3 ") #-march=native -mtune=native")
set(USER_Fortran_FLAGS_DEBUG "-fbacktrace -finit-real=nan -W -Wall -Wuninitialized -g -pg -fcheck=all -fbounds-check -ffpe-trap=invalid,zero,overflow") #")

set(NETCDF_INCLUDE_DIR "/home/opt/cosmo_tica_lib//ompi1.8.1/netcdf-fortran-4.2/install/include/")
set(NETCDF_LIB_1       "/home/opt/cosmo_tica_lib//ompi1.8.1/netcdf-fortran-4.2/install/lib/libnetcdff.a")
set(NETCDF_LIB_2       "/home/opt/cosmo_tica_lib//ompi1.8.1/netcdf-4.3.0/install/lib64/libnetcdf.a")
set(HDF5_LIB_1         "/home/opt/cosmo_tica_lib//ompi1.8.1/hdf5/HDF5-1.8.13-Linux/HDF_Group/HDF5/1.8.13/lib/libhdf5_hl.so")
set(HDF5_LIB_2         "/home/opt/cosmo_tica_lib//ompi1.8.1/hdf5/HDF5-1.8.13-Linux/HDF_Group/HDF5/1.8.13/lib/libhdf5.so")
set(SZIP_LIB           "/home/opt/cosmo_tica_lib//ompi1.8.1/hdf5/HDF5-1.8.13-Linux/HDF_Group/HDF5/1.8.13/lib/libszip.so")

set(LIBS ${NETCDF_LIB_1} ${NETCDF_LIB_2} ${HDF5_LIB_1} ${HDF5_LIB_2} ${SZIP_LIB} m z curl)

# TENSTREAM LIB for 3D Radiation
set(TENSTREAM_DIR       "/home/users/jakub/hdcp2-rt/tenstream/")
set(TENSTREAM_BUILD_DIR "/home/users/jakub/hdcp2-rt/tenstream/build/")

