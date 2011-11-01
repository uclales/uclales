set (ENV{NETCDF_ROOT} /sw/aix61/netcdf-4.1.2-hdf5-threadsafe)
set (ENV{HDF5_ROOT} /sw/aix61/hdf5-1.8.6-threadsafe)
set (ENV{FFTW_ROOT} $ENV{HOME}/fftw)
find_library(LAPACK_LIBRARIES essl)
ADD_DEFINITIONS("-DESSL=TRUE")  # Set precompiler flag