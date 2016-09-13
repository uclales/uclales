# Mistral supercomputer at DKRZ/Hamburg

# TENSTREAM LIB for 3D Radiation
set(TENSTREAM_DIR       "$ENV{HOME}/tenstream/")
set(TENSTREAM_BUILD_DIR "$ENV{HOME}/tenstream/build/")
 
# modules loaded:
#   module load intel intelmpi
#   export MPIMODE=impi
#   export MYmpif90=mpiifort
#   export MYmpicc=mpiicc
#   export MYmpicxx=mpiicpc

# or with mellanox support
#   module load intel mxm fca bullxmpi_mlx
#   export MPIMODE=bullxmpi
#   export MYmpif90=mpif90
#   export MYmpicc=mpicc
#   export MYmpicxx=mpicxx

# Netcdf paths taken from nc-config // nf-config
#
# export NETCDFCROOT = /sw/rhel6-x64/netcdf/netcdf_c-4.3.2-parallel-$MPIMODE-intel14/      
# export NETCDFFROOT = /sw/rhel6-x64/netcdf/netcdf_fortran-4.4.2-parallel-$MPIMODE-intel14/
# export HDF5ROOT    = /sw/rhel6-x64/hdf5/hdf5-1.8.14-parallel-$MPIMODE-intel14/
# export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$NETCDFFROOT/lib:$NETCDFCROOT/lib:

set(NETCDFC_DIR      "$ENV{NETCDFCROOT}")
set(NETCDFF_DIR      "$ENV{NETCDFFROOT}")
set(HDF5_DIR         "$ENV{HDF5ROOT}")

execute_process(COMMAND "${NETCDFF_DIR}/bin/nf-config" "--includedir" OUTPUT_VARIABLE NETCDF_INCLUDE_DIR)
execute_process(COMMAND "${NETCDFF_DIR}/bin/nf-config" "--flibs"      OUTPUT_VARIABLE NETCDFF_LIB)
execute_process(COMMAND "${NETCDFC_DIR}/bin/nc-config" "--libs"       OUTPUT_VARIABLE NETCDFC_LIB)

string(STRIP ${NETCDFF_LIB} NETCDFF_LIB)
string(STRIP ${NETCDFC_LIB} NETCDFC_LIB)

list(APPEND LIBS "${NETCDFF_LIB} ${NETCDFC_LIB}" )

# Petsc installed with:
# export PETSC_ARCH=fast-$MPIMODE
#./configure                           \
#  --with-cc=$(which  $MYmpicc)        \
#  --with-fc=$(which  $MYmpif90)       \
#  --with-cxx=$(which $MYmpicxx)       \
#  --with-shared-libraries=0           \
#  --with-netcdf-dir=$NETCDFCROOT      \
#  --with-hdf5-dir=$HDF5ROOT           \
#  --with-blas-lapack-dir=$MKL         \
#  --with-fortran                      \
#  --with-fortran-interfaces           \
#  --with-cmake=$(which cmake)         \
#  --download-ml                       \
#  --download-hypre                    \
#  --with-debugging=0                  \
#  --with-valgrind-dir=/usr            \
#  COPTFLAGS='-O3 -xAVX -mkl' \
#  FOPTFLAGS='-O3 -xAVX -mkl' \
#  \
#  && make all test install

set(CMAKE_C_COMPILER       "$ENV{MYmpicc}")
set(CMAKE_CXX_COMPILER     "$ENV{MYmpicxx}")
set(CMAKE_Fortran_COMPILER "$ENV{MYmpif90}")

set(USER_C_FLAGS "-std=c99 ")
set(USER_Fortran_FLAGS "-fpp -traceback -r8 -ftz -extend_source -g -no-wrap-margin")
set(USER_Fortran_FLAGS_RELEASE " -O3 -xCORE-AVX2 -fp-model source ")
set(USER_Fortran_FLAGS_DEBUG "-fpe0 -O0 -g -check all -check nopointers -check noarg_temp_created")

