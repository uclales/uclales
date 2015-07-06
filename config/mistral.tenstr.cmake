# Mistral supercomputer at DKRZ/Hamburg

# modules loaded:
#   module load intel intelmpi
#   module load intel mxm fca bullxmpi_mlx

# TENSTREAM LIB for 3D Radiation
set(TENSTREAM_DIR       "$ENV{HOME}/tenstream/")
set(TENSTREAM_BUILD_DIR "$ENV{HOME}/tenstream/build/")
 
# Netcdf paths taken from nc-config // nf-config
#
# export NETCDFCROOT = /sw/rhel6-x64/netcdf/netcdf_c-4.3.2-parallel-bullxmpi-intel14/      
# export NETCDFFROOT = /sw/rhel6-x64/netcdf/netcdf_fortran-4.4.2-parallel-bullxmpi-intel14/
# export HDF5ROOT    = /sw/rhel6-x64/hdf5/hdf5-1.8.14-parallel-bullxmpi-intel14/

set(NETCDFC_DIR      "$ENV{NETCDFCROOT}")
set(NETCDFF_DIR      "$ENV{NETCDFFROOT}")
set(HDF5_DIR         "$ENV{HDF5ROOT}")

# Petsc installed with:
# export PETSC_ARCH=fast
#./configure                           \
#  --prefix="$PETSC_PREFIX/$PETSC_ARCH"\
#  --with-cc=$(which mpicc)            \
#  --with-fc=$(which mpif90)           \
#  --with-cxx=$(which mpicxx)          \
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
#  COPTFLAGS='-O3 -no-prec-div -xAVX -mkl -fno-omit-frame-pointer' \
#  FOPTFLAGS='-O3 -no-prec-div -xAVX -mkl -fno-omit-frame-pointer' \
#  \
#  && make all test install

set(CMAKE_C_COMPILER "mpicc")
set(CMAKE_CXX_COMPILER "mpicxx")
set(CMAKE_Fortran_COMPILER "mpif90")

set(USER_C_FLAGS "-std=c99 ")
set(USER_Fortran_FLAGS "-fpp -traceback -r8 -ftz -extend_source -g ")
set(USER_Fortran_FLAGS_RELEASE " -O3 -no-prec-div -xAVX -fp-model source ")
set(USER_Fortran_FLAGS_DEBUG "-fpe0 -O0 -g -check all -check nopointers -check noarg_temp_created")

execute_process(COMMAND "${NETCDFF_DIR}/bin/nf-config" "--includedir" OUTPUT_VARIABLE NETCDF_INCLUDE_DIR)

execute_process(COMMAND "${NETCDFF_DIR}/bin/nf-config" "--flibs"      OUTPUT_VARIABLE NETCDFF_LIB)
string(STRIP ${NETCDFF_LIB} NETCDFF_LIB )

execute_process(COMMAND "${NETCDFC_DIR}/bin/nc-config" "--libs"       OUTPUT_VARIABLE NETCDFC_LIB)
string(STRIP ${NETCDFC_LIB} NETCDFC_LIB )

list(APPEND LIBS "${NETCDFF_LIB} ${NETCDFC_LIB}" )

