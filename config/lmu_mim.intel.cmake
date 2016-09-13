# Config for setup at LMU Munich using intelmpi
#
# Following modules should be loaded:
#    module load intelmpi icc
#    module load netcdf/4.4.0-icc-16 hdf5
#

set(CMAKE_Fortran_COMPILER   "mpiifort")
set(Fortran_COMPILER_WRAPPER "mpiifort")

set(USER_C_FLAGS "-std=c99 ")
set(USER_Fortran_FLAGS "-fpp -traceback -r8 -ftz -extend_source -g -no-wrap-margin")
set(USER_Fortran_FLAGS_RELEASE " -O3 -xHost -fp-model source ")
set(USER_Fortran_FLAGS_DEBUG "-fpe0 -O0 -g -check all -check nopointers -check noarg_temp_created")

# TENSTREAM LIB for 3D Radiation
set(TENSTREAM_DIR       "~/tenstream/")
set(TENSTREAM_BUILD_DIR "~/tenstream/build/")

execute_process(COMMAND "nf-config" "--includedir" OUTPUT_VARIABLE NETCDF_INCLUDE_DIR)
execute_process(COMMAND "nf-config" "--flibs"      OUTPUT_VARIABLE NETCDFF_LIB)
execute_process(COMMAND "nc-config" "--libs"       OUTPUT_VARIABLE NETCDFC_LIB)

string(STRIP ${NETCDFF_LIB} NETCDFF_LIB)
string(STRIP ${NETCDFC_LIB} NETCDFC_LIB)

list(APPEND LIBS "${NETCDFF_LIB} ${NETCDFC_LIB}" )
