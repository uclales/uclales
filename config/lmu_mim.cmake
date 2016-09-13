# Default GCC
set(CMAKE_Fortran_COMPILER   "mpif90")
set(Fortran_COMPILER_WRAPPER "mpif90")

set(USER_Fortran_FLAGS "-fdefault-real-8 -fdefault-double-8 -ffree-line-length-none") #
set(USER_Fortran_FLAGS_RELEASE "-O3 ") #-march=native -mtune=native")
set(USER_Fortran_FLAGS_DEBUG "-fbacktrace -finit-real=nan -W -Wall -Wuninitialized -g -pg -fcheck=all -fbounds-check -ffpe-trap=invalid,zero,overflow") #")

# TENSTREAM LIB for 3D Radiation
set(TENSTREAM_DIR       "~/tenstream/")
set(TENSTREAM_BUILD_DIR "~/tenstream/build/")

execute_process(COMMAND "nf-config" "--includedir" OUTPUT_VARIABLE NETCDF_INCLUDE_DIR)
execute_process(COMMAND "nf-config" "--flibs"      OUTPUT_VARIABLE NETCDFF_LIB)
execute_process(COMMAND "nc-config" "--libs"       OUTPUT_VARIABLE NETCDFC_LIB)

string(STRIP ${NETCDFF_LIB} NETCDFF_LIB)
string(STRIP ${NETCDFC_LIB} NETCDFC_LIB)

list(APPEND LIBS "${NETCDFF_LIB} ${NETCDFC_LIB}" )
