# ARCH Linux
set(CMAKE_C_COMPILER "mpcc_r")
set(CMAKE_CXX_COMPILER "mpCC_r")
set(CMAKE_Fortran_COMPILER "mpxlf2003_r")
set(Fortran_COMPILER_WRAPPER mpxlf2003_r)

set(USER_Fortran_FLAGS "-qfree=F90 -qrealsize=8  -qwarn64 -qnosave -qinitauto=FFF00000 -qflttrap=en:ov:zero:inv:imp -qflag=w:e -qsuffix=cpp=f90 -qxlf90=autodealloc -q64 ")
set(USER_Fortran_FLAGS_RELEASE "-O4 -qnoipa -qstrict=none:exceptions:zerosigns -qinitauto=ff -qsigtrap")
set(USER_Fortran_FLAGS_DEBUG "-O0 -qfullpath -C -g9 -qflttrp=enable:inexact:invalid:nanq:overflow:zerodivide -qsigtrap -qinitauto")

set(NETCDF_INCLUDE_DIR "/sw/aix61/netcdf-4.2.1.1/include")
set(NETCDF_LIB_1       "/sw/aix61/netcdf-4.2.1.1/lib/libnetcdff.a")
set(NETCDF_LIB_2       "/sw/aix61/netcdf-4.2.1.1/lib/libnetcdf.a")
set(HDF5_LIB_1         "/sw/aix61/hdf5-1.8.8/lib/libhdf5_hl.a")
set(HDF5_LIB_2         "/sw/aix61/hdf5-1.8.8/lib/libhdf5.a")
set(SZIP_LIB           "/sw/aix53/szip-2.1/lib/libsz.a")
set(LAPACK_LIB         "/sw/aix53/lapack-3.2.0/lib/liblapack.a")
include_directories(${NETCDF_INCLUDE_DIR})
set(LIBS ${NETCDF_LIB_1} ${NETCDF_LIB_2} ${HDF5_LIB_1} ${HDF5_LIB_2} ${SZIP_LIB} ${LAPACK_LIB} essl blas mass z dl )

# TENSTREAM LIB for 3D Radiation
set(TENSTREAM_DIR       "/pf/b/b380246/tenstream/")
set(TENSTREAM_BUILD_DIR "/pf/b/b380246/tenstream/build/")

# With optimizations advf.f90 breaks -- need -qstrict=zerosigns:
# set_source_files_properties(${CMAKE_CURRENT_SOURCE_DIR}/src/advf.f90 PROPERTIES COMPILE_FLAGS -qstrict=zerosigns)
