FILE(GLOB sourcefiles "*.?90")
add_executable(${PROJECT_NAME} ${sourcefiles})

include_directories( ${NETCDF_INCLUDE_DIR})
target_link_libraries(${PROJECT_NAME} ${NETCDF_LIBS})
