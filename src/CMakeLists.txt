FILE(GLOB sourcefiles "*.?90")
set_source_files_properties(${sourcefiles} PROPERTIES COMPILE_FLAGS -cpp)
add_executable(${PROJECT_NAME} ${sourcefiles})
include_directories(${INCLUDE_DIRS})
include_directories (${CMAKE_CURRENT_SOURCE_DIR})
target_link_libraries(${PROJECT_NAME} ${LIBS})

set_property(TARGET ${PROJECT_NAME} PROPERTY INSTALL_RPATH_USE_LINK_PATH TRUE)
install(TARGETS ${PROJECT_NAME} RUNTIME DESTINATION bin)
