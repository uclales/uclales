include_directories(${UCLALES_INCLUDE_DIRS})
FILE(GLOB sourcefiles "*.?90")
add_executable(${PROJECT_NAME} ${sourcefiles})
target_link_libraries(${PROJECT_NAME} ${UCLALES_LIBRARIES})
