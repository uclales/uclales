FILE(GLOB sourcefiles "*.?90")
add_executable(${PROJECT_NAME} ${sourcefiles})
include_directories(${INCLUDE_DIRS})
target_link_libraries(${PROJECT_NAME} ${LIBS})

# From MicroHH: define githash, so we can print it in de model output
# retrieve the git hash from the current commit
find_package(Git)
if(GIT_EXECUTABLE)
  exec_program(git ${CMAKE_CURRENT_SOURCE_DIR} ARGS "describe --dirty" OUTPUT_VARIABLE GITHASH RETURN_VALUE GITRETURN)
  if(NOT GITRETURN MATCHES 0)
    set(GITHASH "not available")
  endif()
else()
  set(GITHASH "not available")
endif()

message(STATUS "Git hash " ${GITHASH})
add_definitions(-DGITHASH="${GITHASH}")
