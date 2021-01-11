
find_package(GROMACS 2018.1 REQUIRED)
include_directories(${PROJECT_SOURCE_DIR}/Core/include)
include_directories(${PROJECT_SOURCE_DIR}/Gmx1802/include)
include_directories(${GROMACS_INCLUDE_DIRS})
add_definitions(${GROMACS_DEFINITIONS})

add_executable(${sourcefile} ${sourcefile})
set_target_properties(${sourcefile} PROPERTIES
                      COMPILE_FLAGS "${GROMACS_CXX_FLAGS}")
target_link_libraries(${sourcefile} ${GROMACS_LIBRARIES})

install(TARGETS ${sourcefile} DESTINATION bin)
