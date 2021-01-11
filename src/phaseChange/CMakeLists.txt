
find_package(GROMACS 2020.1 REQUIRED)
include_directories(${GROMACS_INCLUDE_DIRS})
add_definitions(${GROMACS_DEFINITIONS})

add_executable(${sourcefile} ${sourcefile})
set_target_properties(${sourcefile} PROPERTIES
                      COMPILE_FLAGS "${GROMACS_CXX_FLAGS}"
					  LINKER_LANGUAGE CXX)
target_link_libraries(${sourcefile} ${GROMACS_LIBRARIES})

install(TARGETS ${sourcefile} DESTINATION bin)
