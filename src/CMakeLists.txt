# src CMakeLists.txt
# ==================================

set(simulator_sources ${simulator_STAT_SRCS})
add_library(simulatorlib ${simulator_sources})

set(parmesh_sources ${parmesh_STAT_SRCS})
add_library(parmeshlib ${parmesh_sources})


file(GLOB sources *.cpp)
set(SRCS ${sources} ${SRCS})

file(GLOB headers *.h)
set(HDRS ${headers} ${HDRS})


add_executable(PADMEC_EXE main.cpp ${HDRS})
target_link_libraries(PADMEC_EXE simulatorlib parmeshlib)