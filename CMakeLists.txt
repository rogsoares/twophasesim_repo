cmake_minimum_required (VERSION 3.0)
project (PADMEC_EXE)

set (APPLICATIONS /Users/rogerio/programacao/scientific-applications)
set(CMAKE_BUILD_TYPE Release)
set (CMAKE_CXX_FLAGS "-g -Wall -std=c++14 -DMATRIX_DEBUG" )
#set (CMAKE_CXX_FLAGS "-O3 -std=c++14")

set(PETSC_DIR "${APPLICATIONS}/petsc")
set(PETSC_INCLUDE "${PETSC_DIR}/include" "${PETSC_DIR}/c-opt/include" CACHE PATH "Directory where PETSc header files can be found")
set(PETSC_LIB "${PETSC_DIR}/c-opt/lib" CACHE PATH "Directory where the PETSc library can be found" )

#set(PETSC_INCLUDE "${PETSC_DIR}/include" "${PETSC_DIR}/debugging_build/include" CACHE PATH "Directory where PETSc header files can be found")
#set(PETSC_LIB "${PETSC_DIR}/debugging_build/lib" CACHE PATH "Directory where the PETSc library can be found" )

include_directories(${PETSC_INCLUDE})
link_directories(${PETSC_LIB})
link_libraries(petsc)


include_directories(
                    "${PROJECT_SOURCE_DIR}/include"
                    "${PROJECT_SOURCE_DIR}/src/simulator"
                    "${PROJECT_SOURCE_DIR}/src/parmesh"
                    )
                    
add_subdirectory(src/simulator)
add_subdirectory(src/parmesh)
add_subdirectory(src)

