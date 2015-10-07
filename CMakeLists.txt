project(matrix_mult)
cmake_minimum_required(VERSION 2.8)
add_executable(matrix_mult multip_mat.c print_mat.c create_mat.c main.c)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/Modules/")

set(os Linux)
set(mach 64)
#find_package(MKL REQUIRED)
#include_directories(${MKL_INCLUDE_DIRS})
#link_directories(${MKL_LIBRARIES})

set(LIBS ${LIBS} -Wl,--no-as-needed -lmkl_blas95_lp64 -lmkl_lapack95_lp64 -lmkl_intel_lp64 -lmkl_core -lmkl_sequential -lpthread -lm)

target_link_libraries(matrix_mult ${LIBS})
include (${CMAKE_ROOT}/Modules/FindBLAS.cmake)

