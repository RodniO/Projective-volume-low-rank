cmake_minimum_required(VERSION 4.0)

project(pvll VERSION 0.0.0 LANGUAGES Fortran)

option(BUILD_SHARED_LIBS "Build using shared libraries" ON)
option(BUILD_TESTS "Build tests or not." OFF)

option(BLA_STATIC "BLAS/LAPACK static linkage." OFF)
option(BLA_VENDOR "BLAS/LAPACK vendors." "Generic")

set(CMAKE_CONFIGURATION_TYPES "Debug;MinSize;Release;RelWithDebInfo" CACHE
    STRING "Available build configurations" FORCE)

message(STATUS "Available build configurations: ${CMAKE_CONFIGURATION_TYPES}")

set(CMAKE_COLOR_DIAGNOSTICS ON)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

find_package(BLAS REQUIRED)
find_package(LAPACK REQUIRED)

message(STATUS "BLAS_LINKER_FLAGS:   ${BLAS_LINKER_FLAGS}")
message(STATUS "BLAS_LIBRARIES:      ${BLAS_LIBRARIES}")
message(STATUS "LAPACK_LINKER_FLAGS: ${LAPACK_LINKER_FLAGS}")
message(STATUS "LAPACK_LIBRARIES:    ${LAPACK_LIBRARIES}")

add_library(pvll)

target_sources(pvll PRIVATE
    src/MerTwist.f90
    src/ModAppr.f90
    src/ModIntVec.f90
    src/ModMonte.f90
    src/ModMtrx.f90
    src/ModRecon.f90
    src/ModSparse.f90
    src/ModVec.f90
)

target_include_directories(pvll PRIVATE src)

target_link_libraries(pvll PUBLIC BLAS::BLAS LAPACK::LAPACK)

if (BUILD_TESTS)
    add_executable(main src/Main.f90)

    target_link_libraries(main PRIVATE pvll)
endif()
