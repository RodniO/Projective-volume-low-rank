cmake_minimum_required(VERSION 4.0)

project(pvlr VERSION 0.0.0 LANGUAGES CXX Fortran)

option(BUILD_SHARED_LIBS "Build using shared libraries" ON)
option(BUILD_TESTS "Build tests or not." OFF)

option(BLA_STATIC "BLAS/LAPACK static linkage." OFF)
option(BLA_VENDOR "BLAS/LAPACK vendors." "Generic")

set(CMAKE_CONFIGURATION_TYPES "Debug;MinSize;Release;RelWithDebInfo" CACHE
    STRING "Available build configurations" FORCE)

message(STATUS "Available build configurations: ${CMAKE_CONFIGURATION_TYPES}")

set(CMAKE_COLOR_DIAGNOSTICS ON)
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED YES)
set(CMAKE_EXPERIMENTAL_CXX_MODULE_CMAKE_API ON)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

if(CMAKE_CXX_COMPILER_ID MATCHES "Clang|GNU")
    set(CXX_COMPILE_OPTIONS
        "-Wall;-Wextra;-Wpedantic;-Wconversion;-Wshadow;-Wnull-dereference")
elseif(MSVC)
    set(CXX_COMPILE_OPTIONS "/W4 /permissive-")
endif()

find_package(BLAS REQUIRED)
find_package(LAPACK REQUIRED)

message(STATUS "BLAS_LINKER_FLAGS:   ${BLAS_LINKER_FLAGS}")
message(STATUS "BLAS_LIBRARIES:      ${BLAS_LIBRARIES}")
message(STATUS "LAPACK_LINKER_FLAGS: ${LAPACK_LINKER_FLAGS}")
message(STATUS "LAPACK_LIBRARIES:    ${LAPACK_LIBRARIES}")

add_library(pvlr)

target_sources(pvlr PRIVATE
    src/MerTwist.f90
    src/ModAppr.f90
    src/ModIntVec.f90
    src/ModMonte.f90
    src/ModMtrx.f90
    src/ModRecon.f90
    src/ModSparse.f90
    src/ModVec.f90
    src/ModPVLR.f90
)

target_include_directories(pvlr PRIVATE src)

target_link_libraries(pvlr PUBLIC BLAS::BLAS LAPACK::LAPACK)

add_library(pvlr-cc STATIC)

target_sources(pvlr-cc
    PUBLIC FILE_SET CXX_MODULES FILES
        pvlr/pvlr.ccm
    PUBLIC FILE_SET HEADERS FILES
        pvlr/pvlr.h
)

target_link_libraries(pvlr-cc PUBLIC pvlr)

target_compile_options(pvlr-cc PRIVATE "${CXX_COMPILE_OPTIONS}")

if (BUILD_TESTS)
    add_executable(main src/Main.f90)

    target_link_libraries(main PRIVATE pvlr)

    add_executable(pvlr_test)

    target_sources(pvlr_test PRIVATE pvlr/pvlr_test.cc)

    target_compile_options(pvlr_test PRIVATE "${CXX_COMPILE_OPTIONS}")

    target_link_libraries(pvlr_test PRIVATE pvlr-cc)
endif()
