find_package(Python 3.12 COMPONENTS Interpreter Development.Module REQUIRED)

# Detect the installed nanobind package and import it into CMake
execute_process(
    COMMAND "${Python_EXECUTABLE}" -m nanobind --cmake_dir
    OUTPUT_STRIP_TRAILING_WHITESPACE OUTPUT_VARIABLE nanobind_ROOT)
find_package(nanobind CONFIG REQUIRED)

nanobind_add_module(extension NOMINSIZE)

target_sources(extension PRIVATE extension.cc)

target_include_directories(extension PRIVATE ${PROJECT_SOURCE_DIR})

target_link_libraries(extension PRIVATE pvlr-cc)

nanobind_add_stub(extension_stub
    MODULE extension
    OUTPUT extension.pyi
    PYTHON_PATH $<TARGET_FILE_DIR:extension>
    DEPENDS extension
    EXCLUDE_FROM_ALL
)

# Install Python extension if it is build with scikit-build build system. We
# expect that BUILD_SHARED_LIBS=OFF, i.e. `extension.*.so` is self-contained.
if (DEFINED SKBUILD_PROJECT_NAME)
    install(TARGETS extension DESTINATION ${SKBUILD_PROJECT_NAME})
endif()
