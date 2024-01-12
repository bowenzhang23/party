cmake_minimum_required(VERSION 3.20)

project(
  party
  VERSION 0.3
  DESCRIPTION
    "Minimal particle physics toolkit"
  LANGUAGES CXX)

message(STATUS "Building ${CMAKE_PROJECT_NAME} ...")

enable_testing()

# external
add_subdirectory(ext/nanobind)

# internal
add_subdirectory(miniroot)

set(CMAKE_INSTALL_PREFIX
    ${CMAKE_CURRENT_BINARY_DIR}/${CMAKE_CXX_LIBRARY_ARCHITECTURE}-${CMAKE_CXX_COMPILER_VERSION}
)

