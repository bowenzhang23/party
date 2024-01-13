cmake_minimum_required(VERSION 3.20)

project(
  miniroot
  VERSION 0.1
  DESCRIPTION "Root file reader"
  LANGUAGES CXX)

if(NOT DEFINED CMAKE_CXX_STANDARD)
  message(STATUS "Using c++20 standard")
  set(CMAKE_CXX_STANDARD 20)
  set(CMAKE_CXX_STANDARD_REQUIRED ON)
  set(CMAKE_CXX_EXTENSIONS ON)
endif()

find_package(
  Python 3.10
  COMPONENTS Interpreter Development.Module
  REQUIRED)
find_package(ZLIB REQUIRED)

message(STATUS "Building ${CMAKE_PROJECT_NAME}::${PROJECT_NAME} ...")

include_directories(include/)

file(GLOB MINIROOT_SRC src/*.cpp)
file(GLOB MINIROOT_INC src/*.hpp)

add_library(${PROJECT_NAME} SHARED ${MINIROOT_SRC} ${MINIROOT_INC})

nanobind_add_module(pyminiroot bind/module.cpp)
target_link_libraries(pyminiroot PRIVATE ${PROJECT_NAME} ZLIB::ZLIB)

install(
  TARGETS ${PROJECT_NAME} pyminiroot
  LIBRARY DESTINATION lib
  ARCHIVE DESTINATION lib)

install(FILES ${MINIROOT_INC} DESTINATION include/{PROJECT_NAME})

add_executable(Miniroot_Test test/main.cpp)
target_link_libraries(Miniroot_Test PRIVATE ${PROJECT_NAME} ZLIB::ZLIB)

install(TARGETS Miniroot_Test RUNTIME DESTINATION bin)

add_test(NAME Miniroot_Testing
         COMMAND Miniroot_Test ${CMAKE_SOURCE_DIR}/resources/example.root)
