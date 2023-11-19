cmake_minimum_required(VERSION 3.10)

project(layout_routing)

set(CMAKE_CXX_STANDARD 17)

set(CMAKE_BUILD_TYPE Debug)

set(root ${CMAKE_SOURCE_DIR})


include_directories(${root}/include)

file(GLOB SOURCES ${root}/src/*.cpp)



add_executable(layout_routing ${SOURCES})