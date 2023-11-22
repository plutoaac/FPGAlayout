cmake_minimum_required(VERSION 2.8.12.2)

project(layout_routing)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")


set(CMAKE_BUILD_TYPE Debug)

set(root ${CMAKE_SOURCE_DIR})


include_directories(${root}/include)

file(GLOB SOURCES ${root}/src/*.cpp)



add_executable(layout_routing ${SOURCES})