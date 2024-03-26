cmake_minimum_required(VERSION 2.8.12.2)


project(layout_routing)
find_package(OpenMP REQUIRED)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++20")



set(root ${CMAKE_SOURCE_DIR})


include_directories(${root}/include)

file(GLOB SOURCES ${root}/src/*.cpp)



add_executable(layout_routing ${SOURCES})