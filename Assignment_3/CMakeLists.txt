cmake_minimum_required(VERSION 3.9)
PROJECT(lattice CXX)

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -march=native  -Wall -ffast-math -funsafe-loop-optimizations -fstrict-aliasing -fomit-frame-pointer ")


SET(SOURCES src/solvers.cc  src/properties.cc  src/main.cc src/xoroshiro128plus.cc)
INCLUDE_DIRECTORIES(include)
ADD_EXECUTABLE(lattice ${SOURCES})

