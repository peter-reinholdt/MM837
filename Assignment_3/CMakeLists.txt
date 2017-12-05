cmake_minimum_required(VERSION 3.9)
PROJECT(lattice CXX)

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -march=native -Wall")


SET(SOURCES src/solvers.cc  src/properties.cc  src/main.cc )
INCLUDE_DIRECTORIES(include)
ADD_EXECUTABLE(lattice ${SOURCES})

