cmake_minimum_required(VERSION 3.9)
PROJECT(spring CXX)

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O2 -march=native")


SET(SOURCES src/forces.cc  src/initialize.cc  src/integrators.cc  src/main.cc)
INCLUDE_DIRECTORIES(include)
ADD_EXECUTABLE(spring ${SOURCES})

