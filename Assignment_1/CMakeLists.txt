PROJECT(spring CXX)
SET(SOURCES src/forces.cc  src/initialize.cc  src/integrators.cc  src/main.cc)
INCLUDE_DIRECTORIES(include)
ADD_EXECUTABLE(spring ${SOURCES})
