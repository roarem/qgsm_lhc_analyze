cmake_minimum_required(VERSION 2.8)

project("histograms")

set(GCC_COMPILE_FLAGS "--std=c++11")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${GCC_COMPILE_FLAGS}")


set(SOURCES 
  Histogram.cpp
  )

set(HEADERS 
  Histogram.h
  )

include_directories(${PROJECT_SOURCE_DIR}/lib)
add_executable(hist main.cpp ${SOURCES} ${HEADERS})

