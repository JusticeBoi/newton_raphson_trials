cmake_minimum_required (VERSION 3.5)

project(newton_ralphson)


set ( CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++17 -Wall -O3 -Wfloat-equal -Wextra -march=native ")


set(lib_files
        newton_ralphson.cpp
        )

#add_library(lib ${lib_files})


set ( CMAKE_BINARY_DIR ${CMAKE_SOURCE_DIR}/bin)

set ( EXECUTABLE_OUTPUT_PATH ${CMAKE_BINARY_DIR})

set ( source_dir "${PROJECT_SOURCE_DIR}/src/")

include_directories("${PROJECT_SOURCE_DIR}") 




add_executable  (  newton_ralphson  src/newton_ralphson.cpp  )
