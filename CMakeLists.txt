cmake_minimum_required(VERSION 3.15)

project(circle3_exp LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)

find_package(PCL 1.10 REQUIRED COMPONENTS common io visualization registration segmentation filters)
find_package(Eigen3 3 REQUIRED)

if(CMAKE_BUILD_TYPE MATCHES Debug)
    message("Debug mode")
    add_definitions(-DDEBUG)
    add_compile_options(-g -O0)
else()
    message("Release mode")
    add_definitions(-DNDEBUG)
    add_compile_options(-O3 -fopenmp)
endif()

include_directories(${PCL_INCLUDE_DIRS})
include_directories(${EIGEN3_INCLUDE_DIRS})
include_directories(scancontext)
include_directories(${CMAKE_SOURCE_DIR})


add_executable(trace2 trace2.cpp scancontext/Scancontext.cpp  files.cpp registration.cpp datasource.cpp)
target_link_libraries(trace2 -fopenmp ${PCL_LIBRARIES} ${EIGEN3_LIBRARIES})
