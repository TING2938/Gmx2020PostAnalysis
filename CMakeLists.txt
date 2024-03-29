project(Analysis CXX)
cmake_minimum_required(VERSION 2.8.0)

## source file name 
set(sourcefile $ENV{BaseFileName})

## install path
set(CMAKE_INSTALL_PREFIX $ENV{CODE_INSTALL_PATH})

## flags
set(CMAKE_POSITION_INDEPENDENT_CODE ON)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++17 -lstdc++fs -O3")

## options
option(USE_DEBUG "Debug this code." OFF)
option(USE_THREAD "Use thread library." ON)
option(USE_OMP "Use OpenMP." ON)

if (USE_DEBUG)
	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -gdwarf-2")
endif()

if (USE_THREAD)
	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -pthread")
endif()

if (USE_OMP)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fopenmp")
endif()

include_directories(${PROJECT_SOURCE_DIR}/include) 

add_subdirectory($ENV{DirName})


