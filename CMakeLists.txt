project(Analysis CXX)
cmake_minimum_required(VERSION 2.8.0)

## source file name 
set(sourcefile $ENV{BaseFileName})

## install path
set(CMAKE_INSTALL_PREFIX $ENV{CODE_INSTALL_PATH})

## flags
set(CMAKE_POSITION_INDEPENDENT_CODE ON)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++17 -O3")

## options
option(USE_DEBUG "Debug this code." OFF)
option(USE_GTEST "Build test." OFF)
option(USE_THREAD "Use thread library." OFF)
option(USE_OMP "Use OpenMP." OFF)

if (USE_DEBUG)
	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -gdwarf-2")
endif()

if (USE_THREAD)
	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -pthread")
endif()

if (USE_OMP)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fopenmp")
endif()

if (USE_GTEST)
  add_subdirectory(dependencies/googletest EXCLUDE_FROM_ALL)
endif ()

include_directories(${PROJECT_SOURCE_DIR}/Core/include) 

add_subdirectory($ENV{DirName})


