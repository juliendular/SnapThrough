PROJECT(SnapThrough CXX)
cmake_minimum_required(VERSION 2.7)
include(CheckCXXCompilerFlag)

IF(NOT CMAKE_BUILD_TYPE)
    SET( CMAKE_BUILD_TYPE "Release" CACHE STRING "" FORCE)
ENDIF()

IF(MSVC)
    add_definitions(-D_USE_MATH_DEFINES -D_CRT_SECURE_NO_WARNINGS -D_SCL_SECURE_NO_DEPRECATE)
ENDIF(MSVC)

check_cxx_compiler_flag("-std=c++11" STD11CHECK)
IF(STD11CHECK)
    ADD_DEFINITIONS(-std=c++11 )
ENDIF()

FILE(GLOB SRCS1 *.cpp)

ADD_EXECUTABLE(main ${SRCS1} main.cpp)

IF(MINGW)
    TARGET_LINK_LIBRARIES(main psapi) # for "GetProcessMemoryInfo"
ENDIF(MINGW)
