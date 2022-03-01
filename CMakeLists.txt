cmake_minimum_required(VERSION 3.20)

project(crypto CXX)

set(CMAKE_CXX_STANDARD 23)
set(CMAKE_CXX_STANDARD_REQUIRED True)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra -pedantic")

set(SOURCES
    BigInt/BigInt.cpp

    BigInt/Algorithms/Multiplication.cpp
    BigInt/Algorithms/Division.cpp

    EllipticCurve/EllipticCurve.cpp

    Modmath.cpp

    REPL.cpp
)

include_directories(${CMAKE_SOURCE_DIR})

add_executable(crypto ${SOURCES})
