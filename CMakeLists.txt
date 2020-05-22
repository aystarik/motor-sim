cmake_minimum_required(VERSION 3.5)

project(motor-sim LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
include_directories(Ascent/include)

add_executable(motor-sim main.cpp)
