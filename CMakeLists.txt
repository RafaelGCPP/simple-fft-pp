cmake_minimum_required(VERSION 3.16)

# Project name
project(FastFFT_Project LANGUAGES CXX)

# Require C++20 (necessary for constexpr, if constexpr, std::array, etc.)
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# 1. Define an "Interface" library (Header-only)
# This allows you to just "connect" the include/ anywhere
add_library(fft_logic INTERFACE)
target_include_directories(fft_logic INTERFACE 
    ${CMAKE_CURRENT_SOURCE_DIR}/include
)

if(BUILD_TESTSET)
    add_subdirectory(tests)
endif()