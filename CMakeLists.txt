cmake_minimum_required(VERSION 3.16)

# Project name
project(FastFFT_Project LANGUAGES CXX)

# Require C++17 or higher (necessary for advanced constexpr and std::array)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# 1. Define an "Interface" library (Header-only)
# This allows you to just "connect" the include/ anywhere
add_library(fft_logic INTERFACE)
target_include_directories(fft_logic INTERFACE 
    ${CMAKE_CURRENT_SOURCE_DIR}/include
)

# 2. Configure the test executable
add_executable(twiddles_test tests/twiddles_test.cpp)

# 3. Link the test with our logic library
# This automatically adds the 'include' directory to the compiler
target_link_libraries(twiddles_test PRIVATE fft_logic)

add_executable(fixed_point_tests tests/fixed_point_tests.cpp)
target_link_libraries(fixed_point_tests PRIVATE fft_logic)


# Optional: Enable compiler warnings (helps avoid type bugs)
if(MSVC)
    target_compile_options(twiddles_test PRIVATE /W4)
else()
    target_compile_options(twiddles_test PRIVATE -Wall -Wextra -Wpedantic)
endif()