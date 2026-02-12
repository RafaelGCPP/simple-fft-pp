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

add_executable(butterfly_tests tests/butterfly_tests.cpp)
target_link_libraries(butterfly_tests PRIVATE fft_logic)

add_executable(fft_tests tests/fft_tests.cpp)
target_link_libraries(fft_tests PRIVATE fft_logic)

add_executable(roundtrip_tests tests/roundtrip_tests.cpp)
target_link_libraries(roundtrip_tests PRIVATE fft_logic)

add_executable(rfft_tests tests/rfft_tests.cpp)
target_link_libraries(rfft_tests PRIVATE fft_logic)


target_compile_options(twiddles_test PRIVATE -Wall -Wextra -Wpedantic)
target_compile_options(fixed_point_tests PRIVATE -Wall -Wextra -Wpedantic)
target_compile_options(butterfly_tests PRIVATE -Wall -Wextra -Wpedantic)
target_compile_options(fft_tests PRIVATE -Wall -Wextra -Wpedantic)
target_compile_options(roundtrip_tests PRIVATE -Wall -Wextra -Wpedantic)
target_compile_options(rfft_tests PRIVATE -Wall -Wextra -Wpedantic)