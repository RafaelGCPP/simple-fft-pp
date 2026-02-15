add_executable(rfft_benchmark rfft_benchmark.cpp)
target_link_libraries(rfft_benchmark PRIVATE fft_logic)

add_executable(fft_benchmark fft_benchmark.cpp)
target_link_libraries(fft_benchmark PRIVATE fft_logic)


if (CMAKE_CXX_COMPILER_ID STREQUAL "GNU" OR CMAKE_CXX_COMPILER_ID MATCHES "Clang")
    target_compile_options(rfft_benchmark PRIVATE -Wall -Wextra -Wpedantic -O3) 
    target_compile_options(fft_benchmark PRIVATE -Wall -Wextra -Wpedantic -O3)
endif()