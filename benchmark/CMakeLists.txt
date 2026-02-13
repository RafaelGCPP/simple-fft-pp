add_executable(benchmark benchmark.cpp)
target_link_libraries(benchmark PRIVATE fft_logic)


target_compile_options(benchmark PRIVATE -Wall -Wextra -Wpedantic)