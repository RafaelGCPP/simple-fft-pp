add_executable(fixed_point_tests fixed_point_tests.cpp)
target_link_libraries(fixed_point_tests PRIVATE fft_logic)

add_executable(fft_tests fft_tests.cpp)
target_link_libraries(fft_tests PRIVATE fft_logic)

add_executable(roundtrip_tests roundtrip_tests.cpp)
target_link_libraries(roundtrip_tests PRIVATE fft_logic)

add_executable(rfft_tests rfft_tests.cpp)
target_link_libraries(rfft_tests PRIVATE fft_logic)


if (CMAKE_CXX_COMPILER_ID STREQUAL "GNU" OR CMAKE_CXX_COMPILER_ID MATCHES "Clang")
    target_compile_options(fixed_point_tests PRIVATE -Wall -Wextra -Wpedantic) 
    target_compile_options(fft_tests PRIVATE -Wall -Wextra -Wpedantic)
    target_compile_options(roundtrip_tests PRIVATE -Wall -Wextra -Wpedantic)
    target_compile_options(rfft_tests PRIVATE -Wall -Wextra -Wpedantic)
endif()