# Opt-in research build. Does not change the installed production library.
cmake_minimum_required(VERSION 3.16)
project(elips_pairing_research LANGUAGES C)

set(CMAKE_C_STANDARD 11)
set(CMAKE_C_STANDARD_REQUIRED ON)
set(CMAKE_C_EXTENSIONS OFF)
get_filename_component(ELIPS_ROOT "${CMAKE_CURRENT_LIST_DIR}/../../.." ABSOLUTE)
add_subdirectory("${ELIPS_ROOT}" "${CMAKE_CURRENT_BINARY_DIR}/elips" EXCLUDE_FROM_ALL)

add_executable(normalized_line_experiment normalized-line/experiment.c)
add_executable(scaled_square_experiment scaled-square/experiment.c)
add_executable(primitive_ratios norm-and-exponent/ratios.c)
foreach(target normalized_line_experiment scaled_square_experiment primitive_ratios)
  target_link_libraries(${target} PRIVATE elips_arith_BLS12_381)
  target_compile_definitions(${target} PRIVATE _POSIX_C_SOURCE=200809L)
  target_compile_options(${target} PRIVATE -Wall -Wextra -Werror)
endforeach()

enable_testing()
add_test(NAME research.normalized_line COMMAND normalized_line_experiment)
add_test(NAME research.scaled_square COMMAND scaled_square_experiment)
set_tests_properties(research.normalized_line research.scaled_square
                     PROPERTIES TIMEOUT 300)
