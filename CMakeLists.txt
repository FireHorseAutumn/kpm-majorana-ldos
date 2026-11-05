cmake_minimum_required(VERSION 3.20)
project(kpm_majorana_ldos LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()

enable_testing()

# Host-only. Eigen backs the R1 ED oracle (exact_eigenpairs) and nothing
# else; no CUDA target is planned for 1D (DESIGN.md section 5).
find_package(Eigen3 REQUIRED NO_MODULE)

add_library(majorana INTERFACE)
target_include_directories(majorana INTERFACE ${CMAKE_CURRENT_SOURCE_DIR}/include)
target_link_libraries(majorana INTERFACE Eigen3::Eigen)

add_executable(test_stub tests/test_stub.cpp)
target_link_libraries(test_stub PRIVATE majorana)
target_compile_options(test_stub PRIVATE -Wall -Wextra)
add_test(NAME stub COMMAND test_stub)

add_executable(test_r1 tests/test_r1.cpp)
target_link_libraries(test_r1 PRIVATE majorana)
target_compile_options(test_r1 PRIVATE -Wall -Wextra)
add_test(NAME r1 COMMAND test_r1)

add_executable(test_r2 tests/test_r2.cpp)
target_link_libraries(test_r2 PRIVATE majorana)
target_compile_options(test_r2 PRIVATE -Wall -Wextra)
add_test(NAME r2 COMMAND test_r2)

add_executable(test_r3 tests/test_r3.cpp)
target_link_libraries(test_r3 PRIVATE majorana)
target_compile_options(test_r3 PRIVATE -Wall -Wextra)
add_test(NAME r3 COMMAND test_r3)

add_executable(test_r4 tests/test_r4.cpp)
target_link_libraries(test_r4 PRIVATE majorana)
target_compile_options(test_r4 PRIVATE -Wall -Wextra)
add_test(NAME r4 COMMAND test_r4)
