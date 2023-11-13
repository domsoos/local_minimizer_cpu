cmake_minimum_required(VERSION 3.10)
project(MyProject VERSION 1.0)

enable_testing()

# Specify the C++ standard
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED True)

# Include the subdirectories
add_subdirectory(src)
add_subdirectory(lib)
add_subdirectory(tests)
