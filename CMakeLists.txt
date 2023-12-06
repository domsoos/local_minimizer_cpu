cmake_minimum_required(VERSION 3.10)
project(Celsius VERSION 1.0)

enable_testing()

# Specify the C++ standard
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED True)

find_package(dlib REQUIRED)
# Include the subdirectories
add_subdirectory(src)
add_subdirectory(lib)
add_subdirectory(tests)

#include_directories(/opt/homebrew/opt/dlib/include)
find_package(dlib REQUIRED)
#include_directories(${dlib_INCLUDE_DIRS})

#message("Included directories: ${CMAKE_INCLUDE_DIRECTORIES_PROJECT_BEFORE}")