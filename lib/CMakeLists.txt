# lib/CMakeLists.txt

# Create a static library from the utility code and other core functionality
add_library(core_lib STATIC 
    utility.cpp
    ../src/test_functions.cpp
    ../src/optimization.cpp
    ../src/genetic.cpp
    # other source files...
)

# Specify where the library headers are
target_include_directories(core_lib PUBLIC ${CMAKE_CURRENT_SOURCE_DIR} ../include)
