# Add the source files for the main executable
add_executable(main_app main.cpp genetic.cpp optimization.cpp test_functions.cpp)

# Link the main executable with the core_lib
target_link_libraries(main_app core_lib)

# Specify where to find the header files for the core library and test functions
target_include_directories(main_app PUBLIC ${CMAKE_CURRENT_SOURCE_DIR} ../lib)
