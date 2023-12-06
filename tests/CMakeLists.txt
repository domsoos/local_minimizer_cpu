# Enable CTest for running test cases
enable_testing()

# Set the build type if it is not already set
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Debug)
endif()

# Add flags for AddressSanitizer
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -fsanitize=address -fno-omit-frame-pointer -g")
set(CMAKE_LINKER_FLAGS_DEBUG "${CMAKE_LINKER_FLAGS_DEBUG} -fsanitize=address")

add_executable(utility_tests utility_test.cpp)
add_executable(optimization_tests optimization_test.cpp)
#add_executable(genetic_tests genetic_test.cpp ../src/test_functions.cpp)
add_executable(function_tests function_tests.cpp ../src/test_functions.cpp)


# Create test executables
target_link_libraries(utility_tests core_lib)
target_link_libraries(optimization_tests core_lib)
target_link_libraries(function_tests core_lib)
#add_executable(genetic_tests genetic_test.cpp ../src/test_functions.cpp)
#target_link_libraries(genetic_tests core_lib)

# Define the tests
add_test(NAME UtilityTests COMMAND utility_tests)
add_test(NAME OptimizationTests COMMAND optimization_tests)
add_test(NAME FunctionTests COMMAND function_tests)
#add_test(NAME GeneticTests COMMAND genetic_tests)

# Specify include directories for the test files to find the header files
target_include_directories(utility_tests PUBLIC ../lib ../src)
target_include_directories(optimization_tests PUBLIC ../lib ../src)
target_include_directories(function_tests PUBLIC ../lib ../src)
#target_include_directories(genetic_tests PUBLIC ../lib ../src)

