cmake_minimum_required(VERSION 3.15)
project(SingletEWPT VERSION 1.0.0 LANGUAGES CXX)


set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)


# Set output directories
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

## Installation direction
set(INSTALL_DIR ${CMAKE_SOURCE_DIR}/bin)


# Set default built type
if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
	set(CMAKE_BUILD_TYPE Release)
endif()


message(STATUS "Build type: ${CMAKE_BUILD_TYPE}")


set(CMAKE_CXX_FLAGS_DEBUG_INIT "-Wall")

# Compiler-specific settings
if(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
    # Enable warnings
    add_compile_options(-Wall -Wextra)
    # Treat all warnings as errors?
    #add_compile_options(-Werror)
    # Enable more warnings (optional)
    #add_compile_options(-Wpedantic)
    
    ## Add flags for gprof profiling. Debug builds only
    if(CMAKE_BUILD_TYPE STREQUAL "Debug")
    	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -pg")
		add_compile_definitions(DEBUG_MODE)
    endif()
    
    if(CMAKE_BUILD_TYPE STREQUAL "Release")
    	add_compile_options(-O3)
    endif()
    
elseif(CMAKE_CXX_COMPILER_ID STREQUAL "MSVC")
	##
endif()


## Main program target
set(MAIN_EXEC SingletEWPT)

# Main code
add_subdirectory(src)

option(BUILD_TESTS "Build unit tests (fetches GTest)" OFF)

if(CMAKE_BUILD_TYPE STREQUAL "Debug")
	message(STATUS "==== NOTE: Configuring Debug build")
	message(STATUS "")
endif()

if(BUILD_TESTS)
	add_subdirectory(tests)
endif()

