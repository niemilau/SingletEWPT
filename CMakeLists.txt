cmake_minimum_required(VERSION 3.0)

set(PROJECT_NAME SingletEWPT)
project(${PROJECT_NAME} VERSION 0.5.0 LANGUAGES CXX)

# Set C++ standard
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# Set output directories
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)


# Set the path to the installation directories
set(INSTALL_DIR ${CMAKE_SOURCE_DIR}/bin)


# Configure build type if not given by the user

set(CMAKE_CONFIGURATION_TYPES Debug Release)

if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Debug)
    #set(CMAKE_BUILD_TYPE Release)
endif()


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
    if(${CMAKE_BUILD_TYPE} MATCHES Debug)
    	## This flag is needed also during linking so add_compile_options doesn't do the trick apparently
    	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -pg")
    	
		# no optimization here since it can interfere with profiling
    endif()
    
    if(CMAKE_BUILD_TYPE MATCHES "Release")
    	# Optimization
    	add_compile_options(-O3)
    endif()
    
elseif(CMAKE_CXX_COMPILER_ID STREQUAL "MSVC")
    # Enable warnings and treat them as errors
    #add_compile_options(/W4 /WX)
endif()

##---------------------- Unit testing
if(${CMAKE_BUILD_TYPE} MATCHES Debug)
	message("-- Configuring a Debug build")
	add_subdirectory(tests)
endif()

set(TARGET_EXEC ${PROJECT_NAME})

# Add all source files in the src directory
file(GLOB COMMON_SOURCES "src/*.cpp")

# Create an executable target
add_executable(${TARGET_EXEC} ${COMMON_SOURCES})

# Project specific includes (global scope again) 
include_directories(
    ${CMAKE_CURRENT_SOURCE_DIR}/src/include/SingletEWPT
)

## External packages
find_package(dlib REQUIRED)

target_link_libraries(${TARGET_EXEC} PRIVATE dlib::dlib)

#--------------------------


# 'make install'
install(TARGETS ${TARGET_EXEC} DESTINATION ${INSTALL_DIR})


# Enable testing
#enable_testing()
#add_subdirectory(tests)

