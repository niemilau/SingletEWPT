
## This will build everything except for main.cpp (the main function) into a library
## and link that to a 'main' executable (the actual program).
## The library can be used for a separate unit test target without having to compile the sources twice 

set(MAIN_LIB ${MAIN_EXEC}_lib CACHE INTERNAL "")

set(COMMON_SOURCES betafunct.cpp dimred.cpp PassarinoVeltman.cpp renormalization.cpp scanner.cpp selfEnergy.cpp)

add_library(${MAIN_LIB} STATIC ${COMMON_SOURCES})


## ------ Find/fetch external packages
message(STATUS "Loading external libraries")

include(FetchContent)
Set(FETCHCONTENT_QUIET FALSE) ## Show what the fetching is doing
set(FETCHCONTENT_UPDATES_DISCONNECTED ON)

# GSL
find_package(GSL REQUIRED)

# dlib: Fetch if not found on the system
find_package(dlib QUIET)

if (NOT dlib_FOUND)
	message(STATUS "dlib not found, attempting to fetch")
	FetchContent_Declare(dlib
		GIT_REPOSITORY https://github.com/davisking/dlib.git
		GIT_TAG        v19.24.2
		GIT_PROGRESS TRUE
	)
	FetchContent_MakeAvailable(dlib)
endif()

## ------- 

target_include_directories(${MAIN_LIB} PRIVATE 
	${CMAKE_CURRENT_SOURCE_DIR}/include/SingletEWPT
	dlib::dlib # dlib_INCLUDE_DIRS is apparently deprecated
	GSL::gsl
)

## Public linkage to external libs
target_link_libraries(${MAIN_LIB} PUBLIC dlib::dlib GSL::gsl)

## The 'main' program
add_executable(${MAIN_EXEC} main.cpp)


## main exec needs to include headers, just read them from the lib
get_target_property(MAIN_INCLUDES ${MAIN_LIB} INCLUDE_DIRECTORIES)

target_include_directories(${MAIN_EXEC} PRIVATE 
	${MAIN_INCLUDES}
)

target_link_libraries(${MAIN_EXEC} PRIVATE ${MAIN_LIB})


install(TARGETS ${MAIN_EXEC} DESTINATION ${INSTALL_DIR})


