
## This will build everything except for main.cpp (the main function) into a library
## and link that to a 'main' executable (the actual program).
## The library can be used for a separate unit test target without having to compile the sources twice 

set(MAIN_LIB ${MAIN_EXEC}_lib CACHE INTERNAL "")

set(COMMON_SOURCES betafunct.cpp dimred.cpp PassarinoVeltman.cpp renormalization.cpp scanner.cpp selfEnergy.cpp)

add_library(${MAIN_LIB} STATIC ${COMMON_SOURCES})

# Third-party packages
find_package(Boost REQUIRED)
find_package(dlib REQUIRED)


target_include_directories(${MAIN_LIB} PRIVATE 
	${CMAKE_CURRENT_SOURCE_DIR}/include/SingletEWPT
	dlib::dlib # dlib_INCLUDE_DIRS is apparently deprecated
    ${Boost_INCLUDE_DIRS}
)

## Public linkage to external libs
target_link_libraries(${MAIN_LIB} PUBLIC dlib::dlib ${Boost_LIBRARIES})

## The 'main' program
add_executable(${MAIN_EXEC} main.cpp)


## main exec needs to include headers, just read them from the lib
get_target_property(MAIN_INCLUDES ${MAIN_LIB} INCLUDE_DIRECTORIES)

target_include_directories(${MAIN_EXEC} PRIVATE 
	${MAIN_INCLUDES}
)

target_link_libraries(${MAIN_EXEC} PRIVATE ${MAIN_LIB})


install(TARGETS ${MAIN_EXEC} DESTINATION ${INSTALL_DIR})


