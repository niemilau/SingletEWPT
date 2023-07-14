# Get latest release of Google Test from github 
include(FetchContent)

FetchContent_Declare(
  googletest
  GIT_REPOSITORY https://github.com/google/googletest.git
  GIT_TAG        main
)

# Set options for fetching Google Test
set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)
set(gtest_disable_pthreads ON CACHE BOOL "" FORCE)


## Fetch Google Test library
message(STATUS "Fetching Google Test...")
FetchContent_MakeAvailable(googletest)


set(TEST_EXEC ${MAIN_EXEC}_Tests)

set(TEST_SOURCES test_main.cpp test_renormalization.cpp test_Veff3D.cpp)

add_executable(${TEST_EXEC} ${TEST_SOURCES})


# Get includes from the main lib
get_target_property(MAIN_INCLUDES ${MAIN_LIB} INCLUDE_DIRECTORIES)

## Includes for unit tests
target_include_directories(${TEST_EXEC} PRIVATE 
	${CMAKE_CURRENT_SOURCE_DIR}/include
	${MAIN_INCLUDES}
)

# Link against Google Test
target_link_libraries(${TEST_EXEC} PRIVATE 
	gtest 
	gtest_main
	${MAIN_LIB}
)

install(TARGETS ${TEST_EXEC} DESTINATION ${INSTALL_DIR})

