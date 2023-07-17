# Get latest release of Google Test from github 
include(FetchContent)
Set(FETCHCONTENT_QUIET FALSE) ## Show what the fetching is doing
set(FETCHCONTENT_UPDATES_DISCONNECTED ON)

## Fetch Google Test library
include(FetchContent)
FetchContent_Declare(
  googletest
  GIT_REPOSITORY https://github.com/google/googletest.git
  GIT_TAG        v1.13.0
)

# For Windows: Prevent overriding the parent project's compiler/linker settings
set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)

#FetchContent_MakeAvailable(googletest)
## This should be equivalent to FetchContent_MakeAvailable, except that we exclude gtest from install
FetchContent_GetProperties(googletest)
if(NOT googletest_POPULATED)
  FetchContent_Populate(googletest)
  add_subdirectory(${googletest_SOURCE_DIR} ${googletest_BINARY_DIR} EXCLUDE_FROM_ALL)
endif()


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

