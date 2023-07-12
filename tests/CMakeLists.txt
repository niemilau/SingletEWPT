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

# Fetch and configure Google Test
message("Fetching Google Test")
FetchContent_MakeAvailable(googletest)

set(TEST_EXEC UnitTests)

add_executable(${TEST_EXEC} test_main.cpp test_Veff3D.cpp)

# Link against Google Test
target_link_libraries(${TEST_EXEC} gtest gtest_main)
# and the main program
target_link_libraries(${TEST_EXEC} ${EXEC_NAME})
