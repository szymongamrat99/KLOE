# Tests CMakeLists.txt for KLOE project
# Author: Szymon Gamrat
# Description: Configuration file for unit tests of KLOE analysis framework

# Find required packages
find_package(Boost REQUIRED COMPONENTS unit_test_framework system filesystem)
find_package(ROOT REQUIRED COMPONENTS Core RIO Tree)
find_package(CURL REQUIRED)

# Create test executables
add_executable(KLOEClassTest KLOE_class_test.cpp)
add_executable(KaonMassReconstructorTest KaonMassReconstructor_test.cpp)

# Set include directories
target_include_directories(KLOEClassTest PRIVATE
    ${CMAKE_SOURCE_DIR}/Include/Codes/inc  # Main KLOE headers
    ${ROOT_INCLUDE_DIRS}
    ${Boost_INCLUDE_DIRS}
)

# Link required libraries
target_link_libraries(KLOEClassTest PRIVATE
    LibRec                    # Main KLOE library
    ${ROOT_LIBRARIES}
    ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY}
    ${Boost_SYSTEM_LIBRARY}   # Add Boost.System library
    ${Boost_FILESYSTEM_LIBRARY}  # Add Boost.Filesystem library
    CURL::libcurl            # Add libcurl
)

target_link_libraries(KaonMassReconstructorTest PRIVATE
    LibRec                    # Main KLOE library
    ${ROOT_LIBRARIES}
    ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY}
    ${Boost_SYSTEM_LIBRARY}
    ${Boost_FILESYSTEM_LIBRARY}
    CURL::libcurl
)

# Add include directories for KaonMassReconstructorTest
target_include_directories(KaonMassReconstructorTest PRIVATE
    ${CMAKE_SOURCE_DIR}/Include/Codes/inc
    ${ROOT_INCLUDE_DIRS}
    ${Boost_INCLUDE_DIRS}
)

# Enable testing functionality
enable_testing()

# Register tests with CTest
add_test(
    NAME KLOE_Class_Tests
    COMMAND KLOEClassTest
    WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/tests
)

add_test(
    NAME KaonMassReconstructor_Tests
    COMMAND KaonMassReconstructorTest
    WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/tests
)

# Set test properties
set_tests_properties(KLOE_Class_Tests PROPERTIES
    ENVIRONMENT "BOOST_TEST_LOG_LEVEL=all"
    ENVIRONMENT "ROOT_INCLUDE_PATH=${ROOT_INCLUDE_DIRS}"
)

set_tests_properties(KaonMassReconstructor_Tests PROPERTIES
    ENVIRONMENT "BOOST_TEST_LOG_LEVEL=all"
    ENVIRONMENT "ROOT_INCLUDE_PATH=${ROOT_INCLUDE_DIRS}"
)