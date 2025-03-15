project(KLSPM00)

set(CMAKE_CXX_STANDARD 14)

# Inclusion of directories to the project 
add_subdirectory(inc)
add_subdirectory(Subanalysis)

# What should be the result of compilation
add_executable(KLSPM00 CPVAnalysis.cpp)

# Inclusion of libraries and header directories
target_include_directories(KLSPM00 PRIVATE inc)
target_link_libraries(KLSPM00 PRIVATE librec.so)