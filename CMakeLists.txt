cmake_minimum_required(VERSION 3.13)
project(KLSPM00 CXX Fortran)

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED TRUE)
set(CMAKE_CXX_EXTENSIONS OFF)

# Set-up required packages
find_package(ROOT REQUIRED COMPONENTS Tree Physics MathCore Matrix TMVA Hist Graf)
find_package(Boost REQUIRED COMPONENTS system filesystem unit_test_framework)
find_package(CURL REQUIRED)

# Enable tests for the files
enable_testing()

add_subdirectory(tests)


# Inclusion of directories to the project 
add_subdirectory(Subanalysis/CovarianceMatrix)
add_subdirectory(Subanalysis/CPFit)
# add_subdirectory(Subanalysis/EfficiencyAnalysis)
add_subdirectory(Subanalysis/GeneratedVars)
# # add_subdirectory(Subanalysis/InterfFunction)
# # add_subdirectory(Subanalysis/IPRec)
add_subdirectory(Subanalysis/KchRec)
add_subdirectory(Subanalysis/OmegaRec)
add_subdirectory(Subanalysis/Neutrec)
add_subdirectory(../Include/Codes ${CMAKE_BINARY_DIR}/LibRecBuild)
add_subdirectory(../Include/FortranAnalysis ${CMAKE_BINARY_DIR}/FortBuild)
# add_subdirectory(Subanalysis/Pipi)
# add_subdirectory(Subanalysis/Plots)
# add_subdirectory(Subanalysis/RegenerationAnalysis)
# add_subdirectory(Subanalysis/Semileptonic)


# What should be the result of compilation
add_executable(KLSPM00 app/CPVAnalysis.cpp)

# Inclusion of libraries and header directories
target_include_directories(KLSPM00 PRIVATE
  ${Boost_INCLUDE_DIRS}
  ${CURL_INCLUDE_DIRS}
  ${CMAKE_CURRENT_SOURCE_DIR}/inc
  ${CMAKE_CURRENT_SOURCE_DIR}/../Include/Codes/inc
  ${CMAKE_CURRENT_SOURCE_DIR}/../Include/FortranAnalysis/inc
  )
target_link_libraries(KLSPM00 PRIVATE
  CovMatrix
  NeutRec
  KchRec
  GenVars
  CPFit
  OmegaRec
  ${Boost_LIBRARIES}
  ${CURL_LIBRARIES}
  ${ROOT_LIBRARIES}
  LibRec
  LibFort)

# Optional output directories
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)