cmake_minimum_required(VERSION 3.24)
project(OpenMSFLASHTnT VERSION 1.0.0 LANGUAGES CXX)
option(OPENMS4_WARNINGS_AS_ERRORS "Treat package compiler warnings as errors" OFF)
set(CMAKE_COMPILE_WARNING_AS_ERROR ${OPENMS4_WARNINGS_AS_ERRORS})
include(GNUInstallDirs)
include(CTest)
include(cmake/OpenMS4Dependencies.cmake)
openms4_find_core()
openms4_find_package(OpenMSCLI)
openms4_find_package(OpenMSFLASH)
include(cmake/OpenMS4Tools.cmake)
openms4_load_tools(package_tools)
add_library(flashtnt_backend STATIC
  source/ANALYSIS/TOPDOWN/FLASHTnTHelpers.cpp
  source/ANALYSIS/TOPDOWN/FLASHTaggerAlgorithm.cpp
  source/ANALYSIS/TOPDOWN/FLASHExtenderAlgorithm.cpp
  source/ANALYSIS/TOPDOWN/FLASHTnTAlgorithm.cpp
  source/FORMAT/FLASHTnTFile.cpp)
target_compile_features(flashtnt_backend PUBLIC cxx_std_23)
target_include_directories(flashtnt_backend PUBLIC include)
target_link_libraries(flashtnt_backend PUBLIC OpenMS::Core OpenMS::FLASH)
if(OpenMS_WITH_OPENMP)
  target_link_libraries(flashtnt_backend PRIVATE OpenMP::OpenMP_CXX)
endif()
add_executable(FLASHTnT src/FLASHTnT.cpp)
target_link_libraries(FLASHTnT PRIVATE flashtnt_backend OpenMS::CLI)
set_target_properties(FLASHTnT PROPERTIES RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin/$<0:>")
openms4_install_tools(flashtnt bin ${package_tools})
if(BUILD_TESTING)
  add_executable(FLASHTnT_port_contract tests/port_contract.cpp)
  target_link_libraries(FLASHTnT_port_contract PRIVATE flashtnt_backend)
  add_test(NAME FLASHTnT_port_contract COMMAND FLASHTnT_port_contract)
  set_tests_properties(FLASHTnT_port_contract PROPERTIES ENVIRONMENT "OMP_NUM_THREADS=2")
  find_package(Python3 REQUIRED COMPONENTS Interpreter)
  add_test(NAME FLASHTnT_AQPZ_identity COMMAND ${Python3_EXECUTABLE}
    ${CMAKE_CURRENT_SOURCE_DIR}/tests/test_aqpz.py $<TARGET_FILE:FLASHTnT>
    ${CMAKE_CURRENT_BINARY_DIR}/aqpz-results)
  set_tests_properties(FLASHTnT_AQPZ_identity PROPERTIES TIMEOUT 200 ENVIRONMENT
    "OMP_NUM_THREADS=2;OPENMS_TOOL_PREFIX_PATH=${CMAKE_BINARY_DIR};OPENMS_DISABLE_UPDATE_CHECK=ON")
endif()
