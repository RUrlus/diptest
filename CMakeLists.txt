cmake_minimum_required(VERSION 3.16...3.25)
set(STRIPPED_VERSION "")
string(REGEX MATCH "[0-9]+\.[0-9]+\.[0-9]+" STRIPPED_VERSION
             ${SKBUILD_PROJECT_VERSION})
project(
  ${SKBUILD_PROJECT_NAME}
  VERSION ${STRIPPED_VERSION}
  LANGUAGES CXX)
set(CMAKE_MODULE_PATH_SAVED ${CMAKE_MODULE_PATH})
list(INSERT CMAKE_MODULE_PATH 0 "${PROJECT_SOURCE_DIR}/src/diptest-core/cmake")

include(GNUInstallDirs)

# Set build type to Release if not specified
include(BuildType)

set(CMAKE_CXX_STANDARD ${DIPTEST_CPP_STANDARD})
set(DIPTEST_INCLUDE_DIR "${PROJECT_SOURCE_DIR}/src/diptest-core/include")
set(DIPTEST_PCG_INCLUDE_DIR
    "${PROJECT_SOURCE_DIR}/src/external/pcg-cpp/include")
set(DIPTEST_SRC_DIR "${PROJECT_SOURCE_DIR}/src/diptest-core/src/")
set(DIPTEST_SRC_FILES bindings.cpp bootstrap.cpp dipstat.cpp)

if(DIPTEST_ENABLE_EXT_TESTS OR DIPTEST_DEV_MODE)
  add_definitions(-DDIPTEST_BUILD_CPP_TESTS=TRUE)
  list(APPEND DIPTEST_SRC_FILES test_pcg.cpp)
endif()

list(TRANSFORM DIPTEST_SRC_FILES PREPEND ${DIPTEST_SRC_DIR})

if(CMAKE_CXX_COMPILER_ID MATCHES "MSVC")
  # error C2039: '_snprintf': is not a member of 'std'
  add_definitions(-DHAVE_SNPRINTF)
endif()

include(FindDependencies)

pybind11_add_module(_diptest_core MODULE ${DIPTEST_SRC_FILES})
target_include_directories(_diptest_core PUBLIC ${DIPTEST_INCLUDE_DIR})
target_include_directories(_diptest_core PRIVATE ${DIPTEST_PCG_INCLUDE_DIR})
target_link_libraries(_diptest_core PRIVATE pybind11::pybind11)
target_compile_definitions(_diptest_core
                           PRIVATE DIPTEST_VERSION_INFO=${PROJECT_VERSION})

include(ConfigureTarget)

if(DIPTEST_MBUILD)
  set(CMAKE_INSTALL_PREFIX "${PROJECT_SOURCE_DIR}/src")
endif()

install(TARGETS _diptest_core LIBRARY DESTINATION "${PROJECT_NAME}/lib")

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH_SAVED})
