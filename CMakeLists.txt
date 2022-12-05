CMAKE_MINIMUM_REQUIRED(VERSION 3.16)

IF (NOT DEFINED DIPTEST_VERSION_INFO)
    SET(DIPTEST_VERSION_INFO "0.5.2")
ENDIF ()

PROJECT(
    diptest
    VERSION ${DIPTEST_VERSION_INFO}
    DESCRIPTION "Python bindings for Hartigan's diptest."
    LANGUAGES CXX
)

################################################################################
#                                   SETTINGS                                   #
################################################################################
OPTION(DIPTEST_DEV_MODE OFF)
OPTION(DIPTEST_ENABLE_DEBUG OFF)
OPTION(DIPTEST_ENABLE_OPENMP OFF)
OPTION(DIPTEST_ENABLE_EXT_TESTS OFF)
OPTION(DIPTEST_ENABLE_ARCH_FLAGS OFF)
OPTION(DIPTEST_ENABLE_ARCH_FLAGS_SIMPLE OFF)

IF (DIPTEST_ENABLE_EXT_TESTS OR DIPTEST_DEV_MODE)
    ADD_DEFINITIONS(-DDIPTEST_BUILD_CPP_TESTS=TRUE)
ENDIF ()

SET(CMAKE_CXX_STANDARD 11)
SET(CMAKE_POSITION_INDEPENDENT_CODE ON)

FIND_PACKAGE(pybind11 CONFIG REQUIRED)
IF (DIPTEST_ENABLE_OPENMP)
    FIND_PACKAGE(OpenMP REQUIRED)
ENDIF ()

################################################################################
#                                     FLAGS                                    #
################################################################################
INCLUDE(GNUInstallDirs)

SET(CMAKE_MODULE_PATH_SAVED ${CMAKE_MODULE_PATH})
LIST(INSERT CMAKE_MODULE_PATH 0 "${PROJECT_SOURCE_DIR}/cmake")

SET(DEFAULT_BUILD_TYPE "Release")
# Set build type to Release if not specified
INCLUDE(BuildType)

IF (CMAKE_CXX_COMPILER_ID MATCHES "MSVC")
    # error C2039: '_snprintf': is not a member of 'std'
    ADD_DEFINITIONS(-DHAVE_SNPRINTF)
ENDIF ()

################################################################################
#                               PYTHON EXTENSION                               #
################################################################################
SET(MODNAME "_diptest")
ADD_DEFINITIONS(-DEXTENSION_MODULE_NAME=${MODNAME})

pybind11_add_module(${MODNAME}
    MODULE
        ${PROJECT_SOURCE_DIR}/src/bindings.cpp
        ${PROJECT_SOURCE_DIR}/src/wrapper.cpp
        ${PROJECT_SOURCE_DIR}/src/test_pcg.cpp
)

IF (DIPTEST_ENABLE_DEBUG)
    MESSAGE(STATUS "diptest: Building with debug support")
    TARGET_COMPILE_DEFINITIONS(${MODNAME} PRIVATE DIPTEST_DEBUG=TRUE)
ENDIF ()

IF (OpenMP_CXX_FOUND)
    MESSAGE(STATUS "diptest: Building with OpenMP support")
    TARGET_COMPILE_DEFINITIONS(
        ${MODNAME}
        PUBLIC
            DIPTEST_HAS_OPENMP_SUPPORT=TRUE
    )
    TARGET_LINK_LIBRARIES(
        ${MODNAME}
        PRIVATE
            OpenMP::OpenMP_CXX
    )
ENDIF ()

TARGET_INCLUDE_DIRECTORIES(${MODNAME} PUBLIC "${PROJECT_SOURCE_DIR}/include")
TARGET_INCLUDE_DIRECTORIES(${MODNAME} PRIVATE "${PROJECT_SOURCE_DIR}/external/pcg-cpp/include")
TARGET_COMPILE_DEFINITIONS(${MODNAME} PRIVATE VERSION_INFO=${DIPTEST_VERSION_INFO})

IF(DIPTEST_DEV_MODE)
    MESSAGE(STATUS "diptest: Building in development mode")
    TARGET_COMPILE_OPTIONS(${MODNAME}
        PRIVATE
            -Wall
            -Wextra
            -Wunused-variable
            -Wunused-const-variable
    )
ENDIF()

TARGET_COMPILE_OPTIONS( ${MODNAME}
    PRIVATE
        "$<$<CONFIG:RELEASE>:${DIPTEST_ARCHITECTURE_FLAGS}>"
)

# Set required C++11 flags
SET_PROPERTY(TARGET ${MODNAME} PROPERTY CXX_STANDARD_REQUIRED ON)
SET_PROPERTY(TARGET ${MODNAME} PROPERTY CXX_STANDARD 11)
SET_PROPERTY(TARGET ${MODNAME} PROPERTY POSITION_INDEPENDENT_CODE ON)

if(SKBUILD)
    INSTALL(TARGETS ${MODNAME} LIBRARY DESTINATION "${PROJECT_NAME}/lib")
ELSE ()
    INSTALL(
        TARGETS ${MODNAME}
        LIBRARY DESTINATION "${PROJECT_SOURCE_DIR}/diptest/lib"
    )
ENDIF ()
UNSET(DIPTEST_DEV_MODE CACHE)
UNSET(DIPTEST_ENABLE_DEBUG CACHE)
UNSET(DIPTEST_ENABLE_OPENMP CACHE)
UNSET(DIPTEST_ENABLE_EXT_TESTS CACHE)
UNSET(DIPTEST_ENABLE_ARCH_FLAGS CACHE)
UNSET(DIPTEST_ENABLE_ARCH_FLAGS_SIMPLE CACHE)
