IF (NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  MESSAGE(STATUS "Setting build type to '${DEFAULT_BUILD_TYPE}' as none was specified.")
  SET(CMAKE_BUILD_TYPE "${DEFAULT_BUILD_TYPE}" CACHE STRING "Choose the type of build." FORCE)
  # Set the possible values of build type for cmake-gui
  SET_PROPERTY(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS
    "Debug" "Release" "Asan" "MinSizeRel" "RelWithDebInfo")
ENDIF ()

IF (DIPTEST_BUILDING_WHEELS)
    # disable architecture flags
    SET(DIPTEST_RELEASE_FLAGS "")
    MESSAGE(STATUS "diptest: Building for non-native host")
ELSE ()
    SET(DIPTEST_RELEASE_FLAGS "")
    INCLUDE(CheckCXXCompilerFlag)
    CHECK_CXX_COMPILER_FLAG(-march=native COMPILER_SUPPORTS_MARCH_NATIVE)
    IF(COMPILER_SUPPORTS_MARCH_NATIVE)
        SET(DIPTEST_RELEASE_FLAGS ${DIPTEST_RELEASE_FLAGS} -march=native)
    ENDIF()
    CHECK_CXX_COMPILER_FLAG(-mtune=native COMPILER_SUPPORTS_MTUNE_NATIVE)
    IF(COMPILER_SUPPORTS_MTUNE_NATIVE)
        SET(DIPTEST_RELEASE_FLAGS ${DIPTEST_RELEASE_FLAGS} -mtune=native)
    ENDIF()
    CHECK_CXX_COMPILER_FLAG(-ftree-vectorize COMPILER_SUPPORTS_FTREE)
    IF(COMPILER_SUPPORTS_FTREE)
        SET(DIPTEST_RELEASE_FLAGS ${DIPTEST_RELEASE_FLAGS} -ftree-vectorize)
    ENDIF()
    CHECK_CXX_COMPILER_FLAG(-mavx COMPILER_SUPPORTS_MAVX)
    IF(COMPILER_SUPPORTS_MAVX2)
        SET(DIPTEST_RELEASE_FLAGS ${DIPTEST_RELEASE_FLAGS} -mavx)
    ENDIF()
    CHECK_CXX_COMPILER_FLAG(-mavx2 COMPILER_SUPPORTS_MAVX2)
    IF(COMPILER_SUPPORTS_MAVX2)
        SET(DIPTEST_RELEASE_FLAGS ${DIPTEST_RELEASE_FLAGS} -mavx2)
    ENDIF()
ENDIF()


IF (DIPTEST_COVERAGE)
    # --coverage option is used to compile and link code instrumented for coverage analysis.
    # The option is a synonym for
    #    -fprofile-arcs
    #    -ftest-coverage (when compiling)
    #    -lgcov (when linking).
    SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} --coverage")
    SET(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} --coverage")
    SET(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} --coverage")
    SET(CMAKE_SHARED_LINKER_FLAGS "${CMAKE_SHARED_LINKER_FLAGS} --coverage")
ENDIF ()

IF (CMAKE_BUILD_TYPE EQUAL "Asan")
    SET(CMAKE_C_FLAGS_ASAN
        "${CMAKE_C_FLAGS_DEBUG} -fsanitize=address -fno-omit-frame-pointer" CACHE STRING
        "Flags used by the C compiler for Asan build type or configuration." FORCE)
    
    SET(CMAKE_CXX_FLAGS_ASAN
        "${CMAKE_CXX_FLAGS_DEBUG} -fsanitize=address -fno-omit-frame-pointer" CACHE STRING
        "Flags used by the C++ compiler for Asan build type or configuration." FORCE)
    
    SET(CMAKE_EXE_LINKER_FLAGS_ASAN
        "${CMAKE_SHARED_LINKER_FLAGS_DEBUG} -fsanitize=address" CACHE STRING
        "Linker flags to be used to create executables for Asan build type." FORCE)
    
    SET(CMAKE_SHARED_LINKER_FLAGS_ASAN
        "${CMAKE_SHARED_LINKER_FLAGS_DEBUG} -fsanitize=address" CACHE STRING
        "Linker lags to be used to create shared libraries for Asan build type." FORCE)
ENDIF ()
