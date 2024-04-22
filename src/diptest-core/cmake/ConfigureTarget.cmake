# Set required C++ flags
set_property(TARGET _diptest_core PROPERTY CXX_STANDARD ${DIPTEST_CPP_STANDARD})
set_property(TARGET _diptest_core PROPERTY CXX_STANDARD_REQUIRED ON)
set_property(TARGET _diptest_core PROPERTY POSITION_INDEPENDENT_CODE ON)

# -- Optional
if(OpenMP_CXX_FOUND)
    target_link_libraries(_diptest_core PUBLIC OpenMP::OpenMP_CXX)
    target_compile_definitions(_diptest_core PRIVATE DIPTEST_HAS_OPENMP_SUPPORT=TRUE)
endif()

if(DIPTEST_ENABLE_DEBUG)
  message(STATUS "diptest: Building with debug support")
  target_compile_definitions(_diptest_core PRIVATE DIPTEST_DEBUG=TRUE)
endif()

if(DIPTEST_ENABLE_DEVMODE)
  target_compile_options(_diptest_core PRIVATE ${DIPTEST_DEVMODE_OPTIONS})
endif()

if(DIPTEST_ASAN_BUILD)
  target_compile_options(_diptest_core PRIVATE -fsanitize=address
                                               -fno-omit-frame-pointer)
  target_link_options(_diptest_core PRIVATE -fsanitize=address
                      -fno-omit-frame-pointer -shared-libasan)
endif()

include(CheckCXXCompilerFlag)

function(check_cxx_support FLAG DEST)
    string(SUBSTRING ${FLAG} 1 -1 STRIPPED_FLAG)
    string(REGEX REPLACE "=" "_" STRIPPED_FLAG ${STRIPPED_FLAG})
    string(TOUPPER ${STRIPPED_FLAG} STRIPPED_FLAG)
    set(RES_VAR "${STRIPPED_FLAG}_SUPPORTED")
    check_cxx_compiler_flag("${FLAG}" ${RES_VAR})
    if(${RES_VAR})
        set(${DEST} "${${DEST}} ${FLAG}" PARENT_SCOPE)
    endif()
endfunction()

# -- Compiler Flags
if (DIPTEST_ENABLE_ARCH_FLAGS AND "${CMAKE_CXX_FLAGS}" STREQUAL "${CMAKE_CXX_FLAGS_DEFAULT}")
    set(DIPTEST_ARCHITECTURE_FLAGS "")
    if (APPLE AND (${CMAKE_SYSTEM_PROCESSOR} STREQUAL "arm64"))
        # see https://github.com/google/highway/issues/745
        check_cxx_support("-march=native" DIPTEST_ARCHITECTURE_FLAGS)
    else()
        include(FindSse)
        include(FindAvx)
        DIPTEST_CHECK_FOR_SSE()
        add_definitions(${SSE_DEFINITIONS})
        DIPTEST_CHECK_FOR_AVX()
        string(APPEND DIPTEST_ARCHITECTURE_FLAGS "${SSE_FLAGS} ${AVX_FLAGS}")
    endif()

    if (NOT ${DIPTEST_ARCHITECTURE_FLAGS} STREQUAL "")
        string(STRIP ${DIPTEST_ARCHITECTURE_FLAGS} DIPTEST_ARCHITECTURE_FLAGS)
        message(STATUS "diptest | Enabled arch flags: ${DIPTEST_ARCHITECTURE_FLAGS}")
        if (MSVC)
            separate_arguments(DIPTEST_ARCHITECTURE_FLAGS WINDOWS_COMMAND "${DIPTEST_ARCHITECTURE_FLAGS}")
        else()
            separate_arguments(DIPTEST_ARCHITECTURE_FLAGS UNIX_COMMAND "${DIPTEST_ARCHITECTURE_FLAGS}")
        endif()
        target_compile_options(_diptest_core PRIVATE $<$<CONFIG:RELEASE>:${DIPTEST_ARCHITECTURE_FLAGS}>)
    else()
        message(STATUS "diptest | Architecture flags enabled but no valid flags were found")
    endif()
endif()
