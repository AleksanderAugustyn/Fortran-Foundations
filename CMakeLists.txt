cmake_minimum_required(VERSION 3.20)

include(FetchContent)

project(FortranFoundations
        VERSION 2.2.2
        LANGUAGES Fortran
        DESCRIPTION "Foundational Fortran modules: precision, constants, math utilities"
)

# --- Fetch GCC compiler options ---
FetchContent_Declare(
        GCC-Compiler-Options
        GIT_REPOSITORY https://github.com/AleksanderAugustyn/GCC-Compiler-Options.git
        GIT_TAG 1.4.2
)
FetchContent_MakeAvailable(GCC-Compiler-Options)

if (CMAKE_VERSION VERSION_GREATER_EQUAL "3.20")
    set(CMAKE_Fortran_PREPROCESS ON)
endif ()

if (NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
    set(CMAKE_BUILD_TYPE "Release" CACHE STRING "Build type" FORCE)
    set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS "Debug" "RelWithDebInfo" "Release")
endif ()

# --- Compiler flags ---
include(GCCCompilerOptions/FortranCompilerOptions)

create_fortran_library_interface(
        TARGET fortran_foundations_flags
)

# --- Library ---
add_library(fortran_foundations STATIC
        src/precision_utilities_mod.f08
        src/c_bindings_mod.f08
        src/mathematical_and_physical_constants_mod.f08
        src/mathematical_utilities_mod.f08
        src/utility_procedures_mod.f08
)

target_link_libraries(fortran_foundations PRIVATE fortran_foundations_flags)

target_compile_definitions(fortran_foundations PRIVATE
        $<$<CONFIG:Debug>:DEBUG>
        $<$<CONFIG:RelWithDebInfo>:RELWITHDEBINFO>
        $<$<CONFIG:Release>:RELEASE>
)

set_target_properties(fortran_foundations PROPERTIES
        Fortran_MODULE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/mod_files
)

target_include_directories(fortran_foundations
        PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}/mod_files>
)

add_library(FortranFoundations::fortran_foundations ALIAS fortran_foundations)
