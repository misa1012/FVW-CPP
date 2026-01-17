cmake_minimum_required(VERSION 3.10)
project(FVW_Simulation)

# --- 1. Compiler Settings ---
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Profiling Option
option(ENABLE_PROFILING "Enable gprof profiling" OFF)

# Build Type Defaults
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
  message(STATUS "No build type specified, defaulting to 'Release'.")
endif()

message(STATUS "Build type: ${CMAKE_BUILD_TYPE}")

# Compiler Flags
set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -O3 -DNDEBUG")
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -O0 -g")
set(CMAKE_CXX_FLAGS_RELWITHDEBINFO "${CMAKE_CXX_FLAGS_RELWITHDEBINFO} -O2 -g -DNDEBUG")

# --- 2. Dependencies ---
find_package(HDF5 REQUIRED COMPONENTS CXX)
find_package(OpenMP)

# --- 3. Core Library (fvw_core) ---
# Groups: Core, Models, IO, Simulation, Validation
set(CORE_SOURCES
    # Core
    src/core/utils.cpp
    src/core/geometry.cpp
    src/core/position.cpp
    src/core/velocity.cpp
    src/core/wake.cpp
    # Models
    src/models/airfoil.cpp
    src/models/performance.cpp
    src/models/bem.cpp
    # IO
    src/io/config.cpp
    src/io/postprocess.cpp
    src/io/logger.cpp
    # Simulation
    src/simulation/simulation_runner.cpp
)

add_library(fvw_core STATIC ${CORE_SOURCES})

# Optional Validation Library
option(BUILD_VALIDATION "Build validation library" OFF)
if(BUILD_VALIDATION)
    add_library(fvw_validation STATIC
        src/validation/validate.cpp
        src/validation/validate_geometry.cpp
        src/validation/validate_position.cpp
    )
    target_link_libraries(fvw_validation PUBLIC fvw_core)
endif()

# Include Directories
# PUBLIC means any target linking fvw_core will inherently include these paths
target_include_directories(fvw_core PUBLIC 
    ${CMAKE_SOURCE_DIR}/include
    ${CMAKE_SOURCE_DIR}/include/validation
    ${HDF5_INCLUDE_DIRS}
)

# Link Libraries
target_link_libraries(fvw_core PUBLIC ${HDF5_CXX_LIBRARIES})

if(OpenMP_CXX_FOUND)
    target_compile_options(fvw_core PUBLIC ${OpenMP_CXX_FLAGS})
    target_link_libraries(fvw_core PUBLIC ${OpenMP_CXX_LIBRARIES})
endif()

if(ENABLE_PROFILING AND CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
    target_compile_options(fvw_core PUBLIC -pg)
    target_link_options(fvw_core PUBLIC -pg)
endif()

# Compiler Warnings
if(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
    target_compile_options(fvw_core PRIVATE -Wall -Wextra)
elseif(CMAKE_CXX_COMPILER_ID MATCHES "MSVC")
    target_compile_options(fvw_core PRIVATE /W4)
endif()

# --- 4. Executables ---

# Main Simulation
option(BUILD_FVW_CPP "Build main simulation executable" ON)
if(BUILD_FVW_CPP)
    add_executable(fvw_cpp main.cpp)
    target_link_libraries(fvw_cpp PRIVATE fvw_core)
endif()

# Post-processing Tools
option(BUILD_TOOL_GRID "Build VTK grid generator" ON)
option(BUILD_TOOL_PROBES "Build probe sampler" OFF)
option(BUILD_TOOL_DOWNSTREAM "Build downstream velocity calculator" OFF)

if(BUILD_TOOL_GRID)
    add_executable(postprocess_grid tools/postprocess_grid.cpp)
    target_link_libraries(postprocess_grid PRIVATE fvw_core)
endif()

if(BUILD_TOOL_PROBES)
    add_executable(postprocess_probes tools/postprocess_probes.cpp)
    target_link_libraries(postprocess_probes PRIVATE fvw_core)
endif()

if(BUILD_TOOL_DOWNSTREAM)
    add_executable(compute_downstream_velocity tools/compute_downstream.cpp)
    target_link_libraries(compute_downstream_velocity PRIVATE fvw_core)
endif()