cmake_minimum_required(VERSION 3.10)
project(FVW_Simulation)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

option(ENABLE_PROFILING "Enable gprof profiling" OFF)

# Executable options
option(BUILD_FVW_CPP "Build fvw_cpp executable" OFF)
option(BUILD_POSTPROCESS_PROBES "Build postprocess_probes executable" OFF)
option(BUILD_POSTPROCESS_GRID "Build postprocess_grid executable" ON)
option(BUILD_VALIDATE_TSR "Build validate_tsr executable" OFF)
option(BUILD_VALIDATE_NSEGMENTS "Build validate_n_segment executable" OFF)
option(BUILD_COMPUTE_DOWNSTREAM_VELOCITY "Build executable to compute downstream velocity profile at specific x,y=0,z locations" OFF)

# ----------------------------------------
# External Dependencies
find_package(HDF5 REQUIRED COMPONENTS CXX)
find_package(OpenMP)

# ----------------------------------------
include_directories(${CMAKE_SOURCE_DIR}/include)
include_directories(${CMAKE_SOURCE_DIR}/include/validation)

# ----------------------------------------
# Core Library: fvw_core
# ----------------------------------------
set(CORE_SOURCES
  src/core/utils.cpp
  src/models/airfoil.cpp
  src/core/geometry.cpp
  src/core/position.cpp
  src/core/velocity.cpp
  src/models/performance.cpp
  src/models/bem.cpp
  src/core/wake.cpp
  src/io/postprocess.cpp
  src/io/config.cpp
  src/simulation/simulation_runner.cpp
)

set(VALIDATION_SOURCES
  src/validation/validate.cpp
  src/validation/validate_geometry.cpp
  src/validation/validate_position.cpp
)

# 这样所有 executables 都会链接这个静态库，避免重复编译
add_library(fvw_core STATIC ${CORE_SOURCES} ${VALIDATION_SOURCES})

# Target-specific include directories (PUBLIC propagates to consumers)
target_include_directories(fvw_core PUBLIC 
    ${CMAKE_SOURCE_DIR}/include
    ${CMAKE_SOURCE_DIR}/include/validation
    ${HDF5_INCLUDE_DIRS}
)

# Link dependencies to fvw_core
target_link_libraries(fvw_core PUBLIC ${HDF5_CXX_LIBRARIES})

if(OpenMP_CXX_FOUND)
    target_compile_options(fvw_core PUBLIC ${OpenMP_CXX_FLAGS})
    target_link_libraries(fvw_core PUBLIC ${OpenMP_CXX_LIBRARIES})
endif()

if(ENABLE_PROFILING AND CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
    target_compile_options(fvw_core PUBLIC -pg)
    target_link_options(fvw_core PUBLIC -pg)
endif()

# Compiler warnings
if(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
    target_compile_options(fvw_core PRIVATE -Wall -Wextra)
elseif(CMAKE_CXX_COMPILER_ID MATCHES "MSVC")
    target_compile_options(fvw_core PRIVATE /W4)
endif()


# ----------------------------------------
# Executables
# ----------------------------------------

if(BUILD_FVW_CPP)
  add_executable(fvw_cpp main.cpp)
  target_link_libraries(fvw_cpp PRIVATE fvw_core)
endif()

if(BUILD_POSTPROCESS_PROBES)
  add_executable(postprocess_probes tools/postprocess_probes.cpp)
  target_link_libraries(postprocess_probes PRIVATE fvw_core)
endif()

if(BUILD_POSTPROCESS_GRID)
  add_executable(postprocess_grid tools/postprocess_grid.cpp)
  target_link_libraries(postprocess_grid PRIVATE fvw_core)
endif()

if(BUILD_COMPUTE_DOWNSTREAM_VELOCITY)
  add_executable(compute_downstream_velocity tools/compute_downstream.cpp)
  target_link_libraries(compute_downstream_velocity PRIVATE fvw_core)
endif()

# ----------------------------------------
# Build Type Defaults
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
  message(STATUS "No build type specified, defaulting to 'Release'.")
endif()

set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -O3 -DNDEBUG")
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -O0 -g")
set(CMAKE_CXX_FLAGS_RELWITHDEBINFO "${CMAKE_CXX_FLAGS_RELWITHDEBINFO} -O2 -g -DNDEBUG")

message(STATUS "Build type is: ${CMAKE_BUILD_TYPE}")