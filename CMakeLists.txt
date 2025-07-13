cmake_minimum_required(VERSION 3.10)
project(FVW_Simulation)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

option(ENABLE_PROFILING "Enable gprof profiling" OFF)

# Executable options
option(BUILD_FVW_CPP "Build fvw_cpp executable" ON)
option(BUILD_COMPUTE_VELOCITY "Build compute_induced_velocity executable" OFF)
option(BUILD_POSTPROCESS_PROBES "Build postprocess_probes executable" OFF)
option(BUILD_POSTPROCESS_GRID "Build postprocess_grid executable" ON)

# ----------------------------------------
find_package(HDF5 REQUIRED COMPONENTS CXX)
if(HDF5_FOUND)
  include_directories(${HDF5_INCLUDE_DIRS})
else()
  message(FATAL_ERROR "HDF5 not found. Please install libhdf5-dev.")
endif()

include_directories(${CMAKE_SOURCE_DIR}/include)
include_directories(${CMAKE_SOURCE_DIR}/validation/include)

set(SHARED_SOURCES
  src/utils.cpp
  src/airfoil.cpp
  src/geometry.cpp
  src/position.cpp
  src/velocity.cpp
  src/performance.cpp
  src/bem.cpp
  src/wake.cpp
  src/postprocess.cpp
)

set(VALIDATION_SOURCES
  validation/src/validate.cpp
  validation/src/validate_geometry.cpp
  validation/src/validate_position.cpp
)

find_package(OpenMP)

if(BUILD_FVW_CPP)
  add_executable(fvw_cpp main.cpp ${SHARED_SOURCES} ${VALIDATION_SOURCES})
  target_link_libraries(fvw_cpp PRIVATE ${HDF5_CXX_LIBRARIES})
  if(OpenMP_CXX_FOUND)
    target_compile_options(fvw_cpp PRIVATE ${OpenMP_CXX_FLAGS})
    target_link_libraries(fvw_cpp PRIVATE ${OpenMP_CXX_LIBRARIES})
  endif()
  if(ENABLE_PROFILING AND CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
    target_compile_options(fvw_cpp PRIVATE -pg)
    target_link_options(fvw_cpp PRIVATE -pg)
  endif()
  if(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
    target_compile_options(fvw_cpp PRIVATE -Wall -Wextra)
  elseif(CMAKE_CXX_COMPILER_ID MATCHES "MSVC")
    target_compile_options(fvw_cpp PRIVATE /W4)
  endif()
endif()

if(BUILD_COMPUTE_VELOCITY)
  add_executable(compute_induced_velocity compute_induced_velocity.cpp ${SHARED_SOURCES})
  target_link_libraries(compute_induced_velocity PRIVATE ${HDF5_CXX_LIBRARIES})
  if(OpenMP_CXX_FOUND)
    target_compile_options(compute_induced_velocity PRIVATE ${OpenMP_CXX_FLAGS})
    target_link_libraries(compute_induced_velocity PRIVATE ${OpenMP_CXX_LIBRARIES})
  endif()
  if(ENABLE_PROFILING AND CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
    target_compile_options(compute_induced_velocity PRIVATE -pg)
    target_link_options(compute_induced_velocity PRIVATE -pg)
  endif()
  if(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
    target_compile_options(compute_induced_velocity PRIVATE -Wall -Wextra)
  elseif(CMAKE_CXX_COMPILER_ID MATCHES "MSVC")
    target_compile_options(compute_induced_velocity PRIVATE /W4)
  endif()
endif()

if(BUILD_POSTPROCESS_PROBES)
  add_executable(postprocess_probes postprocess_probes.cpp ${SHARED_SOURCES})
  target_link_libraries(postprocess_probes PRIVATE ${HDF5_CXX_LIBRARIES})
  if(OpenMP_CXX_FOUND)
    target_compile_options(postprocess_probes PRIVATE ${OpenMP_CXX_FLAGS})
    target_link_libraries(postprocess_probes PRIVATE ${OpenMP_CXX_LIBRARIES})
  endif()
endif()

if(BUILD_POSTPROCESS_GRID)
  add_executable(postprocess_grid postprocess_grid.cpp ${SHARED_SOURCES})
  target_link_libraries(postprocess_grid PRIVATE ${HDF5_CXX_LIBRARIES})
  if(OpenMP_CXX_FOUND)
    target_compile_options(postprocess_grid PRIVATE ${OpenMP_CXX_FLAGS})
    target_link_libraries(postprocess_grid PRIVATE ${OpenMP_CXX_LIBRARIES})
  endif()
endif()

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
  message(STATUS "No build type specified, defaulting to 'Release'.")
endif()

set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -O3 -DNDEBUG")
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -O0 -g")
set(CMAKE_CXX_FLAGS_RELWITHDEBINFO "${CMAKE_CXX_FLAGS_RELWITHDEBINFO} -O2 -g -DNDEBUG")

message(STATUS "Build type is: ${CMAKE_BUILD_TYPE}")