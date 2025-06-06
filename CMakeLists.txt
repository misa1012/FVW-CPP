cmake_minimum_required(VERSION 3.10)
project(FVW_Simulation)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Option to enable gprof profiling
option(ENABLE_PROFILING "Enable gprof profiling" OFF) # Default, can be changed with -DENABLE_PROFILING=ON

# 查找 HDF5
find_package(HDF5 REQUIRED COMPONENTS CXX)
if(HDF5_FOUND)
    include_directories(${HDF5_INCLUDE_DIRS})
else()
    message(FATAL_ERROR "HDF5 not found. Please install libhdf5-dev.")
endif()

# Include directories
include_directories(${CMAKE_SOURCE_DIR}/include)
include_directories(${CMAKE_SOURCE_DIR}/validation/include)

# Source files
set(SRC_SOURCES
    main.cpp
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

# --- Executable: fvw_cpp ---
add_executable(fvw_cpp ${SRC_SOURCES} ${VALIDATION_SOURCES})

# Link HDF5 libraries for fvw_cpp
target_link_libraries(fvw_cpp PRIVATE ${HDF5_CXX_LIBRARIES})

# --- Executable: compute_induced_velocity ---
add_executable(compute_induced_velocity
    compute_induced_velocity.cpp
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

# Link HDF5 libraries for compute_induced_velocity
target_link_libraries(compute_induced_velocity PRIVATE ${HDF5_CXX_LIBRARIES})

find_package(OpenMP)
if(OpenMP_CXX_FOUND)
    message(STATUS "OpenMP CXX support found. Adding OpenMP flags.")
    # 为目标添加编译和链接选项
    target_compile_options(compute_induced_velocity PRIVATE ${OpenMP_CXX_FLAGS})
    target_link_libraries(compute_induced_velocity PRIVATE ${OpenMP_CXX_LIBRARIES}) # 关键是这个链接库

    # 如果 fvw_cpp 也用 OpenMP
    target_compile_options(fvw_cpp PRIVATE ${OpenMP_CXX_FLAGS})
    target_link_libraries(fvw_cpp PRIVATE ${OpenMP_CXX_LIBRARIES})
endif()

# --- Compiler and Linker Flags ---

# General warning flags (applied to both targets)
if(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
    target_compile_options(fvw_cpp PRIVATE -Wall -Wextra)
    target_compile_options(compute_induced_velocity PRIVATE -Wall -Wextra)
elseif(CMAKE_CXX_COMPILER_ID MATCHES "MSVC")
    target_compile_options(fvw_cpp PRIVATE /W4)
    target_compile_options(compute_induced_velocity PRIVATE /W4)
endif()

# Profiling flags (applied to both targets if ENABLE_PROFILING is ON)
if(ENABLE_PROFILING)
    if(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang") # -pg is for GCC/Clang
        message(STATUS "gprof profiling enabled. Adding -pg flags.")
        target_compile_options(fvw_cpp PRIVATE -pg)
        target_link_options(fvw_cpp PRIVATE -pg)
    else()
        message(WARNING "ENABLE_PROFILING is ON, but the current compiler (${CMAKE_CXX_COMPILER_ID}) might not support gprof (-pg).")
    endif()
endif()

