cmake_minimum_required(VERSION 3.10)
project(FVW_Simulation)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

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

# Executable
add_executable(fvw_cpp ${SRC_SOURCES} ${VALIDATION_SOURCES})

# Link HDF5 libraries
target_link_libraries(fvw_cpp ${HDF5_CXX_LIBRARIES})

# Compiler flags
if(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang|MSVC")
    add_compile_options(-Wall -Wextra)
endif()