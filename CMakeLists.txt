cmake_minimum_required(VERSION 3.10)
project(FVW_Simulation)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Include directories
include_directories(${CMAKE_SOURCE_DIR}/include)
include_directories(${CMAKE_SOURCE_DIR}/validation/include)

# Source files
set(SRC_SOURCES
    src/main.cpp
    src/airfoil.cpp
    src/geometry.cpp
    src/position.cpp
    src/utils.cpp
    src/velocity.cpp
    src/bem.cpp
    src/wake.cpp
)

set(VALIDATION_SOURCES
    validation/src/validate.cpp
    validation/src/validate_geometry.cpp
    validation/src/validate_position.cpp
)

# Executable
add_executable(fvw_cpp ${SRC_SOURCES} ${VALIDATION_SOURCES})

# Compiler flags
if(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang|MSVC")
    add_compile_options(-Wall -Wextra)
endif()