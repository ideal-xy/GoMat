# --- 项目定义 ---
cmake_minimum_required(VERSION 3.15)
set(GOMAT_VERSION "0.0.1")
project(GoMat VERSION ${GOMAT_VERSION} LANGUAGES CXX)

# --- 关键设置：为 macOS 指定目标架构 ---
set(CMAKE_OSX_ARCHITECTURES "x86_64" CACHE STRING "为指定的 macOS 架构构建")
message(STATUS "正在为 macOS 架构构建: ${CMAKE_OSX_ARCHITECTURES}")

# --- C++ 标准 ---
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)
set(CMAKE_CXX_EXTENSIONS OFF)

# --- 库目标：gomat (静态库) ---
add_library(gomat STATIC
    src/linear_solver.cpp
    src/matrix_complex_features.cpp
    src/matrix_basic_features.cpp
    src/matrix_constructor.cpp
    src/matrix_decompositions.cpp
    src/matrix_io.cpp
    src/matrix_operations.cpp
    src/matrix_private_member.cpp
    src/matrix_transformation.cpp
    src/matrix_view.cpp
    src/utils.cpp
    src/vector.cpp
    src/matrix_simdproduct.cpp
)

# --- 头文件目录 ---
target_include_directories(gomat PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
    $<INSTALL_INTERFACE:include>
)

# --- 编译器选项 ---
if(CMAKE_CXX_COMPILER_ID MATCHES "Clang|AppleClang|GNU")
    target_compile_options(gomat PRIVATE -Wall -Wextra -Wpedantic)
endif()

# *** 关键修改：启用 AVX2 指令集支持 (Rosetta 2 兼容) ***
# 从 -mavx512f 降级到 -mavx2，因为 Rosetta 2 不支持 AVX-512。
target_compile_options(gomat PRIVATE -mavx2 -mfma)

# 为不同的构建类型（Debug/Release）分别设置选项
target_compile_options(gomat PRIVATE
    $<$<CONFIG:Debug>:-g>
    $<$<CONFIG:Release>:-O3>
)

# --- 子目录 ---
add_subdirectory(test)

