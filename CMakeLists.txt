cmake_minimum_required(VERSION 3.15)
set(GOMAT_VERSION "0.0.1")
project(GoMat VERSION ${GOMAT_VERSION} LANGUAGES CXX)

#目标架构为x86,一开始在改机器上测试avx指令集，需要面向x86构建（结合Rosetta2）
# set(CMAKE_OSX_ARCHITECTURES "x86-64" CACHE STRING "为指定的 macOS 架构构建")
# message(STATUS "正在为 macOS 架构构建: ${CMAKE_OSX_ARCHITECTURES}")

#
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)
set(CMAKE_CXX_EXTENSIONS OFF)


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

# 头文件目录
target_include_directories(gomat PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
    $<INSTALL_INTERFACE:include>
)


if(CMAKE_CXX_COMPILER_ID MATCHES "Clang|AppleClang|GNU")
    target_compile_options(gomat PRIVATE -Wall -Wextra -Wpedantic)
endif()


if(CMAKE_SYSTEM_PROCESSOR MATCHES "x86_64|amd64")
    target_compile_options(gomat PRIVATE -mavx2 -mfma)
    message(STATUS "启用 x86_64 AVX2/FMA 指令集支持")
elseif(CMAKE_SYSTEM_PROCESSOR MATCHES "arm64|aarch64")
    message(STATUS "启用 ARM64 NEON 指令集支持")
    if(NOT APPLE)
        target_compile_options(gomat PRIVATE -march=armv8-a+simd)  
    endif()
else()
    message(STATUS "未知架构 ${CMAKE_SYSTEM_PROCESSOR}，不启用 SIMD 优化标志")
endif()



target_compile_options(gomat PRIVATE
    $<$<CONFIG:Debug>:-O3 >
    $<$<CONFIG:Release>:-O3>
)

# 在这里写测试代码
add_subdirectory(test)

