add_library(gomat_core STATIC
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

# 指定头文件目录
# CMAKE_CURRENT_SOURCE_DIR 在这里指的是 core/ 目录
target_include_directories(gomat_core PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
    $<INSTALL_INTERFACE:include>
)



if(CMAKE_CXX_COMPILER_ID MATCHES "Clang|AppleClang|GNU")
    target_compile_options(gomat_core PUBLIC -Wall -Wextra -Wpedantic)
endif()

if(CMAKE_SYSTEM_PROCESSOR MATCHES "x86_64|amd64")
    target_compile_options(gomat_core PRIVATE -mavx2 -mfma)
    message(STATUS "为 gomat_core 启用 x86_64 AVX2/FMA 指令集")
elseif(CMAKE_SYSTEM_PROCESSOR MATCHES "arm64|aarch64")
    message(STATUS "为 gomat_core 启用 ARM64 NEON 指令集")
    if(NOT APPLE)
        target_compile_options(gomat_core PRIVATE -march=armv8-a+simd)  
    endif()
else()
    message(STATUS "未知架构 ${CMAKE_SYSTEM_PROCESSOR},gomat_core 不启用 SIMD 优化")
endif()

target_compile_options(gomat_core PUBLIC
    $<$<CONFIG:Debug>:-Og >
    $<$<CONFIG:Release>:-O3>
)