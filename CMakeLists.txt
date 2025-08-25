# 最低 CMake 版本要求
cmake_minimum_required(VERSION 3.10)

# 项目名称
set(GOMAT_VERSION "0.0.1")
project(GoMat VERSION ${GOMAT_VERSION} LANGUAGES CXX)

# 关键添加：生成位置无关代码（兼容Apple Silicon）
set(CMAKE_POSITION_INDEPENDENT_CODE ON)

# 设置 C++ 标准
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)

# 在 add_library 前添加
add_compile_definitions(GOMAT_EXPORT)  # 启用导出
set(CMAKE_CXX_VISIBILITY_PRESET hidden)  # 隐藏非导出符号



# 添加静态库，包含所有源文件
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
)
# 这一行是：为指定的目标（target）设置头文件搜索路径，告诉编译器在哪里查找 #include 的头文件
#可见性修饰符：PUBLIC 表示头文件路径对以下两者有效：目标本身（gomat）：编译 gomat 库（src/*.cpp）时，编译器会搜索指定的头文件路径（include 目录）。
#        依赖目标：任何链接到 gomat 的目标（如 tes，通过 target_link_libraries(tes PRIVATE gomat)）会自动继承这个头文件路径。
#  
#
#
target_include_directories(gomat PUBLIC 
    $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/include>
)

# 设置编译选项
# 修改编译选项（Apple Silicon 必需）
if(APPLE)
    target_compile_options(gomat PRIVATE 
        -fvisibility=default
        -arch arm64  # 明确指定架构
    )
endif()

if(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
    target_compile_options(gomat PRIVATE -Wall -Wextra -Wpedantic)
endif()
set(CMAKE_BUILD_TYPE Debug CACHE STRING "Choose build type: Debug Release")
if(CMAKE_BUILD_TYPE MATCHES Debug)
    target_compile_options(gomat PRIVATE -g)
elseif(CMAKE_BUILD_TYPE MATCHES Release)
    target_compile_options(gomat PRIVATE -O3 -march=native)  # 性能优化
endif()

# 指定头文件目录
# target_include_directories(gomat PUBLIC include)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -g -O0")  # 禁用优化并添加调试符号

# 把test目录作为子目录,这一行确保test目录中的cmake文件也会被解析
add_subdirectory(test)