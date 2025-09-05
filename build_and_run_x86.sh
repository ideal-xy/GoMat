#!/bin/bash

# 此脚本用于配置、构建并运行 GoMat 库的测试。
# 它专门为 Apple Silicon Mac 上的 Rosetta 2 环境指定了 x86_64 架构。

# 如果任何命令失败，则立即停止脚本
set -e

# 定义构建目录的名称
BUILD_DIR="build_x86"

# --- 配置步骤 ---
# 1. 清理旧的构建目录，确保一个干净的开始。
echo "--- 正在清理旧的构建目录... ---"
rm -rf ${BUILD_DIR}

# 2. 使用 CMake 配置项目。
#    -B ${BUILD_DIR}: 指定构建目录。
#    -DCMAKE_BUILD_TYPE=Debug: 设置构建类型为 Debug（包含调试符号，优化较少）。
#    -DCMAKE_OSX_ARCHITECTURES=x86_64: 这是最关键的命令！它告诉 CMake/Clang
#                                     为 Intel x86_64 架构编译代码。
echo "--- 正在为 x86_64 Debug 模式配置 CMake... ---"
cmake -B ${BUILD_DIR} -DCMAKE_BUILD_TYPE=Debug -DCMAKE_OSX_ARCHITECTURES=x86_64

# --- 构建步骤 ---
# 1. 使用上一步的配置来构建项目。
#    --build ${BUILD_DIR}: 指定要构建的目录。
#    --verbose: 显示详细的编译器命令，便于调试。
echo "--- 正在构建项目... ---"
cmake --build ${BUILD_DIR} --verbose

# --- 运行步骤 ---
# 测试程序已在 test/CMakeLists.txt 中通过 POST_BUILD 命令设置为自动运行。
# 构建完成后，测试结果会立刻显示。

