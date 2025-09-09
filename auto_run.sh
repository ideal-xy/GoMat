#!/bin/bash
# 正确的脚本
set -e

BUILD_DIR="build_test" # 建议改名以匹配真实架构

echo "--- 正在清理旧的构建目录... ---"
rm -rf ${BUILD_DIR}

echo "--- 正在为本地 arm64 架构配置 CMake... ---"

cmake -B ${BUILD_DIR} -DCMAKE_BUILD_TYPE=Debug

echo "--- 正在构建项目... ---"
cmake --build ${BUILD_DIR} 

# 假设测试脚本在项目根目录的 test 文件夹下
echo "--- 正在运行测试... ---"
cd ${BUILD_DIR}

cd ../test
python test.py
