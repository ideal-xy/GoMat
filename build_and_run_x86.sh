
# 如果任何命令失败，则立即停止脚本
set -e


BUILD_DIR="build_x86"


echo "--- 正在清理旧的构建目录... ---"
rm -rf ${BUILD_DIR}


echo "--- 正在为 x86_64 Debug 模式配置 CMake... ---"
cmake -B ${BUILD_DIR} -DCMAKE_BUILD_TYPE=Debug -DCMAKE_OSX_ARCHITECTURES=x86_64

echo "--- 正在构建项目... ---"
cmake --build ${BUILD_DIR} --verbose



