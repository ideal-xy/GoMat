
set -e

BUILD_DIR="build_test" 

echo "清理旧的构建目录.."
rm -rf ${BUILD_DIR}

echo "配置 CMake... "

cmake -B ${BUILD_DIR} \
    -DCMAKE_BUILD_TYPE=Release \
    -DCMAKE_CXX_COMPILER=/opt/homebrew/opt/llvm/bin/clang++ \
    -DCMAKE_C_COMPILER=/opt/homebrew/opt/llvm/bin/clang

echo "构建项目... "
cmake --build ${BUILD_DIR} 


echo "运行测试... "
cd ${BUILD_DIR}

cd ../test
python test.py
