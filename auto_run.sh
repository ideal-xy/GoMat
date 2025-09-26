# 这个自动化脚本是用来测试的，而非debug，debug的目录位于gomat/debug

set -e

SCRIPT_DIR=$(cd "$(dirname "$0")" && pwd)
echo "脚本所在的根目录是: ${SCRIPT_DIR}"

BUILD_DIR="${SCRIPT_DIR}/build_test"  

echo "--清理旧的构建目录--"
rm -rf ${BUILD_DIR}

echo "--配置 CMake--(Release Type)"
cmake -B ${BUILD_DIR} \
    -DCMAKE_BUILD_TYPE=Release \
    -DCMAKE_CXX_COMPILER=/opt/homebrew/opt/llvm/bin/clang++ 

echo "--构建项目--" 
cmake --build ${BUILD_DIR}

echo "--运行测试--"
"${BUILD_DIR}/tests/tes"