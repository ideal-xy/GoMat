#!/usr/bin/env bash
set -euo pipefail

# 该脚本用于一键：Release 配置、构建并运行测试。

SCRIPT_DIR=$(cd "$(dirname "$0")" && pwd)
echo "脚本所在的根目录是: ${SCRIPT_DIR}"

BUILD_DIR="${SCRIPT_DIR}/build_release"

echo "--清理旧的构建目录--"
rm -rf "${BUILD_DIR}"
mkdir -p "${BUILD_DIR}"

echo "--配置 CMake--(Release + 优化标志)"

CXXFLAGS="-O3 -ffast-math -fstrict-aliasing -DNDEBUG"
UNAME_M=$(uname -m)
if [[ "$UNAME_M" == "arm64" ]]; then
  CXXFLAGS+=" -mcpu=apple-m3"
fi

cmake -S "${SCRIPT_DIR}" -B "${BUILD_DIR}" \
  -DCMAKE_BUILD_TYPE=Release \
  -DCMAKE_OSX_ARCHITECTURES=arm64 \
  -DCMAKE_CXX_FLAGS="${CXXFLAGS}" \
  -DCMAKE_INTERPROCEDURAL_OPTIMIZATION=ON \
  -DCMAKE_EXPORT_COMPILE_COMMANDS=ON

echo "--并行构建 tes 目标--"
cmake --build "${BUILD_DIR}" --target tes --parallel

# 为 clangd 提供编译数据库，链接到项目根目录
if [ -f "${BUILD_DIR}/compile_commands.json" ]; then
  ln -sf "${BUILD_DIR}/compile_commands.json" "${SCRIPT_DIR}/compile_commands.json"
fi

echo "--运行测试--"

# 为 OpenMP 设置线程绑定与线程数，提升并行一致性
THREADS="$(sysctl -n hw.logicalcpu 2>/dev/null || getconf _NPROCESSORS_ONLN 2>/dev/null || nproc 2>/dev/null || echo 8)"
export OMP_NUM_THREADS="${THREADS}"
export OMP_PROC_BIND=close
export OMP_PLACES=cores

"${BUILD_DIR}/tests/tes"