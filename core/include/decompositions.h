#pragma once
#include "matrix.h"
#include <utility>

namespace gomat {

// 简单 Doolittle LU 分解（不含主元选取），A = L * U
// 返回 (L, U)，其中 L 为单位下三角矩阵，U 为上三角矩阵
std::pair<Matrix<>, Matrix<>> lu(const Matrix<>& A);

}