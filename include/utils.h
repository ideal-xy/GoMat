#ifndef UTILS_H
#define UTILS_H

#include <iostream>
#include <vector>
#include <numeric>
#include <cstddef>
#include <utility>
#include "matrix.h"
namespace gomat {
class Vector;

class Utils
{
public:

    Utils() = default;
    static bool checkVectorSize(std::vector<double> vec1,std::vector<double> vec2);
    static double vecDotProduct(std::vector<double> vec1,std::vector<double> vec2);
    static int gcd(size_t a,size_t b);
    static bool isZero(double value);

    static std::vector<double> mulyiply(double scalar,std::vector<double> vec);
    static double norm_of_stdVector(std::vector<double> vec);
    static std::vector<double> vectorAdd(std::vector<double> vec1,std::vector<double> vec2);
    static std::vector<double> normalize(std::vector<double> vec);

    // 比较特殊的运算符重载
    static std::vector<double> frontVectorMinusRearVector(std::vector<double> vec1,std::vector<double> vec2);

    // householder变换部分的辅助函数
    static Vector householder_normal_vector(const Vector& vec);
    static void householder_left_mulyiply(Matrix& A,Matrix& U,const Vector& vec);
    static void householder_right_multiply(Matrix& A,Matrix& V,const Vector& vec);
    static void bidiagonalize(Matrix& A,Matrix& U,Matrix& V);

    //QR迭代辅助函数
    static void compute_givens(double a,double b,double& cos_theta,double& sin_theta);
    static void apply_givens_left(Matrix& A,size_t i,size_t j,double cos_theta,double sin_theta);
    static void apply_givens_right(Matrix& A,size_t i,size_t j,double cos_theta,double sin_theta);
    static double compute_wilkinson_shift(double a,double b,double c,double d);
    static void QrIteration(Matrix& B,Matrix& U,Matrix& V,double epsilon = 1e-10,int max_iteration = 100);

    // 矩阵的一些辅助检查
    static void checkMultiplicationDims(const Matrix& a, const Matrix& b);
};
} //namespace
#endif