#include <iostream>
#include <chrono>
#include "gomat.h"
#include "decompositions.h"
#include "matrix.h"
#include <random>
#include <thread>
#include <Eigen/Dense>
#include <cmath>

int main()
{
    gomat::Matrix<int> mat(3,3);
    mat << 1,2,3,4,5,6,7,8,9;
    gomat::Matrix<double> mat2(3,3);
    mat2 << 1,2,3,4,5,6,7,8,9;
    
    auto view = mat.transposeView() + mat.scalaredView(2);
    auto view2 = view * view + view;
    auto result = view2.eval();

    gomat::Matrix<double> mat3(2000,5000);
    mat3.random(1,2);

    gomat::Matrix<double> mat4(5000,2000);
    mat4.random(1,2);

    std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();
    auto result2 = mat3 * mat4;
    std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
    std::cout << "Time difference = " << std::chrono::duration_cast<std::chrono::milliseconds>(end - begin).count() << "[ms]" << std::endl;
    
    return 0;

}
