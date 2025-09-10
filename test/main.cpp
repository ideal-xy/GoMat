#include <iostream>
#include <chrono>
#include "gomat.h"
#include <random>
#include <thread>
#include <Eigen/Dense>


int main() 
{
    unsigned int n_threads = std::thread::hardware_concurrency();
    std::cout << "硬件支持 " << n_threads << " 个并发线程。" << std::endl;
    Eigen::setNbThreads(n_threads); 
 
.
    gomat::Matrix mat(1000, 1000, true);
    gomat::Matrix mat2(1000, 1000, true);
    mat.random(1, 2);
    mat2.random(2, 3);


    const int SIZE = 1000;
    Eigen::MatrixXd A = Eigen::MatrixXd::Random(SIZE, SIZE);
    Eigen::MatrixXd B = Eigen::MatrixXd::Random(SIZE, SIZE);
    
    Eigen::MatrixXd D = A * B;
    auto start3 = std::chrono::high_resolution_clock::now();
    Eigen::MatrixXd C = A * B;
    auto end3 = std::chrono::high_resolution_clock::now();
    auto duration3 = std::chrono::duration_cast<std::chrono::microseconds>(end3 - start3);
    std::cout << "Eigen (多线程) 的1000*1000矩阵乘法: " << duration3.count() / 1000000.0 << "s" << std::endl;


    auto start = std::chrono::high_resolution_clock::now();
    mat.multiplyWithMulThread(mat2);
    auto end = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start);
    std::cout << "GoMat (多线程) 的1000*1000矩阵乘法: " << duration.count() / 1'000'000.0 << " s" << std::endl;

    
    return 0;
}