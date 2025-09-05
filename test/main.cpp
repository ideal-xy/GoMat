#include <iostream>
#include <chrono>
#include "gomat.h"

int main() 
{
    
    gomat::Matrix mat(1000,1000);
    gomat::Matrix mat2(1000,1000);

    mat.random(1,2);
    mat2.random(2,3);

    auto start = std::chrono::high_resolution_clock::now();
    mat * mat2;
    auto end = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start);
    std::cout << "SIMD type takes time: " << duration.count() / 1'000'000.0 << " s" << std::endl;

    auto start1 = std::chrono::high_resolution_clock::now();
    mat.multiplyDense(mat2);
    auto end1 = std::chrono::high_resolution_clock::now();
    auto duration1 = std::chrono::duration_cast<std::chrono::microseconds>(end1 - start1);
    std::cout << "Normal type takes time: " << duration1.count() / 1'000'000.0  << " s"<< std::endl;

    std::cout << "Times is : " << (duration1.count()) / (duration.count()) << std::endl;

    return 0;
}