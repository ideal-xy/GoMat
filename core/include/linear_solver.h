#ifndef LINEAR_SOLVER_H
#define LINEAR_SOLVER_H

#include "matrix.h"
#include "vector.h"
#include <cstddef>
namespace gomat {
enum class SolutionType {
    NoSolution,       // 无解
    UniqueSolution,   // 唯一解
    InfiniteSolutions // 无穷多解
};

struct Solution
{
    Solution() = default;

    SolutionType type;

    // 唯一解
    Vector unique_solution;

    // 无穷多解 
    Vector particular_solution;
    size_t nums{0};
    std::vector<std::vector<double>> basic_solution_system;

    // 其他信息
    std::string output_message;
   
    
};

class Linear_Solver
{
private:
    double m_tolerance = 1e-6;
    int m_max_iteration_times = 10000;

public:
    Linear_Solver() = default; 
    Linear_Solver(double tolerance,int max_iteration_times);

    //矩阵的分解算法
    Solution gausslian(Matrix<>& A_mat,Vector b_vec) const;
    Solution Lu(Matrix<>& A_mat,Vector b_vec) const;

    // 说一就是一，说二就是二
    
    
    
    friend std::ostream& operator<<(std::ostream& out,Solution solution);
};
} //. namespace
#endif