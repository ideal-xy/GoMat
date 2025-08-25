#include <iostream>
#include <cmath>
#include <stdexcept>
#include <ranges>
#include <iomanip>
#include <algorithm>
#include "vector.h"
#include "matrix.h"
#include "utils.h"

#define epsilon 1e-10
using namespace gomat;

double Matrix::determiant() const {
    if (m_cols != m_rows) 
        throw std::runtime_error("Only square matrices have determinants");

    auto [swap_times, upTriangle] = this->gaussElimination();
    double det = 1.0;
    
    for (size_t i = 0; i < m_cols; ++i) {
        det *= upTriangle(i, i);  // 主对角线乘积
    }
    
    return det * ((swap_times % 2 == 0) ? 1 : -1);
}



bool Matrix::isInvertible() const
{
    return (*this).determiant() < epsilon;
}

double Matrix::trace() const
{
    if (m_cols != m_rows)
    {
        throw std::domain_error("Only square matrix has trace");
    }

    double traceValue = 0.0;

    for (int i=0;i<m_cols;++i)
    {
        traceValue = traceValue + m_mat[i][i];
    }
    return traceValue;
}


int Matrix::rank() const {
    Matrix temp = std::get<1>(gaussElimination());
    int rank = 0;
    
    for (size_t i = 0; i < m_rows; ++i) {
        bool isNotAllZero = std::any_of(
            temp.m_mat[i].begin(), 
            temp.m_mat[i].end(),
            [](double x) { return std::abs(x) >= 1e-10; }
        );
        if (isNotAllZero) ++rank;
    }
    
    return rank;
}

double Matrix::sum() const
{
    double total = 0.0;
    for (auto& row:m_mat) 
    {
        total += std::accumulate(row.begin(),row.end(),0.0);
    }
    return total;
}

double Matrix::product() const
{
    double total = 1.0;
    // 连续存储类型
    if(m_is_contiguous)
    {
        for(double val:m_data)
        {
            total *= val;
        }
    }
    else
    {
        for(auto row:m_mat)
        {
            for(double val : row)
            {
                total *= val;
            }
        }
    }
    return total;
}

double Matrix::frobeniusNorm() const
{
    double total = 0.0;
    if(m_is_contiguous)
    {
        total = std::accumulate(m_data.begin(),m_data.end(),0.0);
    }
    else{
        for (auto& row:m_mat)
        {
            total += std::accumulate(row.begin(),row.end(),0.0, [](double acc,double ele){return acc + ele*ele;});
        }
    }
    return std::sqrt(total);
}


