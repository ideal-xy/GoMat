#include <iostream>
#include <cmath>
#include <stdexcept>
#include <ranges>
#include <algorithm>
#include "vector.h"
#include "matrix.h"
#include "utils.h"

using namespace gomat;

std::tuple<Matrix,Matrix> Matrix::LuDecomposition() const
{
    if(!this->isSquare())
    {
        throw std::runtime_error("Only square matrix has LuDecomposition!");
    }

    Matrix L(m_rows,m_rows),U(m_rows,m_rows);
    L.m_mat.resize(m_cols,std::vector<double>(m_cols,0.0));
    U.m_mat.resize(m_cols,std::vector<double>(m_cols,0.0));

    for (size_t i=0;i<m_cols;++i)
    {
        L(i,i) = 1.0;
    }

    if((*this)(0,0) == 0.0)
    {
        throw std::runtime_error("Singular matrix is not feasible with Lu decomposition!");
    }

    U(0,0) = (*this)(0,0);

    // 计算出了 L的第一列和U的第一行
    for (size_t i=1;i<m_rows;++i)
    {
        U(0,i) = (*this)(0,i);
        L(i,0) = (*this)(i,0) / U(0,0); 
    }

    for(size_t i=1;i<m_cols;++i)
    {
        for(size_t j=i;j<m_cols;++j)
        {
            double sum = 0.0;
            for (size_t k=1;k<i;++k)
            {
                sum += L(i,k) * U(k,j);
            }

            U(i,j) = (*this)(i,j) - sum;
        }
            
        if (U(i,i) == 0.0)
        {
            throw std::runtime_error("Singular matrix is not feasible with LuDecomposition");  // 后续的运算会使用到U(i,i)
        }

        for(size_t j=i;j<m_cols;++j)
        {
            double sum = 0.0;
            for (size_t k=1;k<i;++k)
            {
                sum += L(j,k)*U(k,i);
            }

            L(j,i) = ((*this)(j,i) - sum) / U(i,i);
        }
    }
    return std::make_tuple(L,U);
}


std::tuple<Matrix,Matrix> Matrix::QrDecompostion() const
{
    Matrix mat = (*this).transpose();
    std::vector<std::vector<double>> vec(m_cols); // 用于存储我们Schmidt正交化后的n个列向量 其数量是确定的
    std::vector<std::vector<double>> coeff(m_cols); // 用于存储Schemidt正交化得到的系数
    Matrix Q,R;

    for(size_t i=0;i<m_cols;++i)
    {
        vec[i].reserve(m_rows);
        coeff[i].reserve(m_cols);
    }

    vec[0] = Utils::normalize(m_mat[0]);
    coeff[0].push_back(1 * Utils::norm_of_stdVector(m_mat[0]));
    coeff[0].insert(coeff[0].end(),coeff[0].capacity()-coeff[0].size(),0.0); // 使用0把剩下的元素填满
   
    
    for(size_t col=1;col<m_cols;++col) // 主循环 遍历被分解的矩阵的每一个列向量
    {
        Vector temp_vec(m_rows); // 创建一个临时的向量 方便我们计算 \beta 也就是施密特正交化后得到的向量

        for(size_t k=0;k<col;++k)
        {
            double scalar = (Utils::vecDotProduct(m_mat[col],vec[col-1-k]))/(Utils::vecDotProduct(vec[col-1-k],vec[col-1-k]));
            temp_vec = temp_vec + Utils::mulyiply(scalar,vec[col-1-k]);
            coeff[col].push_back(scalar * Utils::norm_of_stdVector(vec[col-1-k]));           
        }

        std::vector<double> minus_vec (Utils::frontVectorMinusRearVector(m_mat[col],temp_vec.to_std_vector()));
        vec[col] = Utils::normalize(minus_vec);
        coeff[col].push_back(1 * Utils::norm_of_stdVector(minus_vec)); //

        coeff[col].insert(coeff[col].end(),coeff[col].capacity()-coeff[col].size(),0.0); // 使用0把剩下的元素填满
        
    } //  至此，我们只需要转置操作就可以得到完整的[\beta_1,\beta_2,...,\beta_n]，以及上三角系数矩阵
     
    Matrix mat1= vec;
    Q = mat1.transpose(); //  [\beta_1,\beta_2,...,\beta_n]

    Matrix mat2 = coeff;
    R = mat2.transpose(); // 上三角的系数矩阵

    return std::make_tuple(Q,R);
}

std::tuple<Matrix,Matrix,Matrix> Matrix::SvdDecomposition() const
{
    size_t m = m_rows;
    size_t n = m_cols;

    Matrix U(m,m),V(n,n);
    Matrix A(*this);

    Utils::bidiagonalize(A,U,V); 
    Utils::QrIteration(A,U,V);

    return std::make_tuple(U,A,V);
}

