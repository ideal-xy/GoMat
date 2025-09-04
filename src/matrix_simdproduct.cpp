#include "matrix.h"
#include <immintrin.h>

namespace gomat{

/*
参数的含义：
    对于矩阵乘法 C = A * B，这里第一个参数和第二个参数使用来确定C中需要被计算的分块
    第3个参数，是A的列数，也是B的行数，会用在这个核函数中的主循环之中
    第4，5，6个参数依次是矩阵 C，A，B
*/
void kernel_4x4(size_t start_row,size_t start_col,size_t invariant,Matrix& result,const Matrix& left,const Matrix& right)
{
    // load four simd registers,AVX instructions set
    __m256d res_col_0 = _mm256_loadu_pd(&result(start_row,start_col));
    __m256d res_col_1 = _mm256_loadu_pd(&result(start_row,start_col+1));
    __m256d res_col_2 = _mm256_loadu_pd(&result(start_row,start_col+2));
    __m256d res_col_3 = _mm256_loadu_pd(&result(start_row,start_col+3));

    for (size_t i = 0;i < invariant;i++)
    {
        __m256d left_col_i = _mm256_loadu_pd(&left.at(start_row,i));
    
        __m256d right_row_i_0 = _mm256_set1_pd(right(i,start_col));
        __m256d res_col_0 = _mm256_fmadd_pd(left_col_i,right_row_i_0,res_col_0);

        __m256d right_row_i_1 = _mm256_set1_pd(right(i,start_col+1));
        __m256d res_col_1 = _mm256_fmadd_pd(left_col_i,right_row_i_1,res_col_1);

        __m256d right_row_i_2 = _mm256_set1_pd(right(i,start_col+2));
        __m256d res_col_2 = _mm256_fmadd_pd(left_col_i,right_row_i_2,res_col_2);

         __m256d right_row_i_3 = _mm256_set1_pd(right(i,start_col+3));
         __m256d res_col_3 = _mm256_fmadd_pd(left_col_i,right_row_i_3,res_col_3);
    }

    _mm256_storeu_pd(&result(start_row,start_col),res_col_0);
    _mm256_storeu_pd(&result(start_row,start_col+1),res_col_1);
    _mm256_storeu_pd(&result(start_row,start_col+2),res_col_2);
    _mm256_storeu_pd(&result(start_row,start_col+3),res_col_3);

}


Matrix Matrix::multiplyWithSimd(const Matrix& other) const
{
    Matrix result((*this).getRows(),other.getCols(),true); // contiguous memory with column major

    size_t start_row = 0;
    size_t start_col = 0;
    const size_t rows = getRows();
    const size_t cols = other.getCols();
    const size_t invariant = getCols();

    for (;start_row + 3 < (*this).getRows();start_row+=4)
    {
        for (;start_col + 3 < other.getCols(); start_col+=4)
        {
            kernel_4x4(start_row,start_col,invariant,result,*this,other);
        }
    }

    if (start_row == m_rows)
    {
        for (;start_col < cols;++start_col)
        {
            for (size_t k = 0;k < invariant;++k)
            {
                double val = other(k,start_col);
                for (size_t i = 0;i < rows;++i)
                {
                     result(i,start_col) += val * (*this)(i,k);
                }
            }
        }
    }
    else if (start_col == m_cols)
    {
       for (size_t j =0;j < cols;++j)
       {
           for (size_t k=0;k < invariant;++k)
           {
               double val = other(k,j);
               for (;start_row < rows;++start_row)
               {
                   result(start_row,j) += val * (*this)(start_row,k);
               }
           }
       }
    }
    else
    {
        size_t temp_row = start_row;
        for (size_t j =0;j < cols;++j)
        {
           for (size_t k=0;k < invariant;++k)
           {
               double val = other(k,j);
               for (;start_row < rows;++start_row)
               {
                   result(start_row,j) += val * (*this)(start_row,k);
               }
           }
        }

        for (;start_col < cols;++start_col)
        {
            for (size_t k = 0;k < invariant;++k)
            {
                double val = other(k,start_col);
                for (size_t i = 0;i < temp_row;++i)
                {
                     result(i,start_col) += val * (*this)(i,k);
                }
            }
        }    
    }    
}
    
 
} // namespace 










