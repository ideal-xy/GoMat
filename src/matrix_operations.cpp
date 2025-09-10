#include <iostream>
#include <cassert>
#include <cmath>
#include <algorithm>
#include "vector.h"
#include "matrix.h"
#include "utils.h"
#ifdef __x86_64__
#include <immintrin.h>
#elif defined(__aarch64__)
#include <arm_neon.h>
#endif

#define epsilon 1e-10
using namespace gomat;

double& Matrix::operator()(size_t row,size_t col)
{
    
    return m_is_contiguous ? m_data[col * m_rows + row] : m_mat[row][col];
}


double Matrix::operator()(size_t row,size_t col) const
{
   
    return m_is_contiguous ? m_data[col * m_rows + row] : m_mat[row][col];
}

Matrix Matrix::operator+(const Matrix& mat) const // 加法
{
    // 断言检查：确保两个矩阵的维度相同
    assert(m_cols == mat.m_cols && m_rows == mat.m_rows);
    Matrix sum_mat{m_rows, m_cols};

    // 如果是非连续内存的小矩阵，直接使用双层循环
    if (!m_is_contiguous)
    {
        for (size_t r = 0; r < m_rows; ++r)
        {
            for (size_t c = 0; c < m_cols; ++c)
            {
                sum_mat.m_mat[r][c] = m_mat[r][c] + mat.m_mat[r][c];
            }
        }
    }
    // 如果是使用连续内存的大矩阵，我们使用 SIMD 指令集加速
    else
    {
        const size_t total_size = m_cols * m_rows;
        size_t i = 0;

#if defined(__x86_64__) || defined(_M_X64)
        // --- Intel x86-64 AVX 实现 ---
        // 每次处理 4 个 double
        for (; i + 3 < total_size; i += 4)
        {
            __m256d left_vec = _mm256_loadu_pd(&m_data[i]);
            __m256d right_vec = _mm256_loadu_pd(&(mat.m_data[i]));
            __m256d result_vec = _mm256_add_pd(left_vec, right_vec);
            _mm256_storeu_pd(&(sum_mat.m_data[i]), result_vec);
        }

#elif defined(__aarch64__)
        for (; i + 3 < total_size; i += 4)
        {
            // 加载第一个 2-element vector
            float64x2_t left_vec1 = vld1q_f64(&m_data[i]);
            float64x2_t right_vec1 = vld1q_f64(&(mat.m_data[i]));
            float64x2_t result_vec1 = vaddq_f64(left_vec1, right_vec1);
            vst1q_f64(&(sum_mat.m_data[i]), result_vec1);
            
            // 加载第二个 2-element vector
            float64x2_t left_vec2 = vld1q_f64(&m_data[i + 2]);
            float64x2_t right_vec2 = vld1q_f64(&(mat.m_data[i + 2]));
            float64x2_t result_vec2 = vaddq_f64(left_vec2, right_vec2);
            vst1q_f64(&(sum_mat.m_data[i + 2]), result_vec2);
        }
#endif
        for (; i < total_size; ++i)
        {
            sum_mat.m_data[i] = m_data[i] + mat.m_data[i];
        }
    }

    return sum_mat;
}

Matrix Matrix::operator-(const Matrix& mat) const
{
    assert(m_cols == mat.m_cols || m_rows == mat.m_rows);
    Matrix re_mat{m_rows,m_cols};

    for (size_t r = 0;r < m_rows;++r )
    {
        for (size_t c = 0;c <  m_cols;++c)
        {
            re_mat.m_mat[r][c] = m_mat[r][c] - mat.m_mat[r][c];
        }
    }
    return re_mat;

}

bool Matrix::operator==(const Matrix& other) const // 比较两个矩阵是否相等
{
    if(m_cols != other.m_cols || m_rows != other.m_rows)
    {
        return false;
    }

    for (size_t i=0;i<m_rows;++i)
    {
        if (!std::equal(m_mat[i].begin(),m_mat[i].end(),other.m_mat[i].begin()))
        {
                return false;
        }
    }
    return true;
}

bool Matrix::operator!=(const Matrix& other) const
{
    if(m_cols != other.m_cols || m_rows != other.m_rows)
    {
        return true;
    }

    for (size_t i=0;i<m_rows;++i)
    {
        if (!std::equal(m_mat[i].begin(),m_mat[i].end(),other.m_mat[i].begin()))
        {
            return true;
        }
    }
    return false;
}

Matrix Matrix::operator* (const Matrix& mat) const // 矩阵乘法
{
   return (m_is_contiguous || mat.m_is_contiguous) ? multiplyWithMulThread(mat) : multiplySmall(mat);
}

Matrix Matrix::operator*=(double scalar)
{
    size_t rows = this->getRows();
    size_t cols = this->getCols();
    bool flag = this->isContiguous();
    if (flag)
    {
        Matrix result(rows,cols,true);
        std::transform(m_data.begin(),m_data.end(),result.m_data.begin(),[scalar](double ele){ele *= scalar;});

        return result;
    }
    else
    {
        Matrix result(rows,cols,flag);
        for (size_t i = 0;i < rows;i++)
        {
            std::transform(m_mat[i].begin(),m_mat[i].end(),result.m_mat[i].begin(),[scalar](double ele){ele *= scalar;});
        }

        return result;
    }
}

Matrix& Matrix::operator=(Matrix&& other) noexcept
{
    if (this != &other)
    {
        m_mat = std::move(other.m_mat);
        m_cols = other.m_cols;
        m_rows = other.m_rows;

        other.m_cols = 0;
        other.m_rows = 0;
    }
    return *this;
}

Matrix Matrix::operator,(const Vector& vec) const
{
    if (this->m_rows != vec.getDimension())
    {
        throw std::runtime_error("dismatched rows");
    }

    Matrix result(this->m_rows,this->m_cols+1);
    for(size_t i=0;i<this->m_rows;++i)
    {
        for (size_t j=0;j<this->m_cols;++j)
        {
            result(i,j) = (*this)(i,j);
        }

        result(i,m_cols) = vec[i];
    }
    return result;
}

Matrix Matrix::operator,(const Matrix& mat) const
{
    if (m_rows != mat.m_rows)
    {
        throw std::runtime_error("dismatched rows");
    }

    Matrix result(m_rows,m_cols+mat.m_cols);
    for (size_t i=0;i<m_rows;++i)
    {
        for (size_t j=0;j<m_cols;++j)
        {
            result(i,j) = this->operator()(i,j);
        }
        for(size_t k=m_cols;k<m_cols+mat.m_cols;++k)
        {
            result(i,k) = mat(i,k-m_cols);
        }
    }

    return result;
}

Matrix Matrix::inverse() const
{
    if(m_cols != m_rows) // 检查是否为方阵
    {
        throw std::invalid_argument("Only square matrix has inverse matrix");
    }

    size_t n = m_rows;
    Matrix inverseMat(n,n);
    Matrix augmented_mat(n,2*n); // 空构造
    
    for (size_t i=0;i<n;++i)
    {
        for (size_t j=0;j<n;++j)
        {
            augmented_mat(i,j) = this->operator()(i,j);
        }
        augmented_mat(i,i+n) = 1.0;
    }
    // 我上面写好的gauss消元函数并不能满足使用增广矩阵求逆的需求
    Matrix upper = std::get<1>(augmented_mat.gaussElimination());
    for (size_t i=0;i<n;++i)
    {
        if (std::abs(upper(i,i))< epsilon )
        {
            throw std::runtime_error("This matrix is not invertible");
        }
    }
    Matrix reduced = upper;
    for (size_t col = n;col-->0;) // 从最后一列开始 把增广矩阵的左半部分变为单位矩阵
    {
        double pivot = reduced(col,col); // 定位左边矩阵的主对角线上的元素
        for (size_t j=col;j<2*n;++j) // 主循环，一列一列处理，我们的思路是把主对角线上方的每一个元素都变为0
        {
            reduced(col,j) /= pivot; // 先把主对角线上的元素都变为1
        }
        for(size_t row=col;row-->0;) // 开始处理主对角线上每个元素上面的所有元素，行索引控制
        {
            double factor = reduced(row,col); // 由于该列主对角线元素已经是1
            
            for (size_t j = col;j<2*n;++j) // row行其它的元素也要处理
            {
                reduced(row,j) -= factor * reduced(col,j);
            }
        }
        // 提取右边的矩阵
        // inverseMat = Matrix(n,n);
        for (size_t i =0;i<n;++i)
        {
            for (size_t j=0;j<n;++j)
            {
                inverseMat(i,j) = reduced(i,j+n);
            }
        }
    }
    return inverseMat;
}

Matrix Matrix::multiplyDense(const Matrix& other) const
{
    Matrix mul_mat;
    mul_mat.m_mat.assign(m_rows, std::vector<double>(other.m_cols, 0.0)); // 批量初始化，这比resize好
    mul_mat.m_cols = other.m_cols;
    mul_mat.m_rows = m_rows;

    for (size_t r = 0; r < m_rows; ++r) 
    {
        for (size_t k = 0; k < m_cols; ++k) 
        {
            const double val = (*this)(r,k);
            for (size_t c = 0; c < other.m_cols; ++c) 
            {
                mul_mat(r,c) += val * other(k,c);
            }
        }
    }
    return mul_mat;
}

Matrix Matrix::multiplyDiagonal(const Matrix& other) const  
{
    Matrix mul_mat;
    mul_mat.m_mat.assign(m_rows, std::vector<double>(other.m_cols, 0.0)); // 批量初始化，这比resize好
    mul_mat.m_cols = other.m_cols;
    mul_mat.m_rows = m_rows;

    for(size_t i=0;i<m_rows;++i)
    {
        for(size_t j=0;j<other.m_cols;++j)
        {
            mul_mat(i,j) = (*this)(i,j) * other(j,j);
        }
    }
    return mul_mat;
}

Matrix Matrix::multiplyUpperTriangle(const Matrix& other) const
{
    Matrix mul_mat;
    mul_mat.m_mat.assign(m_rows, std::vector<double>(other.m_cols, 0.0)); // 批量初始化，这比resize好
    mul_mat.m_cols = other.m_cols;
    mul_mat.m_rows = m_rows;

    for (size_t r = 0; r < m_rows; ++r) {
        for (size_t k = 0; k < m_cols; ++k) {
            const double val = (*this)(r,k);
            for (size_t c = 0; c < other.m_cols; ++c) {
                mul_mat(r,c) += val * other(k,c);
            }
        }
    }
   return mul_mat;

}

Matrix Matrix::multiplyLarge(const Matrix& mat) const
{
    constexpr size_t blockSize = 64;// 这个数字取决于机器
    Matrix result((*this).getRows(),mat.getCols(),true); // 我们强制连续存储

    for(size_t total_row = 0;total_row < this->getRows();total_row += blockSize)
    {
        for(size_t total_col = 0;total_col < mat.getCols();total_col += blockSize)
        {
            for(size_t temp_index=0;temp_index < this->getCols();temp_index += blockSize)
            {
                for(size_t i=total_row;i<std::min(total_row+blockSize,this->getRows());++i)
                {
                    for(size_t j=total_col;j<std::min(total_col+blockSize,mat.getCols());++j)
                    {
                        double sum = result(i,j);
                        for(size_t k=temp_index;k<std::min(temp_index+blockSize,this->getCols());++k)
                        {
                            sum += (*this)(i,k) * mat(k,j);
                        }
                        result(i,j) = sum;
                    }
                }      
            }
        }
    }
    return result;
} 

Matrix Matrix::multiplySmall(const Matrix& mat) const
{
    if(mat.isDiagonal(1e-10))
    {
        return multiplyDiagonal(mat);
    }
    return multiplyDense(mat);
}

void Matrix::erase()
{
    this->fillWithOneValue(0.0);
}


Matrix Matrix::subMatrix(size_t start_row,size_t end_row,size_t start_col,size_t end_col)  const
 //子矩阵确定方法是两条横线，两条竖线交叉产生一个矩阵
{
    Matrix mat(end_row - start_row+1,end_col - start_col+1);

    for(size_t i=start_row;i < end_row;++i)
    {
        std::copy(m_mat[i].begin(),m_mat[1].begin()+(end_col-start_col),mat.m_mat[i-start_row].begin());
    }
    
    return mat;
}

double Matrix::blockSum(size_t start_row,size_t end_row,size_t start_col,size_t end_col) const
{
    double sum = 0.0;

    if (this->isContiguous())
    {
        for (size_t i =start_row;i <= end_row;++i )
        {
            for (size_t j = start_col;j <= end_col; ++j)
            {
                sum += this->operator()(i,j);
            }
        }
        return sum;
    }
    else
    {
        for (size_t i=start_row; i <= end_row;++i)
        {
            sum += std::accumulate(m_mat[i].begin()+start_col,m_mat[i].begin()+end_col,0.0, [](double acc,double ele){return acc + ele*ele;});
        }
        
        return sum;
    }
}
