#include <iostream>
#include <cmath>
#include <algorithm>
#include "vector.h"
#include "matrix.h"
#include "utils.h"

using namespace gomat;

void Matrix::toContiguous()
{
    if(m_is_contiguous) return;
    
    m_is_contiguous = true; // 标记更新
    m_data.resize(m_rows * m_cols);
    for(size_t i=1;i<m_rows;++i)
    {
        std::copy(m_mat[i].begin(),m_mat[i].end(),m_data.begin()+i * m_rows);
    }
    m_mat.clear();
}

void Matrix::toNonContiguous()
{
    if(!m_is_contiguous) return;

    m_is_contiguous = false;
    m_mat.resize(m_rows,std::vector<double>(m_cols));
    for(size_t i=0;i<m_rows;++i)
    {
        std::copy(m_data.begin() + i * m_rows,m_data.begin() + ((i+1) * m_rows) -1,m_mat[i].begin());
    }
}

void Matrix::ensureContiguousIfLarge()
{
    if(m_rows * m_cols > 10000 && !m_is_contiguous)
    {
        toContiguous();
    }
}

void Matrix::rowInterchange(size_t i,size_t j)
{
    if(i == j || i > m_rows ||j > m_rows)
        {
            throw std::invalid_argument("INVALID PARAMETERS");
        }
    if(!m_is_contiguous)
    {
        std::swap(m_mat[i],m_mat[j]);
    }
    else
    {
        std::swap_ranges(
            m_data.begin() + i * m_cols,
            m_data.begin() + (i+1) * m_cols,
            m_data.begin() + j * m_cols
        );
    }
}

void Matrix::resize(size_t row,size_t col)
{
    m_rows = row;
    m_cols = col;
    
    if(m_is_contiguous)
    {
        std::vector<double>vec (row * col,0.0);
        m_data.swap(vec);
    }
    if(!m_is_contiguous)
    {
        std::vector<std::vector<double>> vec(row,std::vector<double>(col,0));
        m_mat.swap(vec);
    }
}

void Matrix::rowMultiply(size_t row,double factor)
{
    if (row>m_rows)
    {
        throw std::invalid_argument("INVALID PARAMETERS");
    }

    for (size_t j=0;j<m_cols;++j)
    {
        (*this)(row,j) *= factor;
    }
}

void Matrix::rowAddition(size_t row1,size_t row2,double factor) // 把第row2行的factor倍数 加到row1行
{
    if (row1>m_rows || row2 > m_rows)
    {
        throw std::invalid_argument("INVALID PARAMETERS");
    }

    std::transform(m_mat[row1].begin(),m_mat[row1].begin(),m_mat[row2].begin(),
                m_mat[row1].begin(),[factor](double x,double y){return x + factor * y;});
}

std::tuple<int, Matrix> Matrix::gaussElimination() const 
{
    size_t currentRow = 0;
    int swap_times = 0;
    Matrix copy = *this;

    for (size_t col = 0; col < copy.m_cols && currentRow < copy.m_rows; ++col) {
        // 寻找主元
        size_t pivot_row = currentRow;
        double max_val = std::abs(copy(currentRow, col));
        
        for (size_t r = currentRow + 1; r < copy.m_rows; ++r) {
            if (std::abs(copy(r, col)) > max_val) {
                max_val = std::abs(copy(r, col));
                pivot_row = r;
            }
        }
        // 跳过零列
        if (max_val < 1e-10) continue;
        
        // 行交换
        if (pivot_row != currentRow) {
            copy.rowInterchange(currentRow, pivot_row);
            ++swap_times;
        }
        // 消元
        double pivot = copy(currentRow, col);
        if (std::abs(pivot) < 1e-10) continue;
        for (size_t r = currentRow + 1; r < copy.m_rows; ++r) {
            double factor = copy(r, col) / pivot;
            for (size_t c = col; c < copy.m_cols; ++c) {
                copy(r, c) -= factor * copy(currentRow, c);
            }
        }
        ++currentRow;
    }
    
    return std::make_tuple(swap_times, copy);
}


void Matrix::transposeInPlace()
{
    if (m_cols != m_rows)
    {
        throw std::invalid_argument("INVALID PARAMATERS");
    }

    for (size_t i=0;i<m_rows;++i)
    {
        for (size_t j=i+1;j<m_cols;++j)
        {
            std::swap((*this)(i,j),(*this)(i,j));
        }
    }        
}

Matrix Matrix::transposeByBlock(size_t blockSize) const
{
    Matrix transposed(this->getCols(),this->getRows());
    for (size_t i=0;i<m_rows;i+=blockSize)
    {
        for (size_t j=0;j<m_cols;j+=blockSize)
        {
            // 计算当前块的边界（避免越界）
            size_t i_end = static_cast<size_t>(std::min(i + blockSize, m_rows));
            size_t j_end = static_cast<size_t>(std::min(j + blockSize,m_cols));

            for (size_t ii=i;ii<i_end;++ii)
            {
                for (size_t jj=j;jj<j_end;++jj)
                {
                    transposed(jj,ii) = (*this)(ii,jj);
                }
            }
        }
    }

    return transposed;   
}

Matrix Matrix::transpose() const
{
    if ((*this).isEmpty(1e-10))
    {
        throw std::runtime_error("Empty matrix has no transpose");
    }

    Matrix mat;
    std::vector<std::vector<double>> vec(m_cols,std::vector<double>(m_rows));
    mat.m_cols = m_rows;
    mat.m_rows = m_cols;

    for(size_t j = 0;j<m_cols;++j)
    {
        for(size_t i=0;i<m_rows;++i)
        {
            vec[j][i] = (*this)(i,j);
        }
    }

    mat.m_mat = std::move(vec);
    return mat;
}

void Matrix::fillWithData(const std::vector<std::vector<double>>& data)
{
    if (data.size() != m_rows && data[0].size() != m_cols)
    {
        throw std::invalid_argument("dismatched dimensions");
    }
    m_cols = data[0].size();
    m_rows = data.size();
    m_mat = std::move(data);
}

void Matrix::fillWithOneValue(double value)
{
    if(!m_is_contiguous)
    {
        for (auto& row:m_mat)
        {
            std::fill(row.begin(),row.end(),value);
        }
    }
    else
    {
        std::fill(m_data.begin(),m_data.end(),value);
    }
}

void Matrix::fillScaledIdentity(double k)
{
    if (m_cols != m_rows)
    {
        throw std::invalid_argument("NOT SQUARE MATRIX");
    }
    this->fillWithOneValue(0.0);

    for (size_t i =0;i<m_rows;++i)
    {
        (*this)(i,i) = k;
    }

}

Matrix Matrix::replicate(size_t rowTimes,size_t colTimes) const
{
    if(m_is_contiguous)
    {
        // 
    }

    if(!m_is_contiguous)
    {
        // 
    }
    return Matrix(1,1);
}



    