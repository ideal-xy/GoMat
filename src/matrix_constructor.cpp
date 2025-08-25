#include "matrix.h"


using namespace gomat;
Matrix::Matrix(size_t rows,size_t cols,bool use_contiguous):m_rows(rows),m_cols(cols),m_is_contiguous(use_contiguous)
{
    if(use_contiguous || rows * cols > 10000)
    {
        m_data.reserve(2 * rows * cols);
        m_data.resize(rows * cols);
        m_is_contiguous = true;
    }
    else{
        m_mat.reserve(2 * rows);
        m_mat.resize(rows,std::vector<double>(cols));
    }
}

Matrix::Matrix(const std::vector<std::vector<double>>& data) : m_mat(data),m_rows(data.size()),
                                                m_cols(data[0].size()),m_is_contiguous(false)
{
    
}

Matrix::Matrix(const std::vector<double>& data,size_t rows,size_t cols):m_data(data),m_rows(rows),
                                                 m_cols(cols),m_is_contiguous(true)
{
    if(data.size() != rows * cols)
    {
        throw std::invalid_argument("input error,size dismatched! ");
    }
} 

Matrix::Matrix(Matrix&& other) noexcept
    : m_rows(other.m_rows),
      m_cols(other.m_cols),
      m_is_contiguous(other.m_is_contiguous),
      m_data(std::move(other.m_data)),    // 移动而非拷贝
      m_mat(std::move(other.m_mat))       // 移动而非拷贝
{
    other.m_rows = 0;
    other.m_cols = 0;
    other.m_is_contiguous = true;
}

Matrix::Matrix(std::vector<std::vector<double>>&& other) noexcept : m_mat(std::move(other)),m_rows(other.size()),m_cols(other[0].size())
{

}

Matrix::Matrix(const Matrix& other)
    : m_rows(other.m_rows),
      m_cols(other.m_cols),
      m_is_contiguous(other.m_is_contiguous),
      m_data(other.m_data),    // vector 自动深拷贝
      m_mat(other.m_mat)      // vector<vector> 自动深拷贝
{
    // 确保存储方式与标志一致
    if (m_is_contiguous && !m_mat.empty()) {
        throw std::logic_error("Contiguous matrix should not have m_mat data");
    }
    if (!m_is_contiguous && !m_data.empty()) {
        throw std::logic_error("Non-contiguous matrix should not have m_data");
    }
}
                                                




