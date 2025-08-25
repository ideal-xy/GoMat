#include "matrix.h"

using namespace gomat;
double& Matrix::getElement(size_t row,size_t col)
{
    return m_is_contiguous ? m_data[row * m_cols + col] : m_mat[row][col];
}

double Matrix::getEle(size_t row,size_t col) const
{
    return m_is_contiguous ? m_data[row * m_cols + col] : m_mat[row][col];
}