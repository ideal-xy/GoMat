#include "matrix.h"

using namespace gomat;
inline double& Matrix::getElement(size_t row,size_t col)
{
    return m_is_contiguous ? m_data[col * m_rows + row] : m_mat[row][col];
}

inline double Matrix::getEle(size_t row,size_t col) const
{
    return m_is_contiguous ? m_data[col * m_rows + row] : m_mat[row][col];
}