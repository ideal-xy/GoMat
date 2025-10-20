#ifndef MATRIX_STREAM_H
#define MATRIX_STREAM_H

#include "matrix.h"
namespace gomat {

// 将 MatrixStream 改为模板，使其支持所有 Matrix<T, L, Alloc>
template<typename T, Layout L, typename Alloc>
class MatrixStream
{
private:
    Matrix<T, L, Alloc>& m_subject;
    size_t current_row = 0;
    size_t current_col = 0;

public:
    explicit MatrixStream(Matrix<T, L, Alloc>& mat) : m_subject(mat) {}

    MatrixStream& operator<<(double value)
    {
        ++current_col;
        if(current_row >= m_subject.getRows() && current_col >= 1)
        {
            throw std::out_of_range("Too much input data,matrix overflow!");
        }

        m_subject(current_row,current_col-1) = static_cast<T>(value);
        if(current_col == m_subject.getCols())
        {
            ++current_row;
            current_col = 0;
        }
        return *this;
    }

    MatrixStream& operator,(double value)
    {
        return (*this) << value;
    }
};

} //namespace gomat

// 内联提供 Matrix 的流式输入成员函数模板定义，支持所有模版参数
namespace gomat {
template<typename T, Layout L, typename Alloc>
inline MatrixStream<T, L, Alloc> Matrix<T, L, Alloc>::operator<<(double value)
{
    MatrixStream<T, L, Alloc> stream(static_cast<Matrix<T, L, Alloc>&>(*this));
    stream << value;
    return stream;
}
}
#endif