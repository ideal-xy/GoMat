#ifndef MATRIX_STREAM_H
#define MATRIX_STREAM_H

#include "matrix.h"
namespace gomat {
class MatrixStream
{
private:
    Matrix& m_subject;
    size_t current_row = 0;
    size_t current_col = 0;

public:
    explicit MatrixStream(Matrix& mat) : m_subject(mat) {} // 禁止隐式构造

    MatrixStream& operator<<(double value)
    {
        
        ++current_col;
        if(current_row >= m_subject.getRows() && current_col >= 1)
        {
            throw std::out_of_range("Too much input data,matrix overflow!");
        }

        m_subject(current_row,current_col-1) = value;
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
#endif