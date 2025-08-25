#ifndef VECTOR_STREAM_H
#define VECTOR_STREAM_H
#include "vector.h"

namespace gomat {
class VectorStream
{
private:
    Vector& m_vec;
    size_t current_index = 0;

public:
    explicit VectorStream(Vector& vec) : m_vec(vec) {}

    VectorStream& operator<<(double value)
    {
        ++current_index;
        if(current_index >= m_vec.getDimension() + 1)
        {
            throw std::out_of_range("Input data is out of range,vector overflow!");
        }

        m_vec[current_index-1] = value;
        
        return *this;
    }

    VectorStream& operator,(double value)
    {
        return (*this) << value;
    }  
};
} //. namespace gomat
#endif