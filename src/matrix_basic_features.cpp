#include <iostream>
#include <cassert>
#include <cmath>
#include <stdexcept>
#include <ranges>
#include <iomanip>
#include <optional>
#include <algorithm>
#include <sstream>
#include <fstream>
#include "vector.h"
#include "matrix.h"
#include "utils.h"

using namespace gomat;

size_t Matrix::getCols() const
{
    return m_cols;
}


size_t Matrix::getRows() const
{
    return m_rows;
}

bool Matrix::isSquare() const
{
    return m_rows == m_cols;
}

bool Matrix::isEmpty(double epsilon = 1e-10) const
{
    for (size_t i = 0;i<m_rows;++i)
    {
        for(size_t j=0;j<m_cols;++j)
        {
            if(std::abs((*this)(i,j)) > epsilon)
            {
                return false;
            }
        }
    }

    return true;
}

bool Matrix::isDiagonal(double epsilon) const
{
    for (size_t i = 0;i<m_rows;++i)
    {
        for(size_t j=0;j<m_cols;++j)
        {
            if(i != j && std::abs((*this)(i,j)) < epsilon)
            {
                return false;
            }
        }
    }
    return true;
}

bool Matrix::isUpperTriangular(double epsilon) const
{
    for(size_t i=0;i<m_rows;++i)
    {
        for(size_t j=0;j<i;++j)
        {
            if (std::abs((*this)(i,j)) > epsilon)
            {
                return false;
            }
        }
    }

    return true;
}

bool Matrix::isLowerTriangular(double epsilon) const
{
    for(size_t i=0;i<m_rows;++i)
    {
        for(size_t j=i+1;j<m_cols;++j)
        {
            if(std::abs((*this)(i,j)) > epsilon)
            {
                return false;
            }
        }
    }

    return true;
}

bool Matrix::isLarge() const
{
    return ((m_cols * m_rows) >= 40000);
}