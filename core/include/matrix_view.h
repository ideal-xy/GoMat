#pragma once

#include <functional>
#include "matrix.h"
#include "vector.h"
#include <cstddef>

namespace gomat{

template <typename LeftDerived, typename RightDerived>
class AddExpr;

template <typename LeftDerived, typename RightDerived>
class MulExpr;

template <typename Derived>
class MatrixView
{
public:
    Derived& derived() {return static_cast<Derived&>(*this);}
    const Derived& derived() const {return static_cast<const Derived&>(*this);}

    double operator()(size_t row,size_t col) const
    {
        return derived().operator()(row,col);
    }

    size_t getRows() const {return derived().getRows();}
    size_t getCols() const {return derived().getCols();}

    Matrix eval() const {
    // if (getRows() == 1) {  // 单行矩阵
    //     Matrix mat(1, getCols());
    //     for (size_t j = 0; j < getCols(); ++j) {
    //         std::cout << (*this)(0, j) << std::endl;
    //     }
    //     return mat;
    // }

    // if (getCols() == 1) {  // 单列矩阵
    //     Matrix mat(getRows(), 1);
    //     for (size_t i = 0; i < getRows(); ++i) {
    //         double value = (*this)(i, 0);
    //         mat(i, 0) = value;
    //         std::cout << mat;
    //     }
    //     return mat;
    // }

    // 通用矩阵
    Matrix result(getRows(), getCols(), false);
    for (size_t i = 0; i < getRows(); ++i) {
        for (size_t j = 0; j < getCols(); ++j) {
            result(i, j) = (*this)(i, j);
        }
    }
    return result;
}
    template<typename OtherDerived> // 注意这里的取地址号非常重要😭
    AddExpr<Derived,OtherDerived> operator+(const MatrixView<OtherDerived>& other) const
    {
        return AddExpr<Derived,OtherDerived>(this->derived(),other.derived());
    }
    
    template<typename OtherDerived> // 
    MulExpr<Derived,OtherDerived> operator*(const MatrixView<OtherDerived>& other) const
    {
        return MulExpr<Derived,OtherDerived>(this->derived(),other.derived());
    }
    
}; // class

class RowView : public MatrixView<RowView>
{
private:
    const Matrix& m_source;
    size_t m_row;

public:
    RowView(const Matrix& source,size_t row) : m_source(source),m_row(row)
    {
        if(row >= source.getRows())
        {
            throw std::invalid_argument("row index is out of range");
        }
    }

    double operator()(size_t row,size_t col) const
    {
        if(col < 0 || col >= m_source.getCols())
        {
            throw std::invalid_argument("column index is out of range!");
        }
        return m_source(row+m_row,col);
    }
    
    size_t getRows() const {return 1;}
    size_t getCols() const {return m_source.getCols();}
}; // class

class ColView : public MatrixView<ColView>
{
private:
    const Matrix& m_source;
    size_t m_col;

public:
    ColView(const Matrix& source,size_t col) : m_source(source),m_col(col)
    {
        if(col >= source.getCols() || col < 0)
        {
            throw std::invalid_argument("column index is out of range");
        }
    }

    double operator()(size_t row,size_t col) const
    {
        if(row < 0 || row >= m_source.getRows())
        {
            throw std::invalid_argument("row index is out of range!");
        }

        return m_source(row,m_col+col);
    }

    size_t getRows() const {return m_source.getRows();}
    size_t getCols() const {return 1;}
}; // class

class DiagonalView : public MatrixView<DiagonalView>
{
private:
    const Matrix& m_source;

public:
    DiagonalView(const Matrix& source) : m_source(source) 
    {
        if(!source.isSquare())
        {
            throw std::invalid_argument("Unsquare matrix has no diagonal line!");
        }
    }

    double operator()(size_t row,size_t col) const
    {
        if(row != 0 || col > m_source.getCols())
        {
            throw std::invalid_argument("diagonal view has noly one row !");
        }
        return m_source(col,col);
    }

    size_t getRows() const {return 1;}
    size_t getCols() const {return m_source.getCols();}

}; // class

class ReplicateView : public MatrixView<ReplicateView>
{
private:
    const Matrix& m_source;
    size_t m_rowMuliples;
    size_t m_colMultiples;

public:
    ReplicateView(const Matrix& source,size_t rowMuliples,size_t colMultiples)
                 : m_source(source),m_rowMuliples(rowMuliples),
                   m_colMultiples(colMultiples)
    {
        if(rowMuliples == 0 || colMultiples == 0)
        {
            throw std::invalid_argument("Replication times must be positive numbers!");
        }
    }

    double operator()(size_t row,size_t col) const
    {
        return m_source((row % m_source.getRows()),(col % m_source.getCols()));
    }

    size_t getRows() const { return m_source.getRows() * m_rowMuliples; }
    size_t getCols() const { return m_source.getCols() * m_colMultiples; }
}; // class

class SubMatrixView : public MatrixView<SubMatrixView>
{
private:
    const Matrix& m_source;
    size_t m_start_row,m_start_col;
    size_t m_rows,m_col;

public:
    SubMatrixView(const Matrix& source,size_t start_row,size_t start_col,size_t end_row,size_t end_col)
                 : m_source(source),m_start_row(start_row),m_start_col(start_col),
                   m_rows(end_row-start_row),m_col(end_col-start_col)
    {
        if (start_row >= source.getRows() || end_row > source.getRows() ||
            start_col >= source.getCols() || end_col > source.getCols() ||
            start_row > end_row || start_col > end_col) 
        {
            throw std::out_of_range("Submatrix indices are out of original bounds!");
        }
    }

    size_t getRows() const {return m_rows;}
    size_t getCols() const {return m_col;}

    double operator()(size_t row,size_t col) const
    {
        return m_source(m_start_row + row,m_start_col + col);
    }

}; // class

//几种常用的视图
class TransposeView : public MatrixView<TransposeView>
{
private:
    const Matrix& m_source;

public:
    TransposeView(const Matrix& source) : m_source(source) {}
    
    size_t getRows() const { return m_source.getCols(); }
    size_t getCols() const { return m_source.getRows(); }

    double operator()(size_t row, size_t col) const 
    {  
        return m_source(col, row);
    }
};// class

class ScalaredView : public MatrixView<ScalaredView>
{
private:
    const Matrix& m_source;
    double m_scalar;

public:
    ScalaredView(const Matrix& source,double scalar) : m_source(source),m_scalar(scalar) {}

    double operator()(size_t row,size_t col) const
    {
        if(!(row>=0 && row < m_source.getRows() && col >=0 && col < m_source.getCols()))
        if(m_scalar < 1e-10) 
        {
            return 0.0;
        }

        return m_scalar * m_source(row,col);
    }

    size_t getRows() const { return m_source.getRows(); }
    size_t getCols() const { return m_source.getCols(); }
}; // class

//表达式模板(设计两种，一种是加法表达式，一种是乘法表达式)
template <typename LeftDerived,typename RightDerived>
class AddExpr : public MatrixView<AddExpr<LeftDerived, RightDerived>>
{
private:
    const LeftDerived m_left;
    const RightDerived m_right;

public:
    AddExpr(const LeftDerived left,const RightDerived right) : m_left(std::move(left)),m_right(std::move(right))
    {
        if(left.getRows() != right.getRows() || left.getCols() != right.getCols())
        {
            throw std::invalid_argument("Matrix sizes must match for addition!");
        }
    }

    double operator()(size_t row,size_t col) const
    {
        if (row >= m_left.getRows() || col >= m_left.getCols()) 
        {
        throw std::out_of_range("AddExpr index out of range");
        }
        return m_left(row,col) + m_right(row,col);
    }

    size_t getRows() const {return m_left.getRows();}
    size_t getCols() const {return m_left.getCols();}
};// class


template <typename LeftDerived,typename RightDerived>
class MulExpr : public MatrixView<MulExpr<LeftDerived, RightDerived>>
{
private:
    const LeftDerived& m_left;
    const RightDerived& m_right;

public:
    MulExpr(const LeftDerived& left,const RightDerived& right):
                 m_left(left),m_right(right)
    {
    }
    double operator()(size_t row, size_t col) const 
    {
        double sum = 0.0;
        for (size_t k = 0; k < m_left.getCols(); ++k) 
        {
            sum += m_left(row, k) * m_right(k, col);
        }
        return sum;
}

    size_t getRows() const {return m_left.getRows();}
    size_t getCols() const {return m_right.getCols();}
                
}; // class MulExpr


} //. namespace gomat
