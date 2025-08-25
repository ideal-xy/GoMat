#include "matrix_view.h"

namespace gomat{

TransposeView Matrix::transposeView() const
{
    return TransposeView(*this);
}

SubMatrixView Matrix::subMatrixView(size_t start_row, size_t end_row, size_t start_col, size_t end_col) const
{
    return SubMatrixView(*this,start_row,start_col,end_row,end_col);
}

ReplicateView Matrix::replicateView(size_t rowTimes, size_t colTimes) const
{
    return ReplicateView(*this,rowTimes,colTimes);
}

DiagonalView Matrix::diagonalView() const
{
    return DiagonalView(*this);
}

RowView Matrix::rowView(size_t row) const
{
    return RowView(*this,row);
}

ColView Matrix::colView(size_t col) const
{
    return ColView(*this,col);
}

ScalaredView Matrix::scalaredView(double scalar) const
{
    return ScalaredView(*this,scalar);
}

} //namespace
