#ifndef MATRIX_H
#define MATRIX_H

#include <vector>
#include <cstddef>
#include <string>
#include <tuple>
#include <stdexcept>
#include <algorithm>
#include <ostream>
#include "mem_alloc.h"
#ifdef __x86_64__
#include <immintrin.h>
#elif defined(__aarch64__)
#include <arm_neon.h>
#endif
namespace gomat {
    

enum class Layout { ColMajor, RowMajor };

// 前置声明供 MatrixT 使用视图与流类型
template <typename Derived> class MatrixView;
template<typename T, Layout L, typename Alloc> class TransposeView;
template<typename T, Layout L, typename Alloc> class SubMatrixView;
template<typename T, Layout L, typename Alloc> class ReplicateView;
template<typename T, Layout L, typename Alloc> class DiagonalView;
template<typename T, Layout L, typename Alloc> class RowView;
template<typename T, Layout L, typename Alloc> class ColView;
template<typename T, Layout L, typename Alloc> class ScalaredView;
class Vector;
template<typename T, Layout L, typename Alloc> class MatrixStream;

// 预声明模板：将主模板命名为 Matrix，并为 T 提供默认为 double
template<typename T = double, Layout L = Layout::ColMajor,
         typename Alloc = SpecialAllocator<T, 32>> class Matrix;

template<typename T, Layout L,
         typename Alloc>
class Matrix {
private:
    std::vector<T, Alloc> m_data;
    size_t m_rows{0};
    size_t m_cols{0};
public:
    Matrix() = default;
    Matrix(size_t rows, size_t cols) : m_data(rows * cols), m_rows(rows), m_cols(cols) {}
    Matrix(size_t rows, size_t cols, bool /*use_contiguous*/) : m_data(rows * cols), m_rows(rows), m_cols(cols) {}

    // 兼容旧接口：从二维数组构造（按列主序填充）
    Matrix(const std::vector<std::vector<T>>& data)
    {
        m_rows = data.size();
        m_cols = m_rows ? data[0].size() : 0;
        m_data.resize(m_rows * m_cols);
        for (size_t j = 0; j < m_cols; ++j)
            for (size_t i = 0; i < m_rows; ++i)
                m_data[j * m_rows + i] = static_cast<T>(data[i][j]);
    }

    Matrix(const std::vector<T>& data, size_t rows, size_t cols)
    {
        if (data.size() != rows * cols) throw std::invalid_argument("dismatched dimensions");
        m_rows = rows; m_cols = cols; m_data.resize(rows * cols);
        std::copy(data.begin(), data.end(), m_data.begin());
    }

    inline T* data() const { return m_data.data();}

    inline size_t getRows() const { return m_rows; }
    inline size_t getCols() const { return m_cols; }

    inline T& operator()(size_t row, size_t col) {
        return (L == Layout::ColMajor) ? m_data[col * m_rows + row] : m_data[row * m_cols + col];
    }
    inline const T& operator()(size_t row, size_t col) const {
        return (L == Layout::ColMajor) ? m_data[col * m_rows + row] : m_data[row * m_cols + col];
    }

    inline T* data() { return m_data.data(); }
    inline const T* data() const { return m_data.data(); }

    inline size_t leadingDimensionColMajor() const { return m_rows; }
    inline size_t leadingDimensionRowMajor() const { return m_cols; }

    void resize(size_t rows, size_t cols) 
    {
        m_rows = rows; m_cols = cols; m_data.resize(rows * cols);
    }

    // 兼容旧接口（统一连续存储后为 no-op）
    inline void toContiguous() {}
    inline void toNonContiguous() {}
    inline void ensureContiguousIfLarge() {}

    inline const T& at(size_t row, size_t col) const { return m_data[col * m_rows + row]; }

    // 基本特征
    inline bool isSquare() const { return m_rows == m_cols; }
    inline bool isContiguous() const { return true; }
    inline bool isLarge() const { return (m_rows >= 64 || m_cols >= 64); }

    bool isEmpty(double epsilon) const;

    bool isDiagonal(double epsilon) const;

    bool isUpperTriangular(double epsilon) const;

    bool isLowerTriangular(double epsilon) const;

    // 三类初等行变换
    void rowInterchange(size_t i, size_t j);

    void rowMultiply(size_t row, double factor);

    void rowAddition(size_t row1, size_t row2, double factor);

    // 变换
    void transposeInPlace();

    Matrix transposeByBlock(size_t blockSize) const;

    Matrix transpose() const;

    // 填充与构造
    void fillWithData(const std::vector<std::vector<T>>& data);

    void fillWithOneValue(T value);

    void fillScaledIdentity(T k);

    void random(T lower_bound, T upper_bound);

    // 辅助
    void erase() { m_data.clear(); m_rows = 0; m_cols = 0; }

    Matrix subMatrix(size_t start_row, size_t end_row, size_t start_col, size_t end_col) const;

    T blockSum(size_t start_row, size_t end_row, size_t start_col, size_t end_col) const;

    // 复制扩展
    Matrix replicate(size_t rowTimes, size_t colTimes) const;

    // 算术
    Matrix operator+(const Matrix& other) const;

    Matrix operator-(const Matrix& other) const;

    bool operator==(const Matrix& other) const;

    bool operator!=(const Matrix& other) const;

    // 乘法
    Matrix multiplyDense(const Matrix& other) const;

    Matrix multiplyDiagonal(const Matrix& other) const;

    Matrix multiplyUpperTriangle(const Matrix& other) const;

    Matrix multiplyWithMulThread(const Matrix& other) const;
    Matrix multiplyLarge(const Matrix& other) const;
    Matrix multiplySmall(const Matrix& other) const;

    Matrix operator*(const Matrix& other) const;

    Matrix operator*=(double scalar);

    // 拼接
    Matrix operator,(const Matrix& mat) const;

    Matrix operator,(const Vector& vec) const;

    // 复杂特征与操作
    std::tuple<int, Matrix> gaussElimination() const;

    std::tuple<Matrix, Matrix> LuDecomposition() const;

    std::tuple<Matrix, Matrix> QrDecompostion() const;

    Matrix inverse() const;

    // 复杂特征
    T determiant() const;

    bool isInvertible() const;

    T trace() const;

    int rank() const;

    T sum() const;

    T product() const;

    T frobeniusNorm() const;

    // IO
    void matrixToCsv(const std::string& filename, int precison = 6, char comma = ',');

    void loadFromCsv(const std::string& filename);

    // 流输入（返回一个 MatrixStream 以便连续写入）——模板通用实现
    MatrixStream<T, L, Alloc> operator<<(double value);

    // 视图（模板化，支持所有 Matrix<T, L, Alloc>）
    TransposeView<T, L, Alloc> transposeView() const;
    SubMatrixView<T, L, Alloc> subMatrixView(size_t start_row, size_t end_row, size_t start_col, size_t end_col) const;
    ReplicateView<T, L, Alloc> replicateView(size_t rowTimes, size_t colTimes) const;
    DiagonalView<T, L, Alloc> diagonalView() const;
    RowView<T, L, Alloc> rowView(size_t row) const;
    ColView<T, L, Alloc> colView(size_t col) const;
    ScalaredView<T, L, Alloc> scalaredView(double scalar) const;
};

// 输出到标准流（声明）
template<typename T, Layout L, typename Alloc>
inline std::ostream& operator<<(std::ostream& out, const Matrix<T, L, Alloc>& mat)
{
    const size_t rows = mat.getRows();
    const size_t cols = mat.getCols();
    out << "[";
    for (size_t i = 0; i < rows; ++i) {
        if (i > 0) out << "\n ";
        for (size_t j = 0; j < cols; ++j) {
            if (j > 0) out << ", ";
            out << mat(i, j);
        }
    }
    out << "]";
    return out;
}
// Legacy Matrix class removed completely to meet new requirements.
} // namespace
#endif