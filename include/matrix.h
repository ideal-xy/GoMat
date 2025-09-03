#ifndef MATRIX_H
#define MATRIX_H

#include <vector>
#include <utility>
#include <cstddef>
#include <numeric>
#include <string>
#include <fstream>
#include <tuple>
#include "vector.h"

namespace gomat {
    
template <typename Derived> class MatrixView;
class TransposeView;
class SubMatrixView;
class ReplicateView;
class DiagonalView;
class RowView;
class ColView;
class ScalaredView;
class Vector;
class MatrixStream;

class Matrix
{
private:
    std::vector<std::vector<double>> m_mat;
    std::vector<double> m_data;
    size_t m_rows;
    size_t m_cols;
    bool m_is_contiguous; // 是否使用连续存储以提高大矩阵运算性能

    double& getElement(size_t row,size_t col);
    double getEle(size_t row,size_t col) const;

public:
    
    // 构造函数和析构函数
    Matrix() {};
    Matrix(size_t rows,size_t cols,bool use_contiguous = false);
    // Matrix(std::initializer_list<double> list,size_t rows,size_t cols); // 支持列表初始化
    Matrix(const std::vector<std::vector<double>>& data);
    Matrix(const std::vector<double>& data,size_t rows,size_t cols);
    // Matrix(const std::vector<double>,size_t row,size_t col);
    Matrix(Matrix&& other) noexcept;
    Matrix(const Matrix& other);
    Matrix(std::vector<std::vector<double>>&& other) noexcept;
    ~Matrix() = default;

    // 两种存储形式的相互转化，应付性能不敏感的操作
    void toContiguous(); // 不连续存储转化为连续存储
    void toNonContiguous(); // 连续存储转化为不连续存储
    void ensureContiguousIfLarge(); // 设置存储类型自动检查+转换

    //尺寸变换
    void resize(size_t row,size_t col);

    // 基本特征
    size_t getRows() const;
    size_t getCols() const ;
    bool isEmpty(double epsilon) const;
    bool isSquare() const ;
    bool isDiagonal(double epsilon) const;
    bool isUpperTriangular(double epsilon) const;
    bool isLowerTriangular(double epsilon) const;
    bool isLarge() const;
    bool isContiguous() const;

    //复杂特征
    double determiant() const ; // 行列式
    bool isInvertible() const; // 可逆否
    double trace() const; // 迹是多少
    int rank() const; // 秩是多少
    double sum() const; // 所有元素的和
    double product() const; // 所有元素的乘积
    double frobeniusNorm() const; // Frobenius范数：(所有元素的平方和)^{1/2}

    // 复杂操作
    std::tuple<int,Matrix> gaussElimination() const;// 返回值中的int times 代表的是进行第三类初等变换的次数，也就是交换两行的次数
    std::tuple<Matrix,Matrix> QrDecompostion() const;
    std::tuple<Matrix,Matrix>LuDecomposition() const;
    std::tuple<Matrix,Matrix,Matrix>SvdDecomposition() const; 
    Matrix inverse() const;

    // 元素访问
    double& operator()(size_t row,size_t col); // 允许修改
    double operator()(size_t row,size_t col) const; //不允许修改

    // generally, 分为大小两种矩阵乘法
    Matrix multiplyLarge(const Matrix& other) const;
    Matrix multiplySmall(const Matrix& other) const;

    // concretely,不同的矩阵乘法
    Matrix multiplyDense(const Matrix& other) const;
    Matrix multiplyDiagonal(const Matrix& other) const;
    Matrix multiplyUpperTriangle(const Matrix& other) const;
    Matrix multiplyWithSimd(const Matrix& other) const;
    
    /*
    multiplyBySimd()是使用SIMD指令集进行矩阵乘法计算，在对矩阵进行4 * 4的分块之后，把每一块作为参数
    传递给核函数，核函数内部进行向量与向量的运算
    */
   
    //基本的算数运运算符重载
    Matrix operator+(const Matrix& other) const;
    Matrix operator-(const Matrix& other) const;
    Matrix operator*(const Matrix& mat) const;
    Matrix operator*=(double scalar);
    bool operator==(const Matrix& other) const;
    bool operator!=(const Matrix& other) const;
    Matrix operator,(const Vector& vec) const;
    Matrix operator,(const Matrix& mat) const;
    Matrix& operator=(Matrix&& other) noexcept; // 移动赋值构造函数
   
    // 这个函数单单是为了写矩阵乘法用的,只读不修改
    const double& at(size_t row, size_t col) const 
    {
        return m_data[col * m_cols + row];
    }

    //三类初等变换
    void rowInterchange(size_t i,size_t j);
    void rowMultiply(size_t row,double factor);
    void rowAddition(size_t row1,size_t row2,double factor);

    // 原地转置方阵（修改this）
    void transposeInPlace();
    Matrix transposeByBlock(size_t blockSize) const;
    Matrix transpose() const; // general transpose 适合不是很大的矩阵
     
    // 特殊填充方式
    void fillWithData(const std::vector<std::vector<double>>& data);
    void fillWithOneValue(double value);
    void fillScaledIdentity(double k);
    // void fillRandomValue();

    //特殊矩阵（静态方法）
    // static Matrix ones(size_t row,size_t cols);
    // static Matrix zeros(size_t row,size_t cols);

    // 小矩阵多次复制变为大矩阵
    Matrix replicate(size_t rowTimes,size_t colTimes) const;

    // 辅助操作
    void erase();   //清空矩阵
    Matrix subMatrix(size_t start_row,size_t end_row,size_t start_col,size_t end_col) const; // 取子矩阵
    double blockSum(size_t start_row,size_t end_row,size_t start_col,size_t end_col) const ;  // 计算某一个子矩阵的元素之和


    //导出矩阵到csv文件
    void matrixToCsv(const std::string& filename,int precison = 6,char comma = ',');
    void loadFromCsv(const std::string& filename);

    //流输入
    MatrixStream operator<<(double value);

    //辅助操作
    friend std::istream& operator>> (std::istream& in, Matrix& mat);
    friend std::ostream& operator<< (std::ostream& out,const Matrix& mat);
    

    //返回视图
    TransposeView transposeView() const;
    SubMatrixView subMatrixView(size_t start_row, size_t end_row, size_t start_col, size_t end_col) const;
    ReplicateView replicateView(size_t rowTimes, size_t colTimes) const;
    DiagonalView diagonalView() const;
    RowView rowView(size_t row) const;
    ColView colView(size_t col) const;
    ScalaredView scalaredView(double scalar) const;
};
} // namespace
#endif