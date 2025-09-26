#ifndef VECTOR_H
#define VECTOR_H

#include <cmath>
#include <iostream>
#include <vector>
#include <cstddef>
namespace gomat {
class Matrix;
struct Solution;
class VectorStream;

class Vector
{
private:
    std::vector<double> data;

public:
    Vector() = default; 
    Vector(const Vector& other) : data(other.data) {}
    Vector(size_t dimension): data(std::vector<double>(dimension)) {}
    Vector(size_t dimension,double value) : data(dimension,value) {}
    Vector(std::initializer_list<double> init) : data(init) {}  // 列表初始化 
    Vector(const std::vector<double>& vec) : data(vec)  {}
    Vector(std::vector<double>&& vec) noexcept : data(std::move(vec)) {}
    Vector& operator=(std::vector<double>&& vec) noexcept;
    void operator=(const Vector& other) {data = other.data;}
    
    size_t getDimension() const;
    
    // 基本运算符重载
    double& operator[](size_t index);
    double operator[](size_t index) const;
    Vector operator*(double scalar) const;
    Vector operator+(Vector other);
    Vector operator+(std::vector<double> other);
   
    // 特殊的 
    std::vector<double> to_std_vector() const;
    void push_back(double value);

    //赋值
    void reverseAssignFrom(const std::vector<double>& vec);
    void assignFrom(const std::vector<double>& vec);


    // 数学运算
    Vector crossProduct(Vector vec) const;
    double dotProduct(Vector vec) const;
    Vector normalize();
    double norm() const;

    // 基本特征
    bool isEmpty() const;

    //工厂方法
    static Vector linSpaced(size_t nums,double start,double end);
    static Vector aranged(double start, double end, double step);
    static Vector ones(size_t nums);
    static Vector zeros(size_t nums);

    // replicate方法
    Matrix replicate(size_t rowTimes,size_t colTimes) const;
    
    //流式赋值
    VectorStream operator<<(double value);

    //输出
    friend std::ostream& operator<<(std::ostream& out,const Vector& vec);
    friend std::ostream& operator<<(std::ostream& out,const Solution& sol);
    friend std::ostream& operator<<(std::ostream& out,size_t i);
};
} //namespace gomat
#endif