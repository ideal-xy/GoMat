#include "vector.h"
#include <stdexcept>
#include <cmath>
#include <ranges>
#include <algorithm>
#include "matrix.h"
#include "vector_stream.h"

using namespace gomat;
// 返回dimension
size_t Vector::getDimension() const 
{
    return data.size();
}

// 下标运算符重载（可修改版本）
double& Vector::operator[](size_t index) 
{
    if (index >= data.size()) {
        throw std::out_of_range("Vector index out of range");
    }
    return data[index];
}

// 下标运算符重载（const版本）
double Vector::operator[](size_t index) const 
{
    if (index >= data.size()) {
        throw std::out_of_range("Vector index out of range");
    }
    return data[index];
}

// 向量数乘
Vector Vector::operator*(double scalar) const 
{
    Vector result(data.size());
    for (size_t i = 0; i < data.size(); ++i) {
        result[i] = data[i] * scalar;
    }
    return result;
}

// 向量叉积（仅适用于3维向量）
Vector Vector::crossProduct(Vector vec) const 
{
    if (data.size() != 3 || vec.getDimension() != 3) {
        throw std::runtime_error("Cross product is only defined for 3D vectors");
    }
    
    return Vector{
        data[1] * vec[2] - data[2] * vec[1],
        data[2] * vec[0] - data[0] * vec[2],
        data[0] * vec[1] - data[1] * vec[0]
    };
}

// 向量点积
double Vector::dotProduct(Vector vec) const 
{
    if (data.size() != vec.getDimension()) {
        throw std::runtime_error("Vectors must have same dimension for dot product");
    }
    
    double result = 0.0;
    for (size_t i = 0; i < data.size(); ++i) {
        result += data[i] * vec[i];
    }
    return result;
}

// 向量归一化
Vector Vector::normalize() 
{
    double norm = 0.0;
    for (double val : data) {
        norm += val * val;
    }
    norm = std::sqrt(norm);
    
    if (norm == 0.0) {
        throw std::runtime_error("Cannot normalize zero vector");
    }
    
    for (double& val : data) {
        val /= norm;
    }
    return *this;
}

double Vector::norm() const
{
    double norm = 0.0;
    for (double val : data) {
        norm += val * val;
    }
    return std::sqrt(norm);  
}

Vector& Vector::operator=(std::vector<double>&& vec) noexcept
{
    data = std::move(vec);
    return *this;
}

void Vector::reverseAssignFrom(const std::vector<double>& vec)
{
    data.resize(vec.size());
    data.assign(vec.rbegin(),vec.rend());
}

void Vector::assignFrom(const std::vector<double>& vec)
{
    data.resize(vec.size());
    data = vec;
}

Vector Vector::operator+(Vector other)
{
    std::transform(data.begin(),data.end(),other.data.begin(),data.begin(),[](double a,double b){return a+b;});
    return *this;
}

Vector Vector::operator+(std::vector<double> other)
{
    std::transform(data.begin(),data.end(),other.begin(),data.begin(),[](double a,double b){return a+b;});
    return *this;
}

//基本特征
bool Vector::isEmpty() const
{
    return data.empty();
}

// 工厂方法
Vector Vector::linSpaced(size_t nums,double start,double end)
{
    Vector vec(nums);
    double step = (end-start) / (nums - 1);
    std::generate(
        vec.data.begin(),vec.data.end(),[&start,step](){
            double val = start;
            start += step;
            return val;
        }
    );

    return vec;
}

Vector Vector::aranged(double start, double end, double step)
{
    Vector vec(static_cast<size_t>((end-start)/step));
    std::generate(
        vec.data.begin(),vec.data.end(),[&start,step]{
            double val = start;
            start += step;
            return val;
        }
    );

    return vec;
}

Vector Vector::ones(size_t nums)
{
    Vector vec(nums);
    std::generate(vec.data.begin(),vec.data.end(),[](){return 1.0;});
    return vec;
}

Vector Vector::zeros(size_t nums)
{
    Vector vec(nums);
    std::generate(vec.data.begin(),vec.data.end(),[](){return 0.0;});
    return vec;
}


Matrix Vector::replicate(size_t rowTimes,size_t colTimes) const
{
    if(this->getDimension() > 1000)
    {
        // Matrix(colTimes,rowTimes * this->getDimension(),true);
        std::vector<double> vec;
        for(size_t k=0;k<colTimes;++k)
        {
            for(size_t i=0;i<rowTimes;++i)
            {
                std::copy(this->data.begin(),this->data.end(),std::back_inserter(vec));
            }
        }
        return Matrix(vec,colTimes,rowTimes * this->getDimension());
    }
    else
    {
        std::vector<std::vector<double>> vec;
        std::vector<double> row_vec;
        for(size_t i=0;i<rowTimes;++i)
        {
            std::copy(this->data.begin(),this->data.end(),std::back_inserter(row_vec));
        }

        for(size_t i=0;i<colTimes;++i)
        {
            vec.push_back(row_vec);
        }
        return Matrix(vec);
    }
}

//流式赋值
VectorStream Vector::operator<<(double value)
{
    VectorStream vec_stream(*this);
    return vec_stream << value;
}

// 输出
namespace gomat{
std::ostream& operator<<(std::ostream& out, const Vector& vec) 
{
    out << "[";
    for (size_t i = 0; i < vec.data.size(); ++i) {
        out << vec.data[i];
        if (i != vec.data.size() - 1) {
            out << ", ";
        }
    }
    out << "]";
    return out;
}

std::ostream& operator<<(std::ostream& out,size_t i)
{
    out << static_cast<int>(i);
    return out;
}
}

std::vector<double> Vector::to_std_vector() const
{
    return this->data;
}

void Vector::push_back(double value)
{
    data.push_back(value);
}