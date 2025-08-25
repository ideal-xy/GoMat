#include "utils.h"
#include <iostream>
#include <algorithm>
#include "vector.h"

using namespace gomat;
double Utils::vecDotProduct(std::vector<double> vec1,std::vector<double> vec2)
{
    if(!checkVectorSize(vec1,vec2)) // 检查两个向量维数是否相等
    {
        throw std::invalid_argument("Vectors of different sizes cannot have an inner product");
    }

    double result = std::inner_product(vec1.begin(),vec1.end(),vec2.begin(),0.0); // 我们使用stl来计算这个内积，需要导入numeric库

    return result;

}

bool Utils::checkVectorSize(std::vector<double> vec1,std::vector<double> vec2)
{
    return (vec1.size() == vec2.size());
}


int Utils::gcd(size_t a,size_t b)
{
    int x = static_cast<int>(a);
    int y = static_cast<int>(b);
    
    if(x < y)
    {
        std::swap(x,y);
    }
    if( y == 0)
    {
        return x;
    }
    else
    {
        return gcd(y,x%y);
    }
    
}

bool Utils::isZero(double value)
{
    return std::abs(value) < 1e-8; 
}

std::vector<double> Utils::mulyiply(double scalar,std::vector<double> vec)
{
    std::transform(vec.begin(),vec.end(),vec.begin(),[scalar](double ele){return ele * scalar;});
    return vec;
}


double Utils::norm_of_stdVector(std::vector<double> vec)
{
    double total = 0.0;
    for(double ele:vec)
    {
        total += ele * ele;
    }

    return std::sqrt(total);
}

std::vector<double> Utils::normalize(std::vector<double> vec)
{
    std::transform(vec.begin(),vec.end(),vec.begin(),[vec](double ele){return ele / norm_of_stdVector(vec);});
    return vec;
}


std::vector<double> Utils::frontVectorMinusRearVector(std::vector<double> vec1,std::vector<double> vec2)
{
    if(vec1.size() != vec2.size())
    {
        throw std::invalid_argument("diamatched vcetor dimensions");
    }

    std::transform(vec1.begin(),vec1.end(),vec2.begin(),vec1.begin(),[](double front,double rear){return front - rear;});
    
    return vec1;


}

// householder 辅助函数

Vector Utils::householder_normal_vector(const Vector& vec) // 这个函数是求解 把向量vec变为单位向量 [1,0,0,0,0,0...0]的超平面法向量
{
    if(vec.isEmpty())
    {
        throw std::invalid_argument("empty vector cannot do householder transformation");
    }

    size_t m = vec.getDimension();
    Vector normal_vec = vec;
    
    double sum = 0.0;
    double c = 0.0;
    for(size_t i=1;i<m;++i)
    {
        double y = std::pow(vec[i],2) - c;
        double t = sum + y;
        c = (t-sum) - y;
        sum = t;
    }

    if (sum < 1e-15 * std::fabs(vec[0])) // 如果后面几项的平方和相对于第一项太小的话，返回零向量
    {
        return Vector(m,0.0);
    }

    double vec_norm = std::sqrt(vec[0] * vec[0] + sum);

    double v0 = vec[0] + std::copysign(vec_norm,vec[0]);
    normal_vec[0] = v0;

    double normal_vec_norm = std::sqrt(v0 * v0 + sum);
    if(normal_vec_norm < 1e-15)
    {
        throw std::runtime_error("calculation failure,the norm of normal vecetor limits to 0");
    }

    double factor = 1.0 / normal_vec_norm;
    for(size_t i=0;i<m;++i)
    {
        normal_vec[i] *= factor;
    }

    return normal_vec;
}

void Utils::householder_left_mulyiply(Matrix& A,Matrix& U,const Vector& vec) // 这一步实现 H * U = (I-2v*v^{T}) * U = U - 2 v*(v^{T}*U)
{
    size_t m = A.getRows();
    size_t n = A.getCols();

    std::vector<double> temp(n,0.0);
    for(size_t j=0;j<n;++j) // 计算 v^{T}*A 
    {
        for(size_t i=0;i<m;++i)
        {
            temp[j] += vec[i] * A(i,j); 
        }
    }

    for(size_t i=0;i<m;++i) //  更新矩阵A 
    {
        for(size_t j=0;j<n;++j)
        {
            A(i,j) -= 2 * vec[i] * temp[j];
        }
    }

    for(size_t col=0;col<m;++col)
    {
        double dot = 0.0;
        for(size_t i=0;i<m;++i)
        {
            dot += U(i,col) * vec[i]; // 计算 v^{T} * U 的第col列处的元素
        }
        for(size_t i=0;i<m;++i)
        {
            U(i,col) -= 2 * dot * vec[i];  
        }
    }  
}

void Utils::householder_right_multiply(Matrix& A,Matrix& V,const Vector& vec)
{
    size_t m = A.getRows();
    size_t n = A.getCols();

    std::vector<double> temp(m,0.0);
    for(size_t i=0;i<m;++i)
    {
        for(size_t j=0;j<n;++j)
        {
            temp[i] += A(i,j) * vec[j];
        }
    }

    for(size_t i=0;i<m;++i) // 更新矩阵A
    {
        for(size_t j=0;j<n;++j)
        {
            A(i,j) -= 2 * temp[i] * vec[j];
        }
    }

    for(size_t i=0;i<n;++i) // 更新矩阵V
    {
        double dot = 0.0;
        for(size_t j=0;j<n;++j)
        {
            dot += vec[j] + V(j,i);
        }

        for(size_t j=0;j<n;++j)
        {
            V(i,j) -= 2 * vec[i] * dot;
        }
    }  
}


void Utils::bidiagonalize(Matrix& A,Matrix& U,Matrix& V)
{
    size_t m = A.getRows();
    size_t n = A.getCols();
    size_t min_dim = std::min(m,n);
    bool flag1 = U.isSquare();
    bool flag2 = V.isSquare();
    bool flag3 = (U.getRows() == m);
    bool flag4 = (V.getRows() == n);

    if(!(flag1 && flag2 && flag3 && flag4))
    {
        throw std::invalid_argument("wrong parameters");
    }

    std::cout << "0000" << std::endl;

    U.fillScaledIdentity(1.0);
    V.fillScaledIdentity(1.0);
    
    for(size_t k=0;k<min_dim;++k) 
    {
        Vector col_vec;// 处理第k列，从第k列索引为k的元素开始截断，把他变为[1,0,...,0]
        for(size_t i=k;i<m;++i)
        {
            col_vec.push_back(A(i,k));
        }
        
        Vector normal_vec = householder_normal_vector(col_vec);
        Vector m_dim_vec(m,0.0);

        for(size_t i=0;i<normal_vec.getDimension();++i)
        {
            m_dim_vec[k+i] = normal_vec[i];
        }
        // 左乘变换并更新U
        householder_left_mulyiply(A,U,m_dim_vec);
        std::cout << "5555" << std::endl;

        // 开始处理第k行，但是最后一行除外；
        if(k == n-1)
        {
            break; // 其实continue也可以
        }

        Vector row_vec;
        for(size_t i=k+1;i<n;++i) //开始处理第k行，从k+1开始截取
        {
            row_vec.push_back(A(k,i));
        }

        Vector normal_vec2 = householder_normal_vector(row_vec);

        Vector n_dim_vec(n,0.0);

        for(size_t j=0;j<normal_vec2.getDimension();++j)
        {
            n_dim_vec[k+1+j] = normal_vec2[j];
        }
    
        householder_right_multiply(A,V,n_dim_vec);    
    }
}


// GIvens旋转实现QR迭代辅助函数

void Utils::compute_givens(double a,double b,double& cos_theta,double& sin_theta)
{
    double norm = std::sqrt(a * a + b * b);
    cos_theta = a / norm;
    sin_theta = std::copysign(b / norm,-1.0);
}

void Utils::apply_givens_left(Matrix& A,size_t i,size_t j,double cos_theta,double sin_theta)
{
    for(size_t col=0;col<A.getCols();++col)
    {
        double ai = A(i,col),aj = A(j,col); 

        A(i,col) = cos_theta * ai - sin_theta * aj;
        A(j,col) = sin_theta * ai + cos_theta * aj;
    }
}

void Utils::apply_givens_right(Matrix& A,size_t i,size_t j,double cos_theta,double sin_theta)
{
    for(size_t row=0;row<A.getRows();++row)
    {
        double ai = A(row,i),aj = A(row,j);

        A(row,i) = ai * cos_theta + aj * sin_theta;
        A(row,j) = aj * cos_theta - ai * sin_theta; 
    }
}

double Utils::compute_wilkinson_shift(double a,double b,double c,double d)
{
   return( c - std::copysign(((b*b) / (std::copysign(d,1) + std::sqrt(d*d+b*b))),d));
}

void Utils::QrIteration(Matrix& B,Matrix& U,Matrix& V,double epsilon,int max_iteration)
{
    size_t n = B.getRows();
    for(int iter=0;iter<max_iteration;++iter)
    {
        bool converged = true;
        for(size_t i=0;i<n;++i)
        {
            if((std::abs(B(i,i+1)) > (epsilon * (std::abs(B(i,i)) + std::abs(B(i+1,i+1))))))
            {
                converged = false;
                break;
            }
        }
        if(converged)
        {
            break;
        }

        double mu = compute_wilkinson_shift(B(n-2,n-2),B(n-2,n-1),B(n-1,n-2),B(n-1,n-1));

        for(size_t i=0;i<n-1;++i)
        {
            double cos_theta,sin_theta;
            compute_givens(B(i,i)-mu,B(i,i+1),cos_theta,sin_theta);
            apply_givens_right(B,i,i+1,cos_theta,sin_theta);
            apply_givens_left(V,i,i+1,cos_theta,sin_theta);

            compute_givens(B(i,i)-mu,B(i+1,i),cos_theta,sin_theta);
            apply_givens_left(B,i,i+1,cos_theta,sin_theta);
            apply_givens_right(U,i,i+1,cos_theta,sin_theta);
        }
    }
}

void Utils::checkMultiplicationDims(const Matrix& a, const Matrix& b)
{
    if(a.getCols() != b.getRows())
    {
        throw std::invalid_argument("dismatched matrix cannot have multiplication");
    }
}

