#include <iostream>
#include <cmath>
#include <stdexcept>
#include <iomanip>
#include "vector.h"
#include "matrix.h"
#include "utils.h"
#include "matrix_stream.h"

namespace gomat{

std::vector<std::string> splitCSVLines(const std::string& line)
{
    std::vector<std::string> elements;
    std::string current;
    bool flag = false; // 判断遍历字符是否在引号内 
    

    for (char c:line)
    {
        if (c == '"')
        {
            flag = !flag;
        }
        else if (c == ',' && !flag)
        {
            elements.push_back(current);
            current.clear();
        }
        else{
            current += c;
        }
    }
    elements.push_back(current);
    return elements;
}
}

namespace gomat {
void Matrix::loadFromCsv(const std::string& filename)
{
    size_t dotPos = filename.find_last_of('.');
    std::string posifix = filename.substr(dotPos+1);
    std::transform(posifix.begin(),posifix.end(),posifix.begin(),::tolower);
        
    if (posifix != "csv")
    {
        throw std::invalid_argument("NOT CSV FILE");
    }

    std::ifstream file(filename);
    if (!file.is_open())
    {
        throw std::runtime_error("File open failure"+filename);
    }

    std::vector<std::vector<double>> data; // 存储矩阵
    std::string line;
    size_t expectedCol = 0; // 每一行元素的个数应当是identical的

    while (std::getline(file,line))
    {
        if (line.empty()) // 允许空行的出现
        {
            continue;
        }

        std::vector<std::string> strElement = splitCSVLines(line); // 用于存储文件中某一行的数据
        std::vector<double> row; // 用于存住每一行转化后的数字

        //转换为数值
        for(const std::string& elem:strElement) // 遍历这一行所有的数据
        {
            try
            {
                double val = elem.empty()? 0.0 : std::stod(elem); // 将每一个数据都转化为数值
                row.push_back(val); // 压入存储数字的vector
            }
            catch(const std::exception& e)
            {
                std::cerr << e.what() << '\n';
            }

            if(data.empty()) // 第一次循环，矩阵为空，得到预期的列数
            {
                expectedCol = row.size();
            }
            else if (row.size() != expectedCol) // 如果存在某一行元素个数不等于第一行,那么报错
            {
                throw std::runtime_error("dismatched row sizes in csv file,not suitable for a matrix");
            }

            data.push_back(row); //将存储数字的vector压入矩阵
        }  
        m_cols = data[0].size();
        m_rows = data.size();
        m_mat = data;   
        
        if(m_rows * m_cols > 10000)
        {
            m_is_contiguous = true;
            this->toContiguous();
        }
    }
}

} //  namespace

namespace gomat {
void Matrix::matrixToCsv(const std::string& filename,int precision,char comma)
{
    size_t dotPos = filename.find_last_of('.'); // 检查文件拓展名是否是.csv
    std::string posifix = filename.substr(dotPos+1);
    std::transform(posifix.begin(),posifix.end(),posifix.begin(),::tolower);
        
    if (posifix != "csv")
    {
        throw std::invalid_argument("NOT CSV FILE");
    }

    if(m_cols == 0 || m_rows == 0)
    {
        throw std::runtime_error("Empty matrix cannot be written into csv file");
    }

    std::ofstream file(filename);
    if( !file.is_open())
    {
        throw std::runtime_error("fail to open" + filename);
    }

    file  << std::fixed << std::setprecision(precision);

    // 按行写入矩阵元素
    for(size_t i=0;i<m_rows;++i)
    {
        for (size_t j=0;j<m_cols;++j)
        {
            file << (*this)(i,j);
            if(j != m_cols-1)
            {
                file << comma;
            }
        }

        file << "\n";
    }
    file.close();

} 
}

namespace gomat{
std::istream& operator>> (std::istream& in, Matrix& mat)
{
        for (size_t i = 0; i < mat.m_rows; ++i) 
        {
            std::cout << "please input the element of row" << " " << i << '\n';
            for (size_t j = 0; j < mat.m_cols; ++j) 
            {
                std::cout << "element" << " " << j << ':';
                in >> mat(i,j);
            }
        }
        return in;
    }
} // namespace

namespace gomat{

std::ostream& operator<<(std::ostream& out, const Matrix& mat) {
    const int cell_width = 6;  // 每个元素的打印宽度
    const int last_cell_width = 1;  // 最后一个元素的宽度（不需要对齐）

    const int border_width = (mat.getCols() > 1) 
        ? cell_width * (mat.getCols() - 1) + last_cell_width + 2  
        : last_cell_width + 2;
    const std::string border(border_width, '-');

    out << "The matrix is:\n"
        << "+" << border << "+\n";

    for (size_t i = 0; i < mat.getRows(); ++i) 
    {
        out << "| ";
        for (size_t j = 0; j < mat.getCols(); ++j) 
        {
            if (j == mat.getCols() - 1) 
            {
                out << mat(i, j);
            } else 
            {
                out << std::setw(cell_width) << std::left << mat(i, j);
            }
        }
        out << " |\n";
    }

    out << "+" << border << "+\n";
    return out;
}
} // namespace 

namespace  gomat {
MatrixStream Matrix::operator<<(double value)
{
    MatrixStream matStream(*this);
    return matStream << value;
}

}