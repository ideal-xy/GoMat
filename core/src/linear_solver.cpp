#include "linear_solver.h"
#include "vector.h"
#include "matrix.h"
#include <iostream>

namespace gomat{

Solution Linear_Solver::gausslian(Matrix& A_mat,Vector b_vec) const
{
    Matrix augmented = (A_mat,b_vec);
    size_t rank1 = static_cast<size_t>(augmented.rank());
    size_t rank2 = static_cast<size_t>(A_mat.rank());

    Solution solution;

    if(rank1 != rank2)
    {
        solution.type = SolutionType::NoSolution;
        solution.output_message += "This equation has no solution";
    }

    if (rank1 == A_mat.getCols() && rank2 == A_mat.getCols())
    {
        solution.type = SolutionType::UniqueSolution;
        solution.output_message += "This equation has only one solution";
        double x;
        std::vector<double> vec_of_x;
        vec_of_x.resize(0);
        x = augmented(rank1-1,augmented.getCols()-1) / augmented(rank1-1,augmented.getCols()-2);
        vec_of_x.push_back(x);

        for(size_t row = static_cast<size_t>(rank1-2);row >= 0;--row)
        {
            double total = 0.0;
            for(size_t i=0;i<vec_of_x.size();++i)
            {
                total += vec_of_x[vec_of_x.size()-i-1] * augmented(row,augmented.getCols()-1-vec_of_x.size()+i);
            }
            
            x = (augmented(row,augmented.getCols()-1) - total) / augmented(row,augmented.getCols()-1-vec_of_x.size()-1);
            vec_of_x.push_back(x);        
        }

        solution.unique_solution.reverseAssignFrom(vec_of_x);
    }

    if(rank1 == rank2 && rank1 < A_mat.getCols())
    {
        solution.type = SolutionType::InfiniteSolutions;
        solution.output_message += "This equation have infinite solutions";
        double x;
        
        size_t freeVaribles = A_mat.getCols() - rank1;
        std::vector<double> vec_of_x;
        vec_of_x.resize(0);
        for(size_t i=0;i<freeVaribles;++i)
        {
            vec_of_x.push_back(0.0);
        }

        for (size_t row=rank1;row>=0;--row)
        {
            double total = 0.0;
            for (size_t i=0;i<vec_of_x.size();++i)
            {
                total += augmented(row,augmented.getCols()-1-vec_of_x.size()+i) * vec_of_x[vec_of_x.size()-1-i];
            }
            x = (augmented(row,augmented.getCols()-1) -total) / augmented(row,augmented.getCols()-2-vec_of_x.size());
            vec_of_x.push_back(x);
        }

        solution.particular_solution.reverseAssignFrom(vec_of_x); 
        //接下来求解基础解系 

        double total = 0.0;
        size_t freeVaribles2 = A_mat.getCols()-rank1;
        std::vector<std::vector<double>> vec;
        vec.resize(freeVaribles2);
        for(size_t i=0;i<freeVaribles2;++i)
        {
            vec[i].resize(freeVaribles2);
            vec[i][i] = 1.0;
        }

        for(size_t i = 0;i < vec.size();++i) // 遍历基础解系的每一个尚未填满的解向量
        {
            for (size_t row = rank1;row >= 0;--row) // 填满基础解系中的每一个尚未填满的解向量，
            {
                for(size_t j = 0;j < vec[i].size();++j) 
                {
                    total += augmented(row,augmented.getCols()-1-vec[i].size()) * vec[i][vec[i].size()-1-j];
                }

                x = ((-1) * total) /  augmented(row,augmented.getCols()-2-vec[i].size());
                vec[i].push_back(x);
            }
        }

        solution.basic_solution_system.resize(freeVaribles2); // basic_solution_system 的type为std::vector<std::vector<double>>

        for(size_t i=0;i<freeVaribles2;++i)
        {
            solution.basic_solution_system[i].resize(vec[0].size());
            std::reverse_copy(vec[i].begin(),vec[i].end(),solution.basic_solution_system[i].begin());
        }

        solution.nums = vec.size();
    }
    return solution;

}

std::ostream& operator<<(std::ostream& out,Solution solution)
{
    out << "The result is as follows:" << '\n';
    out << "--------------------------" << '\n';
    out << "The output message is:" << '\n';
    out << solution.output_message << '\n';
    out << "--------------------------" << '\n';

    if (solution.type == SolutionType::NoSolution)
    {
        out<< "This equation has no solution" << '\n';
    }

    if (solution.type == SolutionType::InfiniteSolutions)
    {
        out <<"This is equation has infinite solutions" <<'\n';
        out << "---------------------------------------" << '\n';
        out << solution.particular_solution << " + "<< '\n';
    }

    if (solution.type == SolutionType::UniqueSolution)
    {
        out << "This is equation has infinite solutions" << '\n';
        out << "---------------------------------------" << '\n';
        out <<solution.unique_solution <<'\n';
        for(size_t i=0;i<solution.nums;++i)
        {
            out << "k_"<< i << " * " << "[";
            for(size_t j=0;j<solution.basic_solution_system[0].size();++j)
            {
                out << solution.basic_solution_system[i][j];
                if(j !=solution.basic_solution_system[0].size()-1)
                {
                    out << ",";
                }
            }
            out << "]" << "+" <<'\n';
        }
        out << ",among which k_1,k_2,... is any real number ";
    }

}
}