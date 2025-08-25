#include <iostream>
#include "gomat.h"

int main() {
    
    gomat::Matrix mat(3, 3);
    
    mat << 1, 2, 3, 
           4, 5, 6,  
           1, 3, 6;

    auto view = mat.colView(1);
    auto view2 = mat.rowView(2);
    auto view3 = mat.diagonalView();
    auto view4 = mat.replicateView(1,1);
    auto view5 = mat.scalaredView(2);
    auto view7 = mat.subMatrixView(0,3,0,3);
    auto view6 = view4 * view5 * view4 + view4;

    std::cout << view7.eval() << std::endl;
    

    std::cout << gomat::version();

    return 0;
}