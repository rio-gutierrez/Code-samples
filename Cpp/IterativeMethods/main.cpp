//  main.cpp
//  Iterative Methods (Jacobi, Gauss-Seidel, SOR)
//
//  Created by Mario L Gutierrez on 1/17/21.


#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>

#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <Eigen/SparseQR>
#include <Eigen/QR>
#include <Eigen/Householder>

using namespace std;
using namespace Eigen;

VectorXd Jacobi(const MatrixXd &A, const VectorXd &b, VectorXd &x0, const int dim = 12, const int max_it = 10){
    
    for (int k {0}; k < max_it; k++){
        
        VectorXd vect = VectorXd::Zero(dim);
        VectorXd sumvect = VectorXd::Zero(dim);
        
        for (int i {0}; i < dim; i++){
            for (int j {0}; j < dim; j++){
                if (i !=j){
                vect(i) = A(i,j) * x0(j);
                sumvect(i) += vect(i);
                }
            }
        }
        
        for (int i {0}; i < dim; i++)
            x0(i) = 1.0/A(i,i) * ( b(i) - sumvect(i) );
    }
    return x0;
}



VectorXd GaussSeidel(const MatrixXd &A, const VectorXd &b, VectorXd &x0, const int dim = 12, const int max_it = 10){
    
    for (int k {0}; k < max_it; k++){
        
        VectorXd vect = VectorXd::Zero(dim);
        VectorXd sumvect = VectorXd::Zero(dim);
        VectorXd currentsumvect = VectorXd::Zero(dim);
        
        for (int i {0}; i < dim; i++){
            for (int j {i+1}; j < dim; j++){
                vect(i) = A(i,j) * x0(j);
                sumvect(i) += vect(i);
            }
        }
        
        for (int i {0}; i < dim; i++){
            for (int j {0}; j < i; j++){
                vect(i) = A(i,j) * x0(j);
                currentsumvect(i) += vect(i);
            }
                x0(i) = 1.0/A(i,i) * ( b(i) - currentsumvect(i) - sumvect(i) );
        }
    }
    return x0;
}


VectorXd SOR(const MatrixXd &A, const VectorXd &b, VectorXd &x0, const double omega = 1.1, const int dim = 12, const int max_it = 10){
    
    for (int k {0}; k < max_it; k++){
        
        VectorXd vect = VectorXd::Zero(dim);
        VectorXd sumvect = VectorXd::Zero(dim);
        VectorXd currentsumvect = VectorXd::Zero(dim);
        
        for (int i {0}; i < dim; i++){
            for (int j {i+1}; j < dim; j++){
                vect(i) = A(i,j) * x0(j);
                sumvect(i) += vect(i);
            }
        }
        
        for (int i {0}; i < dim; i++){
            for (int j {0}; j < i; j++){
                vect(i) = A(i,j) * x0(j);
                currentsumvect(i) += vect(i);
            }
                x0(i) = omega/A(i,i) * ( b(i) - currentsumvect(i) - sumvect(i) ) + (1.0 - omega) * x0(i) ;
        }
    }
    return x0;
}



int main(int argc, const char * argv[]) {

    const int n {12};
    
    MatrixXd A(n,n);
    VectorXd b(n);

    //Generate matrix A
    for (int i {0}; i < n; i++) {
        for (int j {0}; j < n; j++){
            if (i == j)
                A(i,j) = 3.0;
            else if ( (i == j+1) || (j == i+1))
               A(i,j) = -1.0;
            else if ( (j == n-i-1) && (  (i != n/2 - 1) || (i != n/2)    )   )
                A(i,j) = 1.0/2.0;
            else
                A(i,j) = 0.0;
        }
    }
    
    //Generate vector b
    for (int j {0}; j < n; j++) {
        if ( (j == 0) || (j == n-1) )
            b(j) = 2.5;
        else if ( (j == n/2 - 1) || (j == n/2)  )
            b(j) = 1.0;
        else
            b(j) = 1.5;
    }
            
//    cout << "A = \n " << "\n \n " << A << "\n \n "  << endl;
//    cout << "b = \n " << "\n \n " << b << "\n \n " << endl;

    
    
    //initial guess x0
    VectorXd x0(n);
    for (int j {0}; j < n; j++)
        x0(j) = 2.0;
    
    //    cout << "x0 = \n " << "\n \n " << x0 << "\n \n "  << endl;
    
    
    cout << setprecision(8);
    
    cout << " Jacobi answer is = \n " << "\n \n " << Jacobi(A, b, x0) << "\n \n " << endl;
    cout << " Gauss Seidel answer is = \n " << "\n \n " << GaussSeidel(A, b, x0) << "\n \n " << endl;
    cout << " SOR answer is = \n " << "\n \n " << SOR(A, b, x0) << "\n \n " << endl;
    
    

//    //Final solution already given in the problem (for error-measuring)
//    VectorXd xf = VectorXd::Ones(n);
//
//    VectorXd errorvec = xf - Jacobi(A, b, x0);
////    VectorXd errorvec = xf - GaussSeidel(A, b, x0);
////    VectorXd errorvec = xf - SOR(A, b, x0);
//
//
////    cout << "xf = \n " << "\n \n " << xf << "\n \n "  << endl;
//////    cout << "xf - SOR = \n " << "\n \n " << xf - SOR(A, b, x0) << "\n \n "  << endl;
////    cout << "xf - GaussSeidel = \n " << "\n \n " << xf - GaussSeidel(A, b, x0) << "\n \n "  << endl;
//    cout << "xf - Jacobi = \n " << "\n \n " << xf - Jacobi(A, b, x0) << "\n \n "  << endl;
////
//    cout << "error norm " << "\n \n " << errorvec.lpNorm<Infinity>() << "\n \n "  << endl;
//

 
}

