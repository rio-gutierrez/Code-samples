//  main.cpp
//  Gram-Schmidt orthogonalization
//
//  Created by Mario L Gutierrez on 1/15/21.

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



int main(int argc, const char * argv[]) {

    const int n {80};
    
    MatrixXd A(n,n);
    MatrixXd Sigma(n,n);
    MatrixXd Vrand = MatrixXd::Random(n,n), V;
    MatrixXd Urand = MatrixXd::Random(n,n), U;

    //Generate matrix Sigma
    for (int i {0}; i < n; i++) {
        for (int j {0}; j < n; j++){
            //construction of matrix Sigma
            if (i == j)
                 Sigma(i,j) = pow(2.0,(-(i+1)));
             else
                 Sigma(i,j) = 0.0;
        }
    }

    //Orthogonalize random matrices Urand and Vrand
    HouseholderQR<MatrixXd> qr1(Vrand);
    HouseholderQR<MatrixXd> qr2(Urand);
    V = qr1.householderQ();
    U = qr2.householderQ();


    A = U * Sigma * V.transpose();

    MatrixXd Q(n,n);
    MatrixXd R = MatrixXd::Zero(n,n);
    VectorXd v(n);

//    Classical Gram-Schmidt
    for (int j {0}; j < n; j++){
        v = A.col(j);
        for (int i {0}; i < j; i++) {
            R(i,j) = Q.col(i).transpose() * A.col(j);
            v = v - (R(i,j) * Q.col(i));
        }
        R(j,j) = v.norm();
        Q.col(j) = v/R(j,j);

        }
    
//    //    Modified Gram-Schmidt
//        for (int i {0}; i < n; i++)
//            Q.col(i) = A.col(i);
//        for (int i {0}; i < n; i++){
//            R(i,i) = Q.col(i).norm();
//            Q.col(i) = Q.col(i)/R(i,i);
//            for (int j {i+1}; j < n; j++) {
//                R(i,j) = Q.col(i).transpose() * Q.col(j);
//                Q.col(j) = Q.col(j) - (R(i,j) * Q.col(i));
//            }
//        }

        cout << "R = \n " << "\n \n " << R << endl;
//        cout << "Q = \n " << "\n \n " << Q << endl;
    

    
    //    //Save r_{jj} TO FILE
        ofstream myyfile ("Rjj.csv");
         for (int j{0}; j < n; j++) {
                 myyfile << R(j,j) << endl;
         }
        myyfile.close();

 
}

