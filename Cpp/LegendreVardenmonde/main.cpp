//  main.cpp
//  LegendreVardenmonde
//
//  Created by Mario L Gutierrez on 1/18/21.

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

    const int n {5};
    const int m {256};
    const int a {-1};
    const int b {1};
    const double h { double(b-a)/double(m)  };
    
    cout << setprecision(8);

    //Generate matrix A
    MatrixXd A(m,n);
    for (int i {0}; i < m; i++) {
        for (int j {0}; j < n; j++)
            A(i,j) = pow( (i-128)*h, j );
    }

    //Reduced QR Factorization of A
    MatrixXd thinQ(MatrixXd::Identity(m,n)), Q;
    HouseholderQR<MatrixXd> qr(A);
    Q = qr.householderQ();
    thinQ = qr.householderQ() * thinQ;
    
    
    //normalize thinQ
    MatrixXd scaleQ(m,n);
    for (int i {0}; i < m; i++) {
        for (int j {0}; j < n; j++){
            if (i == j)
                scaleQ(i,j) = 1.0/thinQ.row(m-1)(i);
            else
                scaleQ(i,j) = 0.0;
        }
    }
    Q = Q * scaleQ;
    
    //Save Q TO FILE
    ofstream myyfile ("Vandermonde.csv");
        for (int i {0}; i < m; i++) {
            for (int j {0}; j < n; j++){
                if (j == n-1)
                    myyfile << Q(i,j) << endl;
                else
                    myyfile << Q(i,j) << ",";
            }
        }
    myyfile.close();

    
    cout << "final Q is:\n" << Q << "\n\n";
}

