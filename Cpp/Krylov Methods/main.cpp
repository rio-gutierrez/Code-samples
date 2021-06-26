//  main.cpp
//  Arnoldi
//  Created by Mario Luis on 6/5/21.

#include <istream>
#include <fstream>
#include <iostream>
#include <cmath>
#include <Eigen/Dense>
#include <Eigen/Sparse>
 
using namespace std;
using namespace Eigen;

typedef Matrix<double, Dynamic, Dynamic, 0, 100, 100> MatrixXd100;    //max number of rows/columns set to 100; avoids dynamic allocation
typedef Matrix<double, Dynamic, 1, 0, 100> VectorXd100;                       //max number of rows set to 100; avoids dynamic allocation

//Function prototypes
void  Arnoldi(const MatrixXd100 &A, const VectorXd100 &b, const size_t &n,  MatrixXd100 &Q, MatrixXd100 &H, const double tol = 1.0e-10);
void Galerkin(const MatrixXd100 &A, const VectorXd100 &b, VectorXd100 &x0, const size_t &n,  MatrixXd100 &Q, MatrixXd100 &H, const double tol = 1.0e-10);
void GMRES(const MatrixXd100 &A, const VectorXd100 &b, VectorXd100 &x0, const size_t &n,  MatrixXd100 &Q, MatrixXd100 &H, const double tol = 1.0e-10);
void Gaussian(MatrixXd100 &A,  VectorXd100 &b, char pivot_switch = 'n');





int main(int argc, const char * argv[]) {
    
    const size_t m {100};
    const size_t n {25};
    
    //Initialize arrays
    MatrixXd100 A = MatrixXd100::Zero(m,m);
    for (size_t i {0}; i < m; ++i){
        A(i,i) = 2.;
        if (i != m-1){
            A(i+1,i) = -1.;
            A(i,i+1) = A(i+1,i);
        }
    }

    VectorXd100 x = VectorXd100::Ones(m);
    VectorXd100 b = A*x;
    VectorXd100 x0  = VectorXd100::Zero(m);
    MatrixXd100 Q = MatrixXd100::Zero(m,n+1);
    MatrixXd100 H = MatrixXd100::Zero(n+1,n);
    
    Galerkin(A, b, x0, n, Q,  H);
//    GMRES(A, b, x0, n, Q, H);
    
    return 0;
}  //End of main function






//Function definitions
void  Arnoldi(const MatrixXd100 &A, const VectorXd100 &b,  const size_t &n,  MatrixXd100 &Q, MatrixXd100 &H, const double tol){
    /*Inputs:
     A = mxm matrix
     b = mx1 initial vector
     n = number of Arnoldi iterations (typically, n << m)
     Q = m x (n+1) matrix (initialized to zero)
     H = (n+1) x n matrix (initialized to zero)
     tol = user-defined tolerance. We need this to avoid dividing by zero (see below)
  */
    /* Purpose:
     - Orthogonalize Q via Modified Gram-Schmidt
     - Put H in upper-Hessenberg form
  */
    
    VectorXd100 v (b.size());
    Q.col(0) = b/b.norm();

    for (size_t k {0}; k < n; ++k){

        v = A*Q.col(k);

        for (size_t j {0}; j <= k; ++j){
            H(j,k) = Q.col(j).transpose() * v;
            v = v - H(j,k)* Q.col(j);
        }

        if (v.norm() <= tol){                                         //Happy breakdown :) No need to continue any further
            H.conservativeResize(k+1,k+1);              //Need to resize H and Q to get rid of unused zero entries
            Q.conservativeResize(A.rows(),k+1);     //using "conservative" to preserve nonzero entries already computed
            break;
        }
            
        H(k+1,k) = v.norm();
        Q.col(k+1) = v/H(k+1,k);
        
        if (k == n-1){
            H.conservativeResize(k+1,k+1);
            Q.conservativeResize(A.rows(),k+1);
        }
    }
    //End of Arnoldi. H is now Hesssenberg and Q is orthonormal
};






void Galerkin(const MatrixXd100 &A, const VectorXd100 &b, VectorXd100 &x0, const size_t &n,  MatrixXd100 &Q, MatrixXd100 &H, const double tol){

    cout << "  Applying FOM to solve Ax=b ... \n "
            << "--------------------------------- \n "<< endl;
    
    vector<size_t> it_vec {};
    vector<double> res_vec {};
    size_t it {0};
    const size_t it_max {100};   //max number of iterations allowed
    double residual {1.0};
    VectorXd100 x;
    VectorXd100 r0 = b - A*x0;
    
    do {
        
        it += 1;
        it_vec.push_back(it);
        
        Arnoldi(A, r0, n,  Q, H);
            
        VectorXd100 e0 = VectorXd100::Zero(H.rows());
        e0(0) = 1.0;

        //Solve the system Hy = ||r_0|| * e_0  for y using Guassian Elimination function
        VectorXd100 y = r0.norm()*e0;
        MatrixXd100 Hk = H;                //make a copy of H; otherwise it'd be destroyed in Gaussian function
        Gaussian(Hk, y, 'Y');
        
        //Compute the correction vector z
        VectorXd100 z = Q*y;
        
        //Add correction vector to original guess to get the new solution
        x = x0 + z;
        
        residual = (b - A*x).norm();
        res_vec.push_back(residual);
        cout << "The residual at iteration " << it << " is " << residual <<  ". \n " << endl;
        
        //update values for next iteration
        x0 = x;
        r0 = b - A*x;
        
        if (it == it_max) {
            cout << "Max number of iterations reached. \n " << endl;
            break;
        }
        
        //resize H and Q for next Arnoldi iteration
        H.conservativeResize(n+1,n);
        Q.conservativeResize(A.rows(),n+1);
        
    } while (residual > tol);
    
    cout << "Convergence obtained at iteration " << it << ". The solution is \n x = \n " << x << " \n \n " << endl;
        
    //Save residuals and iterations to files
    ofstream resfile ("Galekin_residuals.csv");
    ofstream itfile ("Galerkin_iterations.csv");
    for (int i{0}; i < res_vec.size(); ++i) {
        resfile << res_vec.at(i) << endl;
        itfile << it_vec.at(i) << endl;
    }
    resfile.close();
    itfile.close();

};





void GMRES(const MatrixXd100 &A, const VectorXd100 &b,   VectorXd100 &x0,  const size_t &n,  MatrixXd100 &Q, MatrixXd100 &H, const double tol){
    
    cout << "  Applying GMRES to solve Ax=b ... \n "
            << "--------------------------------- \n "<< endl;
    
    vector<size_t> it_vec {};
    vector<double> res_vec {};
    size_t it {0};
    const size_t it_max {100};   //max number of iterations allowed
    double residual {1.0};
    VectorXd100 x;
    VectorXd100 r0 = b - A*x0;
    
    do {
        
        it += 1;
        it_vec.push_back(it);
        
        Arnoldi(A, r0, n,  Q, H);
            
        VectorXd100 e0 = VectorXd100::Zero(H.rows());
        e0(0) = 1.0;
        
        /*Find the vector y that minimizes
         J(y) = || Hy - ||r_0|| * e_0 ||
      using QR least squares
    */
        VectorXd100 temp = r0.norm()*e0;
        VectorXd100 y =  H.colPivHouseholderQr().solve(temp);
        
        //Compute the correction vector z
        VectorXd100 z = Q*y;
        
        //Add correction vector to original guess to get the new solution
        x = x0 + z;
        
        residual = (b - A*x).norm();
        res_vec.push_back(residual);
        cout << "The residual at iteration " << it << " is " << residual <<  ". \n " << endl;
        
        //update values for next iteration
        x0 = x;
        r0 = b - A*x;
        
        if (it == it_max) {
            cout << "Max number of iterations reached. \n " << endl;
            break;
        }
        
        //resize H and Q for next Arnoldi iteration
        H.conservativeResize(n+1,n);
        Q.conservativeResize(A.rows(),n+1);
        
    } while (residual > tol);
    
    cout << "Convergence obtained at iteration " << it << ". The solution is \n x = \n " << x << " \n \n " << endl;
    
    //Save residuals and iterations to files
    ofstream resfile ("GMRES_residuals.csv");
    ofstream itfile ("GMRES_iterations.csv");
    for (int i{0}; i < res_vec.size(); ++i) {
        resfile << res_vec.at(i) << endl;
        itfile << it_vec.at(i) << endl;
    }
    resfile.close();
    itfile.close();
    
};




void Gaussian(MatrixXd100 &A,  VectorXd100 &b, char pivot_switch){
    /* Input A and b to solve system Ax=b for x.
     pivot_switch set to 'y' (yes) applies LU decomposition of A with pivoting;
     pivot_switch set to 'n' (no) or any other char (default case) applies LU decomposition of A without pivoting.
  */
    
    /* In this implementation both A and b are altered, to avoid extra memory allocation.
     An alternate code that leaves A and b untouched would need to allocate memory for
     explicit matrices L, U, P...
  */
    
    double pivot {};
    long pivot_index {};
    long m {A.rows()};       //using long for compatibility with Eigen's Matrix.rows() /.columns() type
    double l_ik {};


    
    switch (pivot_switch) {
            
        case 'Y':
        case 'y':     //yes...pivoting on
            
            for (long k {0}; k < m-1; ++k){
                
                pivot = A.col(k).tail(m-k).cwiseAbs().maxCoeff();      //set pivot as the max coeff (in abs value) of (the truncated) column k
                for (long i {k}; i < m; ++i){
                    if(  (A.col(k)(i) == pivot) || (A.col(k)(i) == -pivot) )
                        pivot_index = i;
                }
      
                A.row(pivot_index).swap(A.row(k));                      //swap rows in both A and b
                b.row(pivot_index).swap(b.row(k));

                for(long i {k+1}; i < m; ++i){

                    l_ik = A(i,k)/A(k,k);
                
                    for(long j {k}; j < m; ++j)
                        A(i,j) = A(i,j) - l_ik * A(k,j);

                    b(i) = b(i) - l_ik * b(k);
                }
            }
            break;
            
        case 'n':
        case 'N':
        default:    //no...pivoting off
            for (long k {0}; k < m-1; ++k){
                for(long i {k+1}; i < m; ++i){
                
                    l_ik = A(i,k)/A(k,k);
                    
                    for(long j {k}; j < m; ++j)
                        A(i,j) = A(i,j) - l_ik * A(k,j);
                    
                    b(i) = b(i) - l_ik * b(k);
                }
            }
            break;
    }           //end of switch statement...moving on to next phase...
    
    /* Final step (backsolving to solve Ux=y)
     Note, again, we are not allocating new memory for y or x;
     Instead, all calculations are done overwriting b
   */
    b(m-1) = b(m-1)/A(m-1,m-1);
    for(long k {m-2}; k >=0; --k){
        for(long j {k+1}; j < m; ++j)
             b(k) = b(k) - A(k,j) * b(j);
        b(k) = b(k)/A(k,k);
    }
    //Gaussian elimination process done, we now have b = x = solution.
};
