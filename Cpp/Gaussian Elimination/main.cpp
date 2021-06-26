//  main.cpp
//  Gaussian Elimination
//  Created by Mario Luis on 6/5/21.

#include <istream>
#include <fstream>
#include <iostream>
#include <cmath>
#include <Eigen/Dense>
#include <Eigen/Sparse>
 
using namespace std;
using namespace Eigen;

//max number of rows/columns set to 100; avoids dynamic allocation
typedef Matrix<double, Dynamic, Dynamic, 0, 100, 100> MatrixXd100;
typedef Matrix<double, Dynamic, 1, 0, 100> VectorXd100;


//Function prototypes
void Gaussian(MatrixXd100 &A,  VectorXd100 &b, char pivot_switch = 'n');



int main(int argc, const char * argv[]) {
    
    const size_t m {100};
    
    //Initialize arrays
    MatrixXd100 A = MatrixXd100::Random(m,m);
    VectorXd100 b =  VectorXd100::Random(m);;

    Gaussian(A,b, 'y');
    
    cout << "The solution, x, to the system Ax=b is \n "  << b << endl;   //b is overwritten in the Gaussian routine to yield the solution x
    
    return 0;
}  //End of main function




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
