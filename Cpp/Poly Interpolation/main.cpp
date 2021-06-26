//  main.cpp
//  Polynomial Interpolation
//
//  Created by Mario L Gutierrez on 10/4/20.



#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
#include <algorithm>
#include <iterator>
#include <numeric>
#include <array>

using namespace std;



// Represent a data point corresponding to x and y = f(x)
struct Data
{
    double x, y;
};

//Define the Runge function
double Runge (double &x){
    return 1.0/(1.0 + 25.0 * pow(x,2));
}

//Define the Lagrange interpolation, to be evaluated at some point xi
double LagrangeInt(Data f[], double xi, const int n){
    double result {};
    double yval {};
  
    for (int i {0}; i <= n; i++){
        yval = f[i].y;      //set yval equal to y value of the ith data point
        for (int j {0}; j <= n; j++)
        {
            if (j!=i)
                yval = yval * (xi - f[j].x)/(f[i].x - f[j].x);
        }
        result += yval;
    }
  
    return result;
}




//Function to find the coefficients at the top of the triangle in Newton's code
void f_coeff(Data f[], int n, vector <vector <double>> &coeff){
   for (int j {0}; j <= n; j++) {
       for (int i {0}; i <= n-j; i++){
           if (j == 0){
               coeff.at(i).at(j) = f[i].y;
               coeff.at(i).push_back(coeff.at(i).at(j));
           }
           else{
               coeff.at(i).at(j) = ( coeff.at(i+1).at(j-1) - coeff.at(i).at(j-1) )/ ( f[i+j].x - f[i].x );
               coeff.at(i).push_back(coeff.at(i).at(j));
           }
       }
   }
}



// Function to find the product term to be used in Newton's code
double product(Data f[], double x, int i){
    double prod {1.0};
    for (int j {0}; j <= i-1; j++) {
        prod = prod * (x - f[j].x);
    }
    return prod;
}





//Define the Newton interpolation, to be evaluated at some point xi
double NewtonInt(Data f[], double x, int n){
    vector <vector <double>> coeff {};
    vector <double> coeff_vec {};
    
    for (int j {0}; j <= n; j++) {         //initialization of vector to size it as an (n-j)xj matrix
        for (int i {0}; i <= n-j; i++){
            coeff_vec.push_back(0.0);
        }
        coeff.push_back(coeff_vec);
    }


    f_coeff(f, n, coeff);

    double result {};

    for (int j {0}; j <= n; j++){
        if (j == 0)
            result = coeff.at(0).at(0);
        else
            result +=  coeff.at(0).at(j) * product(f, x, j);
    }

    return  result;
}




int main(int argc, const char * argv[]) {
   
    

    Data L[21] = {};    //Lagrange data array
    Data N[11] = {};    //Newton data array

    
        double xL{};
        double yL {};
        double dxL {0.01};

        for (int i {0}; i <= 20; i++){
            xL = 5.0 * cos( ( (2.0 * i + 1.0 ) * M_PI)/42.0   );
            /*These nodes go from x= -5 to x =  5;
                   thus the point xi in LagrangeInt must also be between -5 and 5 */
            yL = Runge(xL);
            L[i] = {xL, yL};
        }

    vector <double> Lagrange_vec {};
    for (int i {-500}; i <= 500; i++) {
        Lagrange_vec.push_back(LagrangeInt(L, i * dxL, 20));
    }
    
    
    double xN{};
    double yN {};
    double dxN {0.01};


    for (int i {0}; i <= 10; i++){
        xN = - 5.0 + i;
        /*These nodes go from x= -5 to x =  5;
               thus the point x in NewtonInt must also be between -5 and 5 */
        yN = Runge(xN);
        N[i] = {xN, yN};
    }


    vector <double> Newton_vec {};
    for (int i {-500}; i <= 500; i++) {
        Newton_vec.push_back(NewtonInt(N, i * dxN, 10));
    }

 
    
    //OUTPUT LAGRANGE DATA TO FILE
    ofstream myyfileL ("lagrange_y_data.csv");
     for (int i{0}; i <= 1000; i++) {
         if (i != 1000) {
             myyfileL << Lagrange_vec.at(i) << ",";
         } else {
             myyfileL << Lagrange_vec.at(i) << endl;
         }
     }

    myyfileL.close();
    
    
    
    //OUTPUT NEWTON DATA TO FILE
    ofstream myyfileN ("newton_y_data.csv");
     for (int i{0}; i <= 1000; i++) {
         if (i != 1000) {
             myyfileN << Newton_vec.at(i) << ",";
         } else {
             myyfileN << Newton_vec.at(i) << endl;
         }
     }

    myyfileN.close();



    return 0;
}
