//  main.cpp
//  Secant Method
//
//  Created by Mario L Gutierrez on 9/11/20.
//  Copyright © 2020 Mario L Gutierrez. All rights reserved.


#include <iostream>
#include <fstream>
#include <cmath>


using namespace std;


double f (double x){
    return pow(x,3) - (2.0 * x) - 2.0;
}


void Secant (double (*func)(double), double xo, double x0, double epsilon = 0.5 * 1e-8, int itmax = 100){
    double x {xo};   //initialization of variable x to our initial guess xo
    double prevx  {x0}; //initialization of variable prevx to our initial guess x0
    double newx {};
    
    for (int i {0}; i <= itmax; i++){
        cout << "Iteration # " << i << "." << " x = " << x << "." << endl;
        newx = x - (  func(x) * (x - prevx) )/( func(x) - func(prevx) );
        if ( abs(newx - x) <= epsilon ) {
            cout << "Secant Method has converged within given tolerance. The root found is x = " << newx << "." << endl;
            break;
        }
        if (i == itmax) {
            cout << "Secant Method has failed to converge after the maximum allowed number of iterations. " << endl;
        }
        prevx = x;
        x = newx;
    }
}


int main(int argc, const char * argv[]) {
   
    cout << setprecision(7);
    Secant(f, 0.0, 1.0);

    return 0;
}
