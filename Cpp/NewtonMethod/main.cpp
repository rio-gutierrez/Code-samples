//  main.cpp
//  Newton's Method
//
//  Created by Mario L Gutierrez on 9/10/20.
//  Copyright © 2020 Mario L Gutierrez. All rights reserved.

#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
#include <algorithm>
#include <iterator>
#include <numeric>


using namespace std;


double f (double x){
    return pow(x,5) + x - 1.0;
}

double fder (double x){
    return 5.0 * pow(x,4) + 1.0;
}

double g (double x){
    return sin(x) - (6.0 * x) - 5.0;
}

double gder (double x){
    return cos(x) - 6.0;
}


double h (double x){
    return exp( pow(sin(x), 3) ) + pow(x,6) - ( 2.0 * pow(x,4) ) - pow(x,3) - 1.0;
}

double hder (double x){
    return ( cos(x) * pow(sin(x), 2) * 3.0 * exp( pow(sin(x), 3) ) ) + ( 6.0 * pow(x,5) ) - ( 8.0 * pow(x,3) ) - ( 3.0 * pow(x,2) );
}



void Newton (double (*func)(double), double (*funcder)(double), double x0, double epsilon = 0.5 * 1e-8, int itmax = 100){
    double x {x0};   //initialization of variable x to our initial guess x0
    double newx {};
    
    for (int i {0}; i <= itmax; i++){
        cout << "Iteration # " << i << "." << " x = " << x << "." << endl;
        newx = x - (func(x) )/( funcder(x) );
        if ( abs(newx - x) <= epsilon ) {
            cout << "Newton's Method has converged within given tolerance. The root found is x = " << newx << "." << endl;
            break;
        }
        if (i == itmax) {
            cout << "Newton's Method has failed to converge after the maximum allowed number of iterations. " << endl;
        }
        x = newx;
    }
}


int main(int argc, const char * argv[]) {
   
    cout << setprecision(7);
//    Newton(f, fder, 0.0);
//    Newton(g, gder, 0.0);
    Newton(h, hder, 0.5, 0.5 * 1e-6);
    
 

    return 0;
}
