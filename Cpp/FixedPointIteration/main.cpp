//  main.cpp
//  Fixed Point Iteration
//
//  Created by Mario L Gutierrez on 9/7/20.
//  Copyright © 2020 Mario L Gutierrez. All rights reserved.


#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
#include <algorithm>
#include <iterator>
#include <numeric>

using namespace std;


double g (double x){
    return ( - exp( pow(sin(x), 3) ) - pow(x,6) + ( 2.0 * pow(x,4) ) + 1.0 )/pow(x,2);
}


void FPI (double (*func)(double), double x0, int itmax){
    double x {x0};   //initialization of variable x to our initial guess x0
    int it {0};
    for (int i {0}; i <= itmax; i++){
        cout << "Iteration # " << it << "." << " x = " << x << "." << endl;
        it += 1;
        x = func(x);
        if (abs(func(x) - x) <= 0.5 * 1e-6) {
            break;
        }
    }
    cout << "The solution (that is, the fixed point) is x = " << x << "." << endl;
}


int main(int argc, const char * argv[]) {
   
    cout << setprecision(6);
    FPI(g, 1.0, 100);

    return 0;
}
