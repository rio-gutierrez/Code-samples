//  main.cpp
//  Bisection Method
//
//  Created by Mario L Gutierrez on 9/5/20.
//  Copyright © 2020 Mario L Gutierrez. All rights reserved.

#include <iostream>
#include <fstream>
//#include <math.h>
#include <cmath>
#include <vector>
#include <algorithm>
#include <iterator>
#include <numeric>

using namespace std;

double f (double x){
    return pow(x,5) + x - 1.0;
}

double g (double x){
    return sin(x) - 6.0 * x - 5.0;
}

double f2 (double x){
    return 1.0/x;
}


template <typename T>
int sign (const T &val) {
    return (val > 0) - (val < 0);
}



void Bisect (double (*func)(double), double a, double b, double epsilon, int itmax){  //itmax = max number of bisection iterations  we are willing to test
                                                                        //epsilon = minimum size of interval that we are willing to test
    double c {}; //initialization of variable c
    int it {0};  //iteration variable initialized to 0
    
    if (sign(func(a)) == sign(func(b))) {
        cout << "Signs at the endpoints are equal. There's no root in this interval, puny human." << endl;
    } else {
            while (b-a >= epsilon){
                 c = (a+b)/2.0; //definition of c as the midpoint
                 if (func(c) == 0.0){
                     cout << "Alas! We have found the (exact) root to be c = " << c << " at func = " << func(c) << "." << endl;
                     break;
                 }
                 else if (sign(func(a)) != sign(func(c))) {
                        b = c;
                 } else {
                        a = c;
                 }
                if (it == itmax) {
                    cout << "We have reached the max number of bisection iterations that we are willing to test. No exact root found; our best approximation is c = " << c << " at func = " << func(c) << "." << endl;
                    break;
                }
                cout << "it = " << it << endl;
                it += 1;
            }
        if (b-a < epsilon) {
            cout << "We have reached the minimum size of interval that we are willing to test. No exact root found; our best approximation is c = " << c << " at func = " << func(c) << "." << endl;
        }
    }
}


int main(int argc, const char * argv[]) {

cout << setprecision(8);   //set output precision to 8 significant digits
    
Bisect(f, 0.0, 1.0, 1e-10, 26);
Bisect(g, -1.0, 0.0, 1e-10, 26);
//Bisect(f2, -2.0, 1.0, 0.001, 19);

    return 0;
}
