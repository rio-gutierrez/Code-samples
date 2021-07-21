//  main.cpp
//  VCD3D_lambda1.20
//
//  Created by Mario L Gutierrez on 2/21/20.
//  Copyright © 2020 Mario L Gutierrez. All rights reserved.

#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
#include <algorithm>
#include <iterator>
#include <numeric>

#include <chrono>


using namespace std;
using namespace chrono;


//Function prototypes
double TauAvg(vector<vector<double>> &vect, int taumax);
double RAvg(vector<vector<double>> &vect, int rmax);
double MinDistToAvg(vector<vector<double>> &vect, const double &average, int rmax);
double MinDistToAvgTau(vector<vector<double>> &vect, const double &average, int taumax);
double RadiusFunction(vector<vector<double>> &vect, const double &instanton_value, const double step_size, int rmax);
double RadiusFunctionTau(vector<vector<double>> &vect, const double &instanton_value, const double step_size, int taumax);
double Simpson2D (vector<vector<double>> &Fij,  const double step_size, int taumax, int rmax);
double S_Ansatz (vector <double> &temp, vector <double> &action);


// functor for getting sum of previous result and square inverse of current element
template<typename T>
struct invsqr
{
  T operator()(const T& Left, const T& Right) const
  {
      return (Left + pow(Right*Right, -1) );
  }
};



int main(int argc, const char * argv[]) {


    auto start = high_resolution_clock::now();


    //****************************************************************
    //              Variables
    //****************************************************************
    const double FV {M_PI}, lambda {1.20}, h {0.1}, k{1.0}, R{2.5}, width {1.0};
    const double DeltaS{(1.0 / 8.0) * pow(h,2)};  // e {0.0}, a {1.0 + pow(e,2.0)}, b {1.0 - pow(e,2.0)}
    const int    r_max {120}, it_max {80000};    // r_max = 120 = Phi.size()-1;
    double action {}, total_action {}, s{}, d{}, phi {}, nphi {}, phi_t {}, phi_r {}, phi_tt {}, phi_rr {},
            F_r {}, F_rr {}, F_tt {}, U {}, U_p {}, U_pp {}, f {}, dagF {}, F_max {}, T{}, avg {}, dist {}, radius {}, tau_avg {}, tau_dist {}, tau_radius {};




    //****************************************************************
    //              Vectors
    //****************************************************************
    vector <double> Phi_vec {};
    vector <vector<double>> Phi {};
    vector <vector<double>> pPhi {};
    vector <vector<double>> nPhi {};

    vector <double> evo_vec {};
    vector <vector<double>> evo_grid {};


    vector <double> F_vec {};
    vector <vector<double>> F {};
    vector <vector<double>> evoF {};

    vector <double> action_vec {};

    vector <vector<double>> Fij {};
    vector <double> ACTION {};
    vector <double> TEMP {};
    vector <double> RADIUS {};
    vector <double> RADIUS_tau {};

     vector <double> res_vec {};
     vector <double> alpha_vec {};

    vector <double> ellip {0.0, 0.25, 0.50, 0.75, 0.99};

    for (int e {0}; e <= ellip.size()-1; e++) {

        double a {1.0 - pow(ellip.at(e),2)}, b {1.0 + pow(ellip.at(e),2)};

        //vary the number of grid points in the tau direction to vary the temperature
        for (int tau_max{120}; tau_max >= 20; tau_max -= 10){

            //****************************************************************
            //              Initial data   (Phi)
            //****************************************************************

            for (int i{0}; i <= tau_max; i++) {
                for (int j{0}; j <= r_max; j++) {
                    d = sqrt( (a * pow(i*h,2)) + (b * pow(j*h,2)) );
                    phi = 0.5 * M_PI * ( 2.0 + tanh((d-R)/width) - tanh((d+R)/width) );
                    Phi_vec.push_back(phi);
                }
                Phi.push_back(Phi_vec);
                Phi_vec.clear();
            }


            //****************************************************************
            //              Numerical Solution (Strong Relaxation)
            //****************************************************************


            for (int it{0}; it <= it_max; it++){

                s = DeltaS * it;

                //F values
                //*************************************************************************************//
                for (int i{0}; i <= tau_max; i++) {
                    for (int j{0}; j <= r_max; j++) {
                        //Derivatives of Phi
                        //*****************************//
                        if (i == 0)
                            phi_tt = 2.0 * ( Phi.at(1).at(j) - Phi.at(0).at(j) );
                        else if (i == tau_max)
                            phi_tt = 2.0 * ( Phi.at(tau_max - 1).at(j) - Phi.at(tau_max).at(j) );
                        else
                            phi_tt = Phi.at(i + 1).at(j) - ( 2.0 * Phi.at(i).at(j) ) + Phi.at(i - 1).at(j);

                        if (j == 0){
                            phi_rr = 2.0 * ( Phi.at(i).at(1) - Phi.at(i).at(0) );
                            phi_r = 0.0;    //never actually going to use this value; may omit it without issues
                        }
                        else if (j == r_max){
                            phi_rr = FV - ( 2.0 * Phi.at(i).at(r_max) ) + Phi.at(i).at(r_max - 1);
                            phi_r = FV - Phi.at(i).at(r_max - 1);
                        }
                        else{
                            phi_rr = Phi.at(i).at(j + 1) - ( 2.0 * Phi.at(i).at(j) ) + Phi.at(i).at(j - 1);
                            phi_r = Phi.at(i).at(j + 1) - Phi.at(i).at(j - 1);
                        }
                        //*****************************//


                        U_p = 0.5 * (pow(lambda, 2) * sin(2.0 * Phi.at(i).at(j))) + sin(Phi.at(i).at(j));
                        if (j == 0)
                            f = (1.0/pow(h, 2)) * (phi_tt + (3.0 * phi_rr)) - U_p;
                        else
                            f = (1.0/pow(h, 2)) * (phi_tt + phi_rr) + ( 1.0/(j * pow(h,2)) * phi_r ) - U_p;

                        F_vec.push_back(f);
                    }   //end of 'j' iteration


                    //Calculating residuals
                    //*****************************//
                    auto F_minmax = minmax_element( F_vec.begin(), F_vec.end() ) ;
                    double minval = *F_minmax.first;
                    double maxval = *F_minmax.second;
                    //cout << "The min value is " << minval << " and the max value is " << maxval << endl;

                    if (fabs(minval) > fabs(maxval))
                        F_max = fabs(minval);
                    else
                        F_max = fabs(maxval);

                    res_vec.push_back(F_max);
                    //*****************************//

                    evoF.push_back(F_vec);
                    F_vec.clear();
                }  //end of 'i' iteration


                F = evoF;
                evoF.clear();

                //Display residuals
                //*****************************//
                auto residual = max_element( res_vec.begin(), res_vec.end() ) ;
                res_vec.clear();

                if ((it % (it_max/5)) == 0)
                    cout << "Step # " << it << ". The residual at s = " << s << " is " << *residual << "." << endl;
                //*****************************//



                //time to evolve the grid
                for (int m{0}; m <= tau_max; m++) {
                    for (int n {0}; n <= r_max; n++) {


                        //Second Derivatives of F
                        //*****************************//
                        if (m == 0)
                            F_tt = 2.0 * ( F.at(1).at(n) - F.at(0).at(n) );
                        else if (m == tau_max)
                            F_tt = 2.0 * ( F.at(tau_max-1).at(n) - F.at(tau_max).at(n) );
                        else
                            F_tt = F.at(m+1).at(n) - 2.0 * F.at(m).at(n) + F.at(m-1).at(n);

                        if (n == 0){
                            F_rr = 2.0 * ( F.at(m).at(1) - F.at(m).at(0) );
                            F_r = 0.0;    //never actually going to use this value; may omit it without issues
                        }
                        else if (n == r_max){
                            F_rr = - ( 2.0 * F.at(m).at(r_max) ) + F.at(m).at(r_max-1);
                            F_r = - F.at(m).at(r_max-1);
                        }
                        else{
                            F_rr = F.at(m).at(n+1) - ( 2.0 * F.at(m).at(n) ) + F.at(m).at(n-1);
                            F_r = F.at(m).at(n+1) - F.at(m).at(n-1);
                        }
                        //*****************************//


                        U_pp = pow(lambda,2) * cos(2.0 * Phi.at(m).at(n)) + cos(Phi.at(m).at(n));

                        if (n == 0) {
                            dagF = ( - 1.0/(pow(h,2)) ) * ( F_tt + (3.0 * F_rr) ) + ( U_pp * F.at(m).at(n) );
                        } else {
                            dagF = ( - 1.0/(pow(h,2)) ) * (F_tt + F_rr) - ( 1.0/(n * pow(h,2)) * F_r ) + ( U_pp * F.at(m).at(n) );
                        }


                        // 's' Evolution of Phi
                        //*****************************//
                        if (it == 0)
                            nphi  = ( 1.0/(1.0 + (0.5 * k * DeltaS)) )  * ( ( 2.0 * Phi.at(m).at(n) )  - ( Phi.at(m).at(n) * (1.0 - (0.5 * k * DeltaS)) ) + ( dagF * pow(DeltaS,2) ) );
                        else
                            nphi  = ( 1.0/(1.0 + (0.5 * k * DeltaS)) )  * (  ( 2.0 * Phi.at(m).at(n) ) - ( pPhi.at(m).at(n) * (1.0 - (0.5 * k * DeltaS)) ) + ( dagF * pow(DeltaS,2) ) );

                        evo_vec.push_back(nphi);

                    }   //end of 'n' iteration


                    evo_grid.push_back(evo_vec);
                    evo_vec.clear();

                }   //end of 'm' iteration       //whole Phi grid covered
                nPhi = evo_grid;
                evo_grid.clear();

                pPhi = Phi;
                Phi = nPhi;



                //Calculate final r and tau radiuses of the bubble at the end of last iteration
                  //*****************************//
                  if ( (it == it_max) || (*residual <= 0.01) ) {

                          avg = RAvg(Phi, r_max);
                          tau_avg = TauAvg(Phi, tau_max);

                          dist = MinDistToAvg(Phi, avg, r_max);
                          tau_dist = MinDistToAvgTau(Phi, tau_avg, tau_max);

                          for (double phi_val : Phi.at(0)) {
                              if (abs(phi_val - avg) == dist) {
                                  radius = RadiusFunction(Phi, phi_val, h, r_max);
                              }
                          }

                          for (int i{0}; i <= tau_max; i++) {
                              if (abs(Phi.at(i).at(0) - tau_avg) == tau_dist) {
                                  tau_radius = RadiusFunctionTau(Phi, Phi.at(i).at(0), h, tau_max);
                             }
                         }

                  }

                  if(*residual <= 0.01){
                    cout << "residual reached the threshold at it = " << it << "; no need to reach it_max. Moving on ..." << endl;
                    break;
                  }

            }    //end of it iteration; run is over


            RADIUS.push_back(radius);   //put all the radius values in one vector; there'll be 11 entries in the vector at the end, one for each tau_max value
            RADIUS_tau.push_back(tau_radius);   //put all the tau radius values in one vector


            //****************************************************************
            //              Calculate Action
            //****************************************************************

            for (int i{0}; i <= tau_max; i++) {
                for (int j{0}; j <= r_max; j++) {
                    //Second Derivatives of Phi
                    //*****************************//
                    if (i == 0 || i == tau_max)
                        phi_t = 0.0;
                    else
                        phi_t = Phi.at(i+1).at(j) - Phi.at(i-1).at(j);

                    if (j == 0)
                        phi_r = 0.0;
                    else if (j == r_max)
                        phi_r = FV - Phi.at(i).at(r_max-1);
                    else
                        phi_r = Phi.at(i).at(j+1) - Phi.at(i).at(j-1);
                    //*****************************//

                    U = 0.5 * ( pow(lambda,2) * pow(sin(Phi.at(i).at(j)),2) ) - cos(Phi.at(i).at(j)) - 1.0;
                    action = ( 4.0 * M_PI * (pow(j * h,2)) ) * ( ( 1.0/(8.0 * pow(h,2)) ) * ( pow(phi_t,2) + pow(phi_r,2) ) + U );
                    action_vec.push_back(action);
                }   //end of 'j' iteration

                 Fij.push_back(action_vec);
                 action_vec.clear();
            }  //end of 'i' iteration

            total_action = Simpson2D(Fij, h, tau_max, r_max);
            ACTION.push_back(2.0 * total_action);
            Fij.clear();

            T = 1.0/(2.0 * h * tau_max);
            TEMP.push_back(T);


            cout << "The total action after this run is " << 2.0 * total_action << " at T= " << T
                   << ", with ellipticity = " << ellip.at(e) << ". The radius of the bubble is "
                   << radius << "." << endl;

            Phi.clear();
            pPhi.clear();
            nPhi.clear();

        } //end of taumax iteration; several runs have been performed and we have several values of the action with different temperatures


        //ansatz S = 1/alpha
        if (e == 4) {
            alpha_vec.push_back(S_Ansatz(TEMP, ACTION));
            cout << "alpha = " << alpha_vec.at(0) << endl;
        }

        ACTION.clear();
        TEMP.clear();
        RADIUS.clear();
        RADIUS_tau.clear();


    }  //end of 'e' loop



    auto stop = high_resolution_clock::now();

    auto duration = duration_cast<seconds>(stop - start);

    // To get the value of duration use the count() member function on the duration object
    cout << "Program ended. It took "<< duration.count() << " seconds. (" << duration.count()/60.0 << " minutes)" << endl;


    return 0;

}




/*---------------------------
    FUNCTION DEFINITIONS
 ---------------------------*/

//Function "TauAvg" that determines the average value of vect (Phi) along the tau axis
double TauAvg(vector<vector<double>> &vect, int taumax)
{
//initialize min/max values along tau direction
double tau_minval {M_PI};
double tau_maxval {0.0};
double tau_avg {};

for (int i {0}; i <= taumax; i++){
  if (vect.at(i).at(0) <= tau_minval)
      tau_minval = vect.at(i).at(0);
  else if (vect.at(i).at(0) >= tau_maxval)
      tau_maxval = vect.at(i).at(0);
}

tau_avg = (tau_maxval + tau_minval) / 2.0;

return tau_avg;
}


// Function "RAvg" that determines the average value of vect (Phi) along the r axis
double RAvg(vector<vector<double>> &vect, int rmax)
{
//initialize min/max values along tau direction
double r_minval {M_PI};
double r_maxval {0.0};
double r_avg{};

for (int j {0}; j <= rmax; j++){
  if (vect.at(0).at(j) <= r_minval)
      r_minval = vect.at(0).at(j);
  else if (vect.at(0).at(j) >= r_maxval)
      r_maxval = vect.at(0).at(j);
}

r_avg = (r_maxval + r_minval) / 2.0;

return r_avg;
}



/* Function "MinDistToAvg" that determines minimum difference between any value
along the 'x/r' direction in the 2D grid "vect" and some specified value "average"
*/
double MinDistToAvg(vector<vector<double>> &vect, const double &average, int rmax)
{
// Initialize difference as some (random) small value; 0.5 suffices
double epsilon {0.5};

for (int j {0}; j <= rmax; j++){
if (abs(vect.at(0).at(j) - average) < epsilon)
   epsilon = abs(vect.at(0).at(j) - average);
}

// Return minimum difference
return epsilon;
}



/* Function "MinDistToAvg" that determines minimum difference between any value
 along the 'tau' direction in the 2D grid "vect" and some specified value "average"
*/
double MinDistToAvgTau(vector<vector<double>> &vect, const double &average, int taumax)
{
// Initialize difference as some (random) small value; 0.5 suffices
double epsilon {0.5};

for (int i {0}; i <= taumax; i++){
if (abs(vect.at(i).at(0) - average) < epsilon)
   epsilon = abs(vect.at(i).at(0) - average);
}

// Return minimum difference
return epsilon;
}


// Function that spews out the radius of the bubble   (along the 'r/x' axis)
double RadiusFunction(vector<vector<double>> &vect, const double &instanton_value, const double step_size, int rmax)
{
// Initialize radius
double rad {};

for (int j {0}; j <= rmax; j++){
      if (instanton_value == vect.at(0).at(j))
            rad = j * step_size;
}

// Return radius
return rad;
}



// Function that spews out the radius of the bubble   (along the 'tau' axis)
double RadiusFunctionTau(vector<vector<double>> &vect, const double &instanton_value, const double step_size, int taumax)
{
// Initialize radius
double rad {};

for (int i {0}; i <= taumax; i++){
      if (instanton_value == vect.at(i).at(0))
            rad = i * step_size;
}

// Return radius
return rad;
}



// 2D Simpson's Rule
double Simpson2D (vector<vector<double>> &Fij,  const double step_size, int taumax, int rmax)
{
  // Initialize variables
  vector<vector<double>> Sij (taumax+1, vector<double> ( rmax+1, 0.0 ) );   //initialize 2D  vector Sij to have size 121 x 121
  vector<double> Si {};
  vector<double> Action_Vec {};
  double actionvals {};
  double Action {};
  double Total_Action {};


  for (int i {0}; i <= taumax; i++){
       for (int j {0}; j <= rmax; j++){
           if (i == 0 || i == taumax) {
               if (j == 0 || j == rmax)
                   Sij.at(i).at(j) = 1.0;
               else if (j % 2 != 0)
                   Sij.at(i).at(j) = 4.0;
               else
               Sij.at(i).at(j) = 2.0;
           }
           else if (j == 0 || j == rmax) {
               if (i % 2 != 0)
                   Sij.at(i).at(j) = 4.0;
               else
               Sij.at(i).at(j) = 2.0;
           }
           else {
               if ( (i % 2 != 0) && (j % 2 != 0) )
                   Sij.at(i).at(j) = 16.0;
               else if ( ((i % 2 != 0) && (j % 2 == 0)) || ((i % 2 == 0) && (j % 2 != 0)) )
                   Sij.at(i).at(j) = 8.0;
               else
                   Sij.at(i).at(j) = 4.0;
           }

           actionvals = (Sij.at(i).at(j)) * (Fij.at(i).at(j));
           Si.push_back(actionvals);
       }   //end of main 'j' loop

      Action = accumulate(Si.begin(), Si.end(), 0.0 );
      Action_Vec.push_back(Action);
      Si.clear();
  }    //end of main 'i' loop


  Total_Action = (pow(step_size,2)/9.0) * ( accumulate(Action_Vec.begin(), Action_Vec.end(), 0.0 ) );
  Action_Vec.clear();

  return Total_Action;
}


// S = alpha/T
double S_Ansatz (vector <double> &temp, vector <double> &action)
{
  // Initialize variables
  double xval_invsqr {};
  double yval_over_xval {};
  vector <double> y_over_x {};
  double alpha {};

  for (int i{0}; i <= temp.size()-1; i++) {
      y_over_x.push_back( (action.at(i))/(temp.at(i)));
  }

  yval_over_xval = accumulate(y_over_x.begin(), y_over_x.end(), 0.0 );
  y_over_x.clear();
  xval_invsqr = accumulate(temp.begin(), temp.end(), 0.0, invsqr<double>() );

  alpha = yval_over_xval/xval_invsqr;

  return alpha;
}
