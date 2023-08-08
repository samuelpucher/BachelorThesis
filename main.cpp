/*
    Author: Samuel Pucher
    Description: Time evolution algorithm

    Algo: 
    1. start with initial wavefunction psi(t)
    2. evaluate wavefunction at time psi(t+delta_t) = e^(-delta_t * H)psi(t)
    3. normalize psi
    4. go to 1
*/

//TODO: Normierung für Phi: Rand ->1


//Includes
#include <iostream>
#include <cmath>
#include <iomanip>
#include <fstream>

using namespace std;


//GLobal constants
int GRID = 500;
double TMAX = 10;
int STEPS = 1000000;
double T_STEP = TMAX/STEPS;         //T_STEP should be arround 0.000001 = 1e-6
double HBAR = 1;
double OMEGA = 1;
double MASS = 1;
double LOWER = 0;
double UPPER = 20;
double DENSITY = 1;

//Physical constants
double LAMBDA_IB = 1;
double LAMBDA_BB = 1;

//helper functions
double get_x(int i) {
    double res = (LOWER+(double(i)/GRID)*(UPPER-LOWER));
    return res;
}

double fder(int i, int grid, double* arr) {
    if(i < 2) {
        double eps = (UPPER-LOWER)/double(grid);
        return ((arr[2]-arr[0])/(2*eps));
    }
    else if (i >=(grid-1))
    {
        //double eps = (UPPER-LOWER)/double(grid);
        //return ((arr[grid-3]-arr[grid-1])/(2*eps));
        return 0;
    }
    else {
        double eps = (UPPER-LOWER)/double(grid);
        return ((arr[i+1]-arr[i-1])/(2*eps));
    }
}

double sder(int i, int grid, double* arr) {
    if(i < 2) {
        //std::cout << "Err2: lower bound not differentiable" << std::endl;
        double eps = (UPPER-LOWER)/double(grid);
        return (arr[2]+arr[0]-2*arr[1])/(eps*eps);
    }
    else if (i >=(grid-1))
    {
        //std::cout << "Err2: upper bound not differentiable" << std::endl;
        //double eps = (UPPER-LOWER)/double(grid);
        //return (arr[grid-1]+arr[grid-3]-2*arr[grid-2])/(eps*eps);
        return 0;
    }
    else {
        double eps = (UPPER-LOWER)/double(grid);
        return (arr[i+1]+arr[i-1]-2*arr[i])/(eps*eps);
    }
}

double update_sder(int i, int grid, double* arr) {
    if(i < 2) {
        //std::cout << "Err2: lower bound not differentiable" << std::endl;
        return arr[2];
    }
    else if (i >=(grid-2))
    {
        //std::cout << "Err2: upper bound not differentiable" << std::endl;
        return arr[grid-3];
    }
    else {
        double eps = (UPPER-LOWER)/double(grid);
        return (arr[i+1]+arr[i-1]-2*arr[i])/(eps*eps);
    }
}

void init(double (*func)(double), double* arr1, double* arr2, double* arr3){
    for (int i = 0; i < GRID; i++)
    {
        double x = get_x(i);
        arr1[i] = func(x);
    }
    for (int i = 0; i < GRID; i++)
    {
        arr2[i] = fder(i, GRID, arr1);
        arr3[i] = sder(i, GRID, arr1);
    }
}

double norm(double* arr){

    double sum = 0;
    double eps = (UPPER-LOWER)/double(GRID);

    for(int i= 0; i< GRID;i++) {
        sum+=eps*arr[i]*arr[i]*get_x(i);
    }

    for(int i= 0; i < GRID;i++) {
        arr[i] = arr[i]/(sqrt(sum));
    }
    return sqrt(sum);
}

double norm_phi(double* arr){
    /*
    double sum = 0;
    double eps = (UPPER-LOWER)/double(GRID);
    double volume = M_PI*(UPPER-LOWER)*(UPPER-LOWER);
    double N = volume * DENSITY;
    
    for(int i= 0; i< GRID;i++) {
        sum+=eps*arr[i]*arr[i]*get_x(i);
    }

    for(int i= 0; i < GRID;i++) {
        arr[i] = arr[i]*sqrt(N)/(sqrt(sum));
    }
    return sqrt(sum);
    */
   double lastVal = arr[GRID-1];

   for(int i= 0; i < GRID;i++) {
        arr[i] = arr[i]/lastVal;
    }
    return lastVal;
}

double integral(double* arr){

    double sum = 0;
    double eps = (UPPER-LOWER)/double(GRID);

    for(int i= 0; i< GRID;i++) {
        sum+=eps*arr[i];
    }
    return sum;
}

void update_der(double* arr1, double* arr2, double* arr3){

    for (int i = 0; i < GRID; i++)
    {
        arr2[i] = fder(i, GRID, arr1);
        arr3[i] = update_sder(i, GRID, arr1);
    }

}

double start_func_chi(double x) {
    return (exp(-(x)*(x)));
}

double start_func_phi(double x) {
    return (1);
}

double approx_a(double f_2eps, double f_3eps, double eps){
    double b = (f_3eps-f_2eps)/(5*eps*eps);
    return (f_2eps-4*b*eps*eps);
}


//main function
int main(){

    //Program vars
    double* function_chi = new double[GRID];
    double* function_chi_old = new double[GRID];
    double* arr_fder_chi = new double[GRID];
    double* arr_sder_chi = new double[GRID];
    double* function_phi = new double[GRID];
    double* arr_fder_phi = new double[GRID];
    double* arr_sder_phi = new double[GRID];

    double t0 = 0;
    double time = t0;
    double delta_t = TMAX/STEPS;
    bool printed = false;

    double eps = (UPPER-LOWER)/double(GRID);

    //File mangament
    ofstream file_chi;
    ofstream file_phi;
    ofstream file_chi_final;
    ofstream file_phi_final;
    file_chi.open ("./results/function_chi.txt");
    file_phi.open ("./results/function_phi.txt");
    file_chi_final.open ("./results/function_chi_final.txt");
    file_phi_final.open ("./results/function_phi_final.txt");

    //init Functions
    init(start_func_chi, function_chi, arr_fder_chi, arr_sder_chi);  
    norm(function_chi);
    init(start_func_phi, function_phi, arr_fder_phi, arr_sder_phi);  
    norm_phi(function_phi);



    //Main Loop
    for (int i = 0; i < STEPS; i++)
    {
        //Step 0: save old chi function
        for (int j = 0; j < GRID; j++)
        {
            function_chi_old[j] = function_chi[j];
        }
        
        //Step 1: propagate chi
        for (int j = 2; j < GRID; j++)
        {
            function_chi[j] = exp(-delta_t/2*function_phi[j]*function_phi[j]*LAMBDA_IB)*function_chi[j];
        }
        for (int j = 2; j < GRID; j++)
        {
            function_chi[j] = function_chi[j] + HBAR*HBAR*delta_t/(2*MASS)*(arr_sder_chi[j] + arr_fder_chi[j]/get_x(j));
        }
        for (int j = 2; j < GRID; j++)
        {
            function_chi[j] = exp(-delta_t/2*function_phi[j]*function_phi[j]*LAMBDA_IB)*function_chi[j];
        }

        double a_chi= approx_a(function_chi[2], function_chi[3], eps);
        function_chi[0] = a_chi;
        function_chi[1] = (function_chi[2]+a_chi)/2;

        //Step 2: propagate phi
        for (int j = 2; j < GRID; j++)
        {
            function_phi[j] = exp(-delta_t/2*(function_phi[j]*function_phi[j]*LAMBDA_BB+ LAMBDA_IB*function_chi_old[j]*function_chi_old[j]))*function_phi[j];
        }
        for (int j = 2; j < GRID; j++)
        {
            function_phi[j] = function_phi[j] + HBAR*HBAR*delta_t/(2*MASS)*(arr_sder_phi[j] + arr_fder_phi[j]/get_x(j));
        }
        for (int j = 2; j < GRID; j++)
        {
            function_phi[j] = exp(-delta_t/2*(function_phi[j]*function_phi[j]*LAMBDA_BB+ LAMBDA_IB*function_chi_old[j]*function_chi_old[j]))*function_phi[j];
        }

        double a_phi= approx_a(function_phi[2], function_phi[3], eps);
        function_phi[0] = a_phi;
        function_phi[1] = (function_phi[2]+a_phi)/2;
        function_phi[GRID-1] = function_phi[GRID-2] = function_phi[GRID-3];

        //Step 3: Norm 
        norm(function_chi);
        norm_phi(function_phi);

        //Step 4: increase time
        time+=delta_t;

        //Step 5 update derivative
        update_der(function_chi, arr_fder_chi, arr_sder_chi);
        update_der(function_phi, arr_fder_phi, arr_sder_phi);

        if(i%20000 == 0) {
            for (int j = 0; j < GRID; j++)
            {
                file_chi << get_x(j) << " " << time << " " << function_chi[j] << "\n";    
                file_phi << get_x(j) << " " << time << " " << function_phi[j] << "\n";           
            }
            file_chi << "\n";
            file_phi << "\n";
        }
    }

    //Write to file
    for (int i = 0; i < GRID; i++)
    {
        if(i%1==0){
            file_chi_final << get_x(i) << " " << function_chi[i] << "\n";    
            file_chi_final << "\n";
            file_phi_final << get_x(i) << " "  << function_phi[i] << "\n";           
            file_phi_final << "\n";
        }
    }
    

    //Close file
    file_chi.close();
    file_phi.close();
    file_chi_final.close();
    file_phi_final.close();
    
}