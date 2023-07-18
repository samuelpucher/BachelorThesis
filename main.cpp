/*
    Author: Samuel Pucher
    Description: Time evolution algorithm

    Algo: 
    1. start with initial wavefunction psi(t)
    2. evaluate wavefunction at time psi(t+delta_t) = e^(-delta_t * H)psi(t)
    3. normalize psi
    4. go to 1
*/

//TODO: Normierung für Polarkoordinaten


//Includes
#include <iostream>
#include <cmath>
#include <iomanip>
#include <fstream>

using namespace std;


//GLobal constants
int GRID = 200;
double TMAX = 5;
int STEPS = 150000;
double T_STEP = TMAX/STEPS;
double HBAR = 1;
double OMEGA = 1;
double MASS = 1;
double LOWER = 0;
double UPPER = 20;

//Physical constants
double LAMBDA_IB = 1;
double LAMBDA_BB = 1;

//helper functions
double get_x(int i) {
    double res = (LOWER+(double(i)/GRID)*(UPPER-LOWER));
    return res;
}

double fder(int i, int grid, double* arr) {
    if(i == 0) {
        //std::cout << "Err1: lower bound not differentiable" << std::endl;
        return 0;
    }
    else if (i >(grid-1))
    {
        //std::cout << "Err1: upper bound not differentiable" << std::endl;
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
        return 0;
    }
    else if (i >=(grid-2))
    {
        //std::cout << "Err2: upper bound not differentiable" << std::endl;
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
        sum+=eps*arr[i]*arr[i];
    }

    for(int i= 0; i < GRID;i++) {
        arr[i] = arr[i]/(sqrt(sum));
    }
    return sqrt(sum);
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
    return (1/sqrt(100000));
}


//main function
int main(){

    //Program vars
    double* function_chi = new double[GRID];
    double* arr_fder_chi = new double[GRID];
    double* arr_sder_chi = new double[GRID];
    double* function_phi = new double[GRID];
    double* arr_fder_phi = new double[GRID];
    double* arr_sder_phi = new double[GRID];

    double t0 = 0;
    double time = t0;
    double delta_t = TMAX/STEPS;
    bool printed = false;

    //File mangament
    ofstream file_chi;
    ofstream file_phi;
    file_chi.open ("./results/function_chi.txt");
    file_phi.open ("./results/function_phi.txt");

    //init Functions
    init(start_func_chi, function_chi, arr_fder_chi, arr_sder_chi);  
    norm(function_chi);
    init(start_func_phi, function_phi, arr_fder_phi, arr_sder_phi);  
    norm(function_phi);

    //Main Loop
    for (int i = 0; i < STEPS; i++)
    {
        //Step 1: propagate chi
        for (int j = 0; j < GRID; j++)
        {
            function_chi[j] = exp(-delta_t/2*function_phi[j]*function_phi[j]*LAMBDA_IB)*function_chi[j];
        }
        for (int j = 0; j < GRID; j++)
        {
            function_chi[j] = function_chi[j] + HBAR*HBAR*delta_t/(2*MASS)*arr_sder_chi[j];
        }
        for (int j = 0; j < GRID; j++)
        {
            function_chi[j] = exp(-delta_t/2*function_phi[j]*function_phi[j]*LAMBDA_IB)*function_chi[j];
        }

        //Step 2: propagate phi
        for (int j = 0; j < GRID; j++)
        {
            function_phi[j] = exp(-delta_t/2*(function_phi[j]*function_phi[j]*LAMBDA_BB+ LAMBDA_IB*function_chi[j]*function_chi[j]))*function_phi[j];
        }
        for (int j = 0; j < GRID; j++)
        {
            function_phi[j] = function_phi[j] + HBAR*HBAR*delta_t/(2*MASS)*arr_sder_phi[j];
        }
        for (int j = 0; j < GRID; j++)
        {
            function_phi[j] = exp(-delta_t/2*(function_phi[j]*function_phi[j]*LAMBDA_BB+ LAMBDA_IB*function_chi[j]*function_chi[j]))*function_phi[j];
        }

        //Step 3: Norm 
        norm(function_chi);
        norm(function_phi);

        //Step 4: increase time
        time+=delta_t;

        //Step 5 update second derivative
        update_der(function_chi, arr_fder_chi, arr_sder_chi);
        update_der(function_phi, arr_fder_phi, arr_sder_phi);
    }

    //Write to file
    for (int i = 0; i < GRID; i++)
    {
        if(i%1==0){
            file_chi << get_x(i) << " " << function_chi[i] << "\n";    
            file_chi << "\n";
        }
        if(i%1==0){
            file_phi << get_x(i) << " " << function_phi[i] << "\n";           
            file_phi << "\n";
        }
    }
    

    //Close file
    file_chi.close();
    file_phi.close();
    
}