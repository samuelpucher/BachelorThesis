/*
    Author: Samuel Pucher
    Description: Time evolution algorithm

    GNUplot cmd: splot [*:*] [*:*] [*:*] "res.txt" w l

    Algo: 
    1. start with initial wavefunction psi(t)
    2. evaluate wavefunction at time psi(t+delta_t) = e^(-delta_t * H)psi(t)
    3. normalize psi
    4. go to 1
*/


//Includes
#include <iostream>
#include <cmath>
#include <iomanip>
#include <fstream>

using namespace std;


//GLobal constants
int GRID = 1000;
double TMAX = 8;
int STEPS = 240000;
double T_STEP = TMAX/STEPS;
double HBAR = 1;
double OMEGA = 1;
double MASS = 1;
double LOWER = -5;
double UPPER = 5;

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

double start_func(double x) {
    return (exp(-(x-1)*(x-1)));
}

double potential(double x){
    return (x*x*1/2);
}

double analytic_sol(double x) {
    return (0.7511255*exp(-x*x/2));
}

double approx_a(double f_2eps, double f_3eps, double eps){
    double b = (f_3eps-f_2eps)/(5*eps*eps);
    return (f_2eps-4*b*eps*eps);
}

//main function
int main(){

    //Program vars
    double* function = new double[GRID];
    double* arr_fder = new double[GRID];
    double* arr_sder = new double[GRID];
    double t0 = 0;
    double time = t0;
    double delta_t = TMAX/STEPS;
    bool printed = false;

    //File mangament
    ofstream myfile;
    ofstream omegafile;
    myfile.open ("./results/function.txt");
    omegafile.open ("./results/omega.txt");

    //init Function
    init(start_func, function, arr_fder, arr_sder);  
    norm(function);

    //Main Loop
    for (int i = 0; i < STEPS; i++)
    {
        //Step 1: multiply
        for (int j = 0; j < GRID; j++)
        {
            function[j] = exp(-delta_t/2*potential(get_x(j)))*function[j];
        }

        //Step 2: Kinetic energy operator
        for (int j = 0; j < GRID; j++)
        {
            function[j] = function[j] + HBAR*HBAR*delta_t/(2*MASS)*arr_sder[j];
        }
        
        //Step 3: multiply
        for (int j = 0; j < GRID; j++)
        {
            function[j] = exp(-delta_t/2*potential(get_x(j)))*function[j];
        }

        //Step 4: Norm and Plot omega_0
        double test_sum = norm(function);
        double omega_0 = (-1)*log(test_sum*test_sum)/(2*T_STEP);
        double soll = 0.5;
        double eps = omega_0-soll;
        if(i%100 == 0 ){
            double t = TMAX*i/STEPS;    
            omegafile << time << " "<< omega_0 << "\n";;    
        }

        if(i == STEPS-1) {/*  
            double Eres = 0;
            for (int i = 0; i < GRID; i++)
            {
                double delta_pos = (UPPER-LOWER)/double(GRID);
                Eres +=  (-function[i]*HBAR*HBAR/(2*MASS)*sder(i, GRID, function)+function[i]/2*get_x(i)*get_x(i)*function[i])*delta_pos;
            }*/
        }

        //Step 5: increase time
        time+=delta_t;

        //Step 5 update second derivative
        update_der(function, arr_fder, arr_sder);

        //Write function to Textfile
        if(i%10000==0){
            //Write to file
            for (int j = 0; j < GRID; j++)
            {
                myfile << get_x(j) << " " << time << " " << function[j] << "\n";
            }            
            myfile << "\n";
            myfile << "\n";
        }
    }

    //Close file
    myfile.close();
    omegafile.close();
    
}