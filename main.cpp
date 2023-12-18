/*
    Author: Samuel Pucher
*/

//Includes
#include <iostream>
#include <cmath>
#include <iomanip>
#include <fstream>
#include <string>
#include <math.h>

using namespace std;

//GLobal constants
const int GRID_R = 100;
const int GRID_P = 50;
double TMAX = 1;
double T_STEP = 0;//5*1e-6;
//int STEPS = TMAX/T_STEP;            //T_STEP should be arround 0.000001 = 1e-6     
double LOWER = 0;
double UPPER = 5;
double THRESHOLD = 0.01;            //not used yet, reserved for later

//Physical constants
double GAMMA = 1;
double ALPHA = 1;

//helper functions
double get_r(int i) {
    double res = (LOWER+(double(i)/GRID_R)*(UPPER-LOWER));
    return res;
}

double get_phi(int i) {
    double res = 2*M_PI*(double(i)/double(GRID_P));
    return res;
}

double fder_f_r0(int i, int j, int k, double arr[GRID_R][GRID_R][GRID_P]) {
    double eps = (UPPER-LOWER)/double(GRID_R);
    if(i < 1) {
        return ((arr[2][j][k]-arr[0][j][k])/(2*eps));
    }
    else if (i >=(GRID_R-1))
    {
        return ((arr[GRID_R-1][j][k]-arr[GRID_R-2][j][k])/(eps));
        //return ((arr[GRID_R-1][j][k]-arr[GRID_R-3][j][k])/(2*eps));
    }
    else {
        return ((arr[i+1][j][k]-arr[i-1][j][k])/(2*eps));
    }
}

double fder_f_rj(int i, int j, int k, double arr[GRID_R][GRID_R][GRID_P]) {
    double eps = (UPPER-LOWER)/double(GRID_R);
    if(j < 1) {    
        return ((arr[i][2][k]-arr[i][0][k])/(2*eps));
    }
    else if (j >=(GRID_R-1))
    {
        return ((arr[i][GRID_R-1][k]-arr[i][GRID_R-2][k])/(eps));
        //return ((arr[i][GRID_R-1][k]-arr[i][GRID_R-3][k])/(2*eps));
    }
    else {
        return ((arr[i][j+1][k]-arr[i][j-1][k])/(2*eps));
    } 
}

double fder_f_deltaphi(int i, int j, int k, double arr[GRID_R][GRID_R][GRID_P]) {
    double eps = (UPPER-LOWER)/double(GRID_P);
    if(k < 1) {
        return ((arr[i][j][2]-arr[i][j][0])/(2*eps));
    }
    else if (k >=(GRID_P-1))
    {
        return ((arr[i][j][GRID_P-1]-arr[i][j][GRID_P-3])/(2*eps));
    }
    else {
        return ((arr[i][j][k+1]-arr[i][j][k-1])/(2*eps));
    }
}

double sder_f_r0(int i, int j, int k, double arr[GRID_R][GRID_R][GRID_P]) {
    double eps = (UPPER-LOWER)/double(GRID_R);
    if(i < 1) {
        return ((arr[2][j][k]+arr[0][j][k]-2*arr[1][j][k])/(eps*eps));
    }
    else if (i >=(GRID_R-1))
    {
        return ((arr[GRID_R-1][j][k]+arr[GRID_R-3][j][k]-2*arr[GRID_R-2][j][k])/(eps*eps));
    }
    else {
        return ((arr[i+1][j][k]+arr[i-1][j][k]-2*arr[i][j][k])/(eps*eps));
    }
}

double sder_f_rj(int i, int j, int k, double arr[GRID_R][GRID_R][GRID_P]) {
    double eps = (UPPER-LOWER)/double(GRID_R);
    if(k < 1) {
        return ((arr[i][2][k]+arr[i][0][k]-2*arr[i][1][k])/(eps*eps));
    }
    else if (k >=(GRID_R-1))
    {
        return ((arr[i][GRID_R-1][k]+arr[i][GRID_R-3][k]-2*arr[i][GRID_R-2][k])/(eps*eps));
    }
    else {
        return ((arr[i][j+1][k]+arr[i][j-1][k]-2*arr[i][j][k])/(eps*eps));
    }
}

double sder_f_deltaphi(int i, int j, int k, double arr[GRID_R][GRID_R][GRID_P]) {
    double eps = (UPPER-LOWER)/double(GRID_P);
    if((k < 2 )|| (k>= (GRID_P-2))) {
        return ((arr[i][j][2]+arr[i][j][0]-2*arr[i][j][1])/(eps*eps));
    }
    else {
        return ((arr[i][j][k+1]+arr[i][j][k-1]-2*arr[i][j][k])/(eps*eps));
    }
    
}

double fder_g_r0(int i, double arr[GRID_R]) {
    double eps = (UPPER-LOWER)/double(GRID_R);
    if(i < 1) {
        return ((arr[2]-arr[0])/(2*eps));
    }
    else if (i >=(GRID_R-1))
    {
        return ((arr[GRID_R-1]-arr[GRID_R-3])/(2*eps));
    }
    else {
        return ((arr[i+1]-arr[i-1])/(2*eps));
    }
}

double sder_g_r0(int i, double arr[GRID_R]) {
    double eps = (UPPER-LOWER)/double(GRID_R);
    if(i < 1) {
        return (arr[2]+arr[0]-2.0*arr[1])/(eps*eps);
    }
    else if (i >=(GRID_R-1))
    {
        return (arr[GRID_R-1]+arr[GRID_R-3]-2*arr[GRID_R-2])/(eps*eps);
    }
    else {
        return (arr[i+1]+arr[i-1]-2*arr[i])/(eps*eps);
    }
}

void init_f(double (*func)(double, double, double), double a_func[GRID_R][GRID_R][GRID_P], double a_fder_r0[GRID_R][GRID_R][GRID_P], double a_fder_rj[GRID_R][GRID_R][GRID_P], double a_fder_deltaphi[GRID_R][GRID_R][GRID_P], double a_sder_r0[GRID_R][GRID_R][GRID_P], double a_sder_rj[GRID_R][GRID_R][GRID_P], double a_sder_deltaphi[GRID_R][GRID_R][GRID_P]){
    
    for (int i = 0; i < GRID_R; i++)
    {
        for (int j = 0; j < GRID_R; j++)
        {
            for (int k = 0; k < GRID_P; k++)
            {
                double r0 = get_r(i);
                double rj = get_r(j);
                double deltaphi = get_phi(k);
                a_func[i][j][k] = func(r0,rj,deltaphi);
            }
        }
    }
    for (int i = 0; i < GRID_R; i++)
    {
        for (int j = 0; j < GRID_R; j++)
        {
            for (int k = 0; k < GRID_P; k++)
            {
                a_fder_r0[i][j][k] = fder_f_r0(i,j,k,a_func);
                a_fder_rj[i][j][k] = fder_f_rj(i,j,k,a_func);
                a_fder_deltaphi[i][j][k] = fder_f_deltaphi(i,j,k, a_func);
                a_sder_r0[i][j][k] = sder_f_r0(i,j,k,a_func);
                a_sder_rj[i][j][k] = sder_f_rj(i,j,k,a_func);
                a_sder_deltaphi[i][j][k] = sder_f_deltaphi(i,j,k, a_func);
            }
        }
    }
}

void init_g(double (*func)(double), double a_func[GRID_R], double a_fder_r0[GRID_R], double a_sder_r0[GRID_R]){
    for (int i = 0; i < GRID_R; i++)
    {
        double r0 = get_r(i);
        a_func[i] = func(r0);
    }
    for (int i = 0; i < GRID_R; i++)
    {
        a_fder_r0[i] = fder_g_r0(i, a_func);
        a_sder_r0[i] = sder_g_r0(i, a_func);

    }
}

void init_U(double (*func)(double), double a_potential_U[GRID_R][GRID_R][GRID_P]){
    double r0, rj, deltaphi, dist;
    for (int i = 0; i < GRID_R; i++)
    {
        r0 = get_r(i);
        for (int j = 0; j < GRID_R; j++)
        {
            rj= get_r(j);
            for (int k = 0; k < GRID_P; k++)
            {
                deltaphi = get_phi(k);
                dist = sqrt(r0*r0+rj*rj-2*r0*rj*cos(deltaphi));
                a_potential_U[i][j][k] = func(dist);
            }            
        }        
    }    
}

void update_der_f(double a_func[GRID_R][GRID_R][GRID_P], double a_fder_r0[GRID_R][GRID_R][GRID_P], double a_fder_rj[GRID_R][GRID_R][GRID_P], double a_fder_deltaphi[GRID_R][GRID_R][GRID_P], double a_sder_r0[GRID_R][GRID_R][GRID_P], double a_sder_rj[GRID_R][GRID_R][GRID_P], double a_sder_deltaphi[GRID_R][GRID_R][GRID_P]){
    for (int i = 0; i < GRID_R; i++)
    {
        for (int j = 0; j < GRID_R; j++)
        {
            for (int k = 0; k < GRID_P; k++)
            {
                a_fder_r0[i][j][k] = fder_f_r0(i,j,k,a_func);
                a_fder_rj[i][j][k] = fder_f_rj(i,j,k,a_func);
                a_fder_deltaphi[i][j][k] = fder_f_deltaphi(i,j,k, a_func);
                a_sder_r0[i][j][k] = sder_f_r0(i,j,k,a_func);
                a_sder_rj[i][j][k] = sder_f_rj(i,j,k,a_func);
                a_sder_deltaphi[i][j][k] = sder_f_deltaphi(i,j,k, a_func);
            }
        }
    }
}

void update_der_g(double a_func[GRID_R], double a_fder_r0[GRID_R], double a_sder_r0[GRID_R]){
    for (int i = 0; i < GRID_R; i++)
    {
        a_fder_r0[i] = fder_g_r0(i, a_func);
        a_sder_r0[i] = sder_g_r0(i, a_func);
    }
}

double norm_f(double a_func[GRID_R][GRID_R][GRID_P]){

    double sum = 0;
    double eps_p = 2*M_PI/(double(GRID_P));
    int phi_half = int(GRID_P/2);
    //double f_tilde = a_func[GRID_R-1][GRID_R-1][phi_half];

    //std::cout << "Before normation: " << a_func[50][GRID_R-1][phi_half]; 
    
    /*
    for (int i = 0; i < GRID_R; i++)
    {
        for (int j = 0; j < GRID_R; j++)
        { 
            for (int k = 0; k < GRID_P; k++)
            {
                a_func[i][j][k] = a_func[i][j][k]/f_tilde; 
            }            
        }
    }
    std::cout << " , after normation: " << a_func[50][GRID_R-1][phi_half] << std::endl;
    return f_tilde;

    */

    for (int i = 0; i < GRID_P; i++)
    {
        sum+= a_func[0][GRID_R-1][i];            
    }
    
    double f_tilde = sum*eps_p/(2*M_PI);

    for (int i = 0; i < GRID_R; i++)
    {
        for (int j = 0; j < GRID_R; j++)
        {
            for (int k = 0; k < GRID_P; k++)
            {
                a_func[i][j][k] = a_func[i][j][k]/f_tilde; 
            }
        }
    }
    return f_tilde;

}

double norm_g(double a_func[GRID_R], bool norm_arr){

    double sum = 0;
    double eps_r0 = (UPPER-LOWER)/double(GRID_R);

    for (int i = 0; i < GRID_R; i++)
    {
        sum+= a_func[i]*a_func[i]*get_r(i);
    }
    sum = sum* 2*eps_r0*M_PI;

    if(norm_arr){
        for (int i = 0; i < GRID_R; i++)
            {
                a_func[i] = a_func[i]/sqrt(sum);
            }   
    }
   return sqrt(sum);
}

double start_func_g(double r0) {
    return (exp(-(r0)*(r0)/1.0)+0.2);
}

double start_func_f(double r0, double rj, double deltaphi) {
    return 1.0;
}

double potential_U(double dist){
    return (exp(-(dist)*(dist)/0.3*0.3));
    //return 0.0;
}

double calc_Vg(double a_func_f[GRID_R][GRID_R][GRID_P], double a_potential_U[GRID_R][GRID_R][GRID_P], int i) {    //i corresponds to the position r0
    double sum = 0;
    double r0 = get_r(i);
    double eps_rj = (UPPER-LOWER)/double(GRID_R);
    double eps_p = (2*M_PI)/double(GRID_P);
    //Helping vars
    double fd_f_r0, fd_f_phi, fd_f_rj, part_1, part_2;
    for (int j = 1; j < GRID_R; j++)
    {
        double rj = get_r(j);
        for (int k = 0; k < GRID_P; k++)        
        {       
            fd_f_r0 = fder_f_r0(i,j,k, a_func_f);
            fd_f_rj = fder_f_rj(i,j,k, a_func_f);
            fd_f_phi = fder_f_deltaphi(i,j,k, a_func_f);

            part_1 = ALPHA*pow(fd_f_r0,2.0)/2+pow(fd_f_rj,2.0)/2+pow(fd_f_phi,2.0)*(ALPHA/(2*r0*r0) + 1/(2*rj*rj));                         //responsible for divergence of Vg by r0 = 0!!
            part_2 = pow(a_func_f[i][j][k], 2.0)*a_potential_U[i][j][k]+(pow(a_func_f[i][j][k], 4.0)-2*pow(a_func_f[i][j][k], 2.0)+1.0)/2;      //responsible for (exponential?) time increase of Vg!!

            sum += (part_1 + part_2)*rj;

            //sum = sum + (ALPHA/2*fder_f_r0(i,j,k, a_func_f)*fder_f_r0(i,j,k, a_func_f)+ 0.5*fder_f_rj(i,j,k, a_func_f)*fder_f_rj(i,j,k, a_func_f)+fder_f_deltaphi(i,j,k, a_func_f)*fder_f_deltaphi(i,j,k, a_func_f)*(ALPHA/(2*r0*r0)+ 1/(2*rj*rj))+a_func_f[i][j][k]*a_func_f[i][j][k]*a_potential_U[i][j][k]+0.5*(a_func_f[i][j][k]*a_func_f[i][j][k]*a_func_f[i][j][k]*a_func_f[i][j][k]-2*a_func_f[i][j][k]*a_func_f[i][j][k]+1))*rj;      
        }
    }
    sum = sum*eps_rj*eps_p/GAMMA;
    return sum;
}

double approx_a(double f_1eps, double f_2eps){
    return (4*(f_1eps-f_2eps/4)/3);
}

int main(){
        
    /*
    {
    int argc, char* argv[]
    std::cout << "Beta: " << argv[1] << std::endl;    
    double BETA = (int)argv[1];
    std::string help_path = argv[1];
    std::string path = "./results/final_chi/function_chi_final_" + help_path +".txt";
    std::cout << "Path: " << path << std::endl;   
    }
    */
   
    
    //g-Funtion: initialization g[r0]
    static double arr_function_g[GRID_R];
    static double arr_function_g_old[GRID_R];
    static double arr_fder_g_r0[GRID_R];
    static double arr_sder_g_r0[GRID_R];    

    //f-Function: initialization f[r0][rj][deltaphi]
    static double arr_function_f[GRID_R][GRID_R][GRID_P];
    static double arr_function_f_old[GRID_R][GRID_R][GRID_P];
    static double arr_fder_f_r0[GRID_R][GRID_R][GRID_P];
    static double arr_fder_f_rj[GRID_R][GRID_R][GRID_P];
    static double arr_fder_f_deltaphi[GRID_R][GRID_R][GRID_P];
    static double arr_sder_f_r0[GRID_R][GRID_R][GRID_P];
    static double arr_sder_f_rj[GRID_R][GRID_R][GRID_P];
    static double arr_sder_f_deltaphi[GRID_R][GRID_R][GRID_P];

    //ftilde-Function: initialization f[r0][rj][deltaphi]
    static double arr_function_ftilde[GRID_R][GRID_R][GRID_P];
    static double arr_function_ftilde_old[GRID_R][GRID_R][GRID_P];
    static double arr_fder_ftilde_r0[GRID_R][GRID_R][GRID_P];
    static double arr_fder_ftilde_rj[GRID_R][GRID_R][GRID_P];
    static double arr_fder_ftilde_deltaphi[GRID_R][GRID_R][GRID_P];
    static double arr_sder_ftilde_r0[GRID_R][GRID_R][GRID_P];
    static double arr_sder_ftilde_rj[GRID_R][GRID_R][GRID_P];
    static double arr_sder_ftilde_deltaphi[GRID_R][GRID_R][GRID_P];

    //physical potentials
    static double arr_potential_U[GRID_R][GRID_R][GRID_P];
    static double arr_Vg[GRID_R];
    static double arr_Vf[GRID_R][GRID_R][GRID_P];

    //misc variables
    double time = 0;
    int i_Vg_max = int(GRID_R*1/2);
    int phi_half = int(GRID_P/2);

    //File mangament
    ofstream file_ft;
    ofstream file_g;
    ofstream file_f;
    file_ft.open ("./results/Impurity-BEC/littleBox/");
    file_g.open ("./results/Impurity-BEC/littleBox/function_g.txt");
    file_f.open ("./results/Impurity-BEC/littleBox/function_f.txt");

    //init Functions
    init_g(start_func_g, arr_function_g, arr_fder_g_r0, arr_sder_g_r0); 
    norm_g(arr_function_g, true);
    init_f(start_func_f, arr_function_f, arr_fder_f_r0, arr_fder_f_rj, arr_fder_f_deltaphi, arr_sder_f_r0, arr_sder_f_rj, arr_sder_f_deltaphi);
    norm_f(arr_function_f);
    init_f(start_func_f, arr_function_ftilde, arr_fder_ftilde_r0, arr_fder_ftilde_rj, arr_fder_ftilde_deltaphi, arr_sder_ftilde_r0, arr_sder_ftilde_rj, arr_sder_ftilde_deltaphi);
    norm_f(arr_function_ftilde);
    init_U(potential_U, arr_potential_U);
    
    //Main Loop
    for (int t = 0; t < 1; t++)
    {        
        //Step 1 and Step 2                  
        for (int i = 1; i < GRID_R; i++) 
        {
            //arr_Vg[i] = calc_Vg(arr_function_f, arr_potential_U, i);
            if(i <= i_Vg_max) {
                arr_Vg[i] = calc_Vg(arr_function_f, arr_potential_U, i);
            } else {
                arr_Vg[i] = arr_Vg[i_Vg_max];
            }        
            for (int j = 1; j < GRID_R; j++)
            {
                for (int  k = 0; k < GRID_P; k++)
                {
                    arr_Vf[i][j][k] = arr_Vg[i]+ arr_potential_U[i][j][k]+ arr_function_f[i][j][k]*arr_function_f[i][j][k];
                    arr_function_ftilde[i][j][k] = arr_function_f[i][j][k]* arr_function_g[i];
                }             
            }
        }
        arr_Vg[0] = arr_Vg[1];
        
        //Step 3
        for (int i = 1; i < GRID_R; i++)
        {
            arr_function_g[i] = exp(-T_STEP*arr_Vg[i]/2)*arr_function_g[i];
        }      
        
        //Step 4
        for (int i = 1; i < GRID_R; i++)
        {
            arr_function_g[i] = arr_function_g[i] + ALPHA/2 * T_STEP * (arr_sder_g_r0[i] +arr_fder_g_r0[i]/get_r(i));
        }
        
        //Step 5
        for (int i = 1; i < GRID_R; i++)
        {
            arr_function_g[i] = exp(-T_STEP*arr_Vg[i]/2)*arr_function_g[i];
        }

        arr_function_g[0] = approx_a(arr_function_g[1], arr_function_g[2]);
        
        //Step 6
        norm_g(arr_function_g, true);    

        //Step 7
        for(int i = 1; i < GRID_R; i++)
        {
            for (int j = 1; j < GRID_R; j++)
            {
                for (int  k = 0; k < GRID_P; k++)
                {
                    arr_function_ftilde[i][j][k] = exp(-arr_Vf[i][j][k]*T_STEP/2)*arr_function_ftilde[i][j][k];
                }
            }            
        }       

        //Step 8
        for (int i = 1; i < GRID_R; i++)
        {
            for (int j = 1; j < GRID_R; j++)
            {
                for (int k= 0; k < GRID_P; k++)
                {
                    arr_function_ftilde[i][j][k] = arr_function_ftilde[i][j][k] + T_STEP/2 * (ALPHA*(arr_sder_ftilde_r0[i][j][k] + arr_fder_ftilde_r0[i][j][k]/get_r(i) + arr_sder_ftilde_rj[i][j][k]*1/(get_r(i)*get_r(i))) + arr_sder_ftilde_rj[i][j][k] + arr_fder_ftilde_rj[i][j][k]/get_r(j) + arr_sder_ftilde_deltaphi[i][j][k]*1/(get_r(j)*get_r(j)));
                }
            }
        }                

        //Step 9
        for (int i = 1; i < GRID_R; i++)
        {
            for (int j = 1; j < GRID_R; j++)
            {
                for (int  k = 0; k < GRID_P; k++)
                {
                    arr_function_ftilde[i][j][k] = exp(-arr_Vf[i][j][k]*T_STEP/2)*arr_function_ftilde[i][j][k];
                }
            }   
        }

        for (int k = 0; k < GRID_P; k++)
        {       
            for (int i = 1; i < GRID_R; i++)
            {
                arr_function_ftilde[i][0][k] = approx_a(arr_function_ftilde[i][1][k], arr_function_ftilde[i][2][k]);
                //arr_function_ftilde[i][1][k] = 0.5*(arr_function_ftilde[i][2][k] + arr_function_ftilde[i][0][k]);
            }            
            for (int j = 1; j < GRID_R; j++)
            {
                arr_function_ftilde[0][j][k] = approx_a(arr_function_ftilde[1][j][k], arr_function_ftilde[2][j][k]);
                //arr_function_ftilde[1][j][k] = 0.5*(arr_function_ftilde[0][j][k] + arr_function_ftilde[2][j][k]);                
            }
            arr_function_ftilde[0][0][k] = 0.5*(arr_function_ftilde[0][1][k]+ arr_function_ftilde[1][0][k]);
            //arr_function_ftilde[0][0][k] = arr_function_ftilde[1][1][k] = 0.5*(arr_function_ftilde[2][1][k]+ arr_function_ftilde[1][2][k]);
            //arr_function_ftilde[1][0][k] = 0.5*(arr_function_ftilde[0][0][k] + arr_function_ftilde[2][0][k]); 
            //arr_function_ftilde[0][1][k] = 0.5*(arr_function_ftilde[0][0][k] + arr_function_ftilde[0][2][k]); 
        }

        //Step 10
        for (int i = 0; i < GRID_R; i++)
        {
            for (int j = 0; j < GRID_R; j++)
            {
                for (int  k = 0; k < GRID_P; k++)
                {
                    arr_function_f[i][j][k] = arr_function_ftilde[i][j][k]/arr_function_g[i];
                }
            }   
        }

        //Step 11
        double f_tilde = norm_f(arr_function_f);
        update_der_f(arr_function_f, arr_fder_f_r0, arr_fder_f_rj, arr_fder_f_deltaphi, arr_sder_f_r0, arr_sder_f_rj, arr_sder_f_deltaphi);
        update_der_f(arr_function_ftilde, arr_fder_ftilde_r0, arr_fder_ftilde_rj, arr_fder_ftilde_deltaphi, arr_sder_ftilde_r0, arr_sder_ftilde_rj, arr_sder_ftilde_deltaphi);
        update_der_g(arr_function_g, arr_fder_g_r0, arr_sder_g_r0);

        //increae time
        time+=T_STEP;

        //Write to file
        if(t%1 == 0) {
            std::cout << t << std::endl;
            for (int i = 0; i < GRID_R; i++)
            {
                file_g << get_r(i) << " " << time << " " << arr_function_g[i] << " " << arr_Vg[i] << " " << arr_potential_U[i][0][25] << "\n"; 
                file_f << get_r(i) << " " << time << " " << arr_function_f[i][GRID_R-1][phi_half] << " " << arr_Vf[GRID_R-1][GRID_R-1][i]<<"\n";      //Sollte 1 für große i und klein für kleine i                               
            }            
            file_g << "\n"  << std::flush;
            file_f << "\n"  << std::flush;
        }

        /*Check if function diverges
        if isnan(arr_function_g[10]) {
            std::cout << "Itteration stop, function g is not a number!" 
            return 0; 
        }*/

    }

    //Close files
    file_f.close();
    file_g.close();
    file_ft.close();
   
}