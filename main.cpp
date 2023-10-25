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
const int GRID_P = 100;
double TMAX = 1;
double T_STEP = 1e-6;
int STEPS = TMAX/T_STEP;            //T_STEP should be arround 0.000001 = 1e-6        
double LOWER = 0;
double UPPER = 10;
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
    if(i < 2) {
        return ((arr[2][j][k]-arr[0][j][k])/(2*eps));
    }
    else if (i >=(GRID_R-1))
    {
        return 0;
    }
    else if(i==j) {
        return 0;
    }
    else {
        return ((arr[i+1][j][k]-arr[i-1][j][k])/(2*eps));
    }
}

double fder_f_rj(int i, int j, int k, double arr[GRID_R][GRID_R][GRID_P]) {
    double eps = (UPPER-LOWER)/double(GRID_R);
    if(j < 2) {    
        return ((arr[i][2][k]-arr[i][0][k])/(2*eps));
    }
    else if (j >=(GRID_R-1))
    {
        return 0;
    }
    else if(i==j) {
        return 0;
    }
    else {
        return ((arr[i][j+1][k]-arr[i][j-1][k])/(2*eps));
    } 
}

double fder_f_deltaphi(int i, int j, int k, double arr[GRID_R][GRID_R][GRID_P]) {
    double eps = (UPPER-LOWER)/double(GRID_P);
    if(k < 2) {
        return ((arr[i][j][2]-arr[i][j][0])/(2*eps));
    }
    else if (k >=(GRID_P-1))
    {
        return ((arr[i][j][GRID_P-1]-arr[i][j][GRID_P-3])/(2*eps));
    }
    else {
        return ((arr[i][j][k+1]-arr[i][j][k]-1)/(2*eps));
    }
}

double sder_f_r0(int i, int j, int k, double arr[GRID_R][GRID_R][GRID_P]) {
    double eps = (UPPER-LOWER)/double(GRID_R);
    if(i < 2) {
        return ((arr[2][j][k]+arr[0][j][k]-2*arr[1][j][k])/(eps*eps));
    }
    else if (i >=(GRID_R-1))
    {
        return 0;
    }
    else {
        return ((arr[i+1][j][k]+arr[i-1][j][k]-2*arr[i][j][k])/(eps*eps));
    }
}

double sder_f_rj(int i, int j, int k, double arr[GRID_R][GRID_R][GRID_P]) {
    double eps = (UPPER-LOWER)/double(GRID_R);
    if(k < 2) {
        return ((arr[i][2][k]+arr[i][0][k]-2*arr[i][1][k])/(eps*eps));
    }
    else if (k >=(GRID_R-1))
    {
        return 0;
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
    if(i < 2) {
        return ((arr[2]-arr[0])/(2*eps));
    }
    else if (i >=(GRID_R-1))
    {
        return 0;
    }
    else {
        return ((arr[i+1]-arr[i-1])/(2*eps));
    }
}

double sder_g_r0(int i, double arr[GRID_R]) {
    double eps = (UPPER-LOWER)/double(GRID_R);
    if(i < 2) {
        return (arr[2]+arr[0]-2.0*arr[1])/(eps*eps);
    }
    else if (i >=(GRID_R-1))
    {
        return 0;
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
    
   return sqrt(sum);
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
    return (exp(-(r0)*(r0))+0.2);
}

double start_func_f(double r0, double rj, double deltaphi) {
    return 1.0;
}

double potential_U(double dist){
    return (exp(-(dist)*(dist))); 
}

double calc_Vg(double a_func_f[GRID_R][GRID_R][GRID_P], double a_potential_U[GRID_R][GRID_R][GRID_P], int i) {    //i corresponds to the position r0
    double sum = 0;
    double r0 = get_r(i);
    double eps_rj = (UPPER-LOWER)/double(GRID_R);
    double eps_p = (UPPER-LOWER)/double(GRID_P);
    for (int j = 2; j < GRID_R; j++)
    {
        double rj = get_r(j);
        for (int k = 0; k < GRID_P; k++)        // i corresponds to rj and j to deltaphi 
        {            
            sum = ALPHA/2*fder_f_r0(i,j,k, a_func_f)*fder_f_r0(i,j,k, a_func_f)+ 0.5*fder_f_rj(i,j,k, a_func_f)*fder_f_rj(i,j,k, a_func_f)+fder_f_deltaphi(i,j,k, a_func_f)*fder_f_deltaphi(i,j,k, a_func_f)*(ALPHA/(2*r0*r0)+ 1/(2*rj*rj))+a_func_f[i][j][k]*a_func_f[i][j][k]*a_potential_U[i][j][k]+0.5*(a_func_f[i][j][k]*a_func_f[i][j][k]*a_func_f[i][j][k]*a_func_f[i][j][k]-2*a_func_f[i][j][k]*a_func_f[i][j][k]+1);      
        }
    }
    sum *= eps_rj*eps_p/GAMMA;
    return sum;
}

double approx_a(double f_2eps, double f_3eps){
    double eps = (UPPER-LOWER)/double(GRID_R);
    double b = (f_3eps-f_2eps)/(5*eps*eps);
    return (f_2eps-4*b*eps*eps);
}



int main(){
    
    //NORM_G DOES NOT WORK!!!!!
    
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

    //Physical variables
    double time = 0;

    //File mangament
    ofstream file_ft;
    ofstream file_g;
    ofstream file_f;
    file_ft.open ("./results/Impurity-BEC/21.10/");
    file_g.open ("./results/Impurity-BEC/21.10/function_g.txt");
    file_f.open ("./results/Impurity-BEC/21.10/");

    //init Functions
    init_g(start_func_g, arr_function_g, arr_fder_g_r0, arr_sder_g_r0); 
    norm_g(arr_function_g, true);
    init_f(start_func_f, arr_function_f, arr_fder_f_r0, arr_fder_f_rj, arr_fder_f_deltaphi, arr_sder_f_r0, arr_sder_f_rj, arr_sder_f_deltaphi);
    norm_f(arr_function_f);
    init_f(start_func_f, arr_function_ftilde, arr_fder_ftilde_r0, arr_fder_ftilde_rj, arr_fder_ftilde_deltaphi, arr_sder_ftilde_r0, arr_sder_ftilde_rj, arr_sder_ftilde_deltaphi);
    norm_f(arr_function_ftilde);
    init_U(potential_U, arr_potential_U);
    
          

    //Main Loop
    for (int t = 0; t < STEPS; t++)
    {
        //Step 1 and Step 2                                 TODO: APPROX FUNCTION FOR LOW INDICES!!!!!!!!
        for (int i = 2; i < GRID_R; i++) 
        {
            arr_Vg[i] = calc_Vg(arr_function_f, arr_potential_U, i);        
            for (int j = 2; j < GRID_R; j++)
            {
                for (int  k = 0; k < GRID_P; k++)
                {
                    arr_Vf[i][j][k] = arr_Vg[i]+ arr_potential_U[i][j][k]+ arr_function_f[i][j][k]*arr_function_f[i][j][k];
                    arr_function_ftilde[i][j][k] = arr_function_f[i][j][k]* arr_function_g[i];
                }                    
            }
        }
        
        //Step 3
        for (int i = 2; i < GRID_R; i++)
        {
            arr_function_g[i] = exp(-T_STEP*arr_Vg[i]/2)*arr_function_g[i];
        }      
        
        
        //Step 4
        for (int i = 2; i < GRID_R; i++)
        {
            arr_function_g[i] = arr_function_g[i] + ALPHA/2 * T_STEP*arr_sder_g_r0[i];
        }
        
        //Step 5
        for (int i = 2; i < GRID_R; i++)
        {
            arr_function_g[i] = exp(-T_STEP*arr_Vg[i]/2)*arr_function_g[i];
        }

        arr_function_g[0] = approx_a(arr_function_g[2], arr_function_g[3]);
        arr_function_g[1] = 0.5*(arr_function_g[0]+arr_function_g[2]);
        
        //Step 6
        norm_g(arr_function_g, true);   

        //cout << arr_function_g[2] << endl; 

        //Step 7
        for (int i = 2; i < GRID_R; i++)
        {
            for (int j = 2; j < GRID_R; j++)
            {
                for (int  k = 0; k < GRID_P; k++)
                {
                    arr_function_ftilde[i][j][k] = exp(-arr_Vf[i][j][k]*T_STEP/2)*arr_function_ftilde[i][j][k];
                }
            }            
        }       
        

        //Step 8
        for (int i = 2; i < GRID_R; i++)
        {
            for (int j = 2; j < GRID_R; j++)
            {
                for (int k= 0; k < GRID_P; k++)
                {
                    arr_function_ftilde[i][j][k] = arr_function_ftilde[i][j][k] + ALPHA/2 * arr_sder_ftilde_r0[i][j][k] + 0.5 * arr_sder_ftilde_rj[i][j][k];
                }
            }
        }                

        //Step 9
        for (int i = 2; i < GRID_R; i++)
        {
            for (int j = 2; j < GRID_R; j++)
            {
                for (int  k = 0; k < GRID_P; k++)
                {
                    arr_function_ftilde[i][j][k] = exp(-arr_Vf[i][j][k]*T_STEP/2)*arr_function_ftilde[i][j][k];
                }
            }   
        }

        for (int k = 0; k < GRID_P; k++)
        {       
            for (int i = 2; i < GRID_R; i++)
            {
                arr_function_ftilde[i][0][k] = approx_a(arr_function_ftilde[i][2][k], arr_function_ftilde[i][3][k]);
                arr_function_ftilde[i][1][k] = 0.5*(arr_function_ftilde[i][2][k] + arr_function_ftilde[i][0][k]);
            }            
            for (int j = 2; j < GRID_R; j++)
            {
                arr_function_ftilde[0][j][k] = approx_a(arr_function_ftilde[2][j][k], arr_function_ftilde[3][j][k]);
                arr_function_ftilde[1][j][k] = 0.5*(arr_function_ftilde[0][j][k] + arr_function_ftilde[2][j][k]);                
            }
            arr_function_ftilde[0][0][k] = arr_function_ftilde[1][1][k] = 0.5*(arr_function_ftilde[2][1][k]+ arr_function_ftilde[1][2][k]);
            arr_function_ftilde[1][0][k] = 0.5*(arr_function_ftilde[0][0][k] + arr_function_ftilde[2][0][k]); 
            arr_function_ftilde[0][1][k] = 0.5*(arr_function_ftilde[0][0][k] + arr_function_ftilde[0][2][k]); 
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
        norm_f(arr_function_f);
        update_der_f(arr_function_f, arr_fder_f_r0, arr_fder_f_rj, arr_fder_f_deltaphi, arr_sder_f_r0, arr_sder_f_rj, arr_sder_f_deltaphi);
        update_der_f(arr_function_ftilde, arr_fder_ftilde_r0, arr_fder_ftilde_rj, arr_fder_ftilde_deltaphi, arr_sder_ftilde_r0, arr_sder_ftilde_rj, arr_sder_ftilde_deltaphi);
        update_der_g(arr_function_g, arr_fder_g_r0, arr_sder_g_r0);

        //Write to file
        if(t%10000 == 0) {
            for (int i = 0; i < GRID_R; i++)
            {
                file_g << get_r(i) << " " << time << " " << arr_function_g[i] << "\n";           
            }            
            file_g << "\n";
        }

        //increae time
        time+=T_STEP;

        cout << time << endl;
    }

    //Close files
    file_f.close();
    file_g.close();
    file_g.close();
   
}