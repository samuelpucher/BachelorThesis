
//PARAMS: SIGMA = 1.0 , U_0 = 1

/*
    Author:         Samuel Pucher
    Created:        21.September
    Disc:           Created to see how Box of length +/-10 affects the problem
    Compile with:   g++ main.cpp -lfftw3 -lfftw3_threads -mcmodel=medium -O2 &
    
*/

//Includes
#include <cmath>
#include <iostream>
#include <fstream>
#include <fftw3.h>

using namespace std;

//Spacing Constants
const int GRID = 64;
const double UPPER = 7.5;
const double LOWER = -7.5;
const double EPS = (UPPER-LOWER)/GRID;

//Change
double SIGMA = 2.0;
double U_0 = 0.0625;

//Physical constants
double GAMMA = 0.5;
double ALPHA = 1.0;         //curios: somehow 'double Alpha = 1;' and setting manally to 1 is not the same. Maybe a compiler 'optimization' error?
double chem_I = 0.0;
double chem_B = 0.0;
double T_STEP = 0.01;
double T_MAX = 150;
double STEPS = T_MAX/T_STEP;
double delta_k = 2.0*M_PI/(GRID*EPS);
double current_time = 0;
double TI, TB, k1x, k1y, k0x, k0y;

//Other constants
const int Rank4 = 4;
const int Rank2 = 2;
double pre_norm_f;
double pre_norm_g;
double prop_norm_f;
double prop_norm_g;


//other helper functions
double get_x(int i) {
    double res = (LOWER+(double(i)/GRID)*(UPPER-LOWER));
    return res;
}

//Derivative Functions, they use only the real part of the function!!!
double fder_f_x0(int i, int j, int k, int l, fftw_complex arr[GRID][GRID][GRID][GRID]){
    if (i<1) {
        return ((arr[1][j][k][l][0]-arr[GRID-1][j][k][l][0])/(2.0*EPS));
    } 
    else if(i>=(GRID-1)) {
        return ((arr[0][j][k][l][0]-arr[GRID-2][j][k][l][0])/(2.0*EPS));
    }
    else {
        return ((arr[i+1][j][k][l][0]-arr[i-1][j][k][l][0])/(2.0*EPS));
    }
}

double fder_f_y0(int i, int j, int k, int l, fftw_complex arr[GRID][GRID][GRID][GRID]){
    if (j<1) {
        return ((arr[i][1][k][l][0]-arr[i][GRID-1][k][l][0])/(2.0*EPS));
    } 
    else if(j>=(GRID-1)) {
        return ((arr[i][0][k][l][0]-arr[i][GRID-2][k][l][0])/(2.0*EPS));
    }
    else {
        return ((arr[i][j+1][k][l][0]-arr[i][j-1][k][l][0])/(2.0*EPS));
    }
}

double fder_f_x1(int i, int j, int k, int l, fftw_complex arr[GRID][GRID][GRID][GRID]){
    if (k<1) {
        return ((arr[i][j][1][l][0]-arr[i][j][GRID-1][l][0])/(2.0*EPS));
    } 
    else if(k>=(GRID-1)) {
        return ((arr[i][j][0][l][0]-arr[i][j][GRID-2][l][0])/(2.0*EPS));
    }
    else {
        return ((arr[i][j][k+1][l][0]-arr[i][j][k-1][l][0])/(2.0*EPS));
    }
}

double fder_f_y1(int i, int j, int k, int l, fftw_complex arr[GRID][GRID][GRID][GRID]){
    if (l<1) {
        return ((arr[i][j][k][1][0]-arr[i][j][k][GRID-1][0])/(2.0*EPS));
    } 
    else if(l>=(GRID-1)) {
        return ((arr[i][j][k][0][0]-arr[i][j][k][GRID-2][0])/(2.0*EPS));
    }
    else {
        return ((arr[i][j][k][l+1][0]-arr[i][j][k][l-1][0])/(2.0*EPS));
    }
}

double sder_f_x0(int i, int j, int k, int l, double arr[GRID][GRID][GRID][GRID]){
    if (i<1) {
        return ((arr[1][j][k][l]+arr[GRID-1][j][k][l]-2*arr[0][j][k][l])/(EPS*EPS));
    } 
    else if(i>=(GRID-1)) {
        return ((arr[0][j][k][l]+arr[GRID-2][j][k][l]-2*arr[GRID-1][j][k][l])/(EPS*EPS));
    }
    else {
        return ((arr[i+1][j][k][l]+arr[i-1][j][k][l]-2*arr[i][j][k][l])/(EPS*EPS));
    }
}

double sder_f_y0(int i, int j, int k, int l, double arr[GRID][GRID][GRID][GRID]){
    if (j<1) {
        return ((arr[i][1][k][l]+arr[i][GRID-1][k][l]-2*arr[i][0][k][l])/(EPS*EPS));
    } 
    else if(j>=(GRID-1)) {
        return ((arr[i][0][k][l]+arr[i][GRID-2][k][l]-2*arr[i][GRID-1][k][l])/(EPS*EPS));
    }
    else {
        return ((arr[i][j+1][k][l]+arr[i][j-1][k][l]-2*arr[i][j][k][l])/(EPS*EPS));
    }
}

double sder_f_x1(int i, int j, int k, int l, double arr[GRID][GRID][GRID][GRID]){
    if (k<1) {
        return ((arr[i][j][1][l]+arr[i][j][GRID-1][l]-2*arr[i][j][0][l])/(EPS*EPS));
    } 
    else if(k>=(GRID-1)) {
        return ((arr[i][j][0][l]+arr[i][j][GRID-2][l]-2*arr[i][j][GRID-1][l])/(EPS*EPS));
    }
    else {
        return ((arr[i][j][k+1][l]+arr[i][j][k-1][l]-2*arr[i][j][k][l])/(EPS*EPS));
    }
}

double sder_f_y1(int i, int j, int k, int l, double arr[GRID][GRID][GRID][GRID]){
    if (l<1) {
        return ((arr[i][j][k][1]+arr[i][j][k][GRID-1]-2*arr[i][j][k][0])/(EPS*EPS));
    } 
    else if(l>=(GRID-1)) {
        return ((arr[i][j][k][0]+arr[i][j][k][GRID-2]-2*arr[i][j][k][GRID-1])/(EPS*EPS));
    }
    else {
        return ((arr[i][j][k][l+1]+arr[i][j][k][l-1]-2*arr[i][j][k][l])/(EPS*EPS));
    }
}

double fder_g_x0(int i, int j, double arr[GRID][GRID]){
    if (i<1) {
        return ((arr[1][j]-arr[GRID-1][j])/(2.0*EPS));
    } 
    else if(i>=(GRID-1)) {
        return ((arr[0][j]-arr[GRID-2][j])/(2.0*EPS));
    }
    else {
        return ((arr[i+1][j]-arr[i-1][j])/(2.0*EPS));
    }
}

double fder_g_y0(int i, int j, double arr[GRID][GRID]){
    if (j<1) {
        return ((arr[i][1]-arr[i][GRID-1])/(2.0*EPS));
    } 
    else if(j>=(GRID-1)) {
        return ((arr[i][0]-arr[i][GRID-2])/(2.0*EPS));
    }
    else {
        return ((arr[i][j+1]-arr[i][j-1])/(2.0*EPS));
    }
}

double sder_g_x0(int i, int j, double arr[GRID][GRID]){
    if (i<1) {
        return ((arr[1][j]+arr[GRID-1][j]-2*arr[0][j])/(EPS*EPS));
    } 
    else if(i>=(GRID-1)) {
        return ((arr[0][j]+arr[GRID-2][j]-2*arr[GRID-1][j])/(EPS*EPS));
    }
    else {
        return ((arr[i+1][j]+arr[i-1][j]-2*arr[i][j])/(EPS*EPS));
    }
}

double sder_g_y0(int i, int j, double arr[GRID][GRID]){
    if (j<1) {
        return ((arr[i][1]+arr[i][GRID-1]-2*arr[i][0])/(EPS*EPS));
    } 
    else if(j>=(GRID-1)) {
        return ((arr[i][0]+arr[i][GRID-2]-2*arr[i][GRID-1])/(EPS*EPS));
    }
    else {
        return ((arr[i][j+1]+arr[i][j-1]-2*arr[i][j])/(EPS*EPS));
    }
}

//Special functions
double start_func_g(double x0, double y0) {
    return (exp(-(x0*x0+y0*y0)/1.0)+0.2);
}

double start_func_f(double x0, double y0, double x1, double y1) {
    return 1.0;
}

double potential_U(double dist){
    return (U_0*exp(-dist*dist/(2*SIGMA*SIGMA)));
}

double calc_Vg(fftw_complex a_func_f[GRID][GRID][GRID][GRID], double a_potential_U[GRID][GRID][GRID][GRID], int i, int j) {   

    double x0 = get_x(i);
    double y0 = get_x(j);

    //Helping vars
    double x1, y1, fd_f_x0, fd_f_y0, fd_f_x1, fd_f_y1, part_1, part_2, part_3, part_4;
    x1 = y1 = fd_f_x0 = fd_f_y0 = fd_f_x1 = fd_f_y1 = part_1 = part_2 = part_3 = part_4 = 0;

    for (int k = 0; k < GRID; k++)
    {
        x1 = get_x(k);
        for (int l = 0; l < GRID; l++)        
        {       
            y1 = get_x(l);
            fd_f_x0 = fder_f_x0(i,j,k,l, a_func_f);
            fd_f_y0 = fder_f_y0(i,j,k,l, a_func_f);
            fd_f_x1 = fder_f_x1(i,j,k,l, a_func_f);
            fd_f_y1 = fder_f_y1(i,j,k,l, a_func_f);

            part_1 += ((pow(fd_f_x0,2) + pow(fd_f_y0,2))/2.0);
            part_2 += ((pow(fd_f_x1,2) + pow(fd_f_y1,2))/2.0);
            part_3 += (pow(a_func_f[i][j][k][l][0], 2)*a_potential_U[i][j][k][l]);
            part_4 += (pow(a_func_f[i][j][k][l][0], 4) - 2.0*pow(a_func_f[i][j][k][l][0], 2) + 1.0)/2.0;
        }
    }

    return ((part_1+part_2+part_3+part_4) * pow(EPS,2)/GAMMA);
}


//Initialize Functions
void init_U(double (*func)(double), double a_potential_U[GRID][GRID][GRID][GRID]){
    double x0, y0, x1, y1, d_x, d_y, dist;
    for (int i = 0; i < GRID; i++)
    {
        x0 = get_x(i);
        for (int j = 0; j < GRID; j++)
        {
            y0= get_x(j);
            for (int k = 0; k < GRID; k++)
            {
                x1 = get_x(k);
                d_x = x0 - x1;
                for (int l = 0; l < GRID; l++)
                {
                    y1 = get_x(l);                    
                    d_y = y0 - y1;
                    dist = sqrt(pow(d_x,2)+ pow(d_y,2));
                    a_potential_U[i][j][k][l] = func(dist);
                }                
            }            
        }        
    }    
}

//Normalization Functions
double norm_g(fftw_complex a_func[GRID][GRID]){
    double sum = 0;
    for (int i = 0; i < GRID; i++)
    {
        for (int j = 0; j < GRID; j++)
        {
            sum += (pow(a_func[i][j][0],2) + pow(a_func[i][j][1],2));
        }       
    }
    sum = sum * pow(EPS,2);
    sum = sqrt(sum);
    for (int i = 0; i < GRID; i++)
    {
        for (int j = 0; j < GRID; j++)
        {
            a_func[i][j][0] = a_func[i][j][0]/sum;
            a_func[i][j][1] = a_func[i][j][1]/sum;
        }       
    }
    return sum;
}

double norm_f(fftw_complex a_func[GRID][GRID][GRID][GRID]){
    double f_max = a_func[0][0][GRID/2][GRID/2][0];
    for (int i = 0; i < GRID; i++)
    {
        for (int j = 0; j < GRID; j++)
        {
            for (int k = 0; k < GRID; k++)
            {
                for (int l = 0; l < GRID; l++)
                {
                    a_func[i][j][k][l][0] = a_func[i][j][k][l][0]/f_max;
                    a_func[i][j][k][l][1] = a_func[i][j][k][l][1]/f_max;
                }
            }
        }
    }
    return f_max;
}

/*
double norm_f(fftw_complex a_func[GRID][GRID][GRID][GRID]){
    double sum = 0;
    for (int i = 0; i < GRID; i++)
    {
        for (int j = 0; j < GRID; j++)
        {
            for (int k = 0; k < GRID; k++)
            {
                for (int l = 0; l < GRID; l++)
                {
                    sum += (pow(a_func[i][j][k][l][0],2) + pow(a_func[i][j][k][l][1],2));
                }
            }
        }
    }
    sum = sum * pow(EPS,4);
    sum = sqrt(sum);
    for (int i = 0; i < GRID; i++)
    {
        for (int j = 0; j < GRID; j++)
        {
            for (int k = 0; k < GRID; k++)
            {
                for (int l = 0; l < GRID; l++)
                {
                    a_func[i][j][k][l][0] = a_func[i][j][k][l][0]/sum;
                    a_func[i][j][k][l][1] = a_func[i][j][k][l][1]/sum;
                }
            }
        }
    }
    return sum;
}*/

//IMPORTANT: Declare arrrays outside of main
//g-Funtion: initialization g[x0][y0]
static fftw_complex arr_function_g[GRID][GRID];
//f-Function: initialization f[x0][y0][x1][y1]
static fftw_complex arr_function_f[GRID][GRID][GRID][GRID];
//ftilde-Function: initialization f_tilde[x0][y0][x1][y1]
static fftw_complex arr_function_ftilde[GRID][GRID][GRID][GRID];

int main() {

    //Render Version Number
    cout << "-----------------" << endl;
    cout << "Version 1.0" << endl;
    cout << "-----------------" << endl;

    //Init Threads
    fftw_init_threads();
    fftw_plan_with_nthreads(2);

    //Setup IO
    ofstream file_g;
    ofstream file_f;
    ofstream file_f_final;
    ofstream file_g_final;
    ofstream file_chem;
    file_g.open ("./function_g.txt");
    file_f.open ("./function_f.txt");
    file_f_final.open ("./function_f_final.txt");
    file_g_final.open ("./function_g_final.txt");
    file_chem.open ("./function_chem.txt");

    //Setup FFTW
    const int d4[Rank4] = {GRID, GRID, GRID, GRID};
    const int d2[Rank2] = {GRID, GRID};
    fftw_plan plan_g_forward = fftw_plan_dft(Rank2, d2, &arr_function_g[0][0], &arr_function_g[0][0], FFTW_FORWARD, FFTW_ESTIMATE);
    fftw_plan plan_g_backward = fftw_plan_dft(Rank2, d2, &arr_function_g[0][0], &arr_function_g[0][0], FFTW_BACKWARD, FFTW_ESTIMATE);
    fftw_plan plan_ft_forward = fftw_plan_dft(Rank4, d4, &arr_function_ftilde[0][0][0][0], &arr_function_ftilde[0][0][0][0], FFTW_FORWARD, FFTW_ESTIMATE);
    fftw_plan plan_ft_backward = fftw_plan_dft(Rank4, d4, &arr_function_ftilde[0][0][0][0], &arr_function_ftilde[0][0][0][0], FFTW_BACKWARD, FFTW_ESTIMATE);

    //Physical potentials
    static double arr_potential_U[GRID][GRID][GRID][GRID];
    static double arr_Vg[GRID][GRID];
    static double arr_Vf[GRID][GRID][GRID][GRID];

    //Initialize U
    init_U(potential_U, arr_potential_U);

    //Initialize other functions
    double x0, y0, x1, y1 = 0;
    for (int i = 0; i < GRID; i++)
    {
        x0 = get_x(i);
        for (int j = 0; j < GRID; j++)
        {
            y0 = get_x(j);
            arr_function_g[i][j][0] = start_func_g(x0,y0);
            arr_function_g[i][j][1] = 0;
            for (int k = 0; k < GRID; k++)
            {
                x1 = get_x(k);
                for (int l = 0; l < GRID; l++)
                {
                    y1 = get_x(l);
                    arr_function_f[i][j][k][l][0] = start_func_f(x0,y0,x1,y1);
                    arr_function_f[i][j][k][l][1] = 0;
                }                
            }            
        }        
    }

    //Norm
    norm_g(arr_function_g);
    norm_f(arr_function_f);

    //Calculate Resizes and change k
    double arr_k[GRID];
    double delta_k = 2.0*M_PI/(UPPER-LOWER);
    for (int i = 0; i < GRID; i++)
    {
        if (i >= GRID/2)
        {
            arr_k[i] = delta_k*i- delta_k*GRID;
        } else {
            arr_k[i] = delta_k*i;
        }        
    }
    
    //Main Loop
    for (int t = 0; t < STEPS; t++)
    {

        //Step 1 and 2
        for (int i = 0; i < GRID; i++)
        {
            for (int j = 0; j < GRID; j++)
            {
                arr_Vg[i][j] = calc_Vg(arr_function_f, arr_potential_U, i, j);
                //arr_Vg[i][j] = 0;
                for (int k = 0; k < GRID; k++)
                {
                    for (int l = 0; l < GRID; l++)
                    {
                        arr_Vf[i][j][k][l] = arr_Vg[i][j] + arr_potential_U[i][j][k][l] + pow(arr_function_f[i][j][k][l][0], 2) + pow(arr_function_f[i][j][k][l][1], 2);
                        arr_function_ftilde[i][j][k][l][0] = arr_function_f[i][j][k][l][0] * arr_function_g[i][j][0];
                        arr_function_ftilde[i][j][k][l][1] = arr_function_f[i][j][k][l][1] * arr_function_g[i][j][0];   //Assumes, that imaginary part of g is negligable
                    }                    
                }                
            }            
        }

        //Pre norm functions
        pre_norm_f = arr_function_ftilde[0][0][GRID/2][GRID/2][0];
        pre_norm_g = arr_function_g[0][0][0]; 

        //Step 3
        for (int i = 0; i < GRID; i++)
        {
            for (int j = 0; j < GRID; j++)
            {
                arr_function_g[i][j][0] = exp(-T_STEP*arr_Vg[i][j]/2.0)*arr_function_g[i][j][0];
                arr_function_g[i][j][1] = exp(-T_STEP*arr_Vg[i][j]/2.0)*arr_function_g[i][j][1];
            }   
        }

        //Step 4
        fftw_execute(plan_g_forward);
        for (int i = 0; i < GRID; i++)
        {
            k0x = arr_k[i];
            for (int j = 0; j < GRID; j++)
            {
                k0y = arr_k[j];
                TI = 0.5*(pow(k0x,2)+pow(k0y,2));
                arr_function_g[i][j][0] = exp(-T_STEP*TI)*arr_function_g[i][j][0];
                arr_function_g[i][j][1] = exp(-T_STEP*TI)*arr_function_g[i][j][1];
            }
            
        }
        fftw_execute(plan_g_backward);
        for (int i = 0; i < GRID; i++)
        {
            for (int j = 0; j < GRID; j++)
            {
                arr_function_g[i][j][0] = arr_function_g[i][j][0]/pow(GRID,2);
                arr_function_g[i][j][1] = arr_function_g[i][j][1]/pow(GRID,2);
            }
            
        }

    	//Step 5
        for (int i = 0; i < GRID; i++)
        {
            for (int j = 0; j < GRID; j++)
            {
                arr_function_g[i][j][0] = exp(-T_STEP*arr_Vg[i][j]/2.0)*arr_function_g[i][j][0];
                arr_function_g[i][j][1] = exp(-T_STEP*arr_Vg[i][j]/2.0)*arr_function_g[i][j][1];
            }   
        }

        //Step 6
        prop_norm_g = arr_function_g[0][0][0];
        norm_g(arr_function_g);
        chem_I = log(pre_norm_g/prop_norm_g)/T_STEP;

        //Step 7
        for (int i = 0; i < GRID; i++)
        {
            for (int j = 0; j < GRID; j++)
            {
                for (int k = 0; k < GRID; k++)
                {
                    for (int l = 0; l < GRID; l++)
                    {
                        arr_function_ftilde[i][j][k][l][0] = exp(-arr_Vf[i][j][k][l]*T_STEP/2.0)*arr_function_ftilde[i][j][k][l][0];
                        arr_function_ftilde[i][j][k][l][1] = exp(-arr_Vf[i][j][k][l]*T_STEP/2.0)*arr_function_ftilde[i][j][k][l][1];
                    }                    
                }
            }            
        }        

        //Step 8
        fftw_execute(plan_ft_forward);
        for (int i = 0; i < GRID; i++)
        {
            k0x = arr_k[i];
            for (int j = 0; j < GRID; j++)
            {
                k0y = arr_k[j];
                TI = 0.5*(pow(k0x,2)+pow(k0y,2));
                for (int k = 0; k < GRID; k++)
                {
                    k1x = arr_k[k];
                    for (int l = 0; l < GRID; l++)
                    {
                        k1y = arr_k[l];
                        TB = 0.5*(pow(k1x,2)+pow(k1y,2));
                        arr_function_ftilde[i][j][k][l][0] = exp(-TI*T_STEP-TB*T_STEP)*arr_function_ftilde[i][j][k][l][0];
                        arr_function_ftilde[i][j][k][l][1] = exp(-TI*T_STEP-TB*T_STEP)*arr_function_ftilde[i][j][k][l][1];
                    }                    
                }
            }            
        }
        fftw_execute(plan_ft_backward);
        for (int i = 0; i < GRID; i++)
        {
            for (int j = 0; j < GRID; j++)
            {
                for (int k = 0; k < GRID; k++)
                {
                    for (int l = 0; l < GRID; l++)
                    {
                        arr_function_ftilde[i][j][k][l][0] = arr_function_ftilde[i][j][k][l][0]/pow(GRID,4);
                        arr_function_ftilde[i][j][k][l][1] = arr_function_ftilde[i][j][k][l][1]/pow(GRID,4);
                    }                    
                }
            }            
        }

        //Step 9
        for (int i = 0; i < GRID; i++)
        {
            for (int j = 0; j < GRID; j++)
            {
                for (int k = 0; k < GRID; k++)
                {
                    for (int l = 0; l < GRID; l++)
                    {
                        arr_function_ftilde[i][j][k][l][0] = exp(-arr_Vf[i][j][k][l]*T_STEP/2)*arr_function_ftilde[i][j][k][l][0];
                        arr_function_ftilde[i][j][k][l][1] = exp(-arr_Vf[i][j][k][l]*T_STEP/2)*arr_function_ftilde[i][j][k][l][1];
                    }                    
                }                
            }
        }

        //Calc prop norm
        prop_norm_f = arr_function_ftilde[0][0][GRID/2][GRID/2][0];
        chem_B = log(pre_norm_f/prop_norm_f)/T_STEP-chem_I;

        //Step 10
        for (int i = 0; i < GRID; i++)
        {
            for (int j = 0; j < GRID; j++)
            {
                for (int k = 0; k < GRID; k++)
                {
                    for (int l = 0; l < GRID; l++)
                    {
                        double h = pow(arr_function_g[i][j][0],2) + pow(arr_function_g[i][j][1],2);
                        arr_function_f[i][j][k][l][0] = (arr_function_ftilde[i][j][k][l][0]*arr_function_g[i][j][0] + arr_function_ftilde[i][j][k][l][1]*arr_function_g[i][j][1])/h;
                        arr_function_f[i][j][k][l][1] = (arr_function_ftilde[i][j][k][l][1]*arr_function_g[i][j][0] - arr_function_ftilde[i][j][k][l][0]*arr_function_g[i][j][1])/h;
                    }                    
                }                
            }            
        }

        //Step 11
        prop_norm_f = norm_f(arr_function_f);

        //Write to file
        if(t % 100 == 0){
            cout << "index t: " << t << endl;

            for (int i = 0; i < GRID; i++)
            {
                file_f << current_time << " " << get_x(i) << " " << arr_function_f[i][GRID/2][GRID/2][GRID/2][0] << " " << arr_function_f[GRID/2][i][GRID/2][GRID/2][0] << " " << arr_function_f[0][0][i][0][0] << " " << arr_function_f[0][0][0][i][0] << "\n";
            }
        
            file_f << "\n"  << flush;
            file_chem << current_time << " " << chem_I << " " << chem_B << "\n";

            
            for (int i = 0; i < GRID; i++)
            {
                for (int j = 0; j < GRID; j++)
                {
                    file_g << get_x(i) << " " << get_x(j) << " " << arr_function_g[i][j][0] << "\n";
                }
                file_g << "\n";                
            }
            file_g << "\n" << flush;
            
        }
        current_time += T_STEP;
        
    }//End of main loop   

    for (int i = 0; i < GRID; i++)
    {
        for (int j = 0; j < GRID; j++)
        {
            file_g_final << arr_function_g[i][j][0] << "\n";
            for (int k = 0; k < GRID; k++)
            {
                for (int l = 0; l < GRID; l++)
                {
                    file_f_final << arr_function_f[i][j][k][l][0] << "\n";
                }                
            }            
        }        
    }

    //Close files
    file_f.close();
    file_f_final.close();
    file_g_final.close();
    file_chem.close();
    file_g.close();
}
