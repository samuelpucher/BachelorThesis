#include <iostream>
#include <cmath>
#include <complex>
#include <iomanip>
#include <fstream>

#include <omp.h>
#include <fftw3.h>

using namespace std;

const int N = 64;
const int Rank = 4;
const double UpperBound =2.5;
const double LowerBound = -2.5;

double func(double x1, double x2, double x3, double x4){
    return exp(-(x1*x1+ x2*x2 + x3*x3 + x4*x4)/0.25);
}

double get_x(int i) {
    double j = i;
    double total = UpperBound-LowerBound;
    return LowerBound + total* j/N;
}

int main(){

    //Init threads
    fftw_init_threads();
    fftw_plan_with_nthreads(1);

    //Stream
    ofstream forwardfile;
    forwardfile.open("./result_forward_threads.txt");


    //Declare vars
    static fftw_complex data_in[N][N][N][N], data_out[N][N][N][N], data_buffer[N][N][N][N];
    double x1, x2, x3, x4;
    
    //Set dimensions
    const int d[4] = {N,N,N,N};

    //Init plan    
    fftw_plan plan = fftw_plan_dft(Rank, d, &data_in[0][0][0][0], &data_buffer[0][0][0][0], FFTW_FORWARD, FFTW_ESTIMATE);

    //Init arrays
    for(int i =0; i<N; i++) 
    {
        x1 = get_x(i);
        for (int j = 0; j < N; j++)
        {
            x2 = get_x(j);
            for (int k = 0; k < N; k++)
            {
                x3 = get_x(k);
                for (int l = 0; l < N; l++)
                {
                    x4 = get_x(l);
                    data_in[i][j][k][l][0] = func(x1, x2, x3, x4);
                    data_in[i][j][k][l][1] = 0;
                    data_out[i][j][k][l][0] = 0;
                    data_out[i][j][k][l][1] = 0;
                }
            }
        }
    }

    //Calculate Resizes
    double arr_k[N];
    double delta_k = 2*M_PI/(UpperBound-LowerBound);

    //Change k
    for (int i = 0; i < N; i++)
    {
        if (i >= N/2)
        {
            arr_k[i] = delta_k*i- delta_k*N;
        } else {
            arr_k[i] = delta_k*i;
        }        
    }

    //Transform forward
    fftw_execute(plan);

    //Write file
    for (int i = 0; i < N; i++)
    {
        forwardfile <<  arr_k[i] << " " << data_buffer[i][0][0][0][0] << " " << data_buffer[i][0][0][0][1] << "\n"; 
    }
    forwardfile.close();

    //Performance Testing
    
    for (int i = 0; i < 499; i++)
    {
        fftw_execute(plan);
    }
    

    //Close fftw
    fftw_destroy_plan(plan);
    void fftw_cleanup_threads(void);

    return 0;
}