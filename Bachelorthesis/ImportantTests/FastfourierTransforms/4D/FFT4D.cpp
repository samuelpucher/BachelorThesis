#include <iostream>
#include <fftw3.h>
#include <cmath>
#include <complex>
#include <iomanip>
#include <fstream>


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

//fftw_complex mul_complex_4D(fftw_complex c1, fftw_complex c2, int N1, int N2, int N3, int N4) {
//    fftw_complex c_res[N][N][N][N];
//
//    c_res[][][][][]
//}


int main(){

    //Stream
    ofstream forwardfile;
    forwardfile.open("./result_forward.txt");
    ofstream backwardfile;
    backwardfile.open("./result_backward.txt");

    //Declare vars
    static fftw_complex data_in[N][N][N][N], data_out[N][N][N][N], data_buffer[N][N][N][N];
    double x1, x2, x3, x4;
    
    //Set dimensions
    int d[4];
    d[0] = N;
    d[1] = N;
    d[2] = N;
    d[3] = N;

    //Init plan    
    fftw_plan plan_forward = fftw_plan_dft(Rank, d, &data_in[0][0][0][0], &data_buffer[0][0][0][0], FFTW_FORWARD, FFTW_ESTIMATE);
    fftw_plan plan_backward = fftw_plan_dft(Rank, d, &data_buffer[0][0][0][0], &data_out[0][0][0][0], FFTW_BACKWARD, FFTW_ESTIMATE);

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
    fftw_execute(plan_forward);

    //Write file
    for (int i = 0; i < N; i++)
    {
        forwardfile <<  arr_k[i] << " " << data_buffer[i][0][0][0][0] << " " << data_buffer[i][0][0][0][1] << "\n"; 
    }
    forwardfile.close();


    //Transform backward
    fftw_execute(plan_backward);

    //Note: Due to normalization, divide by N^4 = 64^4

    //Write file
    for (int i = 0; i < N; i++)
    {
        backwardfile << get_x(i) << " " << data_in[i][N/2][N/2][N/2][0] << " " << data_out[i][N/2][N/2][N/2][0]/(pow(64,4)) << " " << data_out[i][N/2][N/2][N/2][1] << "\n"; 
    }
    backwardfile.close();

    //Performance Testing
    for (int i = 0; i < 999; i++)
    {
        fftw_execute(plan_forward);
        fftw_execute(plan_backward);
    }
    

    //Close fftw
    fftw_destroy_plan(plan_forward);
    fftw_destroy_plan(plan_backward);
    fftw_cleanup();

    return 0;
}