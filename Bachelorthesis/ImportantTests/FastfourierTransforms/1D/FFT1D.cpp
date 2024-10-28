//Compilewith: g++ -lfftw3 -O2 

#include <iostream>
#include <fftw3.h>
#include <cmath>
#include <complex>
#include <iomanip>
#include <fstream>

using namespace std;

const int N = 64;
const int Rank = 1;
const double UpperBound =2.5;
const double LowerBound = -2.5;

double get_x(int i) {
    double j = i;
    double total = UpperBound-LowerBound;
    return LowerBound + total* j/N;
}

double func(double x){
    return exp(-(x)*(x)/0.25);
}

int main(){

    //Stream
    ofstream forwardfile;
    forwardfile.open("./result_forward.txt");
    ofstream backwardfile;
    backwardfile.open("./result_backward.txt");

    //Declare vars
    static fftw_complex data_in[N], data_out[N], data_buffer[N];         //FFTW complex type: [0] real part, [1] imaginary part
    double x;
    double arr_k[N];
    double delta_k = 2*M_PI/(UpperBound-LowerBound);

    for (int i = 0; i < N; i++)
    {
        if (i >= N/2)
        {
            arr_k[i] = delta_k*i- delta_k*N;
        } else {
            arr_k[i] = delta_k*i;
        }        
    }    

    //Set dimensions 
    const int d[1] = {N};

    //init array
    for(int i =0; i<N; i++) {
        x = get_x(i);
        data_in[i][0] = func(x);
        data_in[i][1] = 0;
        data_out[i][0] = 0;
        data_out[i][1] = 0;
    }

    cout << data_out[15][0] << endl;

    //Init plan    
    fftw_plan plan_forward = fftw_plan_dft(Rank, d, &data_in[0], &data_buffer[0], FFTW_FORWARD, FFTW_ESTIMATE);
    fftw_plan plan_backward = fftw_plan_dft_1d(N, &data_buffer[0], &data_out[0], FFTW_BACKWARD, FFTW_ESTIMATE);

    //Transform forward
    fftw_execute(plan_forward);

    //Write file
    for (int i = 0; i < N; i++)
    {
        forwardfile <<  arr_k[i] << " " << data_buffer[i][0] << " " << data_buffer[i][1] << "\n"; 
    }
    forwardfile.close();

    //Transform backward
    fftw_execute(plan_backward);

    //Write file
    for (int i = 0; i < N; i++)
    {
        backwardfile << get_x(i) << " " << data_in[i][0] << " " << data_out[i][0] << " " << data_out[i][1] << "\n"; 
    }
    backwardfile.close();

    //Close fftw
    fftw_destroy_plan(plan_forward);
    fftw_destroy_plan(plan_backward);
    fftw_cleanup();

    return 0;

}