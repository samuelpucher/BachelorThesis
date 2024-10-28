
#include <iostream>
#include <fftw3.h>
#include <cmath>
#include <complex>
#include <iomanip>
#include <fstream>
using namespace std;


//DEFINE GRID SIZE (globally to make arrays passable)

const int N=pow(2, 11);


//FUNCTION DEFINITIONS: COMPLEX OPERATIONS

inline void plueq(fftw_complex *a, fftw_complex *b) //addition "from the left": right+=left
{
	(*b)[0]+=(*a)[0];
	(*b)[1]+=(*a)[1];
}
	
inline void muleq(fftw_complex *a, fftw_complex *b) //multiplication with a complex number "from the left": right*=left
{
	double r=(*a)[0]*(*b)[0]-(*a)[1]*(*b)[1];
	double i=(*a)[0]*(*b)[1]+(*a)[1]*(*b)[0];
	(*b)[0]=r;
	(*b)[1]=i;
}

inline void timeq(double c, fftw_complex *b) //multiplication with a real number "from the left": right*=left
{
	(*b)[0]*=c;
	(*b)[1]*=c;
}

inline void expeq(fftw_complex *b) //exponentiation: arg=e^arg
{
	double c=cos((*b)[1]);
	double s=sin((*b)[1]);
	double r=exp((*b)[0]);
	(*b)[0]=r*c;
	(*b)[1]=r*s;
}

inline void coneq(fftw_complex *b) //conjugate: arg=arg*
{
	(*b)[1]*=-1.0;
}


//FUNCTION DEFINITIONS: ARRAY OPERATIONS

inline double maximag(fftw_complex (&a)[N][N])
{
	static double maxim;
	maxim=0;
	double m=0;
	for(int i=0; i<N;i++){for(int j=0; j<N;j++){	//always i<->x , j<->y
		m=(a)[i][j][1];
		if(m<0){m*=-1.0;}
		if(m>maxim){maxim=m;}
		}}
	
	return maxim;
}

inline double maxerror(fftw_complex (&a)[N][N],fftw_complex (&b)[N][N])
{
	static double maxe;
	maxe=0;
	double m=0;
	for(int i=0; i<N;i++){for(int j=0; j<N;j++){
		m=(a)[i][j][0]-(b)[i][j][0];
		if(m<0){m*=-1.0;}
		if(m>maxe){maxe=m;}
		m=(a)[i][j][1]-(b)[i][j][1];
		if(m<0){m*=-1.0;}
		if(m>maxe){maxe=m;}
		}}
	
	return maxe;
}






//MAIN

int main(int argc, char **argv)			//[1]=directory, [2]=hard_cores?, [3]=rho, [4]=alpha, [5]=dr, [6]=dt
{
	
	string mainDirectory = "./";

    if (argc > 1) {
        mainDirectory = argv[1];
    }
	
	//creating the output files
	
	ofstream Efile;
	Efile.open (mainDirectory + "/E_5000.txt");
	
	
	ofstream outputfile;
	outputfile.open (mainDirectory + "/Psi_output.txt");
	
	
	ofstream Psifile;
	Psifile.open (mainDirectory + "/Psi_new.txt");
	
	ofstream testfile;
	testfile.open (mainDirectory + "/test.txt");
	
	
	ofstream gfile;
	gfile.open (mainDirectory + "/g_5000.txt");
	ofstream Sfile;
	Sfile.open (mainDirectory + "/S_5000.txt");
	
	
	
	ofstream gfile2;
	gfile2.open (mainDirectory + "/g_full_of_t.txt");
	ofstream Sfile2;
	Sfile2.open (mainDirectory + "/S_full_of_t.txt");
	ofstream phifile2;
	phifile2.open (mainDirectory + "/phi_full_of_t.txt");
	ofstream Hfile2;
	Hfile2.open (mainDirectory + "/H_full_of_t.txt");
	ofstream Afile2;
	Afile2.open (mainDirectory + "/A_full_of_t.txt");
	
	
	ofstream vfile;
	vfile.open (mainDirectory + "/v.txt");
	
	ofstream wfile;
	wfile.open (mainDirectory + "/FTomegaI.txt");
	
	ofstream wfile2;
	wfile2.open (mainDirectory + "/omegaI.txt");
	
	
	ofstream goftfile;
	goftfile.open (mainDirectory + "/g_of_t_p250.txt");
	
	ofstream Softfile;
	Softfile.open (mainDirectory + "/S_of_t_p250.txt");
	
	ofstream phioftfile;
	phioftfile.open (mainDirectory + "/phi_of_t_p250.txt");
	
	ofstream Hoftfile;
	Hoftfile.open (mainDirectory + "/H_of_t_p250.txt");
	
	ofstream Aoftfile;
	Aoftfile.open (mainDirectory + "/A_of_t_p250.txt");
	
	
	
	
		
	double x=0; //multi purpose storage for real number
	double x2=0; //multi purpose storage for real number
	double y=0; //multi purpose storage for real number
	double y2=0; //multi purpose storage for real number
	int i2=0; //storage for index number2
	int j2=0; //storage for index number2
	
	
	//ADJUSTABLE STUFF
	
	
	//get grid parameters from (already defined) grid size N
	double Nd=N;
	//double dr = 0.2*sqrt(2.0*M_PI/Nd);	//real space step size adjustable
	double dr=0.00125;
	
	if (argc > 5) {
        dr = strtod(argv[5],NULL);
	}
	
	double dk = 2.0*M_PI/(Nd*dr);
	
	//create 1D coordinate "grids" (mirror-at-N/2-version)
	static double r[N];
	for(int i=0;i<N/2;i++){
		r[i]=i*dr;
		}
	for(int i=N/2;i<N;i++){
		r[i]=(i-N)*dr;
		}
	static double k[N];
	for(int i=0;i<N/2;i++){
		k[i]=i*dk;
		}
	for(int i=N/2;i<N;i++){
		k[i]=(i-N)*dk;
		}
		
	
	static double k2[2*N];
	for(int i=0;i<N;i++){
		k2[i]=i*dk*0.5;
		}
	for(int i=N;i<2*N;i++){
		k2[i]=(i-2*N)*dk*0.5;
		}
	
	
	//set physical constants
	//double hm = 1.0; //this is hbar/m -> defines the units i guess???
	double hm=1.0;	//=2 for h^2/2m =1 		????
	
	double rho=40; //density FROM WHERE??
	
	if (argc > 3) {
        rho = strtod(argv[3],NULL);
	}
	
	//set potential at t=0
	static double v[N][N];
	
	double alp=0.6;
	
	if (argc > 4) {
        alp = strtod(argv[4],NULL);
    }
	
	double s2=sin(alp)*sin(alp);
	
	int hc=0;		//use hard core repulsion?? (use as bool)
	
	if (argc > 2) {
        hc = atoi(argv[2]);
    }
    
    
    double spr=100.0/0.2;	//time steps per full rotation of the interaction
    
    
    double aps=2.0*M_PI/spr;	//angle change per timestep
    
    
    
    bool dp=0;  //damped potential
    double dpe=1; //damped potential exponent
    
    
    
    double rot=0;
    double co=cos(rot);
    double si=sin(rot);
    
	
	for(int i=0; i<N;i++){for(int j=0; j<N;j++){
		if((i==0)&&(j==0)){v[i][j]=0;}	//treat seperately!!!
		else{
		x=sqrt(r[i]*r[i]+r[j]*r[j]);
		x2=co*r[i]-si*r[j];
		
		v[i][j]=(1.0-3.0*s2*x2*x2/(x*x))/(x*x*x);
		
		if(dp){v[i][j]*=exp(-1.0*dpe*x);}
		
		if(hc==1){
			y2=0.31/x;
			y=y2*y2*y2*y2*y2*y2*y2*y2*y2*y2*y2*y2;
	
			v[i][j]+=y;
		}
		
		
	}
	}}
	
	
	//choose imag or real(i.e. both) time evolution
	bool imev=0;  			// 0-> only real ; 1 -> include imag
	bool change=0;			// choose whether to change to real time later (if imag is included)
    
    fftw_complex it={1.0,0}; //exponent factor for evolution (imaginery evolution -> I , real evolution -> 1)
    if(imev==0){
		it[0]=0;
		it[1]=1.0;
		}
		
	//choose timestep
	//double dt=0.00125; //imag step  //units???
	double dt=0.0000025*0.2;
	
	if (argc > 6) {
        dt = strtod(argv[6],NULL);
	}
	
	
	double dT=0.0000025*0.2;
	if(imev==0){dt=dT;}
	
	//chosse maximum iteration numbers
	int imagsteps=0;		  //5000 seems ok
	int realsteps=10001;
	
	int maxita=0;
	if(imev==1){ maxita+=imagsteps;}
	if((imev==0)||(change==1)){ maxita+=realsteps;}
	
	//missing: break condition
	
	
	//initialize pseudo-wavefunction
	static fftw_complex Psi[N][N];
	
	/*
	for(int i=0; i<N;i++){for(int j=0; j<N;j++){
		Psi[i][j][0]=1.0; //initial value adjustable
		Psi[i][j][1]=0;
	}}
	*/
	
								//CONTINUE: read in PSI from file. N,dr,dk NOT CHECKED!!!
		ifstream inputfile("Psi_temp.txt");
		double Rpsi=0;
		double Ipsi=0;
		double d1=0;
		double d2=0;
	for(int i=0; i<N;i++){for(int j=0; j<N;j++){
		inputfile>>d1>>d2>>Rpsi>>Ipsi;
		Psi[i][j][0]=Rpsi;
		Psi[i][j][1]=Ipsi;
	}}
	
	

	
	
	
	
	//NOT ADJUSTABLE STUFF
	
	
	//define useful constants and variables
	fftw_complex I={0,1.0};
	
	fftw_complex z={0,0}; //multi purpose storage for complex number
	fftw_complex z1={0,0}; //multi purpose storage for complex number
	fftw_complex z2={0,0}; //multi purpose storage for complex number
	
	
	double pi=M_PI; //pi without call
	double cf=rho*dr*dr; //constant for forward transformation
	double cb=dk*dk/(rho*4.0*pi*pi); //constant for backward transformation
	
	double cf1D=rho*dr; //constant for forward transformation
	double cb1D=dk*0.5/(rho*2.0*pi); //constant for backward transformation
	
	fftw_complex E={0,0};	//Energy per particle	only Re part SHOULD ever be non-zero
	
	
	
	//create the plans
	static fftw_complex inf[N][N]; //in and out for forward plan
	static fftw_complex outf[N][N];
	static fftw_complex inb[N][N]; //in and out for backward plans
	static fftw_complex outb[N][N];
	
	static fftw_complex in1D[2*N]; //in and out for 1D plans for the jump correction due to H
	static fftw_complex out1D[2*N];
	
	
	fftw_plan forwardplan;
	fftw_plan backwardplan;
	
	fftw_plan forwardplan1D;
	fftw_plan backwardplan1D;
	
	
	forwardplan = fftw_plan_dft_2d(N, N, &inf[0][0], &outf[0][0], FFTW_FORWARD, FFTW_MEASURE);
	backwardplan = fftw_plan_dft_2d(N, N, &inb[0][0], &outb[0][0], FFTW_BACKWARD, FFTW_MEASURE);
	
	forwardplan1D = fftw_plan_dft_1d(2*N, &in1D[0], &out1D[0], FFTW_FORWARD, FFTW_MEASURE);
	backwardplan1D = fftw_plan_dft_1d(2*N, &in1D[0], &out1D[0], FFTW_BACKWARD, FFTW_MEASURE);
	
	
	
	//define variable fields
	static fftw_complex FTPsi[N][N]; //FT of Psi
	static fftw_complex S[N][N]; //S(k)
	static fftw_complex FTf1[N][N]; //FT of f1(=g dx phi)
	static fftw_complex FTf2[N][N]; //FT of f2(=g dy phi)
	static fftw_complex H1[N][N];
	static fftw_complex H2[N][N];
	static fftw_complex eV[N][N]; //exp(-it*dt*V/(2*hbar))
	
	static fftw_complex A1[N][N];
	static fftw_complex A2[N][N];
	
	
	
	
	//initialize all fields
	for(int i=0; i<N;i++){for(int j=0; j<N;j++){
		inf[i][j][0]=0;
		inf[i][j][1]=0;
		outf[i][j][0]=0;
		outf[i][j][1]=0;
		inb[i][j][0]=0;
		inb[i][j][1]=0;
		outb[i][j][0]=0;
		outb[i][j][1]=0;
		FTPsi[i][j][0]=0;
		FTPsi[i][j][1]=0;
		S[i][j][0]=0;
		S[i][j][1]=0;
		FTf1[i][j][0]=0;
		FTf1[i][j][1]=0;
		FTf2[i][j][0]=0;
		FTf2[i][j][1]=0;
		H1[i][j][0]=0;
		H1[i][j][1]=0;
		H2[i][j][0]=0;
		H2[i][j][1]=0;
		eV[i][j][0]=0;
		eV[i][j][1]=0;
		
		A1[i][j][0]=0;
		A1[i][j][1]=0;
		A2[i][j][0]=0;
		A2[i][j][1]=0;
	}}	
	
	
	for(int i=0; i<2*N;i++){
		in1D[i][0]=0;
		in1D[i][1]=0;
		out1D[i][0]=0;
		out1D[i][1]=0;
	}
	
	
	
	
	
	//BEGIN ITERATION
	for(int step=0; step<maxita;step++){
		
		timeq(0,&E);					//reset energy
	
	//update potential if necessary AND add v to eV
	if(imev==0){
		
		//update potential
		rot+=aps;
		co=cos(rot);
		si=sin(rot);
    
    
		for(int i=0; i<N;i++){for(int j=0; j<N;j++){
		if((i==0)&&(j==0)){v[i][j]=0;}	//treat seperately!!!
		else{
		x=sqrt(r[i]*r[i]+r[j]*r[j]);
		x2=co*r[i]-si*r[j];
		
		v[i][j]=(1.0-3.0*s2*x2*x2/(x*x))/(x*x*x);
		
		if(dp){v[i][j]*=exp(-1.0*dpe*x);}
		
		if(hc==1){
			y2=0.31/x;
			y=y2*y2*y2*y2*y2*y2*y2*y2*y2*y2*y2*y2;
	
			v[i][j]+=y;
		}
		}
			
		timeq(0,&eV[i][j]);					//eV=0
		plueq(&it,&eV[i][j]);				//eV=I
		timeq(-0.5*dt*v[i][j],&eV[i][j]);	//eV=-I*dt*v/(2*hbar)
		expeq(&eV[i][j]);					//eV=exp(-I*dt*v/(2*hbar))
	}}
	}
	else{
		for(int i=0; i<N;i++){for(int j=0; j<N;j++){
		timeq(0,&eV[i][j]);					//eV=0
		plueq(&it,&eV[i][j]);				//eV=1
		timeq(-0.5*dt*v[i][j],&eV[i][j]);	//eV=-v*dt/(2*hbar)
		expeq(&eV[i][j]);					//eV=exp(-v*dt/(2*hbar))
	}}
	}
	
	cout<<"eV(1)="<<eV[1][0][0]<<"+i*"<<eV[1][0][1]<<endl;
	
	
	//COMPUTE VARIABLE FIELDS
	
	
	
	//compute S AND add omega_I to eV
	for(int i=0; i<N;i++){for(int j=0; j<N;j++){
		timeq(0,&inf[i][j]);			//in=0
		timeq(0,&z);				//z=0
		plueq(&Psi[i][j],&inf[i][j]);//in=Psi
		plueq(&Psi[i][j],&z);		//z=Psi
		coneq(&z);					//z=Psi*
		muleq(&z,&inf[i][j]);		//in=|Psi|^2=g
		inf[i][j][0]-=1.0;			//in=g-1
	}}
	
	
	fftw_execute(forwardplan);		//out=(S-1)*const
	
	for(int i=0; i<N;i++){for(int j=0; j<N;j++){
		timeq(0,&S[i][j]);			//"S"=0
		plueq(&outf[i][j],&S[i][j]);	//"S"=out
		timeq(cf,&S[i][j]);			//"S"=S-1
		S[i][j][0]+=1.0;			//"S"=S
	}}	
	
			
	
			//(rescale Psi s.t. S(0)=0)
			x=Nd*Nd*cf;							//x=V+
			y=sqrt((x-1.0)/(x-1.0+S[0][0][0]));	//y=a for rescale (see notes)
			cout<<"a="<<y<<endl;
			
			for(int i=0; i<N;i++){for(int j=0; j<N;j++){
				timeq(y,&Psi[i][j]);		//rescale
				//recompute S
				timeq(0,&inf[i][j]);			//in=0
				timeq(0,&z);				//z=0
				plueq(&Psi[i][j],&inf[i][j]);//in=Psi
				plueq(&Psi[i][j],&z);		//z=Psi
				coneq(&z);					//z=Psi*
				muleq(&z,&inf[i][j]);		//in=|Psi|^2=g
				
					//add I_v to E
					timeq(0,&z);					//z=0
					plueq(&inf[i][j],&z);			//z=g
					timeq(rho*0.5*v[i][j]*dr*dr,&z);//z=I_v
					plueq(&z,&E);					//zE+=I_v
				
				
				inf[i][j][0]-=1.0;			//in=g-1
			}}
	
			fftw_execute(forwardplan);		//out=(S-1)*const
	
			for(int i=0; i<N;i++){for(int j=0; j<N;j++){
				timeq(0,&S[i][j]);			//"S"=0
				plueq(&outf[i][j],&S[i][j]);	//"S"=out
				timeq(cf,&S[i][j]);			//"S"=S-1
				S[i][j][0]+=1.0;			//"S"=S
			}}	
	
			cout<<"S(0)="<<S[0][0][0]<<endl;
		
	
	/*
	//(shift Psi s.t. S(0)=0)
			x=0;				//I/V see notes
			for(int i=0; i<N;i++){for(int j=0; j<N;j++){
					x+=sqrt(Psi[i][j][0]*Psi[i][j][0]+Psi[i][j][1]*Psi[i][j][1])/(Nd*Nd);
			}}
	
			y=-1.0*x + sqrt(x*x-S[0][0][0]/(rho*dr*dr*Nd*Nd));	//y=d for shift (see notes)
			cout<<"c="<<y<<endl;
			
			for(int i=0; i<N;i++){for(int j=0; j<N;j++){
				Psi[i][j][0]+=y;			//shift
				//recompute S
				timeq(0,&inf[i][j]);		//in=0
				timeq(0,&z);				//z=0
				plueq(&Psi[i][j],&inf[i][j]);//in=Psi
				plueq(&Psi[i][j],&z);		//z=Psi
				coneq(&z);					//z=Psi*
				muleq(&z,&inf[i][j]);		//in=|Psi|^2=g
				inf[i][j][0]-=1.0;			//in=g-1
			}}
	
			fftw_execute(forwardplan);		//out=(S-1)*const
	
			for(int i=0; i<N;i++){for(int j=0; j<N;j++){
				timeq(0,&S[i][j]);			//"S"=0
				plueq(&outf[i][j],&S[i][j]);	//"S"=out
				timeq(cf,&S[i][j]);			//"S"=S-1
				S[i][j][0]+=1.0;			//"S"=S
			}}	
	
			cout<<"S(0)="<<S[0][0][0]<<endl;
	*/
	/*
	S[0][0][0]=0;
	S[0][0][1]=0;
	
	S[0][1][0]=0.5*S[0][2][0];
	S[1][0][0]=0.5*S[2][0][0];
	S[0][N-1][0]=0.5*S[0][N-2][0];
	S[N-1][0][0]=0.5*S[N-2][0][0];
	S[1][1][0]=0.5*S[2][2][0];
	S[1][N-1][0]=0.5*S[2][N-2][0];
	S[N-1][N-1][0]=0.5*S[N-2][N-2][0];
	S[N-1][1][0]=0.5*S[N-2][2][0];
	*/
		//compute omega_I(k)
	for(int i=0; i<N;i++){for(int j=0; j<N;j++){	
		if((i>0)||(j>0)){
			timeq(0,&inb[i][j]);			//in=0
			timeq(0,&z);					//z=0
			x=1.0/(S[i][j][0]*S[i][j][0]+S[i][j][1]*S[i][j][1]);	//x=1/|S|^2				check fordivision by zero???????
			plueq(&S[i][j],&z);										//z=S
			coneq(&z);												//z=S*
			timeq(x,&z);											//z=S* /|S|^2=1/S
			z[0]-=1.0;												//z=(1/S)-1
			plueq(&z,&inb[i][j]);									//in=z=(1/S)-1
			muleq(&z,&inb[i][j]);									//in=z^2=(1-(1/S))^2
			timeq(0,&z);											//z=0
			plueq(&S[i][j],&z);										//z=S
			timeq(2.0,&z);											//z=2S
			z[0]+=1.0;												//z=2S+1
			muleq(&z,&inb[i][j]);									//in=(2S+1)*(1-(1/S))^2
			x=-0.25*(k[i]*k[i]+k[j]*k[j])*hm;						//x=-hbar*k^2/(4m)
			timeq(x,&inb[i][j]);									//in=(-hbar*k^2/(4m))*(2S+1)*(1-(1/S))^2=omega_i(k)/hbar
		}
	}}	
		//compute omega_I(k=0) with de l'hopital
		/*
		x=(S[0][1][0]+S[1][0][0]+S[0][N-1][0]+S[N-1][0][0])/4.0;	//avarage S01
		x=x/dk;														//avarage S'01
		inb[0][0][0]=-0.25*hm/(x*x);
		inb[0][0][1]=0;
		*/
		//compute omega_I(k=0) with 2D quadratic approx
		timeq(0,&z);					//z=0
		timeq(0,&inb[0][0]);			//in=0
		
		plueq(&inb[0][1],&inb[0][0]);	//in=w01
		plueq(&inb[1][0],&inb[0][0]);	//in=w01+w10
		
		plueq(&inb[1][1],&z);			//z=w11
		plueq(&inb[1][N-1],&z);			//z=w11+w1-1
		timeq(-0.5,&z);					//z=-(w11+w1-1)/2
		
		plueq(&z,&inb[0][0]);			//in=w01+w10-(w11+w1-1)/2
		
		
	
	
		//compute omega_I(r)
	fftw_execute(backwardplan);				//out=const*omega_I(r)/hbar
	
	
	if(step==(maxita-1)){
		for(int i=0; i<N;i++){
		if(i==N/2){wfile2<<endl;}
		wfile2<<r[i]<<" "<<cb*outb[i][0][0]<<" "<<cb*outb[i][0][1]<<endl;
	}
	}
	
	x=1.0/(8.0*rho*pi*pi);
	
		//add to eV
	for(int i=0; i<N;i++){for(int j=0; j<N;j++){
						//add I_22 to E
						timeq(0,&z1);						//z1=0
						plueq(&S[i][j],&z1);				//z1=S
						timeq(2.0,&z1);						//z1=2S
						z1[0]+=1.0;							//z1=2S+1
						timeq(0,&z);						//z=0
						plueq(&z1,&z);						//z=2S+1
						coneq(&z);							//z=2S*+1
						timeq(1.0/(z1[0]*z1[0]+z1[1]*z1[1]),&z);//z=(2S*+1)/|2S+1|^2=1/(2S+1)
						muleq(&S[i][j],&z);					//z=S/(2S+1)
						timeq(0,&z1);						//z1=0
						plueq(&S[i][j],&z1);				//z1=S
						z1[0]-=1.0;							//z1=S-1
						muleq(&z1,&z);						//z=S(S-1)/(2S+1)
						muleq(&inb[i][j],&z);				//z=omega_I(k)*S(S-1)/(2S+1)
						timeq(x*dk*dk,&z);					//z=I_22
						plueq(&z,&E);						//zE+=I_22
		
		timeq(0,&z);						//z=0
		plueq(&outb[i][j],&z);				//z=out
		timeq(cb,&z);						//z=omega_I(r)/hbar		
		muleq(&it,&z);						//z=it*omega_I(r)/hbar
		timeq(-0.5*dt,&z);					//z=-it*dt*omega_I(r)/(2*hbar)
		expeq(&z);							//z=exp(-it*dt*omega_I(r)/(2*hbar))
		muleq(&z,&eV[i][j]);				//eV=exp(-it*dt*(v+omega_I)/(2*hbar))
																						//}}
	
	
	
	//compute f (Fourier transformed) and add beta_I to eV
	
		//compute FTPsi and then FTf1
																						//for(int i=0; i<N;i++){for(int j=0; j<N;j++){
		timeq(0,&inf[i][j]);						//in=0
		plueq(&Psi[i][j],&inf[i][j]);			//in=Psi
	}}
	
	fftw_execute(forwardplan);
	
	for(int i=0; i<N;i++){for(int j=0; j<N;j++){
		timeq(0,&FTPsi[i][j]);						//FTPsi=0
		plueq(&outf[i][j],&FTPsi[i][j]);				//FTPsi=out
		timeq(cf,&FTPsi[i][j]);						//"FTPsi"=FTPsi
		//compute d_x Psi
		timeq(0,&inb[i][j]);							//in=0
		plueq(&FTPsi[i][j],&inb[i][j]);				//in=FTPsi
		timeq(k[i],&inb[i][j]);					
		muleq(&I,&inb[i][j]);						//in=i*k_x*FTPsi
	}}
	
	fftw_execute(backwardplan);
	
	for(int i=0; i<N;i++){for(int j=0; j<N;j++){
		timeq(0,&z);						//z=0
		plueq(&outb[i][j],&z);				//z=out
		timeq(cb,&z);						//z=FT^(-1)[i*k_x*FTPsi]=d_x Psi
		
							//add x-parts of I_21 and I_1 to E
							E[0]+=0.5*rho*hm*(z[0]*z[0]+z[1]*z[1])*dr*dr;	//Vorfaktor?????
		
		//compute f1
		timeq(0,&z1);						//z1=0
		plueq(&Psi[i][j],&z1);				//z1=Psi
		coneq(&z1);							//z1=Psi*
		muleq(&z1,&z);						//z=z1*z=Psi* * (d_x Psi)
		inf[i][j][0]=z[1];
		inf[i][j][1]=0;						//in=Im[Psi* * (d_x Psi)]=f1
	}}
	
	fftw_execute(forwardplan);
	
	for(int i=0; i<N;i++){for(int j=0; j<N;j++){
		timeq(0,&FTf1[i][j]);				//FTf1=0
		plueq(&outf[i][j],&FTf1[i][j]);		//FTf1=out
		timeq(cf,&FTf1[i][j]);				//"FTf1"=FTf1
		//compute d_y Psi
		timeq(0,&inb[i][j]);					//in=0
		plueq(&FTPsi[i][j],&inb[i][j]);		//in=FTPsi
		timeq(k[j],&inb[i][j]);				//NOTE: J!!!!!!
		muleq(&I,&inb[i][j]);				//in=i*k_y*FTPsi
	}}
	
	fftw_execute(backwardplan);
	
	for(int i=0; i<N;i++){for(int j=0; j<N;j++){
		timeq(0,&z);						//z=0
		plueq(&outb[i][j],&z);				//z=out
		timeq(cb,&z);						//z=FT^(-1)[i*k_y*FTPsi]=d_y Psi
		
							//add y-parts of I_21 and I_1 to E
							E[0]+=0.5*rho*hm*(z[0]*z[0]+z[1]*z[1])*dr*dr;	//Vorfaktor?????
		
		//compute f2
		timeq(0,&z1);						//z1=0
		plueq(&Psi[i][j],&z1);				//z1=Psi
		coneq(&z1);							//z1=Psi*
		muleq(&z1,&z);						//z=z1*z=Psi* * (d_y Psi)
		inf[i][j][0]=z[1];
		inf[i][j][1]=0;						//in=Im[Psi* * (d_y Psi)]=f2
	}}
	
	fftw_execute(forwardplan);
	
	
	x=-1.0*hm*dk*dk/(8.0*rho*pi*pi);
	
	for(int i=0; i<N;i++){for(int j=0; j<N;j++){
		timeq(0,&FTf2[i][j]);				//FTf2=0
		plueq(&outf[i][j],&FTf2[i][j]);		//FTf2=out
		timeq(cf,&FTf2[i][j]);				//"FTf2"=FTf2
		
								//add I_3 to E
								timeq(0,&z1);						//z1=0
								plueq(&FTf1[i][j],&z1);				//z1=FTf1
								timeq(0,&z2);						//z2=0
								plueq(&z1,&z2);						//z2=FTf1
								muleq(&z1,&z2);						//z2=FTf1^2
								timeq(0,&z);						//z=0
								plueq(&FTf2[i][j],&z);				//z=FTf2
								timeq(0,&z1);						//z1=0
								plueq(&z,&z1);						//z1=FTf2
								muleq(&z,&z1);						//z1=FTf2^2
								plueq(&z1,&z2);						//z2=FTf1^2+FTf2^2
								
								timeq(0,&z);						//z=0
								plueq(&S[i][j],&z);					//z=S
								z[0]-=1.0;							//z=S-1
								muleq(&z2,&z);						//z=(FTf1^2+FTf2^2)*(S-1)
								timeq(x,&z);						//z=I_3
								plueq(&z,&E);						//E+=I_3
									
		
		
		//compute beta_I
		timeq(0,&inb[i][j]);					//in=0
		plueq(&FTf1[i][j],&inb[i][j]);		//in=FTf1
		muleq(&FTf1[i][j],&inb[i][j]);		//in=FTf1^2
		timeq(0,&z);						//in=0
		plueq(&FTf2[i][j],&z);				//in=FTf2
		muleq(&FTf2[i][j],&z);				//in=FTf2^2
		plueq(&z,&inb[i][j]);				//in=FTf1^2+FTf2^2=FTf^2
	}}
	
	fftw_execute(backwardplan);
		
	for(int i=0; i<N;i++){for(int j=0; j<N;j++){
		timeq(0,&z);				//z=0
		plueq(&outb[i][j],&z);		//z=out
		timeq(cb,&z);				//z=FT^(-1)[(FTf)^2]=f (conv) f
		timeq(-1.0*hm,&z);			//z=-hbar/m * f (conv) f = beta_I/hbar
		//add to eV
		muleq(&it,&z);				//z=it*beta_I/hbar
		timeq(-0.5*dt,&z);			//z=-it*dt*beta_I/(2*hbar)
		expeq(&z);					//z=exp(-it*dt*beta_I/(2*hbar))
		muleq(&z,&eV[i][j]);		//eV=exp(-it*dt*(v+omega_I+beta_I)/(2*hbar))
																				//}}
	
	
	
	
//compute H and add A to eV
	
		//compute A1
																				//for(int i=0; i<N;i++){for(int j=0; j<N;j++){
		timeq(0,&inb[i][j]);					//in=0
		plueq(&S[i][j],&inb[i][j]);			//in=S
		inb[i][j][0]-=1.0;					//in=S-1
		muleq(&FTf1[i][j],&inb[i][j]);		//in=(S-1)*FTf1=FTA1	(see notes)
	}}
	
	fftw_execute(backwardplan);				//out=A1/cb
		//add A1 to eV
	for(int i=0; i<N;i++){for(int j=0; j<N;j++){
		timeq(0,&z);						//z=0
		plueq(&outb[i][j],&z);				//z=A1/cb
		timeq(cb,&z);						//z=A1
		
		timeq(0,&A1[i][j]);					//"A1"=0
		plueq(&z,&A1[i][j]);				//"A1"=A1
		
		muleq(&z,&z);						//z=A1^2
		timeq(-1.0*hm,&z);					//z=-hbar/m *A1^2
		muleq(&it,&z);						//z=it*(-hbar/m *A1^2)
		timeq(-0.5*dt,&z);					//z=-it*dt*(-hbar^2/m *A1^2)/(2*hbar)
		expeq(&z);							//z=exp(-it*dt*(-hbar^2/m *A1^2)/(2*hbar))
		muleq(&z,&eV[i][j]);				//eV=exp(-it*dt*(v + omega_I + beta_I -hbar^2/m *A1^2)/(2*hbar))
	}}
		//numerical x-integration A1->H1
	for(int j=0; j<N;j++){
		timeq(0,&z);						//z=0
		timeq(0,&z1);						//z1=0
				//x=0
		timeq(0,&H1[0][j]);					//H1(x=0)=0 necessary for preservation of inversion symmetry
				//x>0
		for(int i=1; i<N/2;i++){
		timeq(0,&z2);						//z2=0
		plueq(&outb[i][j],&z2);				//z2=A1/cb
		plueq(&outb[i-1][j],&z2);			//z2=(A1(x)+A1(x-dx))/cb	Tapeze um größeren negativen Bereich zu kaschieren??
		timeq(0.5*dr*cb,&z2);				//z2=dr*(A1(x)+A1(x-dx))/2
		plueq(&z2,&z);						//z=~sum_i A1(i)*dr
		timeq(0,&H1[i][j]);					//H1=0
		plueq(&z,&H1[i][j]);				//H1=~sum_i A1(i)*dr -> should be OK
		}
				//x=-dr
		timeq(0,&z2);						//z2=0
		plueq(&outb[N-1][j],&z2);			//z2=A1/cb
		plueq(&outb[0][j],&z2);				//z2=(A1(-dx)+A1(0))/cb	Tapeze um größeren negativen Bereich zu kaschieren??
		timeq(-0.5*dr*cb,&z2);				//z2=-dr*(A1(-dx)+A1(0))/2		negative!!!
		plueq(&z2,&z1);						//z1=~-sum_i A1(i)*dr
		timeq(0,&H1[N-1][j]);				//H1=0
		plueq(&z1,&H1[N-1][j]);				//H1=~-sum_i A1(i)*dr -> should be OK
				//x<-dr
		for(int i=2; i<(N/2)+1;i++){
		timeq(0,&z2);						//z2=0
		plueq(&outb[N-i][j],&z2);			//z2=A1/cb
		plueq(&outb[N-i+1][j],&z2);			//z2=(A1(x)+A1(x+dx))/cb	Tapeze um größeren negativen Bereich zu kaschieren??
		timeq(-0.5*dr*cb,&z2);				//z2=dr*(A1(x)+A1(x+dx))/2
		plueq(&z2,&z1);						//z1=~-sum_i A1(i)*dr
		timeq(0,&H1[N-i][j]);					//H1=0
		plueq(&z1,&H1[N-i][j]);				//H1=~-sum_i A1(i)*dr -> should be OK
		}
	}
		
		
		//compute A2
	for(int i=0; i<N;i++){for(int j=0; j<N;j++){																			
		timeq(0,&inb[i][j]);					//in=0
		plueq(&S[i][j],&inb[i][j]);			//in=S
		inb[i][j][0]-=1.0;					//in=S-1
		muleq(&FTf2[i][j],&inb[i][j]);		//in=(S-1)*FTf2=FTA2	(see notes)
	}}
	
	fftw_execute(backwardplan);				//out=A2/cb
		//add A2 to eV
	for(int i=0; i<N;i++){for(int j=0; j<N;j++){
		timeq(0,&z);						//z=0
		plueq(&outb[i][j],&z);				//z=A2/cb
		timeq(cb,&z);						//z=A2
		
		timeq(0,&A2[i][j]);					//"A2"=0
		plueq(&z,&A2[i][j]);				//"A2"=A2
		
		muleq(&z,&z);						//z=A2^2
		timeq(-1.0*hm,&z);					//z=-hbar/m *A2^2
		muleq(&it,&z);						//z=it*(-hbar/m *A2^2)
		timeq(-0.5*dt,&z);					//z=-it*dt*(-hbar^2/m *A2^2)/(2*hbar)
		expeq(&z);							//z=exp(-it*dt*(-hbar^2/m *A2^2)/(2*hbar))
		muleq(&z,&eV[i][j]);				//eV=exp(-it*dt*(v + omega_I + beta_I -hbar^2/m *A1^2 -hbar^2/m *A2^2)/(2*hbar))
											//->eV=exp(-it*dt*(v + omega_I + beta_I -hbar^2/m *|A|^2)/(2*hbar)) -> complete
	}}
		//numerical y-integration A2->H2
	for(int i=0; i<N;i++){
		timeq(0,&z);						//z=0
		timeq(0,&z1);						//z1=0
				//y=0
		timeq(0,&H2[i][0]);					//H2(y=0)=0 necessary for preservation of inversion symmetry
				//y>0
		for(int j=1; j<N/2;j++){
		timeq(0,&z2);						//z2=0
		plueq(&outb[i][j],&z2);				//z2=A2/cb
		plueq(&outb[i][j-1],&z2);			//z2=(A2(y)+A2(y-dy))/cb	Tapeze um größeren negativen Bereich zu kaschieren??
		timeq(0.5*dr*cb,&z2);				//z2=dr*(A2(y)+A2(y-dy))/2
		plueq(&z2,&z);						//z=~sum_j A2(j)*dr
		timeq(0,&H2[i][j]);					//H2=0
		plueq(&z,&H2[i][j]);				//H2=~sum_j A2(j)*dr -> should be OK
		}
				//y=-dr
		timeq(0,&z2);						//z2=0
		plueq(&outb[i][N-1],&z2);			//z2=A2/cb
		plueq(&outb[i][0],&z2);				//z2=(A2(-dy)+A2(0))/cb	Tapeze um größeren negativen Bereich zu kaschieren??
		timeq(-0.5*dr*cb,&z2);				//z2=-dr*(A2(-dy)+A2(0))/2		negative!!!
		plueq(&z2,&z1);						//z1=~-sum_j A2(j)*dr
		timeq(0,&H2[i][N-1]);				//H2=0
		plueq(&z1,&H2[i][N-1]);				//H2=~sum_i A1(i)*dr -> should be OK
				//y<-dr
		for(int j=2; j<(N/2)+1;j++){
		timeq(0,&z2);						//z2=0
		plueq(&outb[i][N-j],&z2);			//z2=A2/cb
		plueq(&outb[i][N-j+1],&z2);			//z2=(A2(y)+A2(y+dy))/cb	Tapeze um größeren negativen Bereich zu kaschieren??
		timeq(-0.5*dr*cb,&z2);				//z2=dr*(A2(y)+A2(y+dy))/2
		plueq(&z2,&z1);						//z=~-sum_j A2(j)*dr
		timeq(0,&H2[i][N-j]);				//H2=0
		plueq(&z1,&H2[i][N-j]);				//H2=~-sum_j A2(j)*dr -> should be OK
		}
	}
		

	
	
	
	//PROPAGATE PSI
	//use in-arrays for storage
	
	
	
	//previous Psi
	for(int i=0; i<N;i++){for(int j=0; j<N;j++){
		timeq(0,&inf[i][j]);						//in=0
		plueq(&Psi[i][j],&inf[i][j]);			//in=Psi_old
	

	
	//first V-factor
		muleq(&eV[i][j],&inf[i][j]);				//in= exp(-it*dt*V/(2*hbar)) Psi_old
		
	
		
	//exp(iH2)-factor
		timeq(0,&z);							//z=0
		plueq(&H2[i][j],&z);					//z=H2
		muleq(&I,&z);							//z=i*H2
		expeq(&z);								//z=exp(i*H2)
		muleq(&z,&inf[i][j]);					//in=exp(i*H2) exp(-it*dt*V/(2*hbar)) Psi_old	
		
		inf[i][j][0]-=1.0;						// try to avoid delta peak, re-add later
	}}	
		
	
	
	
	//first T2-factor							//1D plans!!!!
	for(int i=0; i<N;i++){
		
		for(int j=0; j<N/2;j++){
			timeq(0,&in1D[j]);						//in=0
			plueq(&inf[i][j],&in1D[j]);				//in=ANS
		
			timeq(0,&in1D[N-1-j]);					//see notes
			plueq(&inf[i][j],&in1D[N-1-j]);
			
			timeq(0,&in1D[N+j]);
			plueq(&inf[i][N-1-j],&in1D[N+j]);
			
			timeq(0,&in1D[2*N-1-j]);
			plueq(&inf[i][N-1-j],&in1D[2*N-1-j]);
		}
		
		fftw_execute(forwardplan1D);
		
		for(int j=0; j<2*N;j++){
			timeq(0,&in1D[j]);						//in=0
			plueq(&out1D[j],&in1D[j]);				//in=FT[ans]/cf1D
			timeq(cf1D,&in1D[j]);					//in=FT[ans]
			
			timeq(0,&z);							//z=0
			plueq(&it,&z);							//z=it
			timeq(-0.5*hm*dt*k2[j]*k2[j],&z);		//z=-it*hbar*dt*k_y*k_y/(2m)			use k2 because of different grid!!!!
			expeq(&z);								//z=exp(-it*hbar*dt*k_y*k_y/(2m))
			muleq(&z,&in1D[j]);						//in= exp(-it*hbar*dt*k_y*k_y/(2m)) FT[ans]
		}
		
		fftw_execute(backwardplan1D);
		
		for(int j=0; j<N;j++){
			timeq(0,&inf[i][j]);					//in=0
			if(j<N/2){
				plueq(&out1D[j],&inf[i][j]);		//in=FT^(-1)[...]/cb1D
			}
			else{
				plueq(&out1D[N+j],&inf[i][j]);		//see notes!!
			}
			
			timeq(cb1D,&inf[i][j]);					//in= FT^(-1)[exp(-it*hbar*dt*k_y*k_y/(2m)) FT[ans]]
													//= exp(it*hbar*dt*d_y*d_y/(2m))*ans= exp(-it*dt*T2/(2*hbar))*ans
													//= exp(-it*dt*T2/(2*hbar)) exp(i*H2) exp(-it*dt*V/(2*hbar)) Psi_old
												
			inf[i][j][0]+=1.0;						//readd 1
													
												
	//exp(i*(H1-H2))-factor
			timeq(0,&z);							//z=0
			plueq(&H2[i][j],&z);					//z=H2
			timeq(-1.0,&z);							//z=-H2
			plueq(&H1[i][j],&z);					//z=H1-H2
			muleq(&I,&z);							//z=i*(H1-H2)
			expeq(&z);								//z=exp(i*(H1-H2))
			muleq(&z,&inf[i][j]);					//in=exp(i*(H1-H2)) exp(-it*dt*T2/(2*hbar)) exp(i*H2) exp(-it*dt*V/(2*hbar)) Psi_old
	
			inf[i][j][0]-=1.0;						// try to avoid delta peak, re-add later
		}
	
	}
	
	
	/*
	//first T2-factor
	fftw_execute(forwardplan);
		
	for(int i=0; i<N;i++){for(int j=0; j<N;j++){
		timeq(0,&inb[i][j]);						//in=0
		plueq(&outf[i][j],&inb[i][j]);			//in=FT[ans]/cf
		timeq(cf,&inb[i][j]);					//in=FT[ans]		
		
		timeq(0,&z);							//z=0
		plueq(&it,&z);							//z=it
		timeq(-0.5*hm*dt*k[j]*k[j],&z);			//z=-it*hbar*dt*k_y*k_y/(2m)
		expeq(&z);								//z=exp(-it*hbar*dt*k_y*k_y/(2m))
		muleq(&z,&inb[i][j]);					//in= exp(-it*hbar*dt*k_y*k_y/(2m)) FT[ans]
	}}	
	
	fftw_execute(backwardplan);
		
	for(int i=0; i<N;i++){for(int j=0; j<N;j++){
		timeq(0,&inf[i][j]);						//in=0
		plueq(&outb[i][j],&inf[i][j]);			//in=FT^(-1)[...]/cb
		timeq(cb,&inf[i][j]);					//in= FT^(-1)[exp(-it*hbar*dt*k_y*k_y/(2m)) FT[ans]]
												//= exp(it*hbar*dt*d_y*d_y/(2m))*ans= exp(-it*dt*T2/(2*hbar))*ans
												//= exp(-it*dt*T2/(2*hbar)) exp(i*H2) exp(-it*dt*V/(2*hbar)) Psi_old
												
		inf[i][j][0]+=1.0;						//readd 1
													
												
	//exp(i*(H1-H2))-factor
		timeq(0,&z);							//z=0
		plueq(&H2[i][j],&z);					//z=H2
		timeq(-1.0,&z);							//z=-H2
		plueq(&H1[i][j],&z);					//z=H1-H2
		muleq(&I,&z);							//z=i*(H1-H2)
		expeq(&z);								//z=exp(i*(H1-H2))
		muleq(&z,&inf[i][j]);					//in=exp(i*(H1-H2)) exp(-it*dt*T2/(2*hbar)) exp(i*H2) exp(-it*dt*V/(2*hbar)) Psi_old
	
		inf[i][j][0]-=1.0;						// try to avoid delta peak, re-add later
	}}	
	*/
	
	//T1-factor							//1D plans!!!!
	for(int j=0; j<N;j++){
		
		for(int i=0; i<N/2;i++){
			timeq(0,&in1D[i]);						//in=0
			plueq(&inf[i][j],&in1D[i]);				//in=ANS
		
			timeq(0,&in1D[N-1-i]);					//see notes
			plueq(&inf[i][j],&in1D[N-1-i]);
			
			timeq(0,&in1D[N+i]);
			plueq(&inf[N-1-i][j],&in1D[N+i]);
			
			timeq(0,&in1D[2*N-1-i]);
			plueq(&inf[N-1-i][j],&in1D[2*N-1-i]);
		}
		
		fftw_execute(forwardplan1D);
		
		for(int i=0; i<2*N;i++){
			timeq(0,&in1D[i]);						//in=0
			plueq(&out1D[i],&in1D[i]);				//in=FT[ans]/cf1D
			timeq(cf1D,&in1D[i]);					//in=FT[ans]
			
			timeq(0,&z);							//z=0
			plueq(&it,&z);							//z=it
			timeq(-hm*dt*k2[i]*k2[i],&z);			//z=-it*hbar*dt*k_x*k_x/(m)			use k2 because of different grid!!!!
			expeq(&z);								//z=exp(-it*hbar*dt*k_x*k_x/(m))
			muleq(&z,&in1D[i]);						//in= exp(-it*hbar*dt*k_x*k_x/(m)) FT[ans]
		}
		
		fftw_execute(backwardplan1D);
		
		for(int i=0; i<N;i++){
			timeq(0,&inf[i][j]);					//in=0
			if(i<N/2){
				plueq(&out1D[i],&inf[i][j]);		//in=FT^(-1)[...]/cb1D
			}
			else{
				plueq(&out1D[N+i],&inf[i][j]);		//see notes!!
			}
			
			timeq(cb1D,&inf[i][j]);					//in= FT^(-1)[exp(-it*hbar*dt*k_x*k_x/(m)) FT[ans]]
													//= exp(it*hbar*dt*d_x*d_x/(m))*ans= exp(-it*dt*T1/(hbar))*ans
													//= exp(-it*dt*T1/hbar) exp(i*(H1-H2)) exp(-it*dt*T2/(2*hbar)) exp(i*H2) exp(-it*dt*V/(2*hbar)) Psi_old
												
			inf[i][j][0]+=1.0;						//readd 1
													
												
	//exp(i*(H2-H1))-factor
		timeq(0,&z);							//z=0
		plueq(&H1[i][j],&z);					//z=H1
		timeq(-1.0,&z);							//z=-H1
		plueq(&H2[i][j],&z);					//z=H2-H1
		muleq(&I,&z);							//z=i*(H2-H1)
		expeq(&z);								//z=exp(i*(H2-H1))
		muleq(&z,&inf[i][j]);					//in= exp(i*(H2-H1)) exp(-it*dt*T1/hbar) exp(i*(H1-H2)) exp(-it*dt*T2/(2*hbar)) exp(i*H2) exp(-it*dt*V/(2*hbar)) Psi_old
	
		inf[i][j][0]-=1.0;						// try to avoid delta peak, re-add later
		}
	
	}
	

	/*
	//T1-factor
	fftw_execute(forwardplan);
		
	for(int i=0; i<N;i++){for(int j=0; j<N;j++){
		timeq(0,&inb[i][j]);						//in=0
		plueq(&outf[i][j],&inb[i][j]);			//in=FT[ans]/cf
		timeq(cf,&inb[i][j]);					//in=FT[ans]
		
		timeq(0,&z);							//z=0
		plueq(&it,&z);							//z=it
		timeq(-1.0*hm*dt*k[i]*k[i],&z);			//z=-it*hbar*dt*k_x*k_x/m
		expeq(&z);								//z=exp(-it*hbar*dt*k_x*k_x/m)
		muleq(&z,&inb[i][j]);					//in= exp(-it*hbar*dt*k_x*k_x/m) FT[ans]
	}}	
	
	fftw_execute(backwardplan);
		
	for(int i=0; i<N;i++){for(int j=0; j<N;j++){
		timeq(0,&inf[i][j]);						//in=0
		plueq(&outb[i][j],&inf[i][j]);			//in=FT^(-1)[...]/cb
		timeq(cb,&inf[i][j]);					//in= FT^(-1)[exp(-it*hbar*dt*k_x*k_x/m) FT[ans]]
												//= exp(it*hbar*dt*d_x*d_x/m)*ans= exp(-it*dt*T1/hbar)*ans
												//= exp(-it*dt*T1/hbar) exp(i*(H1-H2)) exp(-it*dt*T2/(2*hbar)) exp(i*H2) exp(-it*dt*V/(2*hbar)) Psi_old
		
		inf[i][j][0]+=1.0;						//readd 1
		
		
		
		//exp(i*(H2-H1))-factor
		timeq(0,&z);							//z=0
		plueq(&H1[i][j],&z);					//z=H1
		timeq(-1.0,&z);							//z=-H1
		plueq(&H2[i][j],&z);					//z=H2-H1
		muleq(&I,&z);							//z=i*(H2-H1)
		expeq(&z);								//z=exp(i*(H2-H1))
		muleq(&z,&inf[i][j]);					//in= exp(i*(H2-H1)) exp(-it*dt*T1/hbar) exp(i*(H1-H2)) exp(-it*dt*T2/(2*hbar)) exp(i*H2) exp(-it*dt*V/(2*hbar)) Psi_old
	
		inf[i][j][0]-=1.0;						// try to avoid delta peak, re-add later
	}}
	*/
	
	//second T2-factor							//1D plans!!!!
	for(int i=0; i<N;i++){
		
		for(int j=0; j<N/2;j++){
			timeq(0,&in1D[j]);						//in=0
			plueq(&inf[i][j],&in1D[j]);				//in=ANS
		
			timeq(0,&in1D[N-1-j]);					//see notes
			plueq(&inf[i][j],&in1D[N-1-j]);
			
			timeq(0,&in1D[N+j]);
			plueq(&inf[i][N-1-j],&in1D[N+j]);
			
			timeq(0,&in1D[2*N-1-j]);
			plueq(&inf[i][N-1-j],&in1D[2*N-1-j]);
		}
		
		fftw_execute(forwardplan1D);
		
		for(int j=0; j<2*N;j++){
			timeq(0,&in1D[j]);						//in=0
			plueq(&out1D[j],&in1D[j]);				//in=FT[ans]/cf1D
			timeq(cf1D,&in1D[j]);					//in=FT[ans]
			
			timeq(0,&z);							//z=0
			plueq(&it,&z);							//z=it
			timeq(-0.5*hm*dt*k2[j]*k2[j],&z);		//z=-it*hbar*dt*k_y*k_y/(2m)			use k2 because of different grid!!!!
			expeq(&z);								//z=exp(-it*hbar*dt*k_y*k_y/(2m))
			muleq(&z,&in1D[j]);						//in= exp(-it*hbar*dt*k_y*k_y/(2m)) FT[ans]
		}
		
		fftw_execute(backwardplan1D);
		
		for(int j=0; j<N;j++){
			timeq(0,&inf[i][j]);					//in=0
			if(j<N/2){
				plueq(&out1D[j],&inf[i][j]);		//in=FT^(-1)[...]/cb1D
			}
			else{
				plueq(&out1D[N+j],&inf[i][j]);		//see notes!!
			}
			
			timeq(cb1D,&inf[i][j]);					//in= FT^(-1)[exp(-it*hbar*dt*k_y*k_y/(2m)) FT[ans]]
													//= exp(it*hbar*dt*d_y*d_y/(2m))*ans= exp(-it*dt*T2/(2*hbar))*ans
													//= exp(-it*dt*T2/(2*hbar)) exp(i*(H2-H1)) exp(-it*dt*T1/hbar) exp(i*(H1-H2)) exp(-it*dt*T2/(2*hbar)) exp(i*H2) exp(-it*dt*V/(2*hbar)) Psi_old 
												
			inf[i][j][0]+=1.0;						//readd 1
													
												
	//exp(-iH2)-factor
		timeq(0,&z);							//z=0
		plueq(&H2[i][j],&z);					//z=H2
		timeq(-1.0,&z);							//z=-H2
		muleq(&I,&z);							//z=-i*H2
		expeq(&z);								//z=exp(-i*H2)
		muleq(&z,&inf[i][j]);					//in=exp(-i*H2) exp(-it*dt*T2/(2*hbar)) exp(i*(H2-H1)) exp(-it*dt*T1/hbar) exp(i*(H1-H2)) exp(-it*dt*T2/(2*hbar)) exp(i*H2) exp(-it*dt*V/(2*hbar)) Psi_old 
	
	
	
	//second V-factor
		muleq(&eV[i][j],&inf[i][j]);			//in= exp(-it*dt*V/(2*hbar)) exp(-i*H2) exp(-it*dt*T2/(2*hbar)) exp(i*(H2-H1)) exp(-it*dt*T1/hbar) exp(i*(H1-H2)) exp(-it*dt*T2/(2*hbar)) exp(i*H2) exp(-it*dt*V/(2*hbar)) Psi_old 
												//= U(dt) Psi_old =Psi_new
		}
	
	}
	
	/*
	//second T2-factor
	fftw_execute(forwardplan);
		
	for(int i=0; i<N;i++){for(int j=0; j<N;j++){
		timeq(0,&inb[i][j]);						//in=0
		plueq(&outf[i][j],&inb[i][j]);			//in=FT[ans]/cf
		timeq(cf,&inb[i][j]);					//in=FT[ans]
		
		timeq(0,&z);							//z=0
		plueq(&it,&z);							//z=it
		timeq(-0.5*hm*dt*k[j]*k[j],&z);			//z=-it*hbar*dt*k_y*k_y/(2m)
		expeq(&z);								//z=exp(-it*hbar*dt*k_y*k_y/(2m))
		muleq(&z,&inb[i][j]);					//in= exp(-it*hbar*dt*k_y*k_y/(2m)) FT[ans]
	}}	
	
	fftw_execute(backwardplan);
		
	for(int i=0; i<N;i++){for(int j=0; j<N;j++){
		timeq(0,&inf[i][j]);						//in=0
		plueq(&outb[i][j],&inf[i][j]);			//in=FT^(-1)[...]/cb
		timeq(cb,&inf[i][j]);					//in= FT^(-1)[exp(-it*hbar*dt*k_y*k_y/(2m)) FT[ans]]
												//= exp(it*hbar*dt*d_y*d_y/(2m))*ans= exp(-it*dt*T2/(2*hbar))*ans
												//= exp(-it*dt*T2/(2*hbar)) exp(i*(H2-H1)) exp(-it*dt*T1/hbar) exp(i*(H1-H2)) exp(-it*dt*T2/(2*hbar)) exp(i*H2) exp(-it*dt*V/(2*hbar)) Psi_old 
	
		inf[i][j][0]+=1.0;						//readd 1
	
	
	
	//exp(-iH2)-factor
		timeq(0,&z);							//z=0
		plueq(&H2[i][j],&z);					//z=H2
		timeq(-1.0,&z);							//z=-H2
		muleq(&I,&z);							//z=-i*H2
		expeq(&z);								//z=exp(-i*H2)
		muleq(&z,&inf[i][j]);					//in=exp(-i*H2) exp(-it*dt*T2/(2*hbar)) exp(i*(H2-H1)) exp(-it*dt*T1/hbar) exp(i*(H1-H2)) exp(-it*dt*T2/(2*hbar)) exp(i*H2) exp(-it*dt*V/(2*hbar)) Psi_old 
	
	
	
	//second V-factor
		muleq(&eV[i][j],&inf[i][j]);				//in= exp(-it*dt*V/(2*hbar)) exp(-i*H2) exp(-it*dt*T2/(2*hbar)) exp(i*(H2-H1)) exp(-it*dt*T1/hbar) exp(i*(H1-H2)) exp(-it*dt*T2/(2*hbar)) exp(i*H2) exp(-it*dt*V/(2*hbar)) Psi_old 
												//= U(dt) Psi_old =Psi_new
	}}
	*/
	
	
	//BRUTE FORCE SYMMETRIZATION
	for(int i=0; i<N;i++){
		i2=N-i;
		if(i2==N){i2=0;}
		
		for(int j=0; j<N/2;j++){
			
		j2=N-j;
		if(j2==N){j2=0;}
			
		timeq(0,&z);							//z=0
		plueq(&Psi[i][j],&z);					//z=Psi(x,y)
		plueq(&Psi[i2][j2],&z);					//z=Psi(x,y)+Psi(-x,-y)
		timeq(0.5,&z);							//z=(Psi(x,y)+Psi(-x,-y))/2
		timeq(0,&Psi[i][j]);					//Psi(x,y)=0
		plueq(&z,&Psi[i][j]);					//Psi(x,y)=z
		timeq(0,&Psi[i2][j2]);					//Psi(-x,-y)=0
		plueq(&z,&Psi[i2][j2]);					//Psi(-x,-y)=z
	}}
	
	
	
	//COMPARE
	x=maxerror(inf,Psi);
	cout<<"max change during step "<<step<<" ="<<x<<endl;
	
	
	
	//WRITE OUT ENERGY
	Efile<<step<<" "<<step*dt<<" "<<E[0]<<" "<<E[1]<<endl;
	
	//WRITE OF-T-FILES
	if(step%250 == 0){
		
		goftfile<<endl;
		Softfile<<endl;
		phioftfile<<endl;
		Hoftfile<<endl;
		Aoftfile<<endl;
		
		goftfile<<endl;
		Softfile<<endl;
		phioftfile<<endl;
		Hoftfile<<endl;
		Aoftfile<<endl;
		
		
		for(int i=N/2; i<3*N/2;i++){
			i2=i;
			if(i2>(N-1)){i2-=N;}
			
		goftfile<<step*dT<<" "<<r[i2]<<" "<<Psi[i2][0][0]*Psi[i2][0][0]+Psi[i2][0][1]*Psi[i2][0][1]<<" "<<Psi[0][i2][0]*Psi[0][i2][0]+Psi[0][i2][1]*Psi[0][i2][1]<<endl;
		Softfile<<step*dT<<" "<<k[i2]<<" "<<S[i2][0][0]<<" "<<S[i2][0][1]<<" "<<S[0][i2][0]<<" "<<S[0][i2][1]<<endl;
		phioftfile<<step*dT<<" "<<r[i2]<<" "<<atan2(Psi[i2][0][1],Psi[i2][0][0])<<" "<<atan2(Psi[0][i2][1],Psi[0][i2][0])<<endl;
		Hoftfile<<step*dT<<" "<<r[i2]<<" "<<H1[i2][0][0]<<" "<<H1[i2][0][1]<<" "<<H2[i2][0][0]<<" "<<H2[i2][0][1]<<" "<<H1[0][i2][0]<<" "<<H1[0][i2][1]<<" "<<H2[0][i2][0]<<" "<<H2[0][i2][1]<<endl;
		Aoftfile<<step*dT<<" "<<r[i2]<<" "<<A1[i2][0][0]<<" "<<A1[i2][0][1]<<" "<<A2[i2][0][0]<<" "<<A2[i2][0][1]<<" "<<A1[0][i2][0]<<" "<<A1[0][i2][1]<<" "<<A2[0][i2][0]<<" "<<A2[0][i2][1]<<endl;
	}
	}
	
	//if(step%2500 == 0){
	if(step%1000 == 0){
		
		gfile2<<endl;
		Sfile2<<endl;
		phifile2<<endl;
		Hfile2<<endl;
		Afile2<<endl;
		
		gfile2<<endl;
		Sfile2<<endl;
		phifile2<<endl;
		Hfile2<<endl;
		Afile2<<endl;
	
		for(int i=N/2; i<3*N/2;i++){
			i2=i;
			if(i2>(N-1)){i2-=N;}
			
			gfile2<<endl;
			Sfile2<<endl;
			phifile2<<endl;
			Hfile2<<endl;
			Afile2<<endl;
			
			
		for(int j=N/2; j<3*N/2;j++){
			j2=j;
			if(j2>(N-1)){j2-=N;}
			
		gfile2<<step*dT<<" "<<r[i2]<<" "<<r[j2]<<" "<<Psi[i2][j2][0]*Psi[i2][j2][0]+Psi[i2][j2][1]*Psi[i2][j2][1]<<endl;
		Sfile2<<step*dT<<" "<<r[i2]<<" "<<r[j2]<<" "<<S[i2][j2][0]<<" "<<S[i2][j2][1]<<endl;
		phifile2<<step*dT<<" "<<r[i2]<<" "<<r[j2]<<" "<<atan2(Psi[i2][j2][1],Psi[i2][j2][0])<<endl;
		Hfile2<<step*dT<<" "<<r[i2]<<" "<<r[j2]<<" "<<H1[i2][j2][0]<<" "<<H1[i2][j2][1]<<" "<<H2[i2][j2][0]<<" "<<H2[i2][j2][1]<<endl;
		Afile2<<step*dT<<" "<<r[i2]<<" "<<r[j2]<<" "<<A1[i2][j2][0]<<" "<<A1[i2][j2][1]<<" "<<A2[i2][j2][0]<<" "<<A2[i2][j2][1]<<endl;
	}}
	
	
	}
	
	
	
	
	//UPDATE PSI
	for(int i=0; i<N;i++){for(int j=0; j<N;j++){
		timeq(0,&Psi[i][j]);					//Psi=0
		plueq(&inf[i][j],&Psi[i][j]);			//Psi=in=Psi_new
	}}
	
	
	
	
	
	
	//CHANGE TO REAL TIME IF NECESSARY
	if(change==1){if(step==(imagsteps-1)){
		imev=0;			//formally change to real evolution (for v update)
		it[0]=1.0;		//change propagation exponent from i to 1
		it[1]=0;
		dt=dT;			//change timestep
		
		cout<<"changing to real evolution!!!!!!!!!!!!!!!!!"<<endl;
		}}											
												
	//missing: break condition	
	}	
	//END ITERATION
	
	
	
	
	//write output files
	for(int i=0; i<N;i++){for(int j=0; j<N;j++){
		//outputfile<<i<<" "<<j<<" "<<Psi[i][j][0]<<" "<<Psi[i][j][1]<<endl;
		outputfile<<Psi[i][j][0]<<" "<<Psi[i][j][1]<<endl;
	}}
	
	
	
	for(int i=0; i<N;i++){for(int j=0; j<N;j++){
		Psifile<<r[i]<<" "<<r[j]<<" "<<Psi[i][j][0]<<" "<<Psi[i][j][1]<<endl;
	}}
	
	for(int i=0; i<N;i++){
		if(i==N/2){gfile<<endl;
			//gfile2<<endl;
			}
		gfile<<r[i]<<" "<<Psi[i][0][0]*Psi[i][0][0]+Psi[i][0][1]*Psi[i][0][1]<<" "<<Psi[0][i][0]*Psi[0][i][0]+Psi[0][i][1]*Psi[0][i][1]<<endl;
		//gfile2<<sqrt(2.0)*r[i]<<" "<<Psi[i][i][0]*Psi[i][i][0]+Psi[i][i][1]*Psi[i][i][1]<<endl;
	}
	for(int i=0; i<N;i++){
		if(i==N/2){Sfile<<endl;
			//Sfile2<<endl;
			}
		Sfile<<k[i]<<" "<<S[i][0][0]<<" "<<S[i][0][1]<<" "<<S[0][i][0]<<" "<<S[0][i][1]<<endl;
		//Sfile2<<sqrt(2.0)*k[i]<<" "<<S[i][i][0]<<" "<<S[i][i][1]<<endl;
	}
	
	
	
	
	for(int i=0; i<N;i++){
		if(i==N/2){vfile<<endl;}
		vfile<<r[i]<<" "<<v[i][0]<<endl;
	}
	
	
	
	//destroy plans
	fftw_destroy_plan(forwardplan);
	fftw_destroy_plan(backwardplan);
	
	fftw_destroy_plan(forwardplan1D);
	fftw_destroy_plan(backwardplan1D);
	
	return 0;
}
