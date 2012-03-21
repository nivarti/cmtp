///////////////////////////////////////////////////////////
// This file contains the 1D Wave Equation Solver	 //
// i.e. the problem:					 //
// 							 //
//         dT/dt + udT/dx = 0				 //
// 	0 < x < 1, 0 < t				 //
// 	given fixed u					 //
// 							 //
// Using a second order upwind flux evaluation		 //
// and a two-stage Runge Kutta time advance	         //
//                                               **gvn   //
///////////////////////////////////////////////////////////
		       
#include<iostream>
#include<cmath>
#include<iomanip>
#include<fstream>

#define MESH_SIZE 40				        	   // Define a mesh size
#define CFL 0.528 						   // Define the cfl number
#define L 1.0							   // Domain size
#define pi 3.14159265358979323846264338327950288419716939937510

const double u = 2.0;						   // Flow field with fixed velocity
const double dx = L/MESH_SIZE;					   // Element size
const double dt = CFL*dx/u;			        	   // Permitted maximum time step

const int Imin = 2;				                   // Helps fix 
const int Imax = MESH_SIZE + 1;   				   // loop lengths

const int SIZE = MESH_SIZE + 4;					   // Array size for definitions
enum Limiter{Unlimited, Superbee};				   // Convenience by using enumerate...can add other schemes later!

using namespace std;

void InitializeField(double T[]);	                           // Initialise the sine wave 
void InitializePulseField(double T[]);			           // Initialise square wave

void PeriodicBC(double T[]);					   // Implement ghost cells

void EvaluateFlux(Limiter Scheme, double T[], double F[]);	   // Bunch of flux evaluation functions
void EvaluateUnlimitedFlux(double T[], double F[]);
void EvaluateLimitedFlux(double T[], double F[]);
double SuperBee(double);

void RK2TimeAdvance(double T[]);                           	   // 
void SolveWaveEquation(double T[], double);			   // Solve equation for given time

void ComputeExactSolution(double Te[], double);			   // Compute exact solution
double CalculateError(double T[], double Te[], double E[]);	   // Error calculating function

void CopyToField(double T1[], double T2[]);			   // Copies T1 data to T2
void CreateArray(double T[]);
void PrintData(double T[]);
void WriteToFile(double T[], string F);

// _______________________________________________________________________________________________________ //

int main(){
  
  double T[SIZE], Te[SIZE], E[SIZE], error;
  ofstream file;
  
  file.open("N3", ios::app);
  InitializeField(T);		                                   // Initialize temperature field   
  ComputeExactSolution(Te, 1.0);				   // Compute exact solution
  //InitializePulseField(Te);
  //InitializePulseField(T);
  
  SolveWaveEquation(T, 2.0);					   // Calculate Wave Solution at t = 1.0
  error = CalculateError(T, Te, E);				   // Calculate error norm, and write to file

  //WriteToFile(T, "Tuld");
  //WriteToFile(Te, "Ted");
  //WriteToFile(E, "N");
  
  cout<<"\nError (L2 Norm): "<<setprecision(15)<<error<<endl;	   // Report error
  file<<dt<<" "<<setprecision(15)<<error<<endl;
  file.close();
  return 0;
}


// _______________________________________________________________________________________________________ //


void CreateArray(double T[]){
  
  int i;
  
  for(i = 0; i < SIZE; i++){    
    T[i] = 0.0;
  }
} // Create Array with Zeros 


void InitializeField(double T[]){
  
  int i;
  double x;

  CreateArray(T);

  for(i = Imin; i <= Imax; i++){    
    x = dx*(2*(i - Imin) + 1)/2.0;	
    T[i] = sin(2*pi*x);						   // Calculate the initial condition T(x,0)
  }

  cout<<"\nTemperature Field Initialized...";
} // End Initialize Fields


void  InitializePulseField(double T[]){
  
  int i;
  double x;
  
  CreateArray(T);

  for(i = Imin; i <= Imax; i++){
    x = dx*(2*(i - Imin) + 1)/2.0;	
    
    if(x >= 0.25 && x <= 0.5){
      
      T[i] = -1.0;
    }
    else
      if(x > 0.5 && x <= 0.75){
	
	T[i] = 1.0;
      }
    
      else
	T[i] = 0.0;
  }
  cout<<"\nPulse field initialized...";


}

void SolveWaveEquation(double T[], double tFinal){
  
  int i, n = 0, N;
  double t = 0.0;
  N = tFinal/dt + 1;
  
  do{
           
    RK2TimeAdvance(T);
    t += dt;
    n++;

    //cout<<"\nRunge-Kutta Solution marched to time, t = "<<t; 
    
    if(n < N && (t - tFinal) > dt/2.0){
    n++;
  }
  
  else if(n == N && (t - tFinal) < -dt/2.0){
    n--;
  }
    
  }while(n < N);				   // Make sure you reach the final time not less not more
        
  cout<<"\nFinal time reached: "<<t<<" in "<<n<<" time steps";
  
} // End Wave Solver


void RK2TimeAdvance(double T[]){
    
  int i;
  double T1[SIZE], F1[SIZE], F[SIZE];				   // Buffer arrays for RK2 intermediate steps
  
  PeriodicBC(T);	                                           // Implement ghost cell values for flux calculation  
  EvaluateFlux(Unlimited, T, F);			         	   // Evaluate Fluxes from solution at time step n. Change to Limited flux as required
 
  for(i = Imin; i <= Imax; i++){
    
    T1[i] = T[i] - (dt/(dx*2.0))*(F[i] - F[i-1]);		   // Calculate intermediate step solution, store in separate array
  }
  
  PeriodicBC(T1);		                                   // Repeat procedure, use fluxes at n+1/2 to calculate solution at n+1 
  EvaluateFlux(Unlimited, T1, F1);				   

  CopyToField(T, T1);

  for(i = Imin; i <= Imax; i ++){
    
    T[i] = T[i] - (dt/dx)*(F1[i] - F1[i-1]);
  }  
}


void PeriodicBC(double T[]){
  
  T[0] = T[MESH_SIZE];			                           // Achtung!
  T[1] = T[MESH_SIZE + 1];        				   // Copy ghost cell values to get periodic values
  T[Imax + 1] = T[2]; 					           // For flux limiting, the wave is copied infront
  T[Imax + 2] = T[3];						   // so as to find values of r for last few cells

}

// _______________________________________Flux Evaluation_________________________________________________________________ //


void EvaluateFlux(Limiter Scheme, double T[], double F[]){

  int i;     
  double r, psi;
  
  for(i = 0; i <= Imax; i++){    
    
    switch(Scheme){
      
    case Superbee: 
      r = (T[i+2] - T[i+1])/(T[i+1] - T[i]);	                   // Evaluate r based on adjacent values of Temperature     
      psi = SuperBee(r);					   // Assign psi from Superbee scheme
      break;
      
    default: psi = 1.0;						   // Do not bother calculating r for unlimited,
      break;							   // Assign psi = 1.0 directly 
    }
    
    F[i + 1] = u*(T[i+1] + psi/2.0*(T[i+1] - T[i]));		   // Generic flux evaluation, valid for all schemes
  }	
}


double SuperBee(double r){
  
  double psi;

  if (r <= 0){ psi = 0.0; }					   // Assign value to psi based on 
  else if(r > 0 && r <= 0.5){ psi = 2*r; }			   // region of r, which changes 
  else if(r > 0.5 && r <= 1.0){ psi = 1.0; }			   // with each control volume
  else if(r > 1.0 && r <= 2.0){ psi = r; }			   // Assign corners to the portion
  else{ psi = 2.0; }						   // less in magnitude 
  
  return psi;  
}

// ____________________________________________Solution Computations__________________________________________________________ //

void ComputeExactSolution(double Te[], double t){
  
  int i;
  double x;
  
  CreateArray(Te);

  for(i = Imin; i <= Imax; i++){
    x = dx*(2*(i - Imin) + 1)/2.0;	
    Te[i] = sin(2*pi*(x - 2*t));    
  }
  
  cout<<"\nExact solution calculated...";
} // End computation of exact solution


double CalculateError(double T[], double Te[], double E[]){

  int i;
  double L2norm = 0.0;

  for(i = Imin; i <= Imax; i++){					   // Copies entire array (including ghosts)
    E[i] = T[i] - Te[i];
    L2norm += E[i]*E[i];
  }

  L2norm = sqrt(L2norm/MESH_SIZE);

  return L2norm;
}


//________________________________________________Other Functions____________________________________________________________//

void PrintData(double T[]){

  int i;
  double x;
  
  cout<<"\nPrinting data...\n";
  for(i = Imin; i <= Imax; i++){     
    x = dx*(2*(i - Imin) + 1)/2.0;
    cout<<"T["<<i<<"] = T(x = "<<x<<") = "<< setw(20) << setprecision(15) <<T[i]<<"\n";		   // Formatted output of 1D array
  } 
  
  cout<<endl;  

} // End Print Solution


void CopyToField(double T1[], double T2[]){

  int i;

  for(i = 0; i <= Imax; i++){					   // Copies entire array (including ghosts)
    T2[i] = T1[i];
  }

}

void WriteToFile(double T[], string F){
  
  int i;
  double x;
  ofstream plotf;
  
  plotf.open(F.c_str());	                                      // Plot in a file with name F 
  
  for(i = Imin; i <= Imax; i++){          
    
    x = (2*(i-Imin) + 1)/2.0*dx;
      
    plotf<<x<<" "<< setprecision(15) << T[i]<<endl;			      // Plot at each point in domain 
  }
  
  plotf.close();						      // Close file
}// End of File Write

