////////////////////////////////////////////////////
// This file contains the solver for 1D Laplacian //
// Stead State Diffusion Problem                  //
//   d2T/dx2 = 0                                  //
// with boundary conditions:		          //
//	1. Constant temperature			  //
//	2. Constant flux			  //
// 	                                          //
// Written by: gvn	                 	  //
////////////////////////////////////////////////////

// Include Relevant Headers
#include<fstream>
#include<iostream>
#include<iomanip>
#include<math.h>

// Define major parameters
#define MESH_SIZE 10						   // Mesh size includes only cells within the domain
#define L 1.0							   // Dimension of domain
#define M 0.000001						   // Tolerance for maximum change of solution per iteration

const int SIZE = MESH_SIZE + 2;	                                   // Declares size of array including ghost cells
const int Xmax = MESH_SIZE;                                        // Iterations go from 1 to Xmax (inside the mesh)
const double dx = L/(MESH_SIZE);				   // Set increment

using namespace std;

void InitializeField(double T[]);                                  // Initialize Fields to 0.0                   
void DirichletBC(double T[]);                                      // Change Coefficient Matrix to Accommodate
void NeumannBC(double T[]);                                        // Change Coefficient Matrix to Accommodate
void SolveLaplace(double T[], double Tol);		           // Laplacian Solver for Specified Tolerance

double GaussSeidel(double T[]);					   // Gauss Seidel iterative solver

void PrintData(double T[]);			                   // Print Solution
void WriteToFile(int iterations, double errorNorm);		   // Write Numerical Scheme data to file

void ComputeExactSolution(double T[]);  			   // Compute Exact Solution
double AssessError(double T[]);					   // Calculates L2Norm of Error


int main(){							   
  double  T[SIZE];	                 			   // Declare T field, Error and Iterations
  
  InitializeField(T);						   // Initialize Temperature field to zero  
  SolveLaplace(T, M);			                	   // Run Laplacian Solver till M decimal Tolerance is achieved     
  PrintData(T); 						   // Print solution to screen

  return 0;							   
} // End of Main


void PrintData(double T[]){

  int i;
  
  for(i = 1; i <= Xmax; i++){
    cout<<setw(10)<<setprecision(6)<<T[i]<<" ";                    // Formated output of array
  } 
  cout<<endl;  
} // End Print Solution


void InitializeField(double T[]){
  
  int i;
  
  for(i = 0; i < SIZE; i++){     
    T[i] = 0.0;							   // Set 0. for Temperature
  }   
  cout<<"Temperature Field Initialized...\n";
} // End of Initialization of Fields


void SolveLaplace(double T[], double Tol){
  
  double dTmax = 0., errorNorm;
  int iter = 0;
  
  do{
    
    iter++;
    DirichletBC(T);		                                   // Impose Dirichlet BC at x = 0
    NeumannBC(T);						   // Impose constant flux BC at x = L
    dTmax = GaussSeidel(T);
    
  }while(dTmax > Tol);
    
  errorNorm = AssessError(T);					   // Calculate L2 Norm of Solution
  WriteToFile(iter, errorNorm);					   // Write iterations and L2 Norm to file
  
} // End of Solve Laplace


void DirichletBC(double T[]){  
  
  double Tw = 100.0;
  T[0] = 2*Tw - T[1];						   // Implement Dirichlet BC
  
} // End of Boundary Condition

void NeumannBC(double T[]){

  double Fw = 100.0;
  T[SIZE - 1] =  Fw*dx + T[SIZE - 2];				   // Implement Neumann BC
  
}

double GaussSeidel(double T[]){
  
  int i;
  double dMax = 0.0, Tnew;
  
  for(i = 1; i <= Xmax; i++){
    
    Tnew = (T[i-1] + T[i+1])/2.0;				   // Assign update value to buffer variable   
    
    if(dMax < fabs(Tnew - T[i])){				   // Assign maximum change to return variable
      dMax = fabs(Tnew - T[i]);
    }
    T[i] = Tnew;						   // Assign update value to array element
  }
  
  return dMax;							   // Return maximum change for every iteration
} // End of Gauss Seidel Scheme


void ComputeExactSolution(double Te[]){

  int i;
  double x;
  
  cout<<"\nExact solution computed...\n";
  for(i = 1; i <= Xmax; i++){      
    x = (2*i-1)/2.*dx;
    Te[i] = 100*x +100 ;					   // Exact solution of 1D Heat Conduction without Source    
  }
  PrintData(Te);

} // End of ComputeExactSolution

double AssessError(double T[SIZE]){

  int i;
  double L2Norm, Error[SIZE], Te[SIZE];

  ComputeExactSolution(Te);
  
  for(i = 1; i <= Xmax; i++){          
    Error[i] = T[i] - Te[i];
    L2Norm += Error[i]*Error[i];
  }
  
  L2Norm = sqrt(L2Norm/SIZE);
  
  return L2Norm;
} // End of Assess Error

void WriteToFile(int iterations, double errorNorm){
  ofstream file;

  file.open("log");
  file<<"Number of Iterations = "<<iterations
      <<", L2Norm = "<<errorNorm<<endl;				   // Print Iterations and L2 Norm
  file.close();

}// End of File Write
