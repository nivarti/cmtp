//////////////////////////////////////////////////
// This file contains 2D Poisson Solver	        //
// for a steady state diffusion problem	        //
//   					        //
//    d2T/dx2 + d2T/dy2 = S		        //
//   					        //
// Following Boundary Conditions can be handled //
//   1. Dirichlet			        //
//   2. von Neumann			        //
//   3. Mixed				        //
// 					        //
//                                    **gvn     //
//////////////////////////////////////////////////

#include<fstream>
#include<iostream>
#include<iomanip>
#include<math.h>
#include<string>

#define MESH_SIZE 10						      // Mesh includes physical domain of the problem                     

#define Lx 1.0							      // Specify size of domain in X direction
#define Ly 1.0							      // Specify size of domain in Y direction
#define M 0.00000001					              // Tolerance for maximum change in solution 
#define pi 3.14159265358979323846264338327950288419716939937510

const int SIZE = MESH_SIZE + 2;	                                      // Set size of matrix to include ghost cells
const int Xmax = MESH_SIZE, Ymax = MESH_SIZE;			      // Iterations go from 1 to Xmax (inside the mesh)
const double dx = Lx/(MESH_SIZE);                                     // Calculate dx
const double dy = Ly/(MESH_SIZE);				      // Calculate dy

using namespace std;						   

void InitializeField(double T[][SIZE]);	                              // Initialize Fields to 0.0
                   
void DirichletBC_Laplace(double T[][SIZE]);			      // Dirichlet Boundary conditions 
void DirichletBC_Poisson(double P[][SIZE]);			      
void NeumannBC_Laplace(double T[][SIZE]);                      	      // Neumann Boundary conditions
void NeumannBC_Poisson(double P[][SIZE]);			      

void EvaluateSource(double S[][SIZE]);				      // Evaluate source term at each control volume center

double GaussSeidel(double T[][SIZE]);				      // Gauss Seidel single iteration
double GaussSeidel(double T[][SIZE], double w);   		      // Gauss Seidel single iteration with SOR
double GaussSeidel(double T[][SIZE], double S[][SIZE], double w);      // Gauss Seidel single iteration for Poisson with SOR

void SolveLaplace(double T[][SIZE], double Tol);	              // Laplacian Solver for Specified Tolerance
void SolvePoisson(double P[][SIZE], double S[][SIZE], double Tol);    // Poisson Solver for Specified Tolerance

void PrintData(double T[][SIZE]);                 		      // Print A|b augmented matrix
void WriteToFile(double T[][SIZE], string F);

double CalculateFieldValue(double P[][SIZE]);			      // Interpolate value to a point
void CalculateExactSolution(double Te[][SIZE]);  		      // Compute Exact Solution
double CalculateError(double T[][SIZE], string F);		      // Calculates L2Norm of Error
void CalculateErrorASME();					      // Calculate an error bound and apparent order 


int main(){							   
  double  T[SIZE][SIZE], P[SIZE][SIZE], S[SIZE][SIZE];                // Declare T field, Error and Iterations
  
  //Solve Laplace equation for Temperature
  SolveLaplace(T, M); 
  
  //Solve Poisson equation for Pressure  
  //SolvePoisson(P, S, M);			               	      // Run Laplacian Solver till M decimal Tolerance is achieved   
  //CalculateErrorASME();		                                      // Run the Poisson solver for 3 different meshes, then comment the line and ...
                                                                      // Run CalculateErrorASME to get error estimates   
  return 0;							   
} // End of Main



/* I) Solution Strategy  <-------------------------------------------------------------------------------------

   1. Initialize fields */
void InitializeField(double T[][SIZE]){
  
  int i,j;
  
  for(i = 0; i < SIZE; i++){     
    for(j = 0; j < SIZE; j++){
      T[i][j] = 0.0;		                                      // Set 0.0 as initial guess for Temperature
    }
  }   
  cout<<"Field Initialized...\n";			      
  
} // End of Initialization of Fields


/* 2. Implement Boundary Conditions */
void DirichletBC_Laplace(double T[][SIZE]){  
    int j;
  double x;
  
  for(j = 1; j <=Xmax; j++){
    x = (2*j-1)/2.0*dx;
    T[0][j] = - T[1][j];
    T[SIZE - 1][j] = 2*cos(pi*x) - T[SIZE - 2][j];		      // Set desired cosine temperature profile for wall
  }

} // End of Dirichlet Boundary Condition

void DirichletBC_Poisson(double T[][SIZE]){  
  
  int i,j;
  double x,y, Tw = 5.0;
  
  for(j = 1; j <= Xmax; j++){
    x = (2*j - 1.0)/2.0*dx;
    T[SIZE - 1][j] = 2*Tw - (1+x*x)*(1+x*x)*(1+x*x) - T[SIZE - 2][j];   
  }

  for(i = 1; i <= Ymax; i++){
    y = (2*i - 1)/2.0*dy;
    T[i][SIZE - 1] = 2*Tw - (1+y*y)*(1+y*y)*(1+y*y) - T[i][SIZE - 2]; // Set desired value by tweaking ghost cells
  }


} // End of Dirichlet Boundary Condition

void NeumannBC_Laplace(double T[][SIZE]){
  int i;
  
  for (i = 1; i <= Ymax; i++){
    T[i][0] = T[i][1];            
    T[i][SIZE-1] = T[i][SIZE-2];				      // Set Conditions for Poisson Solver
  }  
  
}// End of Neumman BC Function

void NeumannBC_Poisson(double T[][SIZE]){

  int i;
  
  for (i = 1; i <= Ymax; i++){
    T[i][0] =  T[i][1];            
    T[0][i] =  T[1][i];	                 			      // Set fluxes for Poisson Solver
  }    
}// End of Neumman BC for Poisson


/* 3. Evaluate Source Terms */
void EvaluateSource(double S[][SIZE]){
  int i,j;
  double x, y;
  
  for(i = 1; i <=Ymax; i++){
    y = (2*i-1.0)/2.0*dy;

    for(j = 1; j <=Xmax; j++){
      x = (2*j-1.0)/2.0*dx;

      S[i][j] = (x*x + y*y);
      S[i][j] = -18*S[i][j]*S[i][j];                                  // Evaluate Source as S = -18*(x^2 + y^2)^2 (sink, really!)
    }
  }
}// End of Evaluate Source for Poisson


/* 4. Solve Equations 
 a) Laplacian Solver */
void SolveLaplace(double T[][SIZE], double Tol){
  
  ofstream file; 
  double T2[SIZE][SIZE], T3[SIZE][SIZE];			      // Declare arrays for convergence tests
  double dTmax = 0.0, dTmax2 = 0.0, dTmax3 = 0.0, errorNorm;	      // Declare variables    
  int iter = 0, flag = 0;
  
  //file.open("SOR");						      // Output SOR Test with different values of omega 
  file.open("Norm", ios::app);					      // Append mesh norm data 
    
  InitializeField(T);						      // Initialize Temperature field to zero to start with   
  
  // Convergence comparison variables
  //InitializeField(T2);
  //InitializeField(T3);
  
  do{
    iter++;
    // // Set Ghost Cell Values
    DirichletBC_Laplace(T);
    NeumannBC_Laplace(T);    

    // // Perform single Gauss Seidel Iteration
    //dTmax = GaussSeidel(T);	                                      // Perform Gauss Seidel without relaxation (Problem 2)
    dTmax = GaussSeidel(T, 1.5);                                      // Perform Gauss Seidel with relaxation 1.5 (Problem 3)
    
    // // Comparisons of convergence characteristics    
    // DirichletBC_Laplace(T2);
    // NeumannBC_Laplace(T2);    
    // dTmax2 = GaussSeidel(T2, 1.3);

    // DirichletBC_Laplace(T3);
    // NeumannBC_Laplace(T3);    
    // dTmax3 = GaussSeidel(T3, 1.5);				      // Used for compilation of errors for meshes 40x40 and 80x80 
    
    // file << iter << " " << dTmax <<" "<< dTmax2 <<" "<< dTmax3 << endl;    

    // Report Early Convergence: Use for SOR Test
    // if(dTmax3 < Tol && flag<1){      
    //   cout<<"Solution with relaxation (1.5) converged after "<<iter<<" iterations"<<endl;
    //   flag++;
    // }
    // else if(dTmax2 < Tol && flag<2){      
    //   cout<<"Solution with relaxation (1.3) converged after "<<iter<<" iterations"<<endl;
    //   flag++;
    // }
    
    // Hand verification toy code
    // cout<<"At the end of "<<iter<<" iterations, solution characteristics:\n";    
    // cout<<"Maximum Change in Solution: "<<dTmax<<endl<<endl;        
    
  }while(dTmax > Tol);						      // Exit when maximum change is less than tolerance 
  
  //cout<<"Finite Volume Solution of Laplacian Computed...\n"<<endl;  
  cout<<"Solution converged after "<<iter<<" iterations\n";
  
  //PrintData(T);
  //WriteToFile(T,"T0");					      // If need be to plot Temperature  
  
  //Calculate Norms and Errors
  errorNorm  = CalculateError(T, "E0");		                      // Solution to Problem 2. a) <-------------------------------------------
  // errorNorm = CalculateError(T2, "E1p3");
  // errorNorm = CalculateError(T3, "E1p5");

  file <<dx<<" "<<errorNorm<<endl; 
  file.close();
} // End of Solve Laplace

/* b) Poisson Solver */
void SolvePoisson(double P[][SIZE], double S[][SIZE], double Tol){
  double dTmax = 0.0;                              		      // Declare variables
  int iter = 0;
  ofstream errf;
  
  errf.open("phi", ios::app);					      // Record phi for 
  
  InitializeField(P);						      // Initialize Pressure Field to 0.0
  EvaluateSource(S);						      // Calculatee Source Terms in the domain

  do{
    iter++;
    // Set Ghost Cell Values
    DirichletBC_Poisson(P);
    NeumannBC_Poisson(P);    
    // Record maximum change per iteration of Gauss Seidel
    dTmax = GaussSeidel(P, S, 1.5);  
    //errf << iter <<" "<< dTmax <<endl;
    
  }while(dTmax > Tol);						      // Exit when maximum change is less than tolerance

  cout<<"Finite Volume Solution Computed for Poisson Equation after "<<iter<<" iterations..."<<endl;  
  
  //PrintData(P);
  WriteToFile(P, "P");					      // If need be to compute Pressure 

  errf <<MESH_SIZE<<" "<< setprecision(10) << CalculateFieldValue(P)<<"\n";  
  errf.close();
 
}// End of Poisson Solver


/* 5. Calculate Errors */
void CalculateExactSolution(double Te[][SIZE]){

  int i,j;
  double x,y;
    
  for(i = 1; i <= Xmax; i++){				               
    y = (2*i-1)/2.0*dy;						      // Sweep in y
    
    for(j = 1; j <= Ymax; j++){
      x = (2*j-1)/2.0*dx;					      // Sweep in x
      Te[i][j] = cos(pi*x)*sinh(pi*y)/sinh(pi);			      // Exact Solution given. Computes to 3 decimal accuracy      
    }
  }
  
  cout<<"Exact Solution Computed...\n"<<endl;
  
} // End of ComputeExactSolution

double CalculateError(double T[][SIZE], string F){

  int i,j;
  double x,y, L2Norm, Te[SIZE][SIZE], Error[SIZE][SIZE];
      
  CalculateExactSolution(Te);					      // Compute the exact solution 
  
  for(i = 1; i <= Xmax; i++){
    y = (2*i-1)/2.0*dy;						      // Sweep in y
    
    for(j = 1; j <= Ymax; j++){
      x = (2*j-1)/2.0*dx;					      // Sweep in x
      Error[i][j] = T[i][j] - Te[i][j];                               // Exact Solution given. Computes to 3 decimal accuracy      
      L2Norm += Error[i][j]*Error[i][j];	    
    }
  }  
  
  L2Norm = sqrt(L2Norm/(Xmax*Ymax));				      //  Finish computing l2 norm Prob 2a)  <------------------------------
  WriteToFile(Error, F);
  cout<<"Errors Calculated...\nError Norm (l2): "<<L2Norm<<"\n";
    
  return L2Norm;
} // End of ComputeExactSolution

double CalculateFieldValue(double P[][SIZE]){

  int i, j;
  double phi;
  
  i = (0.5 + dx)/dx + 0.5;
  j = (0.5 + dy)/dy + 0.5;
  phi = P[i][j] + P[i-1][j] + P[i-1][j-1] + P[i][j-1];  
  phi = phi/4.0;

  return phi;  
}

void CalculateErrorASME(){
  
  int m1,m2,m3;
  double p, e,e32, e21, phi1, phi2, phi3, r21, ea21, gci, phi21;	     // Declare all variables needed 
  ifstream errf;

  errf.open("phi");						      // Open file with error profiles for 3 different meshes 
  
  errf>>m3>>phi3
      >>m2>>phi2
      >>m1>>phi1;
  // Show Mesh Error Data already Recorded
  cout<<"Mesh Size: "<<m3<<", Phi_1 = "<<phi3<<endl;
  cout<<"Mesh Size: "<<m2<<", Phi_2 = "<<phi2<<endl;
  cout<<"Mesh Size: "<<m1<<", Phi_3 = "<<phi1<<endl;
  
  e32 = phi3 - phi2;
  e21 = phi2 - phi1;
  e = e32/e21;
  
  r21 = 2;
  ea21 = fabs((phi1 - phi2)/phi1);

  p = fabs(log(fabs(e))/log(r21));			      // Calculate apparent order using given formula 
  phi21 = (pow(r21,p)*phi1 - phi2)/(pow(r21,p) - 1);		      // Calculate extropalted value of P 
  gci = 1.25*ea21/(pow(r21,p) - 1);				      // Calculate GCI_fine to get error bound on solution 
  
  cout<<"\nApparent Order for solution in P: "<< setprecision(10) << p<<endl;
  cout<<"Solution estimate for P = "<< setprecision(10) << phi21<<" +/- "<< setprecision(10)<<fabs(gci)<<endl;
  
}// End of ASME Estimate of Error bound


/* II) Numerical Schemes ------------------------------------------------------------------------------------->
   1. Gauss Seidel Solvers */
double GaussSeidel(double T[][SIZE]){
  
  int i, j;
  double dTmax = 0.0, dT;    
  
  for(i = Ymax; i >= 1; i--){
    for(j = Xmax; j >= 1; j--){					      // Iterate beginning from highest values first

      dT = T[i-1][j] + T[i+1][j] + T[i][j-1] + T[i][j+1];	   
      dT = dT/4.0 - T[i][j];				              // Update Value of T
    
      if(dTmax < fabs(dT)){
	dTmax = fabs(dT);
      }
      
      T[i][j] += dT;
    }
  }
        
  return dTmax;
} // End of Gauss Seidel without SOR


double GaussSeidel(double T[][SIZE], double w){

  int i, j;
  double dTmax = 0.0, dT;  
  
  for(i = Ymax; i >= 1; i--){
    for(j = Xmax; j >= 1; j--){					      // Iterate beginning from highest values first
       
      dT = T[i-1][j] + T[i+1][j] + T[i][j-1] + T[i][j+1]; 
      dT = dT/4.0 - T[i][j];           
      dT = (double)w*dT;
      
      if(dTmax < fabs(dT)){
	dTmax = fabs(dT);
      }
      
      T[i][j] += dT;
    }
  }  
  
  return dTmax;
} // End of Gauss Seidel with SOR


double GaussSeidel(double T[][SIZE], double S[][SIZE], double w){

  int i, j;
  double dTmax = 0.0, dT = 0.0;  
  
  for(j = Xmax; j >= 1; j--){					      // Iterate beginning from highest values first
    for(i = Ymax; i >= 1; i--){
      
      dT = T[i-1][j] + T[i+1][j] + T[i][j-1] + T[i][j+1];	      // Record delta, and multiply by SOR constant 
      dT = dT/4.0 - S[i][j]/4.0*dx*dx - T[i][j];           
      dT = w*dT;
      
      if(dTmax < fabs(dT)){
	dTmax = fabs(dT);
      }
      
      T[i][j] += dT;
    }
  }  
  
  return dTmax;
} // End of Gauss Seidel with SOR


/* III) Other functions used in the code <-------------------------------------------------------------------------------------
     1. Printing functions */

void PrintData(double T[][SIZE]){

  int i,j;
  
  for(i = Ymax; i >= 1; i--){     
    for(j = 1; j <= Xmax; j++){
      cout<< setw(10) << setprecision(10) <<T[i][j]<<" ";	      // Formatted output of 2D array
    } 
    cout<<endl;							       
  }   
  
  cout<<endl;  
} // End Print Solution

void WriteToFile(double Array[][SIZE], string F){
  
  int i,j;
  double x,y;
  ofstream plotf;
  
  plotf.open(F.c_str());	                                      // Plot in a file with name F 
  
  for(i = 1; i <= Xmax; i++){          
    y = (2*i-1)/2.0*dy;
    
    for(j = 1; j <= Ymax; j++){
      x = (2*j-1)/2.0*dx;
      
      plotf<<x<<" "<<y<<" "<< setprecision(10) << Array[i][j]<<endl;			      // Plot at each point in domain 
    }
    plotf<<endl;
  }
  
  plotf.close();						      // Close file
}// End of File Write


