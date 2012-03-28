//////////////////////////////////////////////////////////
// The following code forms a 2D Energy Equation Solver //
// i.e. Solution of				        //
// dT/dt + udT/dx + vdT/dy = 1/Re*(d2T/dx2 + d2T/dy2)   //
// 						        //
// With given boundary conditions and initial condition //
//////////////////////////////////////////////////////////

#include "header.h"

int main(){  
  
  int MeshX = 25, MeshY = 10;                                      // Set Initial Mesh Size
  double tFinal = 1.0;						   // Set tFinal for transient calculations
  clock_t ti, tf, RunTime;			                   // Set clock variables
  ti = clock();			                                   // Start clock
  
  cout.precision(5);
  cout.width(10);
  
  while(MeshX <= 25){
    
    Grid Domain(MeshX, MeshY, 1);          	                   // Define Grid of required size  
    
    cout<<"\nEvaluating solution at Time: "<<tFinal<<" for Mesh Size: "<<MeshX<<" by "<<MeshY; 
    
    //EvaluateGridParameters(Domain);				   // Set Grid Values for problem 1, 2
    SolveEnergyEquation(Domain, tFinal);			   // Solve energy equation for problem 5
    
    Domain.PrintFieldValues(Temperature);
    //Domain.PrintFieldValues(Temperature, MeshX - 1, -1);        
    
    MeshX *= 2;							   // Double mesh size
    MeshY *= 2;
  }
  
  tf = clock();
  RunTime = (double)(tf - ti)/CLOCKS_PER_SEC;                      // Calculate Run Time in Seconds 
  
  cout<<"\nSolver Run-Time: "<<RunTime<<" seconds\n";  
  return 0;
}
