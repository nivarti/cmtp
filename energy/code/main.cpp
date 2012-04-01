//////////////////////////////////////////////////////////
// The following code forms a 2D Energy Equation Solver //
// i.e. Solution of				        //
// dT/dt + udT/dx + vdT/dy = 1/Re*(d2T/dx2 + d2T/dy2)   //
// 						        //
// With given boundary conditions and initial condition //
//////////////////////////////////////////////////////////

#include "header.h"

int main(){  
  
  int MeshX = 200, MeshY = 80;                                      // Set Initial Mesh Size
  clock_t ti, tf, RunTime;			                   // Set clock variables
  
  setprecision(15);
  setw(15);   
  
  while(MeshX <= 200){
    ti = clock();			                                   // Start clock
   
    Grid Domain(MeshX, MeshY, 1);          	                   // Define Grid of required size  
    
    cout<<"\nEvaluating solution for Mesh Size: "<<MeshX<<" by "<<MeshY; 
    
    //EvaluateGridParameters(Domain);				   // Set Grid Values for problem 1, 2
    SolveEnergyEquation(Domain);         			   // Solve energy equation for problem 5            
    
    MeshX *= 2;							   // Double mesh size
    MeshY *= 2;
    
    tf = clock();
    RunTime = (double)(tf - ti)/CLOCKS_PER_SEC*1000;                      // Calculate Run Time in Seconds 
    
    cout<<"\nSolver Run-Time: "<<RunTime<<" ms\n";     
    
  }
  
  //CalculateErrorBound();
  //CalculateOrder();
  
  return 0;
}
