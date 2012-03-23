//////////////////////////////////////////////////////////
// The following code forms a 2D Energy Equation Solver //
// i.e. Solution of				        //
// dT/dt + udT/dx + vdT/dy = 1/Re*(d2T/dx2 + d2T/dy2)   //
// 						        //
// With given boundary conditions and initial condition //
//////////////////////////////////////////////////////////

#include "solver.h"

int main(){  
  clock_t ti, tf, RunTime;
  ti = clock();
  

  Grid Domain(10,10, 1);	       // Define Grid of required size  
  Domain.EvaluateCellData();
  Domain.PrintFieldValues(3);
  

  
  tf = clock();
  RunTime = (double)(tf - ti)/CLOCKS_PER_SEC;

  cout<<"\nSolver Run-Time: "<<RunTime<<" seconds\n";  
  return 0;
}
