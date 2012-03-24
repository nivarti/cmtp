//////////////////////////////////////////////////////////
// The following code forms a 2D Energy Equation Solver //
// i.e. Solution of				        //
// dT/dt + udT/dx + vdT/dy = 1/Re*(d2T/dx2 + d2T/dy2)   //
// 						        //
// With given boundary conditions and initial condition //
//////////////////////////////////////////////////////////

#include "solver.h"

int main(){  
  
  int x = 10;
  clock_t ti, tf, RunTime;				           //  
  ti = clock();			                                   // Start clock
  
  while(x <= 80){
    
    Grid Domain(x, x, 1);	                                   // Define Grid of required size  
    
    cout<<"\nEvaluating for Mesh Size: "<<x<<" by "<<x;
    EvaluateGridParameters(Domain);
    
    x *= 2;							   // Double mesh size
  }

  tf = clock();
  RunTime = (double)(tf - ti)/CLOCKS_PER_SEC;

  cout<<"\nSolver Run-Time: "<<RunTime<<" seconds\n";  
  return 0;
}
