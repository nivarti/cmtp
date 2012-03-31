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
  clock_t ti, tf, RunTime;			                   // Set clock variables
  
  cout.precision(5);
  cout.width(10);
  
  //CalculateErrorBound();
  
  while(MeshX <= 400){
    ti = clock();			                                   // Start clock
   
    Grid Domain(MeshX, MeshY, 1);          	                   // Define Grid of required size  
    
    cout<<"\nEvaluating solution for Mesh Size: "<<MeshX<<" by "<<MeshY; 
    
    //EvaluateGridParameters(Domain);				   // Set Grid Values for problem 1, 2
    //SolveEnergyEquation(Domain);         			   // Solve energy equation for problem 5        
    //Domain.PrintFieldValues(Temperature);
    //Domain.FieldVerification();
    //Domain.PrintFieldValues(xVelocity);        
    // Domain.PrintFieldValues(yVelocity);        
    
    MeshX *= 2;							   // Double mesh size
    MeshY *= 2;
    
    tf = clock();
    RunTime = (double)(tf - ti)/CLOCKS_PER_SEC*1000;                      // Calculate Run Time in Seconds 
    
    cout<<"\nSolver Run-Time: "<<RunTime<<" ms\n";     
    
  }
  

  // ifstream file;
  
  // file.open("../plots/order/EtEE");
  // double x[3], y[3], slope;
  // int i = 0;
  
  // while(i < 3){
    
  //   file>>x[i]>>y[i];    
  //   i++;

  // };

  // slope = (log(y[2]) - log(y[1]))/(log(x[2]) - log(x[1]));
  // cout<<"\nOrder of method: "<<slope;

  // file.close();

  return 0;
}
