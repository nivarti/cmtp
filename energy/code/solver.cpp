#include "header.h"

void EvaluateGridParameters(Grid &Domain){
  
  Domain.EvaluateCellCoordinates();    
  Domain.EvaluateExactIntegrals();   
  
  //Domain.PrintCellCoordinates();
  // cout<<"\nOutputting Field Values for T...\n";
  // Domain.PrintFieldValues(1);
  // cout<<"\nOutputting Field Values for u...\n";
  // Domain.PrintFieldValues(2);
  // cout<<"\nOutputting Field Values for v...\n";
  // Domain.PrintFieldValues(3);
  
  Domain.EvaluateFluxIntegrals();
  Domain.EvaluateSourceTerms();
  
  // cout<<"\nOutputting Flux Values using Exact Functions...\n";
  // Domain.PrintFluxes();
  // Domain.PrintSources();
  
  Domain.EvaluateL2Norm();
  
}

// void SolveRK4Explicit(Grid &Domain){
      
//   int i,  j;  
  
//   //PeriodicBC(T);	                                           // Implement ghost cell values for flux calculation  
//   EvaluateFluxIntegral(Unlimited, T, F);	         	   // Evaluate Fluxes from solution at time step n. Change to Limited flux as required
 
//   for(i = Imin; i <= Imax; i++){
//     for(j = Jmin; j <= Jmax; j++){
      
//       Mesh[i][j].T = Mesh[i][j].T + dt*Mesh[i][j].FI;		   // Calculate intermediate step solution, store in separate array
//     }
//   }
  
//   //PeriodicBC(T1);		                                   // Repeat procedure, use fluxes at n+1/2 to calculate solution at n+1 
//   EvaluateFluxIntegrals();				   

//   CopyToField(T, T1);

//   for(i = Imin; i <= Imax; i ++){
    
//     T[i] = T[i] - (dt/dx)*(F1[i] - F1[i-1]);
//   }  
// }


// }

// void SolveWaveEquation(double T[], double tFinal){
  
//   int i, n = 0, N;
//   double t = 0.0;
//   N = tFinal/dt + 1;
  
//   do{
           
//     RK2TimeAdvance(T);
//     t += dt;
//     n++;

//     //cout<<"\nRunge-Kutta Solution marched to time, t = "<<t; 
    
//     if(n < N && (t - tFinal) > dt/2.0){
//     n++;
//   }
  
//   else if(n == N && (t - tFinal) < -dt/2.0){
//     n--;
//   }
    
//   }while(n < N);				   // Make sure you reach the final time not less not more
        
//   cout<<"\nFinal time reached: "<<t<<" in "<<n<<" time steps";
  
// } // End Wave Solver


