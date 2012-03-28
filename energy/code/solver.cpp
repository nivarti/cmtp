#include "header.h"

// For problem 1, evaluate exact integrals and compare with numerical values...
void EvaluateGridParameters(Grid &Domain){
  
  Domain.EvaluateCellCoordinates();				   // Evaluate cell coordinates
  Domain.EvaluateExactIntegrals();                                 // Evaluate exact integrals for problem 1, and 2
  
  //Domain.PrintCellCoordinates();
  // cout<<"\nOutputting Field Values for T...\n";
  // Domain.PrintFieldValues(1);
  // cout<<"\nOutputting Field Values for u...\n";
  // Domain.PrintFieldValues(2);
  // cout<<"\nOutputting Field Values for v...\n";
  // Domain.PrintFieldValues(3);  
  
  Domain.EvaluateFluxIntegrals();				   // Evaluate Flux Integrals for problem 1, 2 
  Domain.EvaluateSourceTerms();					   // Evaluate Source Terms for problem 1, 2
  
  cout<<"\nOutputting Flux Values using Exact Functions...\n";
  Domain.PrintFluxes();
  Domain.PrintSources();

  Domain.FluxVerification();
  Domain.SourceVerification();
  
}

// Solve energy equation for problems 5, 6, and 7
void SolveEnergyEquation(Grid &Domain, double tFinal){
  
  int n = 0, N;
  double t = 0.0, dt;
  Field dU;		       
  
  Domain.EvaluateCellCoordinates();				   // Calculate coordinates and store in cells for entire grid  
  Domain.EvaluateInitialFields();				   // Calculate initial fields
  Domain.EvaluateBoundaryConditions();
  // Domain.PrintFieldValues(Temperature);
  // Domain.PrintFieldValues(xVelocity);
  // Domain.PrintFieldValues(yVelocity);

  // Domain.EvaluateBoundaryConditions();
  // Domain.EvaluateFluxIntegrals();
  Domain.EvaluateFluxIntegrals();
  Domain.EvaluateSourceTerms();  
  
  Domain.PrintSources();
  Domain.PrintFieldValues(Temperature);
  //Domain.PrintFieldValues(xVelocity);
  //Domain.PrintFieldValues(yVelocity);
  
  Domain.PrintFluxes();

  //Domain.PrintFieldValues(yVelocity, 1, -1);
  //Domain.PrintFieldValues(yVelocity, 25, -1);
  
  //dt = Domain.EvaluateTimeStep(ExplicitEuler);		   // Calculate stable time step  
  
  //dt = 0.1;
  
  
  Domain.EvaluateTimeStep(ExplicitEuler);
  dU = Domain.EulerExplicitTimeAdvance();
  
  Domain.EvaluateFluxIntegrals();
  //Domain.EvaluateSourceTerms();  
  
  Domain.PrintFluxes();

  Domain.PrintSources();

  
  
  
  // N = tFinal/dt;		                                   // Modify dt to 
  // dt = tFinal/N;		                                   // make sure you march to tFinal
  
  // do{
    
  // dU = Domain.EulerExplicitTimeAdvance();

  // Domain.PrintFieldValues(Temperature);
  // Domain.PrintFieldValues(xVelocity);
  // Domain.PrintFieldValues(yVelocity);

  
  // Domain.PrintFluxes();

    
  // // t += dt;
  // n++;
    
  // //   cout<<"\nSolution marched to time, t = "<<t;         
    
  // }while(n < N);			                   
        
  // cout<<"\nFinal time reached: "<<t<<" in "<<n<<" time steps";
  
}


// void MarchTime(Grid& Domain, TimeScheme TS){
  
//   switch(TS){

//   case 0: Domain.EulerExplicitTimeAdvance();
//     break;
//   // case 1: Domain.EulerImplicitTimeAdvance();
//   //   break;
//   // case 2: Domain.RK2ExplicitTimeAdvance();
//   //   break;
//   default:
//     break;
//   }  

// }

// void SolveRK2Explicit(Grid &Domain){
        
//   double dt = Domain.dt;

//   EvaluateBoundaryConditions();	                                   // Implement ghost cell values for flux calculation  
//   EvaluateFluxIntegrals();	         	                   // Evaluate Fluxes from solution at time step n. Change to Limited flux as required
 
//   for(int i = Imin; i <= Imax; i++){
//     for(int j = Jmin; j <= Jmax; j++){
      
//       Mesh[i][j].Ubuf.T = Mesh[i][j].U.T + dt*Mesh[i][j].FI;	   // Calculate intermediate step solution, store in separate array
//     }
//   }
  
//   EvaluateBoundaryConditions();	                                   // Repeat procedure, use fluxes at n+1/2 to calculate solution at n+1 
//   EvaluateFluxIntegrals();				   
  
//   for(int i = Imin; i <= Imax; i++){
//     for(int j = Jmin; j <= Jmax; j++){
      
//       Mesh[i][j].T = Mesh[i][j].T + dt*Mesh[i][j].FI;		   // Calculate intermediate step solution, store in separate array

//     }  
//   }
  
// }
