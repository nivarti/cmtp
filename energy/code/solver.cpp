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
  
  Domain.EvaluateFluxes();				   // Evaluate Flux Integrals for problem 1, 2 
  Domain.EvaluateSourceTerms();					   // Evaluate Source Terms for problem 1, 2
  
  cout<<"\nOutputting Flux Values using Exact Functions...\n";
  Domain.PrintFluxes();
  Domain.PrintSources();

  Domain.FluxVerification();
  Domain.SourceVerification();
  
}


void SolveEnergyEquation(Grid& Domain, double tFinal){
  
  int n = 0, N;
  double dt, t = 0.0, dT;
  Field dU;
  Domain.EvaluateCellCoordinates();				   // Evaluate cell coordinates
  Domain.EvaluateInitialFields();
  
  Domain.EvaluateSourceTerms();

  //Domain.EvaluateBoundaryConditions();
  //Domain.EvaluateFluxes();
  
  //Domain.PrintFluxes();
  //Domain.PrintSources();
  dt = Domain.EvaluateTimeStep(ImplicitEuler);

  do{    
    
    dU = Domain.EulerImplicitTimeAdvance();
    
    cout<<"Maximum change in solution: "<<dU.T<<endl;
    n++;

    
  }while(dU.T > 0.000001);
    
  //Domain.PrintFluxes();
  //Domain.PrintSources();

  cout<<"\nSolution converged in "<<n<<" steps\n";

}

void SolveThomas(double LHS[NMAX][3], double RHS[NMAX],
		 const int iSize)
{
  int i;
  LHS[0][0] = LHS[iSize+1][2] = 0;
  /* Forward elimination */
  for (i = 0; i < iSize+1; i++) {
    LHS[i][2] /= LHS[i][1];
    RHS[i] /= LHS[i][1];
    LHS[i+1][1] -= LHS[i][2]*LHS[i+1][0];
    RHS[i+1] -= LHS[i+1][0]*RHS[i];
  }
  /* Last line of elimination */
  RHS[iSize+1] /= LHS[iSize+1][1];

  /* Back-substitution */
  for (i = iSize; i >= 0; i--) {
    RHS[i] -= RHS[i+1]*LHS[i][2];
  }
}

void CopyToLHS(double** Dx, double LHS [NMAX][3], const int Size){

  for(int i = 0; i <= Size; i++){
    
    LHS[i][0] = Dx[i][0];
    LHS[i][1] = Dx[i][1];
    LHS[i][2] = Dx[i][2];

  }

}

void CopyToRHS(double** FI, double RHS[NMAX], const int Size, const int J, Direction RC){

  if(RC == Column)
    for(int i = 0; i <= Size; i++){
      
      RHS[i] = FI[i][J];
      RHS[i] = FI[i][J];
      RHS[i] = FI[i][J];
      
    }
  
  else if(RC == Row)
    for(int i = 0; i <= Size; i++){
      
      RHS[i] = FI[J][i];
      RHS[i] = FI[J][i];
      RHS[i] = FI[J][i];
      
    }
}

void CopyFromRHS(double** FI, double RHS[NMAX], const int Size, const int I, Direction RC){

  if(RC == Row)
    for (int i = 0; i <= Size; i++) {
      
      FI[I][i] = RHS[i];
      
    }
  else if(RC == Column)
    for (int i = 0; i <= Size; i++) {
      
      FI[i][I] = RHS[i];
      
    }

}

/// / Solve energy equation for problems 5, 6, and 7
// void SolveEnergyEquation(Grid &Domain, double tFinal){
  
//   int n = 0, N;
//   double t = 0.0, dt = 0.0;
//   Field dU;		       
  
//   Domain.EvaluateCellCoordinates();				   // Calculate coordinates and store in cells for entire grid  
//   Domain.EvaluateInitialFields();				   // Calculate initial fields
  
//   Domain.EvaluateBoundaryConditions();
//   // // Domain.PrintFieldValues(Temperature);
//   // // Domain.PrintFieldValues(xVelocity);
//   // // Domain.PrintFieldValues(yVelocity);

//   // // Domain.EvaluateBoundaryConditions();
//   // // Domain.EvaluateFluxes();
//   // Domain.EvaluateFluxes();
//   // Domain.EvaluateSourceTerms();  
  
//   // Domain.PrintSources();
//   // Domain.PrintFieldValues(Temperature);
//   // //Domain.PrintFieldValues(xVelocity);
//   // //Domain.PrintFieldValues(yVelocity);
  
//   // Domain.PrintFluxes();

//   //Domain.PrintFieldValues(yVelocity, 1, -1);
//   //Domain.PrintFieldValues(yVelocity, 25, -1);
  
//   //dt = Domain.EvaluateTimeStep(ExplicitEuler);		   // Calculate stable time step  
  
//   //dt = 0.1;
  
  
//   dt = Domain.EvaluateTimeStep(ExplicitEuler);
//   // dU = Domain.EulerExplicitTimeAdvance();
  
//   // Domain.EvaluateFluxes();
//   // //Domain.EvaluateSourceTerms();  
  
//   // Domain.PrintFluxes();

//   // Domain.PrintSources();

//   Domain.EvaluateSourceTerms();
  
  
//   N = tFinal/dt;		                                   // Modify dt to 
//   dt = tFinal/N;		                                   // make sure you march to tFinal
  
//   do{
    
//   dU = Domain.EulerExplicitTimeAdvance();

//   Domain.PrintFieldValues(Temperature, 24, -1);
//   // Domain.PrintFieldValues(xVelocity);
//   // Domain.PrintFieldValues(yVelocity);

  
//   Domain.PrintFluxes();
//   Domain.PrintSources();
    
//   // t += dt;
//   n++;
    
//   //   cout<<"\nSolution marched to time, t = "<<t;         
    
//   }while(n < N);			                   
        
//   cout<<"\nFinal time reached: "<<t<<" in "<<n<<" time steps";
  
// }



// // void MarchTime(Grid& Domain, TimeScheme TS){
  
// //   switch(TS){

// //   case 0: Domain.EulerExplicitTimeAdvance();
// //     break;
// //   // case 1: Domain.EulerImplicitTimeAdvance();
// //   //   break;
// //   // case 2: Domain.RK2ExplicitTimeAdvance();
// //   //   break;
// //   default:
// //     break;
// //   }  

// // }

// // void SolveRK2Explicit(Grid &Domain){
        
// //   double dt = Domain.dt;

// //   EvaluateBoundaryConditions();	                                   // Implement ghost cell values for flux calculation  
// //   EvaluateFluxes();	         	                   // Evaluate Fluxes from solution at time step n. Change to Limited flux as required
 
// //   for(int i = Imin; i <= Imax; i++){
// //     for(int j = Jmin; j <= Jmax; j++){
      
// //       Mesh[i][j].Ubuf.T = Mesh[i][j].U.T + dt*Mesh[i][j].FI;	   // Calculate intermediate step solution, store in separate array
// //     }
// //   }
  
// //   EvaluateBoundaryConditions();	                                   // Repeat procedure, use fluxes at n+1/2 to calculate solution at n+1 
// //   EvaluateFluxes();				   
  
// //   for(int i = Imin; i <= Imax; i++){
// //     for(int j = Jmin; j <= Jmax; j++){
      
// //       Mesh[i][j].T = Mesh[i][j].T + dt*Mesh[i][j].FI;		   // Calculate intermediate step solution, store in separate array

// //     }  
// //   }
  
// // }
