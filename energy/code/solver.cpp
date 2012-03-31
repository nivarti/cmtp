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


void SolveEnergyEquation(Grid& Domain){
  
  int n = 0;
  Field dU;
  Domain.EvaluateCellCoordinates();				   // Evaluate cell coordinates
  Domain.EvaluateInitialFields();
  //Domain.EvaluateExactFields();  
  //Domain.EvaluateBoundaryConditions();
  Domain.EvaluateSourceTerms();
  
  //Domain.EvaluateBoundaryConditions();
  //Domain.EvaluateFluxes();
  
  //Domain.PrintFluxes();
  //Domain.PrintSources();
  
  Domain.EvaluateTimeStep(ImplicitEuler);
  
  do{    
    
    dU = Domain.EulerImplicitTimeAdvance();
    n++;
    
  }while(dU.T > 0.00000001);
  
  Domain.EvaluateBoundaryConditions();
  Domain.EvaluateGradients();
  
  //Domain.PrintFluxes();
  //Domain.PrintSources();
  cout<<"\nMaximum change in solution: "<<dU.T<<endl;      
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


void CalculateErrorBound (){
  
  double p, e, e32, e21, phi1, phi2, phi3, r21, ea21, gci, phi21;	     // Declare all variables needed 
  ifstream errf;
  
  errf.open("phi_v");						      // Open file with error profiles for 3 different meshes 
  
  errf>>phi3
      >>phi2
      >>phi1;
  // Show Mesh Error Data already Recorded
  cout<<"Mesh Size: 0.00625"<<", Phi_1 = "<<phi3<<endl;
  cout<<"Mesh Size: 0.0125"<<", Phi_2 = "<<phi2<<endl;
  cout<<"Mesh Size: 0.025"<<", Phi_3 = "<<phi1<<endl;
  
  e32 = phi3 - phi2;
  e21 = phi2 - phi1;
  e = e32/e21;
  
  r21 = 2.0;
  ea21 = fabs((phi1 - phi2)/phi1);
  
  p = fabs(log(fabs(e))/log(r21));			      // Calculate apparent order using given formula 
  phi21 = (pow(r21,p)*phi1 - phi2)/(pow(r21,p) - 1);		      // Calculate extropalted value of P 
  gci = 1.25*ea21/(pow(r21,p) - 1);				      // Calculate GCI_fine to get error bound on solution 
  
  cout<<"\nApparent Order for solution in Position: "<< setprecision(10) << p<<endl;
  cout<<"Solution estimate for Position = "<< setprecision(10) <<phi21<<" +/- "<< setprecision(10)<<fabs(gci)<<endl;
  
}// End of ASME Estimate of Error bound
