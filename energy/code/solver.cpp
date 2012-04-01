#include "header.h"

// For problem 1,2 evaluate exact integrals and 
// compare with numerical values...
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
  
  Domain.EvaluateFluxes();			         	   // Evaluate Flux Integrals for problem 1, 2 
  Domain.EvaluateSourceTerms();			 		   // Evaluate Source Terms for problem 1, 2
  
  cout<<"\nOutputting Flux Values using Exact Functions...\n";
  Domain.PrintFluxes();
  Domain.PrintSources();

  Domain.FluxVerification();
  Domain.SourceVerification();
  
}


/* 
 * Solve the energy equation using appropriate time marching
 * scheme. 
 */

void SolveEnergyEquation(Grid& Domain){
  
  int n = 0;
  Field dU;

  Domain.EvaluateCellCoordinates();				   // Evaluate cell coordinates
  Domain.EvaluateInitialFields();
  Domain.EvaluateSourceTerms();
  
  Domain.EvaluateTimeStep(ImplicitEuler);
  
  do{    
    
    dU = Domain.EulerImplicitTimeAdvance();
    n++;    
    
  }while(fabs(dU.T) > 0.00000001);

  //Domain.EvaluateExactFields();      
  //Domain.FieldVerification();  
  Domain.PrintFieldValues(Temperature);

  cout<<"\nMaximum change in solution: "<<dU.T;          
  cout<<"\nSolution converged in "<<n<<" steps";
  
  Domain.EvaluateBoundaryConditions();
  Domain.EvaluateGradients();
  
  //Domain.PrintFluxes();
  //Domain.PrintSources();

}

/*
 * Calculate order of method from the slope of log plot
 * Calculate Error Bounds using ASME Method
 * use function from Laplace code
 */

void CalculateOrder(){
  
  ifstream file;
  
  file.open("ErrorNorm");
  double x[3], y[3], slope;
  int i = 0;
  
  while(i < 3){
    
    file>>x[i]>>y[i];    
    i++;

  };

  slope = (log(y[2]) - log(y[1]))/(log(x[2]) - log(x[1]));
  cout<<"\nOrder of method: "<<slope;

  file.close();

}


void CalculateErrorBound (){
  
  double p, e, e32, e21, phi1, phi2, phi3, r21, ea21, gci, phi21;     // Declare all variables needed 
  ifstream errf;
  
  errf.open("phi_p");						      // Open file with error profiles for 3 different meshes 
  
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
  
  p = fabs(log(fabs(e))/log(r21));			              // Calculate apparent order using given formula 
  phi21 = (pow(r21,p)*phi1 - phi2)/(pow(r21,p) - 1);		      // Calculate extropalted value of P 
  gci = 1.25*ea21/(pow(r21,p) - 1);				      // Calculate GCI_fine to get error bound on solution 
  
  cout<<"\nApparent Order for solution in Position: "<< setprecision(10) << p<<endl;
  cout<<"Solution estimate for Position = "<< setprecision(10) <<phi21<<" +/- "<< setprecision(10)<<fabs(gci)<<endl;
  
}// End of ASME Estimate of Error bound
