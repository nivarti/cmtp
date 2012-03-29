#include "header.h"


Field::Field(){
  
  T = 0.0;
  u = 0.0;
  v = 0.0;

}

Field& Field::operator=(const Field &RHS){

  T = RHS.T;
  u = RHS.u;
  v = RHS.v;

  return *this;

}

Field& Field::operator+=(const Field &RHS){
  
  T += RHS.T;
  u += RHS.u;
  v += RHS.v;
  
  return *this;

}

Cell::Cell(){
  
  x = 0.0;
  y = 0.0;

  eF = 0.0;
  eS = 0.0;

  S = 0.0;
  F = 0.0;
  
}


void Cell::SetCellCoordinates(double X, double Y){
  
  x = X;
  y = Y;
  
}


void Cell::SetCellField(double U1, double U2, double U3){

  U.T = U1;
  U.u = U2;
  U.v = U3;

}

// Compute Exact Values of Fields from Formulae
void Cell::ComputeExactField(){
  
  U.T = T0*cos(pi*x)*sin(pi*y);  
  U.u = u0*y*sin(pi*x);
  U.v = v0*x*cos(pi*y);

}

void Cell::SetInitialField(){

  U.T = y;
  U.u = 6.0*Uav*y*(1 - y);	                                   // Exact fully developed profile 
  U.v = 0.0;			                                   // Pertinent to our real problem 

}
  
double Cell::SetFDTemperatureField(){
  
  double a = 0.0;

  a = 1 - 2.0*y;
  a = a*a*a*a;
  
  U.T = y + 0.75*Pr*Ec*Uav*Uav*(1 - a);
  
  return U.T;
}

void Cell::ComputeExactFluxIntegral(){

  eF = u0*T0*pi*cos(2*pi*x)*y*sin(pi*y) 
    + v0*T0*pi*x*cos(pi*x)*cos(2*pi*y) 
    + 1/(Re*Pr)*2*T0*pi*pi*cos(pi*x)*sin(pi*y);	 	 
  eF = -eF;			/* flip sign while taking integral to RHS */
  
}


void Cell::ComputeExactSourceTerm(){
      
  eS  = 2*(u0*pi*cos(pi*x)*y)*(u0*pi*cos(pi*x)*y);
  eS += 2*(v0*pi*sin(pi*y)*x)*(v0*pi*sin(pi*y)*x);
  eS += (u0*sin(pi*x) + v0*cos(pi*y))*(u0*sin(pi*x) + v0*cos(pi*y));
  eS = (Ec/Re)*eS;
  
}


void Cell::PrintCoordinates(){

  cout<<x<<",";
  cout<<y; 

}

/* Print Required Field of Cell */
void Cell::PrintField(FieldName FN){
  
  switch(FN){
  
  case 0:
    cout<<U.T;
    break;
  case 1:
    cout<<U.u;
    break;
  case 2:
    cout<<U.v;
    break;
  default:
    cout<<"\nError: Field Unavailable";
    exit(0);
    break;
  }
  
}


//________________________________________________________________________________________//

// Construct Grid Object 
Grid::Grid(int Nx, int Ny, int nGC){

  int SIZEx;
  int SIZEy;

  dx = Lx/Nx;
  dy = Ly/Ny;

  SIZEx = Nx + 2*nGC;
  SIZEy = Ny + 2*nGC;
    
  Imin = nGC;
  Jmin = nGC;
  
  Imax = SIZEx - nGC - 1;
  Jmax = SIZEy - nGC - 1;

  Mesh = new Cell*[SIZEx];
  for(int i = 0; i < SIZEx; i++)
    Mesh[i] = new Cell[SIZEy];
  
}

// Delete Grid Object
Grid::~Grid(){

  for(int i = 0; i <= Imax + 1 ; i++)
    delete Mesh[i];

}

// Evaluate Values of Coordinates for each Cell & Store them 
void Grid::EvaluateCellCoordinates(){
 
  double x, y;
  
  /* Store values in array, that reflect i,j similar to math */
  /* Move in x first then in y, and so on */
  for(int j = 0; j <= Jmax + Jmin; j++){
    for(int i = 0; i <= Imax + Imin; i++){     
      
      x = dx*(2.0*(double)i - (double)Imin)/2.0; /* Give values to x, that change with i */
      y = dy*(2.0*(double)j - (double)Jmin)/2.0; /* Give values to y, that change with j */
      
      Mesh[i][j].SetCellCoordinates(x, y);
            
    }
  }  
  
  cout<<"\nEvaluated Cell Coordinates for Grid...";    
  
}

// Function to evaluate integrals/fields using exact solution
void Grid::EvaluateExactIntegrals(){

    for(int i = 0; i <= Imax + Imin; i++){
      for(int j = 0; j <= Jmax + Jmin; j++){ 	

	Mesh[i][j].ComputeExactField();        /* Compute Fields using given Exact Functions */
	Mesh[i][j].ComputeExactFluxIntegral(); /* Compute Flux Integrals using Exact Functions */
	Mesh[i][j].ComputeExactSourceTerm();   /* Similarly with Source Terms */
	
      }
    }
    
    cout<<"\nEvaluated Exact Field Values for Grid...";
  
}

// Evaluate Initial Values of all Fields
void Grid::EvaluateInitialFields(){

  for(int j = 0; j <= Jmax + Jmin; j++){
    for(int i = 0; i <= Imax + Imin; i++){
      
      Mesh[i][j].SetInitialField(); // Set given Initial Field for the problem (t = 0)
    }
  }

  cout<<"\nEvaluating Initial Condition on Fields...";
}

// Evaluate ghost cell values using BCs
void Grid::EvaluateBoundaryConditions(){

  double Tfd = 0.0;
  
  for(int j = Jmin; j <= Jmax; j++){
        
    Tfd = Mesh[0][j].SetFDTemperatureField();
    Mesh[0][j].U.T = 2.0*Tfd - Mesh[1][j].U.T;   
 
    Mesh[Imin + Imax][j].U.T = Mesh[Imax][j].U.T;
  }  
  
  for(int i = Imin; i <= Imax; i++){
    
    Mesh[i][0].U.T = - Mesh[i][1].U.T; 
    Mesh[i][Jmin + Jmax].U.T = 2.0 - Mesh[i][Jmax].U.T;
    
    Mesh[i][0].U.u = -Mesh[i][0].U.u;
    Mesh[i][Jmin + Jmax].U.u = -Mesh[i][Jmax].U.u;
  }
  
}

// Evaluate Flux Integrals using CV Averages
void Grid::EvaluateFluxes(){
  
  double a = 0.0, b = 0.0, d = 0.0;
  
  cout<<endl;
  for(int j = Jmin; j <= Jmax; j++){
    for(int i = Imin; i <= Imax; i++){
      
      a = -1.0/(2.0*dx)*(Mesh[i+1][j].U.u*Mesh[i+1][j].U.T - Mesh[i-1][j].U.u*Mesh[i-1][j].U.T);		     
      b = -1.0/(2.0*dy)*(Mesh[i][j+1].U.v*Mesh[i][j+1].U.T - Mesh[i][j-1].U.v*Mesh[i][j-1].U.T);
      
      d = 1.0/(Re*Pr*dx*dx)*(Mesh[i+1][j].U.T - 2.0*Mesh[i][j].U.T + Mesh[i-1][j].U.T);
      d += 1.0/(Re*Pr*dy*dy)*(Mesh[i][j+1].U.T - 2.0*Mesh[i][j].U.T + Mesh[i][j-1].U.T);
      
      Mesh[i][j].F = a + b + d;

    }
  } 
  
}

// Evaluate Source Terms using CV Averages
void Grid::EvaluateSourceTerms(){
  
  double a = 0.0, b = 0.0, c = 0.0, d = 0.0;
  
  for(int j = Jmin; j <= Jmax; j++){
    for(int i = Imin; i <= Imax; i++){
                  
      a = (Mesh[i+1][j].U.u - Mesh[i-1][j].U.u)/(2.0*dx);
      a = 2.0*a*a;
      
      b = (Mesh[i][j+1].U.v - Mesh[i][j-1].U.v)/(2.0*dy);
      b = 2.0*b*b;
      
      c = (Mesh[i+1][j].U.v - Mesh[i-1][j].U.v)/(2.0*dx);
      d = (Mesh[i][j+1].U.u - Mesh[i][j-1].U.u)/(2.0*dy);
    
      Mesh[i][j].S = (Ec/Re)*(a + b + (c + d)*(c + d));     
      
    }
  }     

  cout<<"\nEvaluating Source Terms on Grid...";

}


void Grid::InitialiseDxDyFI(){

  Dx = new double*[Imax + 2];
  for(int i = 0; i < Imax + 2; i++)
    Dx[i] = new double[3];

  Dy = new double*[Jmax + 2];
  for(int i = 0; i < Jmax + 2; i++)
    Dy[i] = new double[3];

  FI = new double*[Imax + 2];
  for(int i = 0; i < Imax + 2; i++)
    FI[i] = new double[Jmax + 2];
    
  
  for(int i = 0; i <= Imax + Imin; i++){  
    for(int j = 0; j < 3; j++){
            
      Dx[i][j] = 0.0;

    }
  }

  for(int j = 0; j <= Jmax + Jmin; j++){  
    for(int i = 0; i < 3; i++){
      
      Dy[j][i] = 0.0;

    }

  }

  for(int i = 0; i <= Imax + Imin; i++){  
    for(int j = 0; j <= Jmax + Jmin; j++){
            
      FI[i][j] = 0.0;

    }
  }

}


void Grid::EvaluateDx(int J){
   
  Dx[0][0] = 0.0;
  Dx[0][1] = 1.0;
  Dx[0][2] = 1.0;
  
  for(int i = Imin; i <= Imax; i++){
    
    Dx[i][0] = - (Mesh[i - 1][J].U.u*dt)/(2.0*dx) - dt/(Re*Pr*dx*dx);
    Dx[i][1] = 1.0 + (2.0*dt)/(Re*Pr*dx*dx);
    Dx[i][2] = (Mesh[i + 1][J].U.u*dt)/(2.0*dx) - dt/(Re*Pr*dx*dx);
    
  }  
  
  Dx[Imax + 1][0] = -1.0;
  Dx[Imax + 1][1] = 1.0;
  Dx[Imax + 1][2] = 0.0;

}

void Grid::EvaluateDy(int I){    

  Dy[0][0] = 0.0;
  Dy[0][1] = 1.0;
  Dy[0][2] = 1.0;

  for(int j = Jmin; j <= Jmax; j++){
    
    Dy[j][0] = -(Mesh[I][j - 1].U.v*dt)/(2.0*dy) - dt/(Re*Pr*dy*dy);
    Dy[j][1] = 1.0 + (2.0*dt)/(Re*Pr*dy*dy);
    Dy[j][2] = (Mesh[I][j + 1].U.v*dt)/(2.0*dy) - dt/(Re*Pr*dy*dy);
    
  }  
  
  Dy[Jmax + 1][0] = 1.0;
  Dy[Jmax + 1][1] = 1.0;
  Dy[Jmax + 1][2] = 0.0;
    
}

void Grid::EvaluateFI(int J){
  
  for(int i = Imin; i <= Imax; i++){
    
    FI[i][J] = dt*(Mesh[i][J].F + Mesh[i][J].S);

  }

  FI[0][J] = FI[Imax + Imin][J] = 0.0;

}


double Grid::EvaluateTimeStep(TimeScheme TS){

  double Umax;

  switch(TS){

  case 0: 
    Umax = 9.0/2.0;
    dt = 0.6*dx/Umax;   
    dt = 0.01;
    cout<<"\nCalculating time step for Explicit Euler...";
    break;
  case 1:
    dt = 0.01;
    cout<<"\nCalculating time step for Implicit Euler...";
    break;
  default:
    break;
    
  }

  return dt;
}

// void Grid::EulerExplicitTimeAdvance(){
  
//   EvaluateBoundaryConditions();	                           // Implement ghost cell values for flux calculation  
//   EvaluateFluxes();	         	                   // Evaluate Fluxes from solution at time step n. Change to Limited flux as required
//   //  EvaluateSourceTerms();					   // Source terms could change with time for NSE

//   for(int i = Imin; i <= Imax; i++){
//     for(int j = Jmin; j <= Jmax; j++){
      
//       Mesh[i][j].U.T = Mesh[i][j].U.T
// 	+ dt*Mesh[i][j].F 
// 	+ dt*Mesh[i][j].S;	                                   // Calculate intermediate step solution, store in separate array
//     }
//   }
  
// }

Field Grid::EulerExplicitTimeAdvance(){

  Field dU, dUmax;

  EvaluateBoundaryConditions();	                                   // Implement ghost cell values for flux calculation  
  EvaluateFluxes();	         	                           // Evaluate Fluxes from solution at time step n. Change to Limited flux as required

  for(int i = Imin; i <= Imax; i++){
    for(int j = Jmin; j <= Jmax; j++){      
      
      dU.T = dt*(Mesh[i][j].S + Mesh[i][j].F);
      Mesh[i][j].U.T += dU.T;
      
      if(dU.T >= dUmax.T)
	dUmax = dU;						   // Find maximum change in solution
      
    }
  }

  return dUmax;
}

//Compute solution using Implicit scheme
Field Grid::EulerImplicitTimeAdvance(){

  Field dU, dUmax;
  double RHS[NMAX], LHS[NMAX][3];

  EvaluateBoundaryConditions();
  EvaluateFluxes();
  
  //PrintFluxes();
  
  InitialiseDxDyFI();
  
    for(int j = Jmin; j <= Jmax; j++){      
            
      EvaluateDx(j);				          	   // Evaluate tridiagonal on x
      EvaluateFI(j);				                   // evaluate the RHS which is (FI + S)*dt
      cout<<endl;
      
      CopyToLHS(Dx, LHS, Imax + 1);
                  
      CopyToRHS(FI, RHS, Imax + 1, j);
      for(int i = 0; i <= Imax + Imin; i++){
	
	cout<<LHS[i][0]<<" "<<LHS[i][1]<<" "<<LHS[i][2]<<endl;
	//cout<<RHS[i]<<endl;
	  }      

      /*
       * Solve Thomas seems to be giving nan
       * Copying to RHS and LHS works fine
       * Must add functionality to copy to RHS & copy from RHS
       * add copy from RHS after every solve thomas routine
       */
      SolveThomas(LHS, RHS, Imax);      
      
      


      //CopyFromRHS(FI, RHS, Imax + 1, j);
    }

    for(int i = Imin; i <= Imax; i++){

      EvaluateDy(i);  
      CopyToLHS(Dy, LHS, Jmax + 1);
      //CopyToRHS(FI, RHS, Jmax + 1, i);
      SolveThomas(LHS, RHS, Jmax);  
      CopyFromRHS(FI, RHS, Jmax + 1, i);
      
    }
    
    
    for(int i = Imin; i <= Imax; i++){
      for(int j = Jmin; j <= Jmax; j++){      
	
      dU.T = FI[i][j];
      Mesh[i][j].U.T += dU.T;
      
      if(dU.T >= dUmax.T)
	dUmax = dU;						   // Find maximum change in solution
      
      }
    }
    
    return dU;    
}

//Verify fluxes when exact values are provided
void Grid::FluxVerification(){
  
  double ErrorF = 0.0, L2NormF = 0.0;
  ofstream fileF;

  for(int i = Imin; i <= Imax; i++){     
    for(int j = Jmin; j <= Jmax; j++){
     
      ErrorF = Mesh[i][j].F - Mesh[i][j].eF;                                 
      L2NormF += ErrorF*ErrorF;	    

    }
  }  
  
  L2NormF = sqrt(L2NormF/(Imax*Jmax));				      
  
  fileF.open("Ef", ios::app);	                                      // Plot in a file with name F   
  fileF<<setprecision(15)<<dx<<" "<<L2NormF<<endl;		      // Plot at each point in domain 
  
  cout<<"Flux Errors Calculated and Written to File...\n";
  cout<<"Error Norm (l2): "<<L2NormF<<"\n";
  
  fileF.close();						      // Close file

}// End of File Write


// 
void Grid::SourceVerification(){
  
  double ErrorS = 0.0, L2NormS = 0.0;
  ofstream fileS;

  for(int i = Imin; i <= Imax; i++){     
    for(int j = Jmin; j <= Jmax; j++){
     
      ErrorS = Mesh[i][j].S - Mesh[i][j].eS;
      L2NormS += ErrorS*ErrorS;

    }
  }  
  
  L2NormS = sqrt(L2NormS/(Imax*Jmax));				      
  
  fileS.open("Es", ios::app);
  fileS<<setprecision(15)<<dx<<" "<<L2NormS<<endl;
  
  cout<<"Source Errors Calculated and Written to File...\n";
  cout<<"Error Norm (l2): "<<L2NormS<<"\n";      
  
  fileS.close();

}// End of File Write


void Grid::FieldVerification(){

  double ErrorT = 0.0, L2NormT = 0.0;
  ofstream fileT;

  for(int i = Imin; i <= Imax; i++){     
    for(int j = Jmin; j <= Jmax; j++){
      
      ErrorT = Mesh[i][j].U.T - Mesh[i][j].eU.T;
      L2NormT += ErrorT*ErrorT;
      
    }
  }  
  
  L2NormT = sqrt(L2NormT/(Imax*Jmax));				      
  
  fileT.open("Et", ios::app);
  fileT<<setprecision(15)<<dx<<" "<<L2NormT<<endl;
  
  cout<<"Temperature Field Errors Calculated and Written to File...\n";
  cout<<"Error Norm (l2): "<<L2NormT<<"\n";      
  
  fileT.close();

}

//Print Cell coordinates for entire mesh
void Grid::PrintCellCoordinates(){

  int i, j;
  
  cout<<endl;
	
  
  for(j = Jmax; j >= Jmin; j--){
    for(i = Imin; i <= Imax; i++){           
      
      Mesh[i][j].PrintCoordinates();
      cout<<" ";
      
    }
    cout<<endl;
  }
  
}

// Printing field values over entire mesh
void Grid::PrintFieldValues(FieldName FN){      
  
  cout<<"\nPrinting Field Values...\n";

  for(int j = Jmax + Jmin; j >= 0; j--){
    for(int i = 0; i <= Imin + Imax; i++){     
      
      Mesh[i][j].PrintField(FN);      
      cout<<" ";

    } 
    cout<<endl;							       
  }     
  
}

// Printing field values over specific regions of mesh
void Grid::PrintFieldValues(FieldName FN, int Row, int Column){

  cout<<"\nPrinting Field Values...\n";
  if(Row == -1){
    for(int i = Imax; i >= Imax; i--)
      {
	Mesh[i][Column].PrintField(FN);
	cout<<endl;
      }
  }
  else if(Column == -1){
    for(int j = Jmax; j >= Jmin; j--)
      {
	Mesh[Row][j].PrintField(FN); 
	cout<<endl;
      }
  }
  else
    Mesh[Row][Column].PrintField(FN);
  
}

void Grid::PrintFluxes(){
  
  int i, j;
  
  cout<<"\nPrinting Numerical Fluxes...\n";
  
  for(j = Jmax; j >= Jmin; j--){
    for(i = Imin; i <=  Imax; i++){     
      
      cout<<Mesh[i][j].F<<" ";
      
    } 
    cout<<endl;							       
  }     

  // cout<<"\nPrinting Exact Fluxes...\n";
  
  // for(j = Jmax; j >= Jmin; j--){
  //   for(i = Imin; i <= Imax; i++){     
      
  //     cout<<Mesh[i][j].eF;
      
  //   } 
  //   cout<<endl;			
  // }				       
  
}

void Grid::PrintSources(){
  
  int i, j;
  
  cout<<"\nPrinting Numerical Sources...\n";
  
  for(j = Jmax + Jmin; j >= 0; j--){
    for(i = 0; i <= Imax + Imin; i++){     
      
      cout<<Mesh[i][j].S;
      
    } 
    cout<<endl;							       
  }     
  
  // cout<<"\nPrinting Exact Sources...\n";  

  // for(j = Jmax; j >= Jmin; j--){
  //   for(i = Imin; i <= Imax; i++){     
      
  //     cout<<Mesh[i][j].eS;
      
  //   } 
  //   cout<<endl;			
  // }				       
  
}
