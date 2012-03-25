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


Cell::Cell(){
  
  x = 0.0;
  y = 0.0;

  eFI = 0.0;
  eS = 0.0;

  S = 0.0;
  FI = 0.0;
  
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
  /* U.u = 6.0*Uav*y*(1 - y); */
  /* U.v = 0.0; */
}
  
void Cell::ComputeExactFluxIntegral(){

  eFI = u0*T0*pi*cos(2*pi*x)*y*sin(pi*y) 
    + v0*T0*pi*x*cos(pi*x)*cos(2*pi*y) 
    + 1/(Re*Pr)*2*T0*pi*pi*cos(pi*x)*sin(pi*y);	 	 
  eFI = -eFI;			/* flip sign while taking integral to RHS */
  
}


void Cell::ComputeExactSourceTerm(){
      
  eS  = 2*(u0*pi*cos(pi*x)*y)*(u0*pi*cos(pi*x)*y);
  eS += 2*(v0*pi*sin(pi*y)*x)*(v0*pi*sin(pi*y)*x);
  eS += (u0*sin(pi*x) + v0*cos(pi*y))*(u0*sin(pi*x) + v0*cos(pi*y));
  eS = (Ec/Re)*eS;
  
}


void Cell::PrintCoordinates(){

  cout<<setw(5)<<setprecision(4)<<x<<",";
  cout<<setw(5)<<setprecision(4)<<y; 

}

/* Print Required Field of Cell */
void Cell::PrintField(int f){
  
  switch(f){
  
  case 1:
    cout<<setw(10)<<setprecision(5)<<U.T;
    break;
  case 2:
    cout<<setw(10)<<setprecision(5)<<U.u;
    break;
  case 3:
    cout<<setw(10)<<setprecision(5)<<U.v;
    break;
  default:
    cout<<"\nError: Field Unavailable";
    exit(0);
    break;
  }
  
}


//________________________________________________________________________________________//

/* Construct Grid Object */
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

  /* Mesh = (Cell**)malloc(SIZEx*sizeof(int)); */
  
  /* for(int i = 0; i < SIZEx; i++) */
  /*   Mesh[i] = (Cell*)malloc(SIZEy*sizeof(Cell)); */
  
}

/* Delete Grid Object */
Grid::~Grid(){

  for(int i = 0; i <= Imax + 1 ; i++)
    delete Mesh[i];

}

/* Evaluate Values of Coordinates for each Cell & Store them */
void Grid::EvaluateCellCoordinates(){
 
  double x, y;
  
  /* Store values in array, that reflect i,j similar to math */
  /* Move in x first then in y, and so on */
  for(int j = 0; j <= Jmax + Imin; j++){
    for(int i = 0; i <= Imax + Imin; i++){     
      
      x = dx*(2.0*(double)i - (double)Imin)/2.0; /* Give values to x, that change with i */
      y = dy*(2.0*(double)j - (double)Jmin)/2.0; /* Give values to y, that change with j */
      
      Mesh[i][j].SetCellCoordinates(x, y);
      
      /* //cout<<"i = "<<i<<", j = "<<j<<": "; */
      /* Mesh[i][j].PrintCoordinates(); */
      /* cout<<" "; */
      
    }
  }  
  
  cout<<"\nEvaluated Cell Coordinates for Grid...";    
  
}


/* Function to evaluate integrals/fields using exact solution */
void Grid::EvaluateExactIntegrals(){

    for(int i = 0; i <= Imax + Imin; i++){
      for(int j = 0; j <= Jmax + Imin; j++){ 	

	Mesh[i][j].ComputeExactField(); /* Compute Fields using given Exact Functions */
	Mesh[i][j].ComputeExactFluxIntegral(); /* Compute Flux Integrals using Exact Functions */
	Mesh[i][j].ComputeExactSourceTerm();   /* Similarly with Source Terms */
	
      }
    }
    
    cout<<"\nEvaluated Exact Field Values for Grid...";
  
}


void Grid::EvaluateFluxIntegrals(){
  
  double a = 0.0, b = 0.0, d = 0.0;
  
  cout<<endl;
  for(int j = Jmin; j <= Jmax; j++){
    for(int i = Imin; i<= Imax;i++){
      
      a = -1/(2.0*dx)*(Mesh[i+1][j].U.u*Mesh[i+1][j].U.T - Mesh[i-1][j].U.u*Mesh[i-1][j].U.T);		     
      b = -1/(2.0*dy)*(Mesh[i][j+1].U.v*Mesh[i][j+1].U.T - Mesh[i][j-1].U.v*Mesh[i][j-1].U.T);
      
      d = 1/(Re*Pr*dx*dx)*(Mesh[i+1][j].U.T - 2*Mesh[i][j].U.T + Mesh[i-1][j].U.T);
      d += 1/(Re*Pr*dy*dy)*(Mesh[i][j+1].U.T - 2*Mesh[i][j].U.T + Mesh[i][j-1].U.T);
      
      Mesh[i][j].FI = a + b + d;
      /* cout<<i-1<<","<<j<<" : "; */
      /* Mesh[i-1][j].PrintCoordinates(); */
      /* cout<<" : T = "; */
      /* Mesh[i-1][j].PrintField(1); */
      /* cout<<", u = "; */
      /* Mesh[i-1][j].PrintField(2); */
      /* cout<<", v = "; */
      /* Mesh[i-1][j].PrintField(3); */
      /* cout<<", FI = "<<Mesh[i-1][j].FI;       */
      /* cout<<endl; */
      
      /* cout<<i<<","<<j<<" : "; */
      /* Mesh[i][j].PrintCoordinates();       */
      /* cout<<" : T = "; */
      /* Mesh[i][j].PrintField(1); */
      /* cout<<", u = "; */
      /* Mesh[i][j].PrintField(2); */
      /* cout<<", v = "; */
      /* Mesh[i][j].PrintField(3); */
      /* cout<<", FI = "<<Mesh[i][j].FI;       */
      /* cout<<endl; */

      /* cout<<i+1<<","<<j<<" : "; */
      /* Mesh[i+1][j].PrintCoordinates();       */
      /* cout<<" : T = "; */
      /* Mesh[i+1][j].PrintField(1); */
      /* cout<<", u = "; */
      /* Mesh[i+1][j].PrintField(2); */
      /* cout<<", v = "; */
      /* Mesh[i+1][j].PrintField(3); */
      /* cout<<", FI = "<<Mesh[i+1][j].FI;       */
      /* cout<<endl<<endl; */
      
    }
  } 
  
}

void Grid::EvaluateSourceTerms(){
  
  double a = 0.0, b = 0.0, c = 0.0, d = 0.0;
  
  for(int j = Jmin; j <= Jmax; j++){
    for(int i = Imin; i<= Imax;i++){
            
      
      a = (Mesh[i+1][j].U.u - Mesh[i-1][j].U.u)/(2.0*dx);
      a = 2.0*a*a;
      
      b = (Mesh[i][j+1].U.v - Mesh[i][j-1].U.v)/(2.0*dy);
      b = 2.0*b*b;

      

      c = (Mesh[i+1][j].U.v - Mesh[i-1][j].U.v)/(2.0*dx);
      d = (Mesh[i][j+1].U.u - Mesh[i][j-1].U.u)/(2.0*dy);
    
      Mesh[i][j].S = (Ec/Re)*(a + b + (c + d)*(c + d));     
      
    }
  }   
  
}

// void Grid::EvaluateBoundaryConditions(){


// }

void Grid::EvaluateL2Norm(){
  
  double ErrorS = 0.0, ErrorF = 0.0, L2NormF = 0.0, L2NormS = 0.0;
  ofstream fileF;
  ofstream fileS;
  int i, j;

  for(i = Imin; i <= Imax; i++){     
    for(j = Jmin; j <= Jmax; j++){
     
      ErrorF = Mesh[i][j].FI - Mesh[i][j].eFI;                                 
      ErrorS = Mesh[i][j].S - Mesh[i][j].eS;
      L2NormF += ErrorF*ErrorF;	    
      L2NormS += ErrorS*ErrorS;

    }
  }  
  
  L2NormF = sqrt(L2NormF/(Lx*Ly));				      
  L2NormS = sqrt(L2NormS/(Lx*Ly));				      
  
  fileF.open("Ef", ios::app);	                                      // Plot in a file with name F   
  fileS.open("Es", ios::app);
  fileF<<setprecision(15)<<dx<<" "<<L2NormF<<endl;		      // Plot at each point in domain 
  fileS<<setprecision(15)<<dx<<" "<<L2NormS<<endl;
  
  cout<<"Errors Calculated and Written to File...\n";
  cout<<"Flux Error Norm (l2): "<<L2NormF<<"\n";
  cout<<"Source Error Norm (l2): "<<L2NormS<<"\n";      
  
  fileF.close();						      // Close file
  fileS.close();
}// End of File Write


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


void Grid::PrintFieldValues(int f){
    
  int i,j;   
  
  /* cout<<"\nImin: "<<Imin<<", Imax: "<<Imax; */
  /* cout<<"\nJmin: "<<Jmin<<", Jmax: "<<Jmax; */

  cout<<endl;

  for(j = Jmax; j >= Jmin; j--){
    for(i = Imin; i <= Imax; i++){     
      
      Mesh[i][j].PrintField(f);
      //cout<<setw(10)<<setprecision(5)<<Mesh[i][j].FI;

    } 
    cout<<endl;							       
  }     
  
  /* cout<<endl;   */

  /* for(i = Imin - 1; i <= Imax + 1; i++){      */
  /*   for(j = Jmin - 1; j <= Jmax + 1; j++){ */
      
  /*     cout<<setw(10)<<setprecision(5)<<Mesh[i][j].eFI; */


  /*   }  */
  /*   cout<<endl;							        */
  /* }      */
  /* cout<<endl;   */

  
}

void Grid::PrintFluxes(){
  
  int i, j;
  
  cout<<"\nPrinting Numerical Fluxes...\n";
  

  for(j = Jmax; j >= Jmin; j--){
    for(i = Imin; i <= Imax; i++){     
      
      //Mesh[i][j].PrintField(f);
      cout<<setw(10)<<setprecision(5)<<Mesh[i][j].FI;
      
    } 
    cout<<endl;							       
  }     


  cout<<"\nPrinting Exact Fluxes...\n";
  

  for(j = Jmax; j >= Jmin; j--){
    for(i = Imin; i <= Imax; i++){     
      
      
      //Mesh[i][j].PrintField(f);
      cout<<setw(10)<<setprecision(5)<<Mesh[i][j].eFI;
      
    } 
    cout<<endl;			
  }				       
  
}

void Grid::PrintSources(){
  
  int i, j;
  
  cout<<"\nPrinting Numerical Sources...\n";
  

  for(j = Jmax; j >= Jmin; j--){
    for(i = Imin; i <= Imax; i++){     
      
      //Mesh[i][j].PrintField(f);
      cout<<setw(10)<<setprecision(5)<<Mesh[i][j].S;
      
    } 
    cout<<endl;							       
  }     
  
  cout<<"\nPrinting Exact Sources...\n";
  

  for(j = Jmax; j >= Jmin; j--){
    for(i = Imin; i <= Imax; i++){     
      
      
      //Mesh[i][j].PrintField(f);
      cout<<setw(10)<<setprecision(5)<<Mesh[i][j].eS;
      
    } 
    cout<<endl;			
  }				       
  
}
