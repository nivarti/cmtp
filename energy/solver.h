/*********************************************/
/* This file contains header file includes,  */
/* and class definitions for relevant data   */
/* structure for a 2D energy equation solver */
/*********************************************/

#include<iostream>		
#include<iomanip>

#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<time.h>

#include"constant.h"

using namespace std;


//______________________________________________________________________________________//

struct Field{
  
  double T;							   // Control Volume Average Temperature
  double u;							   // Control Volume Average Velocity (x)
  double v;							   // Control Volume Average Velocity (y)
public:
  Field();

};							

struct Cell{
  
  Field U;							   // Field contains solution variables T, u, v
  double x;
  double y;
    
  double S;			                                   // Source within Cell
  double FI;							   // Flux Integral for given control volume

public:
  // Computations of Source can be done within CV
  Cell();
  
  void SetCellCoordinates(double, double);
  void ComputeExactField();					   // Computing Exact Values is easy within CV
  void ComputeExactFluxIntegral();  
  void ComputeExactSourceTerm();
  
  void PrintField(int f);
    
};// End of Class Definition


class Grid{
  
  Cell** Mesh;							   // 
  
  int Imin, Imax;						   // 
  int Jmin, Jmax;    

  double dx;		 					   // 
  double dy;							   // In essence cell size should stay with cell
  

public:
  
  Grid();
  Grid(int Nx, int Ny, int nGC);
  ~Grid();			
  
  void PrintFieldValues(int f);
  void WriteFieldValues(int f);
  
  void EvaluateCellData();
  void EvaluateFluxIntegrals();
  void EvaluateSourceTerms();
  
  void FluxVerification();
  void SourceVerification();
  
};


//________________________________________________________________________________________//

Field::Field(){
  
  T = 0.0;
  u = 0.0;
  v = 0.0;

}

Cell::Cell(){

  x = 0.0;
  y = 0.0;
  
}


void Cell::SetCellCoordinates(double X, double Y){
  
  x = X;
  y = Y;
  
}

// Compute Exact Values of Fields from Formulae
void Cell::ComputeExactField(){
  
  U.T = T0*cos(pi*x)*sin(pi*y);
  U.u = u0*y*sin(pi*x);
  U.v = v0*x*cos(pi*y);
  
}


void Cell::ComputeExactSourceTerm(){
  
  S = (Ec/Re)*(2*(u0*pi*cos(pi*x)*y)*(u0*pi*cos(pi*x)*y)
	       + 2*(v0*pi*sin(pi*y)*x)*(v0*pi*sin(pi*y)*x)
	       + (u0*sin(pi*x) + v0*cos(pi*y))*(u0*sin(pi*x) + v0*cos(pi*y)));
  
}

  
void Cell::ComputeExactFluxIntegral(){

  FI = u0*T0*pi*cos(2*pi*x)*y*sin(pi*y) 
    + v0*T0*pi*x*cos(pi*x)*cos(2*pi*x) 
    + 1/(Re*Pr)*2*T0*pi*pi*cos(pi*x)*sin(pi*y);	 	 
  
}

// Output Required Field on Screen
void Cell::PrintField(int f){
  
  switch(f){
  
  case 1:
    cout<<setw(10)<<setprecision(5)<<U.T;
    break;
  case 2:
    cout<<setw(15)<<setprecision(10)<<U.u;
    break;
  case 3:
    cout<<setw(15)<<setprecision(10)<<U.v;
    break;
  default:
    cout<<"\nError: Field Unavailable";
    exit(0);
  }
  
}


//________________________________________________________________________________________//

void Grid::EvaluateFluxIntegrals(){

  double FI;
  
  for(int i = Imin; i<= Imax;i++){
    for(int j = Imin; j <= Jmax; j++)
      {	 
	
	//FI =  -1/dx*( Mesh[i][j].U.u)
	
      }
  }
  
}

//void Grid::FluxVerification(){

//   double dx, dy, x, y, u, v, T, FI;
//   double u0, v0, T0;
  
//   u0 = 1.0;
//   v0 = 1.0;
//   T0 = 1.0;
  
//   for(int i = Imin; i<= Imax;i++){
//     for(int j = Imin; j <= Jmax; j++)
//       {	 
	
//       }
//   }
  
// }



// void Grid::SourceVerification(){

//   double dx, dy, x, y, u, v, T, S;
//   double u0, v0, T0;
  
//   u0 = 1.0;
//   v0 = 1.0;
//   T0 = 1.0;
  
//   for(int i = Imin; i<= Imax;i++){
//     for(int j = Imin; j <= Jmax; j++)
//        {

	 
// 	 Mesh[i][j].S;
//        }
//   }

  
// }

void Grid::EvaluateCellData(){
 
  double x, y;

  for(int i = 0; i <= Imax + Imin;i++){
    for(int j = 0; j <= Jmax + Imin; j++)
      { 	
	
	x = dx*(2.0*(double)i - (double)Imin)/2.0;
	y = dy*(2.0*(double)j - (double)Jmin)/2.0;	
	
	Mesh[i][j].SetCellCoordinates(x, y);
  
	Mesh[i][j].ComputeExactField();
	Mesh[i][j].ComputeExactFluxIntegral();
	Mesh[i][j].ComputeExactSourceTerm();
      }
  }

  cout<<"\nEvaluated Exact Field Values for Grid";
}

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

Grid::~Grid(){

  for(int i = 0; i <= Imax + 1 ; i++)
    delete Mesh[i];

}

 void Grid::PrintFieldValues(int f){
    
  int i,j;   
  
  cout<<"\nImin: "<<Imin<<", Imax: "<<Imax;
  cout<<"\nJmin: "<<Jmin<<", Jmax: "<<Jmax;

  cout<<endl;

  for(i = Imin - 1; i <= Imax + 1; i++){     
    for(j = Jmin - 1; j <= Jmax + 1; j++){
      
      Mesh[i][j].PrintField(f);

    } 
    cout<<endl;							       
  }     
  cout<<endl;  
  
}
