//////////////////////////////////////////////////////////
// The following code forms a 2D Energy Equation Solver //
// i.e. Solution of				        //
// dT/dt + udT/dx + vdT/dy = 1/Re*(d2T/dx2 + d2T/dy2)   //
// 						        //
// With given boundary conditions and initial condition //
//////////////////////////////////////////////////////////

#include"header.h"

int main(){  
  clock_t ti, tf;
  ti = clock();
  
  Grid Domain(10,10, 1);	       // Define Grid of required size
  
  //SolveGoverningEquation(&Domain);
  Domain.PrintCVData();
  
  tf = clock();

  cout<<"\nSolver Run-Time: "<<(double)(tf - ti)/CLOCKS_PER_SEC<<" seconds"<<endl;
  return 0;
}


//________________________________________________________________________________________//


Grid::Grid(int Nx, int Ny, int nGC){
    
  SIZEx = Nx;
  SIZEy = Ny;
  
  Imin = nGC;
  Imax = SIZEx - nGC - 1;					   // 
  
  Jmin = nGC;
  Jmax = SIZEy - nGC - 1;					   // 
  
  GhostCells = nGC;

  Mesh = (Cell**)malloc(Nx*sizeof(int));
  
  for(int i = 0; i < Nx; i++)
    Mesh[i] = (Cell*)malloc(Ny*sizeof(Cell));

  for(int i = 0; i < Nx; i++)
    for(int j = 0; j < Ny; j++)
      Mesh[i][j].SetFieldValues(2.*i/j, 0.0, 0.0);
  
}

Grid::~Grid(){

  for(int i = 0; i <SIZEx ; i++)
    delete Mesh[i];

}

void Grid::PrintCVData(){
    
  int i,j;
  
  for(i = Imin; i <= Imax; i++){     
    for(j = Jmin; j <= Jmax; j++){
      Mesh[i][j].PrintFieldValue();
    } 
    cout<<endl;							       
  }   
  
  cout<<endl;  

  
}


void Cell::SetFieldValues(double Temp, double Ux, double Uy){

  T = Temp;
  u = Ux;
  v = Uy;
  
}

void Cell::PrintFieldValue(){
  
  cout<<setw(15)<<setprecision(10)<<T;

}

// void ExperimentArray(Cell **Mesh){

//   for(int i = 0; i<10;i++){
//     for(int j = 0; j < 10; j++){
//       Mesh[i][j].setValue(i*j);
//     }
//   }
  
//   for(int i = 0; i<10;i++){
//     for(int j = 0; j < 10; j++)
//       cout<<Mesh[i][j].showValue()<<" ";
//     cout<<endl;
//   }
//   cout<<endl;

// }



void Grid::ComputeExactFluxIntegral(){

  double dx, dy, x, y, u, v, T, FI;
  double u0, v0, T0;
  
  u0 = 1.0;
  v0 = 1.0;
  T0 = 1.0;
  
  for(int i = Imin; i<= Imax;i++){
    for(int j = Imin; j <= Jmax; j++)
       {

	 x = dx*(2.0*(double)i -1)/2;
	 y = dy*(2.0*(double)j - 1)/2;

	 T = T0*cos(pi*x)*sin(pi*y);
	 u = u0*y*sin(pi*x);
	 v = v0*x*cos(pi*y);
	 
	 FI = u0*T0*pi*cos(2*pi*x)*y*sin(pi*y) 
	   + v0*T0*pi*x*cos(pi*x)*cos(2*pi*x) 
	   + 2*T0*pi*pi*cos(pi*x)*sin(pi*y)/(Re*Pr);	 	 
	 
	 Mesh[i][j].SetFluxIntegral(FI);	 	 
       }
  }

}
