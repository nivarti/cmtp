//////////////////////////////////////////////////////////
// The following code forms a 2D Energy Equation Solver //
// i.e. Solution of				        //
// dT/dt + udT/dx + vdT/dy = 1/Re*(d2T/dx2 + d2T/dy2)   //
// 						        //
// With given boundary conditions and initial condition //
//////////////////////////////////////////////////////////


#include<iostream>
#include<stdio.h>
#include<stdlib.h>

using namespace std;

#define pi 3.14159265358979323846264338327950288419716939937510    // Value of pi

const double CFL = 0.4;					           // CFL number 
const double Lx  = 5.0;						   // Domain size x
const double Ly = 1.0;						   // Domain size y

const double Re = 50;						   // Reynolds Number
const double Pr = 0.7;						   // Prandtl Number
const double Ec = 0.1;						   // Eckert Number

//______________________________________________________________________________________//


class Cell{
  
  double T;							   // Control Volume Average Temperature
  double u;							   // Control Volume Average Velocity (x)
  double v;							   // Control Volume Average Velocity (y)
  
  double dx;							   // Cell Size (x direction)
  double dy;							   // Cell Size (y direction)
  int I;
  int J;
     
public:
  
  Cell(){
    T = 0.0;
    
  }
  
  Cell(double Temp){
    T = Temp;
  }


  void setFields(double Temp, double Velx, double Vely){
    T = Temp;
    u = Velx;
    v = Vely;
  }
  
  double showFields(){
    return T;}

  // Computations of Source can be done within CV
  // Printing of fields can be done within CV
  // Flux Integrals have to be computed outside CV

};// End of Class Definition


Cell** buildArray(int width, int height);
void computeExactFluxes(Cell** Mesh);


int main(){
  
  
  int x = 10;
  Cell **Mesh = buildArray(x,x);
  
  

  free(Mesh); 
  return 0;
}


//________________________________________________________________________________________//


Cell** buildArray(int width, int height){
  
  Cell** array = (Cell**)malloc(width*sizeof(int));
  
  for(int i = 0; i < width; i++)
    array[i] = (Cell*)malloc(height*sizeof(Cell));
  
  return array;
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


// void initializeCell(Cell **Mesh){
  
//   int i,j;
  
//   for(i = 0; i < SIZE; i++){     
//     for(j = 0; j < SIZE; j++){
//       Mesh[i][j].setvalue(0.0);		                                      // Set 0.0 as initial guess for Temperature
//     }
//   }   
//   cout<<"Field Initialized...\n";			      
  
// } // End of Initialization of Fields
