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

#include"constants.h"

using namespace std;


//______________________________________________________________________________________//


class Cell{
  
  double T;							   // Control Volume Average Temperature
  double u;							   // Control Volume Average Velocity (x)
  double v;							   // Control Volume Average Velocity (y)

  double dx;							   // Cell Size (x direction)
  double dy;							   // Cell Size (y direction)
          
  double S;			                                   // Source within Cell
  double FI;							   // Flux Integral for given control volume

public:
  
  void SetFieldValues(double, double, double);
  void SetFluxIntegral(double){}
  void PrintFieldValue();
  
  // Computations of Source can be done within CV
  
};// End of Class Definition

class Grid{
  
  Cell** Mesh;
  
  int SIZEx;
  int SIZEy;

  int Imin;
  int Imax;

  int Jmin;
  int Jmax;

  int GhostCells;

  // Printing of fields can be done within CV
  // Flux Integrals have to be computed outside CV
public:
  
  Grid();
  Grid(int Nx, int Ny, int nGC);
  ~Grid();			
  
  void PrintCVData();
  void WriteFieldValues();
  void EvaluateFluxIntegrals();
  void ComputeExactFluxIntegral();  
  

};
