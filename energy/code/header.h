/*********************************************/
/* This file contains header file includes,  */
/* and class definitions for relevant data   */
/* structure for a 2D energy equation solver */
/*********************************************/

#include<iostream>
#include<fstream>		
#include<iomanip>

#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<time.h>

#include "constant.h"

using namespace std;

//__________________________________Field Definition____________________________________//

struct Field{
  
  double T;							   // Control Volume Average Temperature
  double u;							   // Control Volume Average Velocity (x)
  double v;							   // Control Volume Average Velocity (y)

  Field();
  Field& operator=(const Field &RHS);

};							

//__________________________________Cell Definition_____________________________________//


struct Cell{
  
  Field U;							   // Field contains solution variables T, u, v
  Field Ubuf;		                                       	/* Stores intermediate time step information */

  double x;
  double y;
    
  double eFI;			                                /*  */
  double eS;

  double S;			                                   // Source within Cell
  double FI;							   // Flux Integral for given control volume

public:
  // Computations of Source can be done within CV
  Cell();
  
  void SetCellCoordinates(double, double);
  void SetCellField(double, double, double);
  
  void ComputeExactField();					   // Computing Exact Values is easy within CV
  void ComputeExactFluxIntegral();  
  void ComputeExactSourceTerm();
  
  void PrintCoordinates();
  void PrintField(int f);
    
};// End of Class Definition

//_________________________________Grid Definition________________________________________//

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
  void PrintCellCoordinates();
  void PrintFluxes();
  void PrintSources();
  void WriteFieldValues(int f);
  
  void EvaluateCellCoordinates();
  void EvaluateCellData();

  void EvaluateExactIntegrals();
  void EvaluateFluxIntegrals();
  void EvaluateSourceTerms();
  
  void EvaluateL2Norm();
  void FluxVerification();
  void SourceVerification();
  
};

//_________________________________Function Declarations__________________________________//

void EvaluateGridParameters(Grid &Domain);
