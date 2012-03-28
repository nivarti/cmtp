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

enum FieldName{Temperature, xVelocity, yVelocity};
enum TimeScheme{ExplicitEuler, ImplicitEuler, ExplicitRK2, ExplicitRK4};

//__________________________________Field Definition____________________________________//

struct Field{
  
  double T;							   // Control Volume Average Temperature
  double u;							   // Control Volume Average Velocity (x)
  double v;							   // Control Volume Average Velocity (y)

  Field();
  Field& operator=(const Field &RHS);
  Field& operator+=(const Field &RHS);
  //friend ostream& operator<<(ostream &out, const Field U); 

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

  void SetInitialField();
  double SetFDTemperatureField();

  void ComputeExactField();					   // Computing Exact Values is easy within CV
  void ComputeExactFluxIntegral();  
  void ComputeExactSourceTerm();
  
  void PrintCoordinates();
  void PrintField(FieldName);
    
};// End of Class Definition

//_________________________________Grid Definition________________________________________//

class Grid{
  
  Cell** Mesh;							   // 
  
  int Imin, Imax;						   // 
  int Jmin, Jmax;    

  double dx;		 					   // 
  double dy;							   // In essence cell size should stay with cell

  double dt;

public:
  
  Grid();
  Grid(int, int, int);
  ~Grid();			
  
  void PrintCellCoordinates();
  void PrintFluxes();
  void PrintSources();
  
  void PrintFieldValues(FieldName);
  void PrintFieldValues(FieldName, int, int);
  void WriteFieldValues(FieldName);

  void EvaluateExactIntegrals();

  void EvaluateCellCoordinates();
  void EvaluateInitialFields();   
  void EvaluateFluxIntegrals();
  void EvaluateSourceTerms();
  void EvaluateBoundaryConditions();
  
  double EvaluateTimeStep(TimeScheme);
  
  //void EulerExplicitTimeAdvance();
  Field EulerExplicitTimeAdvance();
  void EulerImplicitTimeAdvance();  
  void RK2ExplicitTimeAdvance();

  void FluxVerification();
  void SourceVerification();
  
};

//_________________________________Function Declarations__________________________________//

void EvaluateGridParameters(Grid&);
void MarchTime(Grid&, TimeScheme);
void SolveEnergyEquation(Grid&, double);
