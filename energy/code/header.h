/*********************************************/
/* This file contains header file includes,  */
/* and class definitions for relevant data   */
/* structure for a 2D energy equation solver */
/*********************************************/

/* Include important header files */
#include<iostream>
#include<fstream>		
#include<iomanip>

#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<time.h>

/* Include header file containing relevant constants */
#include "constant.h"

using namespace std;

/* Create enumerations to make life easy */
enum FieldName{Temperature, xVelocity, yVelocity};
enum TimeScheme{ExplicitEuler, ImplicitEuler, ExplicitRK2, ExplicitRK4};
enum Direction{Row, Column};
//__________________________________Field Definition____________________________________//

struct Field{
  
  double T;							   // Control Volume Average Temperature
  double u;							   // Control Volume Average Velocity (x)
  double v;							   // Control Volume Average Velocity (y)

  Field();
  Field& operator=(const Field &RHS);				/* Overload equation operator to equate fields directly */
  Field& operator+=(const Field &RHS);				/* Overload add to self operator */
  //friend ostream& operator<<(ostream &out, const Field U);    /* Overload cout operator later */

};							

//__________________________________Cell Definition_____________________________________//


struct Cell{
  
  Field U;							   // Field contains solution variables T, u, v
  Field eU;		                                       	   // Field containing exact field solution values

  double x;							/* Create coordinate system within Cell */
  double y;
    
  double eF;			                                /* Store exact flux integral per control volume */
  double eS;							/* Store exact Source term */

  double S;			                                   // Source within Cell
  double F;							   // Flux Integral for given control volume

public:
  // Computations of Source can be done within CV
  Cell();
  
  void SetCellCoordinates(double, double);			/* Set coordinates for a cell within a cell */
  void SetCellField(double, double, double);			/* Set field within a cell */

  void SetInitialField();					/* Set initial field for pertinent problem */
  double SetFDTemperatureField();				/* Set fully developed temperature field */

  void ComputeExactField();					   // Computing Exact Values is easy within CV
  void ComputeExactFluxIntegral();				/* Compute exact flux integral */
  void ComputeExactSourceTerm();				/* Compute exact source term */
  
  void PrintCoordinates();					/* Print coordinates of cell */
  void PrintField(FieldName);					/* Print Field values of given field */
    
};// End of Class Definition

//_________________________________Grid Definition________________________________________//

class Grid{
  
  Cell** Mesh;							   // Construct array of cells
  
  int Imin, Imax;						   // Set array bounds in x direction
  int Jmin, Jmax;						/* Set array bounds in y directoin */

  double dx;		 					   // This is possible for unifrom mesh, why not do it?
  double dy;							   // In essence cell size should stay with cell

  double dt;							/* The mesh needs to know the time step */

  double** Dx;							/* This is the tri diagonal matrix for X */
  double** Dy;							/* This is the tri diagonal matrix in Y */

  double** FI;

public:
  
  Grid();							/* construct grid object */
  Grid(int, int, int);						/* construct grid given sizes */
  ~Grid();							/* destroy grid object once used */
  
  void PrintCellCoordinates();					/* obvious functions for printing purposes */
  void PrintFluxes();
  void PrintSources();
  
  void PrintFieldValues(FieldName);				/* Printing fields for entire grid */
  void PrintFieldValues(FieldName, int, int);
  void WriteFieldValues(FieldName);

  void EvaluateExactIntegrals();				/* Evaluate exact functions, fields, etc */
  void EvaluateExactFieldValues();

  void EvaluateCellCoordinates();				/* Stepwise in sequence desired */
  void EvaluateInitialFields();   
  void EvaluateFluxes();
  void EvaluateSourceTerms();
  void EvaluateBoundaryConditions();
  
  void InitialiseDxDyFI();
  void EvaluateDx(int);						/* Evaluate coefficient matrices */
  void EvaluateDy(int);
  void EvaluateFI(int);

  double EvaluateTimeStep(TimeScheme);				/* take a different time step based on scheme */
  
  //void EulerExplicitTimeAdvance();
  //void RK2ExplicitTimeAdvance();
  Field EulerExplicitTimeAdvance();				/* key time advance schemes */
  Field EulerImplicitTimeAdvance();  

  void FieldVerification();					/* verify stuff once calculated */
  void FluxVerification();
  void SourceVerification();
  
};

//_________________________________Function Declarations__________________________________//

void EvaluateGridParameters(Grid&);				/*  */
void MarchTime(Grid&, TimeScheme);				/* march in time with given scheme */
void SolveEnergyEquation(Grid&, double);			/* Solve Energy equation for  */
void SolveThomas(double[NMAX][3], double[], const int);		/* Use Carl's Thomas Algorithm */
void CopyToRHS(double**, double RHS[], const int, const int, Direction RC);
void CopyToLHS(double**, double LHS [NMAX][3], const int);
void CopyFromRHS(double**, double RHS[], const int, const int, Direction RC);
