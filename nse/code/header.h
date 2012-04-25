/**
 * \brief This header file defines the class structure.
 *
 * \struct Field
 *			Provides a base for operations on fields
 *
 * \struct Cell
 *
 * \struct Rectangle
 *
 *
 * \class Grid
 */

#include <iostream>
#include <iomanip>
#include<fstream>		
#include<sstream>

#include<string.h>
#include<stdio.h>
#include<stdlib.h>

#include <math.h>
#include <time.h>

#include "constant.h"

using namespace std;

enum Direction{Row, Column};
enum FieldName{Pressure, xVelocity, yVelocity, StreamFunction, Vorticity, All};

struct Field
{
	double C[3];
	
	Field();
	
	// Overload operators for field algebra
	Field& operator=(const double &RHS);
	Field& operator=(const Field &RHS);
	Field& operator+=(const Field &RHS);
	Field operator/(const double &RHS);
	Field& operator/=(const double &RHS);
	Field& operator/=(const Field &RHS);
	Field& operator*=(const double &RHS);
	Field& operator*=(const Field &RHS);
	Field operator*(const double A[][3]);
	Field operator*(const Field& F);	
	bool operator>=(const Field &RHS);
	Field operator-(const Field &RHS);

	friend ostream& operator<<(ostream&, Field&);
};

struct Cell
{
	Field U, eU, dU, Un;
	Field FI, eFI, FIn;
	
	// Define streamfunction, and vorticity
	double psi, omega;
	
	double x, y;

	double Ax[3][3], Ay[3][3];
	double Bx[3][3], By[3][3];
	double Cx[3][3], Cy[3][3];
	double IDx[3][3][3], IDy[3][3][3];
	
	Cell();
	
	// Basic Calculations of Exact Variables
	void set_XY(double, double);
	void set_eU();
	void set_eFI();
	void set_eQ();
};

struct Rectangle
{
	int Imin, Imax;
	int Jmin, Jmax;

	double dx, dy;
		
	// Constructors for rectangle
	Rectangle();
	Rectangle(int, int, int);
	Rectangle& operator=(const Rectangle &RHS);
};

struct Grid
{
	Cell** mesh;
	Rectangle domain;
	
public:

	Grid(Rectangle);
	~Grid();

	// Calculate exact values
	void calc_eU();
	void calc_eFI();
	
	// Calculate control volume based quantities
	void calc_Q();
	void calc_BC();
	void calc_IBC();
	
	void calc_FI(int, int);
	void calc_J(int, int);
	void calc_EJ();	
	void calc_IDxDy(double);
	void calc_PSI();
	void calc_OMEGA();
	
	void add_J(int, int);
	void add_FI(int, int);
	
	void spew_mesh(int);
	void spew_field(FieldName);
	void spew_field(FieldName, int, int);
	void mirror_U();
	void find_center();
	// Functions for Thomas solution
	void calc_LHS(double LHS[][3][3][3], Direction, int);
	void calc_RHS(double RHS[][3], Direction, int);
	void solve_Thomas(double LHS[][3][3][3], double RHS[][3], Direction, int);
	
	// Implicit Euler time march with SOR
	Field march_IE(double, double);	
	
	Field calc_Uav();
	Field ver_FI();	
	void plot_U();
	void plot_uSL();
	void plot_uSL(string);
	void slice_U(int, int);
};

// Field operations
Field sqrt(Field);
Field fabs(Field);
double max(Field);

// Matrix operations
void spew_matrix(double M[][3]);
void mult_matrix(double M[][3], double K);
void init_matrix(double M[][3], double K);
void copy_matrix(double S[][3], double T[][3]);

// Operations on Navier Stokes Equations
void setup(Grid&);
void tune(Grid&, double);
void solve(Grid&);
void verify(Grid&);

// Carl's Thomasian Solver
void SolveBlockTri(double LHS[MAXSIZE][3][3][3], double RHS[MAXSIZE][3], int);

// tabulating and plotting
void tab_L2N(int, int, Field);
void tab_EJ(int, int, Field);
void plot_CH(int, Field, string);
void est_GCI();
void calc_GCI();
