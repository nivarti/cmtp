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

struct Field
{
	double C[3];
	
	Field();
	void sqroot();
	void abs();

	Field& operator=(const double &RHS);
	Field& operator=(const Field &RHS);
	Field& operator+=(const Field &RHS);
	Field& operator/=(const double &RHS);
	Field& operator/=(const Field &RHS);
	Field& operator*=(const double &RHS);
	Field& operator*=(const Field &RHS);
	Field operator*(const double A[][3]);
	bool operator>=(const Field &RHS);

	Field operator-(const Field &RHS);
	friend ostream& operator<<(ostream&, Field&);
};

struct Cell
{
	Field U, eU, dU, Un;
	Field FI, eFI, FIn;

	double x, y;

	double Ax[3][3], Ay[3][3];
	double Bx[3][3], By[3][3];
	double Cx[3][3], Cy[3][3];
	double IDx[3][3], IDy[3][3];
	
	Cell();
	
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
		
	Rectangle();
	Rectangle(int, int, int);
	Rectangle& operator=(const Rectangle &RHS);
};

class Grid
{
	Cell mesh[MAXSIZE][MAXSIZE];
	Rectangle domain;
	
public:

	Grid(Rectangle);
	
	void calc_Q();
	void calc_eU();
	void calc_eFI();
	void calc_FI(int, int);
	void calc_BC();
	void calc_IBC();
	void calc_J(int, int);
	void calc_EJ();
	
	void calc_LHS(double LHS[][3][3][3], Direction, int);
	void calc_RHS(double RHS[][3], Direction, int);
	void solve_Thomas(double LHS[][3][3][3], double RHS[][3], Direction, int);

	Field march_IE(double);

	Field ver_FI();
};

void spew_matrix(double M[][3]);
void mult_matrix(double M[][3], double K);
void init_matrix(double M[][3], double K);
void copy_matrix(double S[][3], double T[][3]);

//void plot_l2norm(Field, double);
void tab_L2N(int, int, Field);
void tab_EJ(int, int, Field);

//void solve(Grid&);
void SolveBlockTri(double LHS[MAXSIZE][3][3][3], double RHS[MAXSIZE][3], int);
