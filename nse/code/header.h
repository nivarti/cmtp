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

struct Field
{
	double C[3];
	
	Field();
	void sqroot();

	Field& operator=(const double &RHS);
	Field& operator=(const Field &RHS);
	Field& operator+=(const Field &RHS);
	Field& operator/=(const double &RHS);
	Field& operator/=(const Field &RHS);
	Field& operator*=(const double &RHS);
	Field& operator*=(const Field &RHS);
	Field operator*(const double A[][3]);

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
	void calc_J(int, int);
	void calc_EJ();
	
	Field ver_FI();
};


void plot_l2norm(Field, double);
void tab_l2norm(Field, double);

void SolveBlockTri(double LHS[MAXSIZE][3][3][3], double RHS[MAXSIZE][3], int);
