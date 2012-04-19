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
	Field operator-(const Field &RHS);

	friend ostream& operator<<(ostream&, Field&);
};

struct Cell
{
	Field U, eU;
	Field FI, eFI;
	
	double x, y;

	double jP[3][3];
	double jE[3][3];
	double jN[3][3];
	double jW[3][3];
	double jS[3][3];
	
	Cell();
	void init_cell(double, double);
	void set_e_field();
	void set_e_fi();
};

struct Rectangle
{
	int Imin, Imax;
	int Jmin, Jmax;

	double dx;
	double dy;

	Rectangle();
	Rectangle(int, int, int);
	void init_rect(int, int, int);
	Rectangle& operator=(const Rectangle &RHS);
};

class Grid
{
	Cell mesh[MAXSIZE][MAXSIZE];
	Rectangle domain;

public:

	//Grid();
	Grid(int, int, int);
	Grid(Rectangle);
	
	void init_grid(Rectangle);
	void calc_cell_coord();
	void calc_e_field();
	void calc_e_fi();
	void calc_fi();
	Field verif_fi();
	void calc_Jacobian(int, int);
};


void plot_l2norm(Field, double);
void tab_l2norm(Field, double);

void SolveBlockTri(double LHS[MAXSIZE][3][3][3], double RHS[MAXSIZE][3], int);
