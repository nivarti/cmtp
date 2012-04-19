#include "header.h"


Cell::Cell()
{
	x = y = 0.0;
	U = eU = 0.0;
	FI = eFI = 0.0;

	for(int i = 0; i < 3; i++){
		for(int j = 0; j < 3; j++){
			
			jP[i][j] = 0.0;
			jE[i][j] = 0.0;
			jW[i][j] = 0.0;
			jN[i][j] = 0.0;
			jS[i][j] = 0.0;					
			
		}
	}
}

void Cell::init_cell(double X, double Y)
{
	x = X; y = Y;
	set_e_field();
	set_e_fi();
	U = eU;
	
}

void Cell::set_e_field()
{
	eU.C[0] = P0*cos(pi*x)*cos(pi*y);
	eU.C[1] = u0*sin(pi*x)*sin(2.0*pi*y);
	eU.C[2] = v0*sin(2.0*pi*x)*sin(pi*y);
}

void Cell::set_e_fi()
{
	double Cx, Cy, Sx, Sy, C2x, C2y, S2x, S2y;

	Cx = cos(pi*x);
	Cy = cos(pi*y);
	Sx = sin(pi*x);
	Sy = sin(pi*y);

	C2x = cos(2.0*pi*x);
	C2y = cos(2.0*pi*y);
	S2x = sin(2.0*pi*x);
	S2y = sin(2.0*pi*y);

	eFI.C[0] = -(pi/beta)*(u0*Cx*S2y + v0*S2x*Cy);

	eFI.C[1] = P0*pi*Sx*Cy - u0*u0*pi*S2x*S2y*S2y
		- u0*v0*pi*Sx*S2x*(Cy*S2y + 2.0*C2y*Sy)
		- (u0/Re)*(5.0*pi*pi*Sx*S2y);

	eFI.C[2] = P0*pi*Cx*Sy - v0*v0*pi*S2x*S2x*S2y
		- u0*v0*pi*Sy*S2y*(Cx*S2x + 2.0*C2x*Sx)
		- (v0/Re)*(5.0*pi*pi*S2x*Sy);
}
