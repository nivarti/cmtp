#include "header.h"


Cell::Cell()
{
	x = y = 0.0;
	U = eU = dU = Un = 0.0;
	FI = eFI = FIn = 0.0;

	for(int i = 0; i < 3; i++){
		for(int j = 0; j < 3; j++){
			
			Ax[i][j] = Ay[i][j] = 0.0;
			Bx[i][j] = By[i][j] = 0.0;
			Cx[i][j] = Cy[i][j] = 0.0;	

			for(int k = 0; k < 3; k++){
				
				IDx[i][j][k] = 0.0;
				IDy[i][j][k] = 0.0;
			}		
		}
	}
}

void Cell::set_XY(double X, double Y)
{
	x = X; y = Y;
}

void Cell::set_eU()
{
	eU.C[0] = P0*cos(pi*x)*cos(pi*y);
	eU.C[1] = u0*sin(pi*x)*sin(2.0*pi*y);
	eU.C[2] = v0*sin(2.0*pi*x)*sin(pi*y);
}

void Cell::set_eFI()
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

void Cell::set_eQ()
{
	set_eU();
	set_eFI();
	U = eU;
}
