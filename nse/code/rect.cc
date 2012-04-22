#include "header.h"

Rectangle::Rectangle()
{
	dx = dy = 1.0;
	Imin = Imax = 1;
	Jmin = Jmax = 1;
}

Rectangle::Rectangle(int Nx, int Ny, int Ngc)
{
	// Calculate cell sizes in either direction
	dx = Lx/Nx;
	dy = Ly/Ny;
	
	// Set bounds for array indices
	Imin = Jmin = Ngc;
	Imax = Nx + Ngc - 1;
	Jmax = Ny + Ngc - 1;
	
	// Output rectangle information
	cout<<"\nSetting up geometry...";
	cout<<"\n Cavity width: "<<Lx<<", cavity height: "<<Ly;
}

Rectangle& Rectangle::operator=(const Rectangle &RHS)
{
	Imin = RHS.Imin;
	Jmin = RHS.Jmin;

	Imax = RHS.Imax;
	Jmax = RHS.Jmax;

	dx = RHS.dx;
	dy = RHS.dy;

	return *this;
}
