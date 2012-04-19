#include "header.h"


Rectangle::Rectangle()
{
	Imin = Imax = 1;
	Jmin = Jmax = 1;

	dx = dy = 1.0;
}
Rectangle::Rectangle(int Nx, int Ny, int Ngc)
{
	  dx = Lx/Nx;
	  dy = Ly/Ny;

	  Imin = Jmin = Ngc;
	  Imax = Nx + Ngc - 1;
	  Jmax = Ny + Ngc - 1;
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


void Rectangle::init_rect(int Nx, int Ny, int Ngc)
{
	  dx = Lx/Nx;
	  dy = Ly/Ny;
	  Imin = Jmin = Ngc;
	  Imax = Nx + Ngc - 1;
	  Jmax = Ny + Ngc - 1;
}
