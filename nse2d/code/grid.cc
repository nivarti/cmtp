#include "header.h"

// Grid::Grid()
// {
// 	;
// }

Grid::Grid(Rectangle rect)
{
	double x, y;

	domain = rect;
	cout<<"\nInitialising mesh of size: "<<domain.Imax<<" x "<<domain.Jmax;
	for (int j = 0; j <= domain.Jmax + domain.Jmin; j++)
		{
			for (int i = 0; i <= domain.Imax + domain.Imin; i++)
				{
					x = domain.dx*(2.0*(double)i - (double)domain.Imin)/2.0;
					y = domain.dy*(2.0*(double)j - (double)domain.Jmin)/2.0;

					mesh[i][j].init_cell(x, y);
				}
		}    
	cout<<"\nEvaluated cell fields...";
}

// Grid::Grid(int Nx, int Ny, int Ngc)
// {
// 	domain.init_rect(Nx, Ny, Ngc);
// 	mesh = new Cell*[domain.Imax + 2];
// 	for(int i = 0; i < domain.Imax + 2; i++)
// 		mesh[i] = new Cell[domain.Jmax + 2];
// }

// // Grid::~Grid()
// // {
// // 	for(int i = 0; i < domain.Imax + 2; i++)
// // 		delete mesh[i];
	
// // 	delete mesh;
// // }


// void Grid::init_grid(Rectangle rect)
// {
// 	domain = rect;
// 	for (int j = 0; j <= domain.Jmax + domain.Jmin; j++)
// 	{
// 		for (int i = 0; i <= domain.Imax + domain.Imin; i++)
// 		{
// 			mesh[i][j] = Cell();
// 			mesh[i][j].init_cell();
// 		}
// 	}
// }

// void Grid::calc_cell_coord()
// {
// 	double x, y;

// 	for (int j = 0; j <= domain.Jmax + domain.Jmin; j++)
// 	{
// 		for (int i = 0; i <= domain.Imax + domain.Imin; i++)
// 		{
// 			x = domain.dx * (2.0 * (double) i - (double) domain.Imin) / 2.0;
// 			y = domain.dy * (2.0 * (double) j - (double) domain.Jmin) / 2.0;

// 			mesh[i][j].set_coord(x, y);
// 		}
// 	}

// 	cout << "\nEvaluated cell coordinates for entire grid...";
// }

// void Grid::calc_e_field()
// {
// 	for(int j = 0; j <= domain.Jmax + domain.Jmin; j++)
// 	{
// 		for(int i = 0; i <= domain.Imax + domain.Imin; i++)
// 		{
// 			mesh[i][j].set_e_field();
// 		}
// 	}
// 	cout<<"\nEvaluated exact fields...";
// }

// void Grid::calc_e_fi()
// {
// 	for(int j = 0; j <= domain.Jmax + domain.Jmin; j++)
// 	{
// 		for(int i = 0; i <= domain.Imax + domain.Imin; i++)
// 		{
// 			mesh[i][j].set_e_fi();
// 		}
// 	}
// 	cout<<"\nEvaluated exact flux integrals...";
// }

void Grid::calc_fi()
{
	Field Fx, Gy;
	double pN, pS, pE, pW, uN, uS, uE, uW, uP, vN, vS, vW, vE, vP;

	for(int j = domain.Jmin; j <= domain.Jmax; j++)
		{
			for(int i = domain.Imin; i <= domain.Imax; i++)
				{
					pN = mesh[i][j+1].U.C[0];
					pS = mesh[i][j-1].U.C[0];
					pE = mesh[i+1][j].U.C[0];
					pW = mesh[i-1][j].U.C[0];
					
					uP = mesh[i][j].U.C[1];
					uE = mesh[i+1][j].U.C[1];
					uW = mesh[i-1][j].U.C[1];
					uN = mesh[i][j+1].U.C[1];
					uS = mesh[i][j-1].U.C[1];
					
					vP = mesh[i][j].U.C[2];
					vN = mesh[i][j+1].U.C[2];
					vS = mesh[i][j-1].U.C[2];
					vE = mesh[i+1][j].U.C[2];
					vW = mesh[i-1][j].U.C[2];
					
					// Store, dF/dx, add dG/dx, invert sign
					Fx.C[0] = (uE - uW)/(2.0*beta);
					Fx.C[1] = ((uP + uE)*(uP + uE) - (uP + uW)*(uP + uW))/4.0
						+ (pE - pW)/2.0
						- (1/(Re*domain.dx))*(uE - 2.0*uP + uW);
					Fx.C[2] = ((uP + uE)*(vP + vE) - (uP + uW)*(vP + vW))/4.0
						- (1/(Re*domain.dx))*(vE - 2.0*vP + vW);
					Fx /= -domain.dx;
					
					Gy.C[0] = (vN - vS)/(2.0*beta);
					Gy.C[1] = ((uP + uN)*(vP + vN) - (uP + uS)*(vP + vS))/4.0
						- (1/(Re*domain.dy))*(uN - 2.0*uP + uS);
					Gy.C[2] = ((vP + vN)*(vP + vN) - (vP + vS)*(vP + vS))/4.0
						+ (pN - pS)/2.0
						- (1/(Re*domain.dy))*(vN - 2.0*vP + vS);
					Gy /= -domain.dy;
					
					mesh[i][j].FI += Fx;
					mesh[i][j].FI += Gy;
				}
		}

	cout<<"\nInformation from cell 1,1";
	cout<<"\nCoordinates: "<<mesh[1][1].x<<", "<<mesh[1][1].y;
	cout<<"\nValue of eFI: "<<mesh[1][1].eFI;
	cout<<"\nValue of FI: "<<mesh[1][1].FI;
	
	cout<<"\nEvaluated control volume average flux integrals...";
}

Field Grid::verif_fi()
{
	Field error, l2norm;

	error = l2norm = 0.0;

	for(int j = domain.Jmin; j <= domain.Jmax; j++)
		{
			for(int i = domain.Imin; i <= domain.Imax; i++)
				{
					error = mesh[i][j].eFI - mesh[i][j].FI;
					error *= error;
					l2norm += error;
				}
		}

	l2norm.sqroot();
	l2norm /= sqrt(domain.Imax*domain.Jmax);

	cout<<"\nL2 Norm for flux integral:"<<l2norm;
	return l2norm;
}
