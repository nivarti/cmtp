#include "header.h"


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
					
					calc_Jacobian(i, j);
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

void Grid::calc_Jacobian(int i, int j)
{
	Field Fx, Gy;
	double uN, uS, uE, uW, uP, vN, vS, vW, vE, vP;
	double dx = domain.dx, dy = domain.dy;

	double dt = 0.1;

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
						
	// I + dt*Bx + dt*By
	mesh[i][j].jP[0][0] = 1.0;
	mesh[i][j].jP[0][1] = 0.0;
	mesh[i][j].jP[0][2] = 0.0;
	
	mesh[i][j].jP[1][0] = 0.0;
	mesh[i][j].jP[1][1] = 1.0 + dt*(((uE - uW)/2.0 + 2.0/(Re*dx))/dx
					+ ((vE - vW)/2.0 + 2.0/(Re*dy))/dy);
	mesh[i][j].jP[1][2] = dt/dx*((uE - uW)/4.0);
	
	mesh[i][j].jP[2][0] = 0.0;
	mesh[i][j].jP[2][1] = dt/dx*((vE - vW)/4.0);
	mesh[i][j].jP[2][2] = 1.0 + dt*(((uE - uW)/2.0 + 2.0/(Re*dx))/dx
					+ ((vE - vW)/2.0 + 2.0/(Re*dy))/dy);
	
	// Cx
	mesh[i][j].jE[0][0] = 0.0;
	mesh[i][j].jE[0][1] = dt/dx*1.0/(2.0*beta);
	mesh[i][j].jE[0][2] = 0.0;
	
	mesh[i][j].jE[1][0] = dt/dx*0.5;
	mesh[i][j].jE[1][1] = dt/dx*((uE + uP)/2.0 - 1.0/(Re*dx));
	mesh[i][j].jE[1][2] = 0.0;
	
	mesh[i][j].jE[2][0] = 0.0;
	mesh[i][j].jE[2][1] = dt/dx*((vE + vP)/4.0);
	mesh[i][j].jE[2][2] = dt/dx*((uE + uP)/4.0 - 1.0/(Re*dx));
	
	// Ax
	mesh[i][j].jW[0][0] = 0.0;
	mesh[i][j].jW[0][1] = -dt/dx*1.0/(2.0*beta);
	mesh[i][j].jW[0][2] = 0.0;
	
	mesh[i][j].jW[1][0] = -dt/dx*0.5;
	mesh[i][j].jW[1][1] = -dt/dx*((uW + uP)/2.0 + 1.0/(Re*dx));
	mesh[i][j].jW[1][2] = 0.0;
	
	mesh[i][j].jW[2][0] = 0.0;
	mesh[i][j].jW[2][1] = -dt/dx*((vW + vP)/4.0);
	mesh[i][j].jW[2][2] = -dt/dx*((uW + uP)/4.0 + 1.0/(Re*dx));
	
	// Cy 
	mesh[i][j].jN[0][0] = 0.0;
	mesh[i][j].jN[0][1] = dt/dy*1.0/(2.0*beta);
	mesh[i][j].jN[0][2] = 0.0;

	mesh[i][j].jN[1][0] = dt/dy*0.5;
	mesh[i][j].jN[1][1] = dt/dy*((vN + vP)/2.0 - 1.0/(Re*dy));
	mesh[i][j].jN[1][2] = dt/dy*((uN + uP)/4.0);
	
	mesh[i][j].jN[2][0] = 0.0;
	mesh[i][j].jN[2][1] = 0.0;
	mesh[i][j].jN[2][2] = dt/dy*((vN + vP)/4.0 - 1.0/(Re*dy));
	
	// Ay
	mesh[i][j].jS[0][0] = 0.0;
	mesh[i][j].jS[0][1] = -dt/dy*1.0/(2.0*beta);
	mesh[i][j].jS[0][2] = 0.0;
	
	mesh[i][j].jS[1][0] = -dt/dy*0.5;
	mesh[i][j].jS[1][1] = -dt/dy*((vS + vP)/2.0 + 1.0/(Re*dy));
	mesh[i][j].jS[1][2] = -dt/dy*((uS + uP)/4.0);
	
	mesh[i][j].jS[2][0] = 0.0;
	mesh[i][j].jS[2][1] = 0.0;
	mesh[i][j].jS[2][2] = -dt/dy*((vS + vP)/4.0 + 1.0/(Re*dy));	
}
