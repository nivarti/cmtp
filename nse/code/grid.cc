#include "header.h"


Grid::Grid(Rectangle rect)
{
	double x, y;

	domain = rect;	
	cout<<"\nInitialising mesh for cavity... Size:  "<<domain.Imax<<" x "<<domain.Jmax;
	
	// Allocate space for mesh
	mesh = new Cell*[domain.Imax + domain.Imin + 1];
	for(int i = 0; i < domain.Imax + domain.Imin + 1; i++)
		mesh[i] = new Cell[domain.Jmax + domain.Jmin + 1];
	
	for (int j = 0; j <= domain.Jmax + domain.Jmin; j++){
		for (int i = 0; i <= domain.Imax + domain.Imin; i++){
			
			x = domain.dx*(2.0*(double)i - (double)domain.Imin)/2.0;
			y = domain.dy*(2.0*(double)j - (double)domain.Jmin)/2.0;
			
			// Store coordinates, and exact values of fields
			mesh[i][j].set_XY(x, y);
			mesh[i][j].set_eQ();		
		}
	}    
	cout<<"\nEvaluated cell coordinates...";
	cout<<"\nEvaluated exact fields, and exact flux integrals...";
}

Grid::~Grid(){
	
	cout<<"\nDeleting mesh...";
	// delete mesh
	for(int i = 0; i <= domain.Imax + domain.Imin; i++)
		delete mesh[i];
	delete mesh;	
}

void Grid::calc_FI(int i, int j)
{
	Field Fx, Gy;
	double pN, pS, pE, pW, uN, uS, uE, uW, uP, vN, vS, vW, vE, vP;
	
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
	
	// Set flux integral as -(Fx + Gy)
	mesh[i][j].FI = Fx;
	mesh[i][j].FI += Gy;
}

void Grid::add_FI(int i, int j)
{
	double dx, dy, pP, pN, pS, pE, pW, dFx, dGy;
	
	dx = domain.dx;
	dy = domain.dy;
	
	pP = mesh[i][j].U.C[0];
	pN = mesh[i][j+1].U.C[0];
	pS = mesh[i][j-1].U.C[0];
	pE = mesh[i+1][j].U.C[0];
	pW = mesh[i-1][j].U.C[0];
	
	// Calculate addition required in flux integral
	dFx = ((A*dy)/(beta*dx))*(pE - 2.0*pP + pW);
	dGy = ((A*dx)/(beta*dy))*(pN - 2.0*pP + pS);
	
	mesh[i][j].FI.C[0] += dFx;
	mesh[i][j].FI.C[0] += dGy;
}

void Grid::calc_J(int i, int j)
{
	double uN, uS, uE, uW, uP, vN, vS, vW, vE, vP, dx, dy;

	dx = domain.dx;
	dy = domain.dy;

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
	
	// Ax
	mesh[i][j].Ax[0][0] = 0.0;
	mesh[i][j].Ax[0][1] = -1.0/(2.0*beta*dx);
	mesh[i][j].Ax[0][2] = 0.0;
	
	mesh[i][j].Ax[1][0] = -0.5/dx;
	mesh[i][j].Ax[1][1] = -((uW + uP)/2.0 + 1.0/(Re*dx))/dx;
	mesh[i][j].Ax[1][2] = 0.0;
	
	mesh[i][j].Ax[2][0] = 0.0;
	mesh[i][j].Ax[2][1] = -((vW + vP)/4.0)/dx;
	mesh[i][j].Ax[2][2] = -((uW + uP)/4.0 + 1.0/(Re*dx))/dx;
					
	// Bx
	mesh[i][j].Bx[0][0] = 0.0;
	mesh[i][j].Bx[0][1] = 0.0;
	mesh[i][j].Bx[0][2] = 0.0;
	
	mesh[i][j].Bx[1][0] = 0.0;
	mesh[i][j].Bx[1][1] = ((uE - uW)/2.0 + 2.0/(Re*dx))/dx;
	mesh[i][j].Bx[1][2] = 0.0;
	
	mesh[i][j].Bx[2][0] = 0.0;
	mesh[i][j].Bx[2][1] = ((vE - vW)/4.0)/dx;
	mesh[i][j].Bx[2][2] = ((uE - uW)/4.0 + 2.0/(Re*dx))/dx;
	
	// Cx
	mesh[i][j].Cx[0][0] = 0.0;
	mesh[i][j].Cx[0][1] = 1.0/(2.0*beta*dx);
	mesh[i][j].Cx[0][2] = 0.0;
	
	mesh[i][j].Cx[1][0] = 0.5/dx;
	mesh[i][j].Cx[1][1] = ((uE + uP)/2.0 - 1.0/(Re*dx))/dx;
	mesh[i][j].Cx[1][2] = 0.0;
	
	mesh[i][j].Cx[2][0] = 0.0;
	mesh[i][j].Cx[2][1] = ((vE + vP)/(4.0*dx));
	mesh[i][j].Cx[2][2] = ((uE + uP)/4.0 - 1.0/(Re*dx))/dx;
	
	
	// Ay
	mesh[i][j].Ay[0][0] = 0.0;
	mesh[i][j].Ay[0][1] = 0.0;
	mesh[i][j].Ay[0][2] = -1.0/(2.0*beta*dy);
	
	mesh[i][j].Ay[1][0] = 0.0;
	mesh[i][j].Ay[1][1] = -((vS + vP)/4.0 + 1.0/(Re*dy))/dy;
	mesh[i][j].Ay[1][2] = -((uS + uP)/4.0)/dy;
	
	mesh[i][j].Ay[2][0] = -0.5/dy;
	mesh[i][j].Ay[2][1] = 0.0;
	mesh[i][j].Ay[2][2] = -((vS + vP)/2.0 + 1.0/(Re*dy))/dy;	
	
	// By
	mesh[i][j].By[0][0] = 0.0;
	mesh[i][j].By[0][1] = 0.0;
	mesh[i][j].By[0][2] = 0.0;
	
	mesh[i][j].By[1][0] = 0.0;
	mesh[i][j].By[1][1] = ((vN - vS)/4.0 + 2.0/(Re*dy))/dy;
	mesh[i][j].By[1][2] = ((uN - uS)/4.0)/dy;

	mesh[i][j].By[2][0] = 0.0;
	mesh[i][j].By[2][1] = 0.0;
	mesh[i][j].By[2][2] = ((vN - vS)/2.0 + 2.0/(Re*dy))/dy;

	// Cy 
	mesh[i][j].Cy[0][0] = 0.0;
	mesh[i][j].Cy[0][1] = 0.0;
	mesh[i][j].Cy[0][2] = 1.0/(2.0*beta*dy);

	mesh[i][j].Cy[1][0] = 0.0;
	mesh[i][j].Cy[1][1] = ((vN + vP)/4.0 - 1.0/(Re*dy))/dy;
	mesh[i][j].Cy[1][2] = ((uN + uP)/(4.0*dy));
	
	mesh[i][j].Cy[2][0] = 0.5/dy;
	mesh[i][j].Cy[2][1] = 0.0;
	mesh[i][j].Cy[2][2] = ((vN + vP)/2.0 - 1.0/(Re*dy))/dy;

}

void Grid::add_J(int i, int j)
{
	double dx, dy;

	dx = domain.dx;
	dy = domain.dy;
	
	// Set Jacobians to accommodate smoothing Pressure
	mesh[i][j].Ax[0][0] += - (A*dy)/(beta*dx);
	mesh[i][j].Bx[0][0] += (2.0*A*dy)/(beta*dx);
	mesh[i][j].Cx[0][0] += - (A*dy)/(beta*dx);
	
	mesh[i][j].Ay[0][0] += - (A*dx)/(beta*dy);
	mesh[i][j].By[0][0] += (2.0*A*dx)/(beta*dy);
	mesh[i][j].Cy[0][0] += - (A*dx)/(beta*dy);
}

void Grid::calc_BC()
{		
	for(int j = domain.Jmin; j <= domain.Jmax; j++){
		
		// Neumann on Pressure
		mesh[0][j].U.C[0] = mesh[1][j].U.C[0];
		mesh[domain.Imin + domain.Imax][j].U.C[0] = mesh[domain.Imax][j].U.C[0];
		
		// No Slip on velocities
		mesh[0][j].U.C[1] =  - mesh[1][j].U.C[1];   		
		mesh[domain.Imin + domain.Imax][j].U.C[1] = - mesh[domain.Imax][j].U.C[1];
		mesh[0][j].U.C[2] = - mesh[1][j].U.C[2];   		
		mesh[domain.Imin + domain.Imax][j].U.C[2] = - mesh[domain.Imax][j].U.C[2];
	}  
	
	for(int i = domain.Imin; i <= domain.Imax; i++){
		
		// Neumann on Pressure
		mesh[i][0].U.C[0] = mesh[i][1].U.C[0];   		
		mesh[i][domain.Jmin + domain.Jmax].U.C[0] = mesh[i][domain.Jmax].U.C[0];
		
		// No Slip on velocities
		mesh[i][0].U.C[1] = - mesh[i][1].U.C[1];
		mesh[i][domain.Jmin + domain.Jmax].U.C[1] = 2.0*U_TOP - mesh[i][domain.Jmax].U.C[1];

		mesh[i][0].U.C[2] = -mesh[i][1].U.C[2];   		
		mesh[i][domain.Jmin + domain.Jmax].U.C[2] = - mesh[i][domain.Jmax].U.C[2];
  }  
}

void Grid::calc_IBC()
{
	for(int j = domain.Jmin; j <= domain.Jmax; j++){
		
		// Neumann on Pressure, No Slip of velocities
		mesh[0][j].IDx[1][0][0] = 1.0;
		mesh[0][j].IDx[1][1][1] = 1.0;
		mesh[0][j].IDx[1][2][2] = 1.0;		

		mesh[0][j].IDx[2][0][0] = -1.0;
		mesh[0][j].IDx[2][1][1] = 1.0;
		mesh[0][j].IDx[2][2][2] = 1.0;		
		
		mesh[domain.Imax + 1][j].IDx[0][0][0] = 1.0;
		mesh[domain.Imax + 1][j].IDx[0][1][1] = 1.0;
		mesh[domain.Imax + 1][j].IDx[0][2][2] = 1.0;		

		mesh[domain.Imax + 1][j].IDx[1][0][0] = -1.0;
		mesh[domain.Imax + 1][j].IDx[1][1][1] = 1.0;
		mesh[domain.Imax + 1][j].IDx[1][2][2] = 1.0;		
	}

	for(int i = domain.Imin; i <= domain.Imax; i++){
		
		// Neumann on Pressure, No Slip of velocities
		mesh[i][0].IDy[1][0][0] = 1.0;
		mesh[i][0].IDy[1][1][1] = 1.0;
		mesh[i][0].IDy[1][2][2] = 1.0;

		mesh[i][0].IDy[2][0][0] = -1.0;
		mesh[i][0].IDy[2][1][1] = 1.0;
		mesh[i][0].IDy[2][2][2] = 1.0;		
		
		mesh[i][domain.Jmax + 1].IDy[0][0][0] = 1.0;
		mesh[i][domain.Jmax + 1].IDy[0][1][1] = 1.0;
		mesh[i][domain.Jmax + 1].IDy[0][2][2] = 1.0;		

		mesh[i][domain.Jmax + 1].IDy[1][0][0] = -1.0;
		mesh[i][domain.Jmax + 1].IDy[1][1][1] = 1.0;
		mesh[i][domain.Jmax + 1].IDy[1][2][2] = 1.0;		
	}
}

void Grid::calc_IDxDy(double dt)
{	
	for(int j = domain.Jmin; j <= domain.Jmax; j++){
		for(int i = domain.Imin; i <= domain.Imax; i++){
			
			mesh[i][j].FI *= dt;
			
			for(int k = 0; k < 3; k++){
				for(int l = 0; l < 3; l++){
					
					mesh[i][j].IDx[0][k][l] = dt*mesh[i][j].Ax[k][l];
					mesh[i][j].IDx[1][k][l] = dt*mesh[i][j].Bx[k][l];					
					mesh[i][j].IDx[2][k][l] = dt*mesh[i][j].Cx[k][l];					
					
					mesh[i][j].IDy[0][k][l] = dt*mesh[i][j].Ay[k][l];
					mesh[i][j].IDy[1][k][l] = dt*mesh[i][j].By[k][l];					
					mesh[i][j].IDy[2][k][l] = dt*mesh[i][j].Cy[k][l];
				}

				mesh[i][j].IDx[1][k][k] += 1.0;
				mesh[i][j].IDy[1][k][k] += 1.0;
			}						
		}
	}
}

void Grid::calc_Q()
{	
	for(int j = domain.Jmin; j <= domain.Jmax; j++){
		for(int i = domain.Imin; i <= domain.Imax; i++){
			
			// Calculate flux integrals, and jacobians
			calc_FI(i, j);
			calc_J(i, j);
			
			// Smoothen by adding pressure laplacian
			add_FI(i, j);
			add_J(i, j);
		}
	}
	
	//calc_IBC();
}

// Verify Flux Integral calculations w.r.t exact Flux Integrals
Field Grid::ver_FI()
{
	Field error, l2norm;

	error = l2norm = 0.0;

	for(int j = domain.Jmin; j <= domain.Jmax; j++){
		for(int i = domain.Imin; i <= domain.Imax; i++){
			
			error = mesh[i][j].eFI - mesh[i][j].FI;
			error *= error;
			l2norm += error;
		}
	}

	l2norm = sqrt(l2norm/(domain.Imax*domain.Jmax));

	return l2norm;
}

// Calculate changes in Jacobians
void Grid::calc_EJ()
{		
	Field LHS, RHS, E;
	mesh[10][10].dU = 0.000001;		
	
	for(int j = domain.Jmin; j <= domain.Jmax; j++){
		for(int i = domain.Imin; i <= domain.Imax; i++){			
			
			mesh[i][j].FIn = mesh[i][j].FI;
			mesh[i][j].Un = mesh[i][j].U;
			mesh[i][j].U += mesh[i][j].dU;
			
			// if (i==10 && j ==10){
				
			// 	setprecision(30);
			// 	setw(40);
			// 	cout<<"\nInformation from cell: "<<i<<", "<<j;
			// 	cout<<"\n Coordinates: "<<mesh[i][j].x<<", "<<mesh[i][j].y;
			// 	cout<<"\n\n Value of eFI: "<<mesh[i][j].eFI;
			// 	cout<<"\n Value of FI(n): "<<mesh[i][j].FIn;				
				
			// 	cout<<"\n\n Exact Solution: "<<mesh[i][j].eU;
			// 	cout<<"\n Numerical Solution: "<<mesh[i][j].Un;
			// 	cout<<"\n Change in Solution: "<<mesh[i][j].dU;
			// 	cout<<"\n New Solution: "<<mesh[i][j].U;				
				
			// }
		}
	}

	for(int j = domain.Jmin; j <= domain.Jmax; j++){
		for(int i = domain.Imin; i <= domain.Imax; i++){			
			
			calc_FI(i, j);
			LHS = mesh[i][j].FIn - mesh[i][j].FI;
			
			RHS = mesh[i-1][j].dU*mesh[i][j].Ax;
			RHS += mesh[i][j].dU*mesh[i][j].Bx;
			RHS += mesh[i+1][j].dU*mesh[i][j].Cx;
						
			RHS += mesh[i][j-1].dU*mesh[i][j].Ay;
			RHS += mesh[i][j].dU*mesh[i][j].By;
			RHS += mesh[i][j+1].dU*mesh[i][j].Cy;
						
			E = LHS - RHS;
			
			if (i>=9 && i <= 11 && j >= 9 && j <= 11)
				{
					cout<<"\n Value of FI(n + 1): "<<mesh[i][i].FI;				
					cout<<"\n\n Information from cell: "<<i<<", "<<j;
					cout<<"\n LHS: "<<LHS;
					cout<<"\n RHS: "<<RHS;
					cout<<"\n Error  = LHS - RHS = "<<E;
					tab_EJ(i, j, E);
				}			
		}
	}
}

Field Grid::march_IE(double dt, double SOR)
{	
	static int n = 1;
	Field dU, dUmax, L2N;
	double RHS[MAXSIZE][3], LHS[MAXSIZE][3][3][3];

	dU = dUmax = L2N = 0.0;

	// Construct flux integrals, and boundary conditions	
	calc_BC();
	calc_Q();
	//add_J();
	
	// Flux integrals have been found to change in ghost cells, please change code!
	calc_IDxDy(dt);
	calc_IBC();
	
	//cout<<"\n Flux integral...\n";
	//spew_mesh(1);

	// Start approximate factorisation solution
	for(int j = domain.Jmin; j <= domain.Jmax; j++){
		
		calc_LHS(LHS, Column, j);
		calc_RHS(RHS, Column, j);
		solve_Thomas(LHS, RHS, Column, j);
	}
	//cout<<"\n Lines of constant j...\n";
	//spew_mesh(1);

	for(int i = domain.Imin; i <= domain.Imax; i++){
		
		calc_LHS(LHS, Row, i);
		calc_RHS(RHS, Row, i);
		solve_Thomas(LHS, RHS, Row, i);
	}
	// End approximate factorisation solution
	
	//cout<<"\n Lines of constant i...\n";
	//spew_mesh(1);
	
	// Update solution, use successive-over-relaxation
	for(int i = domain.Imin; i <= domain.Imax; i++){
		for(int j = domain.Jmin; j <= domain.Jmax; j++){      
			
			mesh[i][j].dU = mesh[i][j].FI;
			mesh[i][j].dU *= SOR;
			mesh[i][j].U += mesh[i][j].dU;
			
			dU = mesh[i][j].dU;
			L2N += (dU*dU);			
		}
	}
	
	L2N = sqrt(L2N/(domain.Imax*domain.Jmax));
		
	n++;	
	return L2N;	
}

// Mean pressure, velocities over domain
Field Grid::calc_Uav()
{
	Field Uav;
	
	Uav = 0.0;
		
	for(int i = domain.Imin; i <= domain.Imax; i++){
		for(int j = domain.Jmin; j <= domain.Jmax; j++){
			
			// Average fields over entire domain
			Uav += mesh[i][j].U;
		}
	}
	
	return Uav;
}

// Plot fields
void Grid::plot_U()
{
	ofstream pfile, ufile, vfile, sfile;
	
	pfile.open("../plot/basic/pressure/P3");
	ufile.open("../plot/basic/velocity/u");
	vfile.open("../plot/basic/velocity/v");
	//sfile.open("../plot/vortex/streamline/psi");

	for(int i = domain.Imin; i <= domain.Imax; i++){
		for(int j = domain.Jmin; j <= domain.Jmax; j++){
			
			pfile<<mesh[i][j].x<<" "<<mesh[i][j].y<<" "<<mesh[i][j].U.C[0]<<endl;
			ufile<<mesh[i][j].x<<" "<<mesh[i][j].y<<" "<<mesh[i][j].U.C[1]<<endl;
			vfile<<mesh[i][j].x<<" "<<mesh[i][j].y<<" "<<mesh[i][j].U.C[2]<<endl;
			//sfile<<mesh[i][j].x<<" "<<mesh[i][j].y<<" "<<mesh[i][j].psi<<endl;
		}
		
		pfile<<endl;
		ufile<<endl;
		vfile<<endl;
		//sfile<<endl;
	}
		
	pfile.close();
	ufile.close();
	vfile.close();
	//sfile.close();
}

// Velocity along line of symmetry
void Grid::plot_uSL()
{
	double uM;
	int I = domain.Imax/2;
	ofstream file;
	
	file.open("../plot/basic/symmetry/uvel2");
	
	for(int j = domain.Jmin; j <= domain.Jmax; j++){			
		
		uM = (mesh[I][j].U.C[1] + mesh[I+1][j].U.C[1])/2.0;
		file<<mesh[I][j].y<<" "<<uM<<endl;
	}

	file.close();				
}


void Grid::plot_uSL(string F)
{
	double uM;
	static int i = 4;
	int I = domain.Imax/2;
	stringstream myFN, myFN2;
	string name, name2;
	ofstream file, file2;
	
	myFN << F << i;
	myFN2 << F << i << i;
	name = myFN.str();
	name2 = myFN2.str();

	file.open(name.c_str());
	file2.open(name2.c_str());
	
	for(int j = domain.Jmin; j <= domain.Jmax; j++){			
		
		uM = (mesh[I][j].U.C[1] + mesh[I+1][j].U.C[1])/2.0;
		file<<mesh[I][j].y<<" "<<uM<<endl;
	}
	
	// double yC = 0.05, yM;
	// for(int j = domain.Jmin; j < domain.Jmax; j++){
		
	// 	yM = (mesh[I][j].y + mesh[I][j+1].y)/2.0;
		
	// 	if(yM >= yC -0.00001 && yM <= yC + 0.00001){
			
	// 		uM = (mesh[I][j].U.C[1] + mesh[I+1][j].U.C[1] + mesh[I][j+1].U.C[1] + mesh[I+1][j+1].U.C[1])/4.0;		
	// 		file2<<yM<<" "<<uM<<endl;				
	// 		yC += 0.1;
			
	// 	}
	// }
	
	i++;
	file.close();
	file2.close();
}

void Grid::mirror_U()
{	
	for(int i = domain.Imin; i <= domain.Imax; i++){
		for(int j = domain.Jmin; j <= domain.Jmax; j++){
			
			mesh[i][j].Un = mesh[i][j].U;
		}
	}
	for(int i = 0; i <= domain.Imax - 1; i++){
		for(int j = domain.Jmin; j <= domain.Jmax; j++){
						
			mesh[i + 1][j].U = mesh[domain.Imax - i][j].Un;
		}
	}	
}

void Grid::slice_U(int I, int J)
{
	double pM;
	//int I = domain.Imax/2;
	ofstream file1, file2;
	
	file1.open("../plot/oscillation/P6i");
	file2.open("../plot/oscillation/P6j");
	
	for(int j = domain.Jmin; j <= domain.Jmax; j++){			
		
		pM = (mesh[I][j].U.C[0] + mesh[I][j].U.C[0])/2.0;
		file1<<mesh[I][j].y<<" "<<pM<<endl;
	}

	for(int i = domain.Imin; i <= domain.Imax; i++){			
		
		pM = (mesh[i][J].U.C[0] + mesh[i][J].U.C[0])/2.0;
		file2<<mesh[i][J].x<<" "<<pM<<endl;
	}

	file1.close();					
	file2.close();
}

void Grid::find_center()
{
	double Tol = 0.0001, x, y, uP, uE, uSE, uS, vP, vE, vSE, vS, uC1, uC2, vC1, vC2;
	int N = 0, uFlag = 0, vFlag = 0;

	for(int i = 2; i < domain.Imax; i++){
		for(int j = 2; j < domain.Jmax; j++){
			
			x = (mesh[i][j].x + mesh[i+1][j].x + mesh[i][j-1].x + mesh[i+1][j-1].x)/4.0;
			y = (mesh[i][j].y + mesh[i+1][j].y + mesh[i][j-1].y + mesh[i+1][j-1].y)/4.0;
			
			uP = mesh[i][j].U.C[1];
			uE = mesh[i+1][j].U.C[1];
			uSE = mesh[i+1][j-1].U.C[1];
			uS = mesh[i][j-1].U.C[1];
			
			vP = mesh[i][j].U.C[2];
			vE = mesh[i+1][j].U.C[2];
			vSE = mesh[i+1][j-1].U.C[2];
			vS = mesh[i][j-1].U.C[2];
						
			uC1 = (uP + uSE)/2.0;
			uC2 = (uE + uS)/2.0;
			
			vC1 = (vP + vSE)/2.0;
			vC2 = (vE + vS)/2.0;
			
			if(fabs(uC1) <= Tol)
				if(fabs(uC2) <= Tol)
					uFlag = 1;
			
			if(fabs(vC1) <= Tol)
				if(fabs(vC2) <= Tol)
					if(uFlag == 1)
						vFlag = 1;

			if(uFlag == 1 && vFlag == 1){				
				cout<<"\n Vortex number: "<<N<<", centered at: "<<x<<", "<<y;
				N++;
			}
			uFlag = vFlag = 0;
		}
	}
}

void Grid::calc_PSI()
{
	double sN = 0.0, sN1 = 0.0;
	
	for(int j = domain.Imin; j <= domain.Jmax; j++){
		for(int i = domain.Imin; i <= domain.Imax; i++){
			
			// Stream Function, psi = int{v*dx
			sN1 += domain.dx*mesh[i][j].U.C[2];
			if(i > domain.Imin)
				sN += domain.dx*mesh[i-1][j].U.C[2];
			
			mesh[i][j].psi = (sN1 + sN)/2.0; 
		}
		sN = sN1 = 0.0;
	}

	// for(int j = domain.Imin; j <= domain.Jmax; j++){
	// 	for(int i = domain.Imin; i <= domain.Imax; i++){
			
	// 		// Stream Function, psi = int{v*dx
	// 		sN1 += domain.dx*mesh[i][j].U.C[2];
	// 		if(i > domain.Imin)
	// 			sN += domain.dx*mesh[i-1][j].U.C[2];
			
	// 		mesh[i][j].psi = (sN1 + sN)/2.0; 
	// 	}
	// 	sN = sN1 = 0.0;
	// }

}

void Grid::calc_OMEGA()
{
	double dx, dy, uN, uS, vE, vW;
	
	for(int j = domain.Imin; j <= domain.Jmax; j++){
		for(int i = domain.Imin; i <= domain.Imax; i++){
			
			dx = domain.dx;
			dy = domain.dy;
			
			uN = mesh[i][j+1].U.C[1];
			uS = mesh[i][j-1].U.C[1];
			
			vE = mesh[i+1][j].U.C[2];
			vW = mesh[i-1][j].U.C[2];
			
			// vorticity function						
			mesh[i][j].omega = fabs((vE-vW)/dx - (uN-uS)/dy); 
		}
	}
}
