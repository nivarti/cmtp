#include "header.h"


Grid::Grid(Rectangle rect)
{
	double x, y;

	domain = rect;	
	cout<<"\nInitialising mesh of size: "<<domain.Imax<<" x "<<domain.Jmax;

	for (int j = 0; j <= domain.Jmax + domain.Jmin; j++){
		for (int i = 0; i <= domain.Imax + domain.Imin; i++){
			
			x = domain.dx*(2.0*(double)i - (double)domain.Imin)/2.0;
			y = domain.dy*(2.0*(double)j - (double)domain.Jmin)/2.0;
			
			mesh[i][j].set_XY(x, y);
			mesh[i][j].set_eQ();								
		}
	}    
	cout<<"\nEvaluated cell coordinates...";
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
	
	mesh[i][j].FI = Fx;
	mesh[i][j].FI += Gy;
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
	
	if (i == 10 && j == 10){
		cout<<"\n U: " <<setprecision(20)<< uP <<"  "<< uN <<"  "<< uE <<"  "<< uW <<"  "<< uS;
		cout<<"\n V: " << vP <<"  "<< vN <<"  "<< vE <<"  "<< vW <<"  "<< vS;		
	}

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

void Grid::calc_BC(){

	
	for(int j = domain.Jmin; j <= domain.Jmax; j++){
		
		// Neumann on Pressure
		mesh[0][j].U.C[0] = mesh[1][j].U.C[0];   		
		mesh[domain.Imin + domain.Imax][j].U.C[0] = mesh[domain.Imax][j].U.C[0];
		
		// No Slip on velocities
		mesh[0][j].U.C[1] = -mesh[1][j].U.C[1];   		
		mesh[domain.Imin + domain.Imax][j].U.C[1] = -mesh[domain.Imax][j].U.C[1];
		mesh[0][j].U.C[2] = -mesh[1][j].U.C[2];   		
		mesh[domain.Imin + domain.Imax][j].U.C[2] = -mesh[domain.Imax][j].U.C[2];

	}  
  
	for(int i = domain.Imin; i <= domain.Imax; i++){
		
		// Neumann on Pressure
		mesh[i][0].U.C[0] = mesh[i][1].U.C[0];   		
		mesh[i][domain.Jmin + domain.Jmax].U.C[0] = mesh[i][domain.Jmax].U.C[0];

		// No Slip on velocities
		mesh[i][0].U.C[1] = -mesh[i][1].U.C[1];   		
		mesh[i][domain.Jmin + domain.Jmax].U.C[1] = -mesh[i][domain.Jmax].U.C[1];

		mesh[i][0].U.C[2] = -mesh[i][1].U.C[2];   		
		mesh[i][domain.Jmin + domain.Jmax].U.C[2] = -mesh[i][domain.Jmax].U.C[2];

  }  
}

void Grid::calc_IBC()
{
	for(int j = domain.Jmin; j <= domain.Jmax; j++){
		
		// Neumann on Pressure, No Slip of velocities
		mesh[0][j].Cx[0][0] = -1.0;
		mesh[0][j].Cx[1][1] = 1.0;
		mesh[0][j].Cx[2][2] = 1.0;		
		
		mesh[domain.Imax + domain.Imin][j].Ax[0][0] = -1.0;
		mesh[domain.Imax + domain.Imin][j].Ax[1][1] = 1.0;
		mesh[domain.Imax + domain.Imin][j].Ax[2][2] = 1.0;
		
	}

	for(int i = domain.Imin; i <= domain.Imax; i++){
		
		// Neumann on Pressure, No Slip of velocities
		mesh[i][0].Cy[0][0] = -1.0;
		mesh[i][0].Cy[1][1] = 1.0;
		mesh[i][0].Cy[2][2] = 1.0;		
		
		mesh[i][domain.Jmax + domain.Jmin].Ay[0][0] = -1.0;
		mesh[i][domain.Jmax + domain.Jmin].Ay[1][1] = 1.0;
		mesh[i][domain.Jmax + domain.Jmin].Ay[2][2] = 1.0;
	}
}

void Grid::calc_Q()
{
	
	for(int j = domain.Jmin; j <= domain.Jmax; j++){
		for(int i = domain.Imin; i <= domain.Imax; i++){
						
			calc_FI(i, j);
			calc_J(i, j);	

			if(i == 10 && j == 10){
				cout<<"\n Ax:";
				spew_matrix(mesh[i][j].Ax);
				
				cout<<"\n Ay:";
				spew_matrix(mesh[i][j].Ay);
				
				cout<<"\n Bx:";
				spew_matrix(mesh[i][j].Bx);
				
				cout<<"\n By:";
				spew_matrix(mesh[i][j].By);
				
				cout<<"\n Cx:";
				spew_matrix(mesh[i][j].Cx);
				
				cout<<"\n Cy:";
				spew_matrix(mesh[i][j].Cy);
				
			}
		}
	}
		
	ver_FI();	
	calc_EJ();
}

void Grid::calc_EJ()
{		
	Field LHS, RHS, E;
	mesh[10][10].dU = 0.000001;		
	
	for(int j = domain.Jmin; j <= domain.Jmax; j++){
		for(int i = domain.Imin; i <= domain.Imax; i++){			
			
			mesh[i][j].FIn = mesh[i][j].FI;
			mesh[i][j].Un = mesh[i][j].U;
			mesh[i][j].U += mesh[i][j].dU;
			
			if (i==10 && j ==10){
				
				setprecision(30);
				setw(40);
				cout<<"\nInformation from cell: "<<i<<", "<<j;
				cout<<"\n Coordinates: "<<mesh[i][j].x<<", "<<mesh[i][j].y;
				cout<<"\n\n Value of eFI: "<<mesh[i][j].eFI;
				cout<<"\n Value of FI(n): "<<mesh[i][j].FIn;				
				
				cout<<"\n\n Exact Solution: "<<mesh[i][j].eU;
				cout<<"\n Numerical Solution: "<<mesh[i][j].Un;
				cout<<"\n Change in Solution: "<<mesh[i][j].dU;
				cout<<"\n New Solution: "<<mesh[i][j].U;				
				
			}
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
			
			if (i>=9 && i <= 11 && j >= 9 && j <= 11){
			//if (i == 10 && j == 10){				
				//cout<<"\n Value of FI(n + 1): "<<mesh[i][i].FI;				
				cout<<"\n\n Information from cell: "<<i<<", "<<j;
				cout<<"\n LHS: "<<LHS;
				cout<<"\n RHS: "<<RHS;
				cout<<"\n Error  = LHS - RHS = "<<E;
				tab_EJ(i, j, E);
			}
		}
	}			
}

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

	l2norm.sqroot();
	l2norm /= sqrt(domain.Imax*domain.Jmax);

	cout<<"\n\nL2 Norm for flux integral:"<<l2norm;
	return l2norm;
}

// double Grid::march_IE(double dt)
// {
	
// 	Field dU = dUmax = 0.0;
// 	double RHS[MAXSIZE], LHS[MAXSIZE][3][3][3];
	
// 	calc_BC();
// 	calc_FI();
	
// 	// Start approximate factorisation solution
// 	for(int j = Jmin; j <= Jmax; j++){
		
// 		LHS = calc_LHS(Column, j);
// 		RHS = calc_RHS(Column, j);
		
// 		SolveBlockTri(LHS, RHS, N);	
		
// 	}
	
// 	for(int j = Jmin; j <= Jmax; j++){
		
// 		LHS = calc_LHS();
// 		RHS = calc_RHS(FIn);
		
// 		SolveBlockTri(LHS, RHS, N);
// 	}
// 	// End approximate factorisation solution
			
//   // 		CopyToLHS(Dx, LHS, Imax + 1);
		
//   // 		CopyToRHS(FI, RHS, Imax + 1, j, Column);
		
//   // 		//   if (i == 1){
//   // 		//     for(int i = Imin; i <= Imax; i++){
		
//   // 		// 	cout<<RHS[i]<<"\n";
	
//   // 		//     }
//   // 		//     cout<<endl;
//   // 		// }
		
//   // 		SolveBlockTri(LHS, RHS, Imax);
//   // 		CopyFromRHS(FI, RHS, Imax + 1, j, Column);
		
//   // 	}
	
//   // 	for(int i = Imin; i <= Imax; i++){
		
//   //   EvaluateDy(i);  
//   //   CopyToLHS(Dy, LHS, Jmax + 1);
    
    
//   //   CopyToRHS(FI, RHS, Jmax + 1, i, Row);
//   //   // if (i == 1){
//   //   //   cout<<"\nbringing up dT tildas for i = 1\n";
    
//   //   //   for(int j = Jmin; j <= Jmax; j++){	
//   //   // 	cout<<RHS[j]<<"\n";	
//   //   //   }
//   //   //   cout<<endl;
//   //   // }
    
//   //   SolveThomas(LHS, RHS, Jmax);  
//   //   CopyFromRHS(FI, RHS, Jmax + 1, i, Row);
    
//   // }  
//   // End of Approximate factorisation solution
	
// 	// Update solution
// 	for(int i = domain.Imin; i <= domain.Imax; i++){
// 		for(int j = domain.Jmin; j <= domain.Jmax; j++){      
			
// 			dU = mesh[i][j].dU = FI[i][j];
// 			mesh[i][j].U += dU;
			
// 			dU.abs();
// 			if(dU >= dUmax)
// 				dUmax = dU;
// 		}
// 	}
	
// 	return dUmax;	
// }
	
// void Grid::calc_LHS(Direction RC, int IJ)
// {
// 	double LHS[MAXSIZE][3][3][3];
	
// 	if(RC == Column)
// 		{
// 			for(int k = 0; k <= domain.Imax + domain.Imin; k++)
// 				{
// 					LHS[k][0] = mesh[k][IJ].Ay;
// 					LHS[k][1] = mesh[k][IJ].IDy;
// 					LHS[k][2] = mesh[k][IJ].Cy;
// 				}
// 		}

	

// 	if(RC == Row)
// 		{
// 			for(int k = 0; k <= domain.Jmax + domain.Jmin; k++)
// 				{
// 					LHS[k][0] = mesh[IJ][k].Ax;
// 					LHS[k][1] = mesh[IJ][k].IDx;
// 					LHS[k][2] = mesh[IJ][k].Cx;
// 				}
// 		}
	
// 	return LHS;
// }

// void Grid::calc_RHS(Direction RC, int IJ)
// {
// 	double RHS[MAXSIZE][3];
	
// 	if(RC == Column)
// 		{
// 			for(int k = 0; k <= domain.Imax + domain.Imin; k++)
// 				{
// 					RHS[k][0] = mesh[k][IJ].FI.C[0];
// 					RHS[k][1] = mesh[k][IJ].FI.C[1];
// 					RHS[k][2] = mesh[k][IJ].FI.C[2];
// 				}
// 		}
	
// 	if(RC == Row)
// 		{
// 			for(int k = 0; k <= domain.Jmax + domain.Jmin; k++)
// 				{
// 					RHS[k][0] = mesh[IJ][k].FI.C[0];
// 					RHS[k][1] = mesh[IJ][k].FI.C[1];
// 					RHS[k][2] = mesh[IJ][k].FI.C[2];
// 				}
// 		}

// 	return RHS;
// }

// void Grid::solve_T(LHS, RHS)
// {



// }
