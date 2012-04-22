#include "header.h"

Field sqrt(Field F)
{
	Field sF;
	
	sF.C[0] = sqrt(F.C[0]);
	sF.C[1] = sqrt(F.C[1]);
	sF.C[2] = sqrt(F.C[2]);
	
	return sF;	
}

Field fabs(Field F)
{
	Field aF;
	
	aF.C[0] = fabs(F.C[0]);
	aF.C[1] = fabs(F.C[2]);
	aF.C[2] = fabs(F.C[2]);
	
	return aF;	
}

double max(Field F)
{
	double Fmax = F.C[0];
	
	if(F.C[1] > Fmax)
		Fmax = F.C[1];
	if(F.C[2] > Fmax)
		Fmax = F.C[2];
	
	return Fmax;
}

void spew_matrix(double M[][3])
{	
	for(int i = 0; i < 3; i++){
		cout<<endl;
		for(int j = 0; j < 3; j++){
			
			cout<<setw(10)<<M[i][j]<<" ";
		}
	}			
}

void mult_matrix(double M[][3], double K)
{
	for(int i = 0; i < 3; i++){
		for(int j = 0; j < 3; j++){
			
			M[i][j] *= K;
		}
	}			
}

void copy_matrix(double S[][3], double T[][3])
{
	for(int i = 0; i < 3; i++){
		for(int j = 0; j < 3; j++){
			
			T[i][j] = S[i][j];
		}
	}			
}

void init_matrix(double M[][3], double K)
{
	for(int i = 0; i < 3; i++){
		for(int j = 0; j < 3; j++){
			
			if(i == j)
				M[i][j] = K;
			else
				M[i][j] = 0.0;
		}
	}			
}

// Tabulate change in l2 norms
void tab_L2N(int Nx, int Ny, Field L2N)
{
	static int I = 1;
	ofstream file;
	file.open("./norm/L2N.tex", ios::app);

	if(I == 1){
		file<<"\n\\begin{table}\n\\begin{center}\n\\begin{tabular}{|l | r | r | r |}\n\\hline";				
		file<<"\nMesh Size & Error (U1) & Error (U2) & Error (U3) \\\\\n\\hline";
	}
	
	file<<"\n"<<Nx<<", "<<Ny<<" & "<<L2N.C[0]<<" & "<<L2N.C[1]<<" & "<<L2N.C[2]<<" \\\\";	
	
	if(I == 2){
		file<<"\n\\hline\n\\end{tabular}\n\\caption{Tabulation of $L_2$ Norms for different mesh sizes.}";
		file<<"\n\\label{approxE}\n\\end{center}\n\\end{table}";	
	}
	I++;
	file.close();	
}

// Tabulate errors in Jacobian
void tab_EJ(int i, int j, Field E)
{	
	static int I = 1;
	ofstream file;
	file.open("./jacobian/EJ.tex", ios::app);

	if(I == 1){
		file<<"\n\\begin{table}\n\\begin{center}\n\\begin{tabular}{|l | r | r | r |}\n\\hline";				
		file<<"\nCell Index & Error (U1) & Error (U2) & Error (U3) \\\\\n\\hline";
	}
	
	file<<"\n"<<i<<", "<<j<<" & "<<E.C[0]<<" & "<<E.C[1]<<" & "<<E.C[2]<<" \\\\";	
	
	if(I == 9){
		file<<"\n\\hline\n\\end{tabular}\n\\caption{Approximation error in change in flux integral with time step.}";
		file<<"\n\\label{approxE}\n\\end{center}\n\\end{table}";	
	}
	I++;
	file.close();
}

// Mean pressure
void plot_avP(int n, Field avP)
{
	ofstream file;
	file.open("../plot/oscillation/avP", ios::app);
	
	for(int i = domain.Imin; i <= domain.Imax; i++){
		for(int j = domain.Jmin; j <= domain.Jmax; j++){
			
	file<<n<<" "<<avP<<endl;
	
	file.close();
}



// Plot convergence history
void plot_CH(int n, Field L2N)
{
	ofstream file;
	file.open("../plot/stability/conv", ios::app);
	
	file<<n<<" "<<L2N<<endl;
	
	file.close();
}

// Symmetry check by comparison of solutions
void plot_U(Grid& grid1, Grid& grid2)
{
	ofstream file;
	double x, y, d;
	
	file.open("../plot/symmetry/uvel");
	
	for(int i = grid1.domain.Imin; i <= grid1.domain.Imax; i++){
		for(int j = grid1.domain.Jmin; j <= grid1.domain.Jmax; j++){
			
			x = grid1.mesh[i][j].x;
			y = grid2.mesh[i][j].y;
			d = grid1.mesh[i][j].U.C[1] - grid1.mesh[grid1.domain.Imax - i][j].U.C[1];

			file<<x<<" "<<y<<" "<<d<<endl;
		}
	}
		
	file.close();				
}
