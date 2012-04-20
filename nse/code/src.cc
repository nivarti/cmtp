#include "header.h"

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

// void plot_l2norm(Field F, double dx)
// {
// 	static int I = 1;
// 	ostringstream counter;
// 	string name;
// 	ofstream file;

// 	counter<<"./norm/mesh"<<I;
// 	name = counter.str();	
// 	file.open(name.c_str(), ios::app);

// 	file<<dx<<" "<<F<<endl;	
	
// 	file.close();
// }

void tab_L2N(int Nx, int Ny, Field L2N)
{
	static int I = 1;
	ofstream file;
	file.open("./norm/L2N", ios::app);

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

void tab_EJ(int i, int j, Field E)
{	
	static int I = 1;
	ofstream file;
	file.open("./jacobian/EJ", ios::app);

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
