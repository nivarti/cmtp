#include "header.h"

void setup(Grid& grid)
{
	grid.calc_Q();
}

void tune(Grid& grid)
{



}

void solve(Grid& grid)
{
	int n = 0;
	double dt = 0.05;       
	Field dU, dUtol;
	dUtol = 0.00000001;
	
	do{
		//dU = grid.march_EE(dt);
		dU = grid.march_IE(dt);
		n++;
		
	}while(dU >= dUtol);

	cout<<"\n Solution converged to "<<dUtol.C[0]<<" in "<<n<<" steps";
}

void verify(Grid& grid)
{
	Field L2norm;
	L2norm = grid.ver_FI();
	cout<<"\nL2 norm = "<<L2norm;
	//tab_L2N(Nx, Ny, L2norm);

}
