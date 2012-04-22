#include "header.h"

void setup(Grid& grid)
{
	grid.calc_Q();
}

void tune(Grid& grid)
{
	//grid.add_FI();
	grid.add_J();

}

void solve(Grid& grid)
{
	int n = 0;
	double dt = 0.05;
	Field dU, dUtol;
	dUtol = 0.000001;
		
	do{
		dU = grid.march_IE(dt);
		n++;

		cout<<"\n Maximum change in solution at the end of "<<n<<" steps = "<<dU;
		
	}while(n < 600);
	
	cout<<"\n\n Solution converged to "<<max(dU)<<" in "<<n<<" steps\n\n";

	grid.spew_field(Pressure);cout<<endl<<endl;
	grid.spew_field(xVelocity);cout<<endl<<endl;
	grid.spew_field(yVelocity);
}

void verify(Grid& grid)
{
	Field L2norm;
	L2norm = grid.ver_FI();
	cout<<"\nL2 norm = "<<L2norm;
	//tab_L2N(Nx, Ny, L2norm);

}
