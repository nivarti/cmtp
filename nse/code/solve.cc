#include "header.h"


void setup(Grid& grid)
{	
	grid.calc_Q();
}

void tune(Grid& grid, double T)
{
	int n = 0;
	double dt = 0.05, SOR = 1.0;
	Field dU, dUtol;
	ofstream file;	
	
	dUtol = 0.000001;
	
	//SOR = T;
	//dt = T;

	do{
		// March in time with implicit euler
		dU = grid.march_IE(dt, SOR);
		n++;		
		
		//cout<<"\n Maximum change in solution at the end of "<<n<<" steps = "<<dU;
				
		//plot_CH(n, dU, "../plot/relaxation/dU6");
		
	}while(dU >= dUtol);
	
	cout<<"\n Solution with dt: " <<dt<<" converged to "<<dU<<" in "<<n<<" steps";
	//cout<<"\n\n Solution with SOR: " <<SOR<<" converged to "<<dU<<" in "<<n<<" steps";				
	
	//grid.slice_U(11, 11);
	//grid.spew_field(Pressure, 11, -1);
	//grid.spew_field(Pressure, -1, 11);

	// for(int i = grid.domain.Imin; i <= grid.domain.Imax; i++){
	// 	grid.spew_field(Pressure, i, -1);
	// 	cout<<endl;
	// }		
	// for(int j = grid.domain.Jmin; j <= grid.domain.Jmax; j++){
		
	// 	grid.spew_field(Pressure, -1, j);
	// 	cout<<endl;
	// }
	
	
	//cout<<"\n Pressure solution:\n";
	//grid.spew_field(Pressure);
	
	//cout<<"\n Pressure solution:\n";
	//grid.spew_field(Pressure);
	
	//cout<<"\n Pressure solution:\n";
	//grid.spew_field(Pressure);
	
	//file.close();		
}

void solve(Grid& grid)
{
	int n = 0;
	double dt = 0.05, SOR = 1.0;
	Field dU, dUtol;
	
	// Specify tolerance for convergence
	dUtol = 0.000001;
	
	do{
		// March in time with implicit euler
		dU = grid.march_IE(dt, SOR);
		n++;
		
		cout<<"\n Maximum change in solution at the end of "<<n<<" steps = "<<dU;
		
		//Plot L2 norm, average field, etc for every step
		//plot_CH(n, dU, "../plot/oscillation/dU");
		//plot_CH(n, grid.calc_Uav(), "../plot/convergence/Uav");
		//if(n <= 200)
		
		//plot_CH(n, dU, "../plot/vortex/conv");

	}while(dU >= dUtol);
	
	cout<<"\n Solution with SOR: " <<SOR<<" converged to "<<dU<<" in "<<n<<" steps";
	//cout<<"\n Solution with beta: " <<beta<<" converged to "<<dU<<" in "<<n<<" steps";
	// dU = grid.calc_Uav();
	// cout<<"\n Average fields: "<<dU;
	
	//grid.mirror_U();
	//grid.find_center();
	//grid.calc_PSI();
	//grid.plot_U();
	
	//cout<<"\n Pressure solution:\n";
	//grid.spew_field(xVelocity);
	//cout<<endl<<endl;
	//cout<<"\n Pressure solution:\n";
	//grid.spew_field(yVelocity);
	//grid.plot_uSL();
	grid.plot_uSL("../plot/basic/gci/mesh");
}

void verify(Grid& grid)
{
	Field L2norm;
	
	L2norm = grid.ver_FI();
	tab_L2N(grid.domain.Imax, grid.domain.Jmax, L2norm);
	
	cout<<"\nL2 norm for solution = "<<L2norm;
	
	grid.calc_EJ();
}
