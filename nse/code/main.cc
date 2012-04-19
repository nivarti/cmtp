/**
 * \brief 
 *
 *
 */


#include "header.h"

int main()
{
	int Nx = 20, Ny = 20;
	clock_t ti, tf;
	Field norm;

	ti = clock();
	while(Nx <= 20)
	{
		Rectangle cavity(Nx, Ny, 1);
		Grid nse(cavity);
		
		nse.calc_Q();
		//norm = nse.ver_FI;
		//plot_l2norm(norm, Lx/Nx);
		

		cout<<"\nAll correct";
		//solve_nse(cavity);
		Nx *= 2;
		Ny *= 2;

	}
	tf = clock();

	cout<<"\nCode run time = "<<(tf - ti)/CLOCKS_PER_SEC*1000<<" ms\n";
	return 0;
}
