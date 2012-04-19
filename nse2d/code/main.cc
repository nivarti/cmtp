/**
 * \brief 
 *
 *
 */


#include "header.h"

int main()
{
	int Nx = 10, Ny = 10;
	clock_t ti, tf;
	Field norm;

	ti = clock();
	while(Nx <= 20)
	{
		Rectangle cavity(Nx, Ny, 1);
		Grid nse2d(cavity);
		
		nse2d.calc_fi();
		norm = nse2d.verif_fi();
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
