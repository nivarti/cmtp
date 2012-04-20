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
	Field L2norm;

	ti = clock();
	while(Nx <= 40)
	{
		Rectangle cavity(Nx, Ny, 1);
		Grid nse(cavity);
		
		nse.calc_Q();
		L2norm = nse.ver_FI();
		tab_L2N(Nx, Ny, L2norm);		

		
		Nx *= 2;
		Ny *= 2;

	}
	tf = clock();

	cout<<"\n\nCode run time = "<<(tf - ti)/CLOCKS_PER_SEC*1000<<" ms\n\n";
	return 0;
}
