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
	Field L2norm;
	
	// Start run-time calculation
	ti = clock();
	
	while(Nx <= 10)
	{
		Rectangle cavity(Nx, Ny, 1);
		Grid nse(cavity);				
		
		setup(nse);
		solve(nse);
		
		Nx *= 2;
		Ny *= 2;
	};
	
	// End run-time calculation
	tf = clock();

	cout<<"\n\nCode run time = "<<(tf - ti)/CLOCKS_PER_SEC*1000<<" ms\n\n";
	return 0;
}
