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
	//double T = 0.05;  // tunable parameter for nse
		
	// Start run-time calculation
	ti = clock();
	
	do{
		// Initialise geometry, and mesh
		Rectangle cavity(Nx, Ny, 1);
		Grid nse(cavity);												
		//setup(nse);
		//tune(nse, T);
		
		solve(nse);
		//verify(nse);
		
		// Increase mesh size
		Nx *= 2;
		Ny *= 2;

		//T += 0.05;
		
	}while(Nx <= 20);	
	// End run-time calculation
	tf = clock();	
	cout<<"\n\nCode run time = "<<(tf - ti)/CLOCKS_PER_SEC*1000<<" ms\n\n";
	
	return 0;
}
