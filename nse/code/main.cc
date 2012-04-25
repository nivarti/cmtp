/**
 * \brief 
 *
 *
 */


#include "header.h"

int main()
{
	int Nx = 10, Ny = 30;
	clock_t ti, tf;
	
	// tuning parameter for nse
	double T = 0.05;
	
	// Start run-time calculation
	ti = clock();
	//est_GCI();
	
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
		
	}while(Nx <= 10);
	
	// End run-time calculation
	tf = clock();	
	cout<<"\n\nCode run time = "<<(tf - ti)/CLOCKS_PER_SEC*1000<<" ms\n\n";
	
	return 0;
}
