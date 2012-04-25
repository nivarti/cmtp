/**
 * \brief 
 *
 *
 */


#include "header.h"

int main()
{
	int Nx = 80, Ny = 240;
	clock_t ti, tf;
	
	// tuning parameter for nse
	double T = 0.05;
	
	// Start run-time calculation
	ti = clock();
	//est_GCI();

	cout<<"\nCaution! Solving the Navier Stokes Equations...";
	//cout<<"Paramters used:\n Mesh = ";
	
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
