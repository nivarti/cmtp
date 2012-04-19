#include "header.h"

void plot_l2norm(Field F, double dx)
{
	static int I = 1;
	ostringstream counter;
	string name;
	ofstream file;

	counter<<"./norm/mesh"<<I;
	name = counter.str();	
	file.open(name.c_str(), ios::app);

	file<<dx<<" "<<F<<endl;	
	
	file.close();
}
