#include "header.h"


Field::Field()
{
	C[0] = C[1] = C[2] = 0.0;
}

Field& Field::operator=(const double &RHS)
{
	C[0] = RHS;
	C[1] = RHS;
	C[2] = RHS;

	return *this;
}

Field& Field::operator=(const Field &RHS)
{
	C[0] = RHS.C[0];
	C[1] = RHS.C[1];
	C[2] = RHS.C[2];

	return *this;
}

Field& Field::operator+=(const Field &RHS)
{
	C[0] += RHS.C[0];
	C[1] += RHS.C[1];
	C[2] += RHS.C[2];

	return *this;
}

Field Field::operator-(const Field &RHS)
{
	Field temp;

	temp.C[0] = this->C[0] - RHS.C[0];
	temp.C[1] = this->C[1] - RHS.C[1];
	temp.C[2] = this->C[2] - RHS.C[2];

	return (temp);
}

Field& Field::operator/=(const double &RHS)
{
	C[0] /= RHS;
	C[1] /= RHS;
	C[2] /= RHS;

	return *this;
}

Field& Field::operator/=(const Field &RHS)
{
	C[0] /= RHS.C[0];
	C[1] /= RHS.C[1];
	C[2] /= RHS.C[2];

	return *this;
}

Field& Field::operator*=(const double &RHS)
{
	C[0] *= RHS;
	C[1] *= RHS;
	C[2] *= RHS;

	return *this;
}

Field& Field::operator*=(const Field &RHS)
{
	C[0] *= RHS.C[0];
	C[1] *= RHS.C[1];
	C[2] *= RHS.C[2];

	return *this;
}


ostream& operator<<(ostream &out, Field &f)
{
	// out<<"\nU1 =  "<<f.C[0];
	// out<<"\nU2 =  "<<f.C[1];
	// out<<"\nU3 =  "<<f.C[2];
	out<<setprecision(10)<<f.C[0]<<" "<<f.C[1]<<" "<<f.C[2];

	return out;
}

void Field::sqroot()
{
	C[0] = sqrt(C[0]);
	C[1] = sqrt(C[1]);
	C[2] = sqrt(C[2]);     
}
