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

Field Field::operator/(const double &RHS)
{
	Field temp;

	temp.C[0] = this->C[0]/RHS;
	temp.C[1] = this->C[1]/RHS;
	temp.C[2] = this->C[2]/RHS;

	return (temp);
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

Field Field::operator*(const Field& F)
{
	Field Result;
	
	Result.C[0] = F.C[0]*this->C[0]; 
	Result.C[1] = F.C[1]*this->C[1]; 
	Result.C[2] = F.C[2]*this->C[2]; 

	return Result;	
}

Field Field::operator*(const double A[][3])
{
	Field Result;
	
	Result.C[0] = A[0][0]*C[0] + A[0][1]*C[1] + A[0][2]*C[2]; 
	Result.C[1] = A[1][0]*C[0] + A[1][1]*C[1] + A[1][2]*C[2]; 
	Result.C[2] = A[2][0]*C[0] + A[2][1]*C[1] + A[2][2]*C[2]; 

	return Result;	
}

bool Field::operator>=(const Field &RHS)
{
	double F1, F2;
	
	F1 = max(*this);
	F2 = max(RHS);

	if(F1 > F2)
		return true;
	else
		return false;
}

ostream& operator<<(ostream &out, Field &f)
{
	// out<<"\nU1 =  "<<f.C[0];
	// out<<"\nU2 =  "<<f.C[1];
	// out<<"\nU3 =  "<<f.C[2];
	//out<<setprecision(5)<<"U1: "<<setw(10)<<f.C[0]<<"  U2: "<<setw(10)<<f.C[1]<<"  U3: "<<setw(10)<<f.C[2];
	out<<setprecision(15)<<f.C[0]<<"  "<<f.C[1]<<"  "<<f.C[2];
	return out;
}
