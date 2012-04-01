///////////////////////////////////////////////////////
// This file contains thomas algorithm function	     //
// and side functions that are needed to mesh thomas //
// into the rest of the code			     //
///////////////////////////////////////////////////////

#include "header.h"

void SolveThomas(double LHS[NMAX][3], double RHS[NMAX],
		 const int iSize)
{
  int i;
  LHS[0][0] = LHS[iSize+1][2] = 0;
  /* Forward elimination */
  for (i = 0; i < iSize+1; i++) {
    LHS[i][2] /= LHS[i][1];
    RHS[i] /= LHS[i][1];
    LHS[i+1][1] -= LHS[i][2]*LHS[i+1][0];
    RHS[i+1] -= LHS[i+1][0]*RHS[i];
  }
  /* Last line of elimination */
  RHS[iSize+1] /= LHS[iSize+1][1];

  /* Back-substitution */
  for (i = iSize; i >= 0; i--) {
    RHS[i] -= RHS[i+1]*LHS[i][2];
  }
}

void CopyToLHS(double** Dx, double LHS [NMAX][3], const int Size){

  for (int i = 0; i < NMAX; i++) {
    
    LHS[i][0] = 0.0;
    LHS[i][1] = 0.0;
    LHS[i][2] = 0.0;

  }

  for(int i = 0; i <= Size; i++){
    
    LHS[i][0] = Dx[i][0];
    LHS[i][1] = Dx[i][1];
    LHS[i][2] = Dx[i][2];

  }

}

void CopyToRHS(double** FI, double RHS[NMAX], const int Size, const int J, Direction RC){

  for (int i = 0; i < NMAX; i++) {
    
    RHS[i] = 0.0;
    
  }
    
  if(RC == Column)
    for(int i = 0; i <= Size; i++){
      
      RHS[i] = FI[i][J];
      RHS[i] = FI[i][J];
      RHS[i] = FI[i][J];
      
    }
  
  else if(RC == Row)
    for(int i = 0; i <= Size; i++){
      
      RHS[i] = FI[J][i];
      RHS[i] = FI[J][i];
      RHS[i] = FI[J][i];
      
    }
}

void CopyFromRHS(double** FI, double RHS[NMAX], const int Size, const int I, Direction RC){

  if(RC == Row)
    for (int i = 0; i <= Size; i++) {
      
      FI[I][i] = RHS[i];
      
    }
  else if(RC == Column)
    for (int i = 0; i <= Size; i++) {
      
      FI[i][I] = RHS[i];
      
    }

}
