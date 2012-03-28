#define NMAX 100

static void SolveThomas(double LHS[NMAX][3], double RHS[NMAX],
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

int main()
{
  double LHS[12][3], RHS[12];
  int i;
  
  for (i = 0; i <= 11; i++) {
    LHS[i][0] = 1;
    LHS[i][1] = 2+i;
    LHS[i][2] = 3;
    RHS[i]   = i;
  }

  SolveThomas(LHS, RHS, 10);
      
  for (i = 0; i <= 11; i++) 
    printf("%f\n", RHS[i]);

  return(0);
}
