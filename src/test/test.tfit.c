#include "../sources/tfit.h"
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

int main(int argc, char** argv)
{
  // Prepare inputs:
  srand(time(NULL));
  int N = atoi(argv[1]); //number of vertices
  double alpha = atof(argv[2]); //density parameter
  int* P = malloc(N*sizeof(int));
  int** A = malloc(N*sizeof(int*));
  for (int i=0; i<N; i++)
  {
    A[i] = calloc(N, sizeof(int));
    for (int j=0; j<i; j++)
    {
      if ((double)rand() / RAND_MAX < alpha)
      {
        A[i][j] = 1;
        A[j][i] = 1;
      }
    }
  }

  // Main call
  tfit_core(A, N, P);

  // Show result:
  printf("Partition:\n");
  for (int i=0; i<N; i++)
  {
    printf("%i", P[i]);
    if (i < N-1)
      printf(" ");
  }
  printf("\n");

  return 0;
}
