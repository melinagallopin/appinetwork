#include "../sources/tfit.h"
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

// Call without aguments to have a small help
int main(int argc, char** argv)
{
  if (argc <= 1)
  {
    printf("Usage: ./test.exe verticesPerCluster nbOfClusters [density [noise]]\n");
    printf("'density' and 'noise' are the respective probabilities for an ");
    printf("edge to appear in or outside a cluster (default 1, 0)\n");
    return 0;
  }

  int nbPerCluster = atoi(argv[1]);
  int K = atoi(argv[2]);
  double density = (argc >= 4 ? atof(argv[3]) : 1.);
  double noise = (argc >= 5 ? atof(argv[4]) : 0.);
  int N = nbPerCluster * K;

  // Prepare inputs:
  srand(time(NULL));
  int* P = malloc(N*sizeof(int));
  int** A = malloc(N*sizeof(int*));
  for (int i=0; i<N; i++)
    A[i] = calloc(N, sizeof(int));
  for (int k=0; k<K; k++)
  {
    // Cluster k:
    int baseIdx = k * nbPerCluster;
    for (int i=0; i<nbPerCluster; i++)
    {
      // i-th vertex of cluster k (fill row and column symmetrically)
      for (int j=0; j<N; j++)
      {
        if (j == baseIdx + i)
          continue; //zeroes on the diagonal (no self-link)
        double threshold = noise;
        if (baseIdx <= j && j < baseIdx + nbPerCluster)
          threshold = density;
        if ((double)rand() / RAND_MAX < threshold)
        {
          A[baseIdx+i][j] = 1;
          A[j][baseIdx+i] = 1;
        }
      }
    }
  }
  if (N <= 100)
  {
    printf("Graph:\n");
    for (int i=0; i<N; i++)
    {
      for (int j=0; j<N; j++)
        printf("%i ", A[i][j]);
      printf("\n");
    }
    printf("\n");
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

  // Release memory
  free(P);
  for (int i=0; i<N; i++)
    free(A[i]);
  free(A);

  return 0;
}
