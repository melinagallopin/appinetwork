#include <stdlib.h>
#include <math.h>

// Méthode de Transfert-Fusion itérée
// Lexique:
//   N = nbVertices
//   K = nbClasses
//   A = adjacency matrix
//   B = modularity matrix
//   P = partition

// For debug:
//#include <stdio.h>
//void printMat(double** mat, int n, int m, char* name)
//{
//  printf("print %s\n", name);
//  for (int i=0; i<n; i++)
//  {
//    for (int j=0; j<m; j++)
//      printf("%f ", mat[i][j]);
//    printf("\n");
//  }
//}
//void printPart(int* p, int N, char* msg)
//{
//  printf("partition: [ %s ]\n", msg);
//  for (int i=0; i<N; i++)
//    printf("%i ", p[i]);
//  printf("\n");
//}

// Compute the total score (sum of intra-clusters weights)
double GetScore(double** B, int N, int* P)
{
  double score = 0.;
  for (int i=0; i<N; i++)
  {
    for (int j=0; j<i; j++)
    {
      if (P[i] == P[j])
        score += B[i][j];
    }
  }
  return score;
}

// Compute the modularity matrix B (weighting of edges)
double** GetModularityMatrix(int N, int** A)
{
  // NOTE: in current version B is an integer matrix
  double** B = malloc(N*sizeof(double*));
  for (int i=0; i<N; i++)
    B[i] = malloc (N*sizeof(double));

  int nbEdges = 0; //in fact 2 x number of edges
  int* degrees = malloc(N*sizeof(int));
  for (int i=0; i<N; i++)
  {
    degrees[i] = 0;
    B[i][i] = 0.;
    for (int j=0; j<N; j++)
      degrees[i] += A[i][j];
    nbEdges += degrees[i];
  }

  for (int i=1; i<N; i++)
  {
    for (int j=0; j<i; j++)
    {
      B[i][j] = nbEdges*A[i][j] - degrees[i]*degrees[j];
      B[j][i] = B[i][j];
    }
  }

  free(degrees);
  return B;
}

// Renumbering of clusters, potentially changing K
void Renumber(int N, int* P, int* K)
{
  int* nbElems = calloc(N, sizeof(int)); //number of items per class
  for (int i=0; i<N; i++)
    nbElems[P[i]]++;
  int classIdx = 0;
  for (int k=0; k<N; k++)
  {
    if (nbElems[k] == 0)
      continue;
    for (int i=0; i<N; i++)
    {
      if (P[i] == k)
        P[i] = classIdx; //classIdx <= k, so it's OK
    }
    classIdx++;
  }
  free(nbElems);
  *K = classIdx; //next class index
}

// Iteratively merge clusters while the resulting score increases
int Fusion(double** B, int N, int* P, int K)
{
  // Precompute the inter-classes scores (weights):
  double** W = malloc(K*sizeof(double*));
  for (int k=0; k<K; k++)
    W[k] = calloc(K, sizeof(double));
  for (int i=0; i<N; i++)
  {
    for (int j=0; j<i; j++)
    {
      W[P[i]][P[j]] += B[i][j];
      W[P[j]][P[i]] += B[i][j];
    }
  }
  int* emptyClass = calloc(K, sizeof(int));

  int atLeastOneFusion = 0;
  while (1)
  {
    double maxDeltaScore = -1.;
    int k1 = 0;
    int k2 = 0;
    for (int k=0; k<K; k++)
    {
      if (emptyClass[k])
        continue;
      for (int kk=0; kk<k; kk++)
      {
        if (emptyClass[kk])
          continue;
        // Compute the score variation if classes k and kk are merged
        if (W[k][kk] > maxDeltaScore)
        {
          maxDeltaScore = W[k][kk];
          k1 = k;
          k2 = kk;
        }
      }
    }

    if (maxDeltaScore <= 0.)
      break; //no more fusions to achieve
    atLeastOneFusion = 1;

    // Apply fusion between classes k1 and k2 (discard k2 --> empty)
    for (int k=0; k<K; k++)
    {
      if (k == k1 || k == k2)
        continue;
      W[k][k1] += W[k][k2];
      W[k1][k] += W[k][k2];
    }
    for (int i=0; i<N; i++)
    {
      if (P[i] == k2)
        P[i] = k1;
    }
    emptyClass[k2] = 1;
  }

  for (int k=0; k<K; k++)
    free(W[k]);
  free(W);
  free(emptyClass);

  return atLeastOneFusion;
}

// Iteratively transfer vertices from one cluster to another,
// until the resulting score cannot increase.
int Transfer(double** B, int N, int* P, int K)
{
  // Compute the contributions of each vertex to all classes
  double** contribs = malloc(N*sizeof(double*));
  for (int i=0; i<N; i++)
    contribs[i] = calloc(K, sizeof(double));
  for (int i=0; i<N; i++)
  {
    for (int j=0; j<i; j++)
    {
      contribs[i][P[j]] += B[i][j];
      contribs[j][P[i]] += B[i][j];
    }
  }

  int atLeastOneTransfer = 0;
  while (1)
  {
    // Try to increase score:
    double maxScoreIncrease = -1.;
    int index = 0;
    int cluster = 0;
    for (int i=0; i<N; i++)
    {
      // If the contribution of vertex i to its own class is negative,
      // then it's better out (potentially creating a new class): case k == K
      for (int k=0; k <= K; k++)
      {
        if (k == P[i]) //moving to its own class is no move
          continue;
        double score = (k < K ? contribs[i][k] : 0.) - contribs[i][P[i]];
        if (score > maxScoreIncrease)
        {
          maxScoreIncrease = score;
          index = i;
          cluster = k;
        }
      }
    }

    if (maxScoreIncrease <= 0.)
      break;
    atLeastOneTransfer = 1;

    // Apply transfert:
    for (int i=0; i<N; i++)
    {
      contribs[i][P[index]] -= B[i][index];
      contribs[i][cluster] += B[i][index];
    }
    P[index] = cluster;
    if (cluster == K)
      K++;
  }

  for (int i=0; i<N; i++)
    free(contribs[i]);
  free(contribs);

  return atLeastOneTransfer;
}

// Perturb the vertices assignation (randomly changing classes of some of
// them), in the hope to find a better partition.
void StochasticOptimization(double** B, int N, int* P, int K)
{
  if (K == 1) //only one partition, nothing to do
    return;
  int* bestP = malloc(N*sizeof(int)); //best partition found so far
  for (int i=0; i<N; i++)
    bestP[i] = P[i];
  int bestK = K;
  double bestScore = GetScore(B, N, P);
  int maxNbAttempts = fmin(N, 500);
  for (int attempt = 0; attempt < maxNbAttempts; attempt++)
  {
    // Start from best partition found so far:
    for (int i=0; i<N; i++)
      P[i] = bestP[i];
    K = bestK;
    // TODO: next nbSwaps value is arbitrary ("magic expression").
    // The idea is to decrease nbSwaps over time ("decrease temperature")
    int nbSwaps = floor(pow(N, (maxNbAttempts-attempt)/maxNbAttempts));
    int swapIdx = 0;
    while (swapIdx < nbSwaps)
    {
      // Draw random indexes between 0 and N-1 (included):
      int i = floor(((double)rand() / RAND_MAX) * N);
      int j = floor(((double)rand() / RAND_MAX) * N);
      if (P[i] != P[j])
      {
        int tmp = P[i];
        P[i] = P[j];
        P[j] = tmp;
        swapIdx++;
      }
    }

    if (!Transfer(B, N, P, K))
      continue; //no transfert achieved
    Renumber(N, P, &K);

    double score = GetScore(B, N, P);
    if (score > bestScore)
    {
      for (int i=0; i<N; i++)
        bestP[i] = P[i];
      bestScore = score;
      bestK = K;
    }
  }

  for (int i=0; i<N; i++)
    P[i] = bestP[i];
  free(bestP);
}

// Find a partition which optimize a modularity criterion:
// Iterated Fusion-Transfer algorithm.
void tfit_core(int** A, int N, int* P)
{
  // Compute the modularity matrix (edges' weights)
  double** B = GetModularityMatrix(N, A);
  // Initialize partition: each vertex in its own class
  for (int i=0; i<N; i++)
    P[i] = i;
  int K = N;

  while (1)
  {
    if (!Fusion(B, N, P, K))
      break;
    Renumber(N, P, &K);
    if (!Transfer(B, N, P, K))
      break;
    Renumber(N, P, &K);
  }
  StochasticOptimization(B, N, P, K);

  for (int i=0; i<N; i++)
    free(B[i]);
  free(B);

  // R expect integers indexes to start at 1:
  for (int i=0; i<N; i++)
    P[i]++;
}
