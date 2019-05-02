#include <R.h>
#include <Rdefines.h>
#include "tfit.h"

SEXP tfit(SEXP adjacencyMatrix_)
{
	// Get number of vertices (size of input matrix)
	int nbVertices = INTEGER(getAttrib(adjacencyMatrix_, R_DimSymbol))[0];

  // Prepare input
  int** adjacencyMatrix = malloc(nbVertices*sizeof(double*));
  const int* A = INTEGER(adjacencyMatrix_);
  for (int i=0; i<nbVertices; i++)
  {
    adjacencyMatrix[i] = malloc(nbVertices*sizeof(double));
    for (int j=0; j<nbVertices; j++)
      adjacencyMatrix[i][j] = A[j*nbVertices+i];
  }

  // Prepare output
  SEXP partition;
	PROTECT(partition = allocVector(INTSXP, nbVertices));
	int* pPartition = INTEGER(partition);

  // The adjacency matrix is given by rows:
	tfit_core(adjacencyMatrix, nbVertices, pPartition);

  // Release memory
	UNPROTECT(1);
  for (int i=0; i<nbVertices; i++)
    free(adjacencyMatrix[i]);
  free(adjacencyMatrix);

	return partition;
}
