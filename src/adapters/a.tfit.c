#include <R.h>
#include <Rdefines.h>
#include "tfit.h"

SEXP tfit(SEXP adjacencyMatrix_)
{
	// Get number of vertices (size of input matrix)
	int n = INTEGER(getAttrib(adjacencyMatrix_, R_DimSymbol))[0];
	const int* phiInit = LOGICAL(phiInit_);
	const char* out = CHAR(STRING_PTR(out_)[0]);

  // Prepare output:
  SEXP partition;
	PROTECT(partition = allocVector(INTSXP, n));
	int* pPartition = INTEGER(partition);

  // Matrix is given by columns (as in R)
	tfit_core(adjacencyMatrix, n, pPartition);

	UNPROTECT(1);
	return partition;
}
