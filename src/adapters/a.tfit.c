#include <R.h>
#include <Rdefines.h>
#include "tfit.h"

// See comments in src/sources/EMGLLF.c and R/EMGLLF.R (wrapper)
SEXP tfit(
	SEXP X_
) {
	// TODO: io...
	const char* X = CHAR(STRING_PTR(X_)[0]);

	////////////////////
	// Call to function //
	////////////////////
	char* out = calloc(64, sizeof(char));
	tfit_core(X,out);

	// TODO: build output and return it
	SEXP output = PROTECT(allocVector(STRSXP, 1));
  SET_STRING_ELT(output, 0, mkChar(out));
	UNPROTECT(1);

	return output;
}
