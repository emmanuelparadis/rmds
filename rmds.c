#include <R.h>
#include <Rinternals.h>
#include <Rmath.h>

/* returns 1 if both bases are different surely, 0 otherwise */
#define DifferentBase(a, b) ((a & b) < 16)

void mysample2(unsigned char *taken, int *n, int *i)
{
    int k;
    for (;;) {
	GetRNGstate();
	k = floor(n[0] * unif_rand());
	PutRNGstate();
	if (!taken[k]) {
	    i[0] = k;
	    break;
	}
    }
}

int mysample(unsigned char *taken, int n)
{
    int i;
    for (;;) {
	GetRNGstate();
	i = floor(n * unif_rand());
	PutRNGstate();
	if (!taken[i]) return i;
    }
}

SEXP select_Euclidean(SEXP X, SEXP M)
{
    int n, p, m, i, ii, j, k, o, *r;
    unsigned char *taken;
    double *x, dmax, d;
    SEXP res;

    PROTECT(X = coerceVector(X, REALSXP));
    PROTECT(M = coerceVector(M, INTSXP));
    m = INTEGER(M)[0];
    PROTECT(res = allocVector(INTSXP, m));
    r = INTEGER(res);

    n = nrows(X);
    p = ncols(X);
    x = REAL(X);

    taken = (unsigned char*)R_alloc(n, sizeof(unsigned char));
    memset(taken, 0, n*sizeof(unsigned char));

    //i = 0; // not random
    i = mysample(taken, n);
    o = 0;
    r[o] = i;
    o++;
    taken[i] = 1;

    for (;;) {
	dmax = 0;
	k = 0;
	for (ii = 0; ii < n; ii++) {
	    if (taken[ii] || ii == i) continue;
	    d = 0;
	    for (j = 0; j < p; j++)
		d += R_pow_di(x[i + n * j] - x[ii + n * j], 2);
	    if (d > dmax) {
		k = ii;
		dmax = d;
	    }
	}
	r[o] = k;
	taken[k] = 1;
	o++;
	if (o == m) break;
	i = mysample(taken, n);
	r[o] = i;
	taken[i] = 1;
	o++;
	if (o == m) break;
    }

    for (i = 0; i < m; i++) (r[i])++;

    UNPROTECT(3);
    return res;
}

SEXP select_DNAbin(SEXP X, SEXP M)
{
    int n, p, m, i, ii, j, k, o, *r;
    unsigned char *x, *taken;
    double dmax, d;
    SEXP res;

    PROTECT(X = coerceVector(X, RAWSXP));
    PROTECT(M = coerceVector(M, INTSXP));
    m = INTEGER(M)[0];
    PROTECT(res = allocVector(INTSXP, m));
    r = INTEGER(res);

    n = nrows(X);
    p = ncols(X);
    x = RAW(X);

    taken = (unsigned char*)R_alloc(n, sizeof(unsigned char));
    memset(taken, 0, n*sizeof(unsigned char));

    //i = 0; // not random
    i = mysample(taken, n);
    o = 0;
    r[o] = i;
    o++;
    taken[i] = 1;

    for (;;) {
	dmax = 0;
	k = 0;
	for (ii = 0; ii < n; ii++) {
	    if (taken[ii] || ii == i) continue;
	    d = 0;
	    for (j = 0; j < p; j++) {
		if (x[i + n * j] == 4 || x[ii + n * j] == 4) continue;
		d += DifferentBase(x[i + n * j], x[ii + n * j]);
	    }
	    if (d > dmax) {
		k = ii;
		dmax = d;
	    }
	}
	r[o] = k;
	o++;
	if (o == m) break;
	i = mysample(taken, n);
	r[o] = i;
	o++;
	if (o == m) break;
    }

    for (i = 0; i < m; i++) (r[i])++;

    UNPROTECT(3);
    return res;
}

SEXP dist_dna_vector_to_matrix(SEXP vec, SEXP mat)
{
    int i, j, n, p;
    unsigned char *x, *v;
    double *r, d;
    SEXP res;

    PROTECT(vec = coerceVector(vec, RAWSXP));
    PROTECT(mat = coerceVector(mat, RAWSXP));
    n = nrows(mat);
    p = ncols(mat);
    x = RAW(mat);
    v = RAW(vec);

    PROTECT(res = allocVector(REALSXP, n));
    r = REAL(res);

    for (i = 0; i < n; i++) {
	d = 0;
	for (j = 0; j < p; j++) {
	    if (x[i + n * j] == 4 || v[j] == 4) continue;
	    d += DifferentBase(x[i + n * j], v[j]);
	}
	r[i] = d;
    }

    UNPROTECT(3);
    return res;
}
