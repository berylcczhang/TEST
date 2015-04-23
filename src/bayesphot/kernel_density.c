/*********************************************************************
Copyright (C) 2014 Robert da Silva, Michele Fumagalli, Mark Krumholz
This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.
*********************************************************************/

#include <assert.h>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_sf_gamma.h>
#include "geometry.h"
#include "kernel_density.h"

/*********************************************************************/
/* Useful macros                                                     */
/*********************************************************************/

#define SQR(x) ((x)*(x))
#define NODEBLOCKSIZE 1024

/*********************************************************************/
/* Static functions                                                  */
/*********************************************************************/

/* Functions to compute PDFs and integrals thereof on single tree nodes */
static inline
double kd_pdf_node(const kernel_density *kd, const double *x, 
		   const unsigned int curnode);

static inline
double kd_pdf_node_int(const kernel_density *kd, const double *x, 
		       const unsigned int *dims, const unsigned int ndim,
		       const unsigned int ndim_int, const double fac,
		       const unsigned int curnode);


/*********************************************************************/
/* Function to evaluate the PDF approximately using a kernel_density */
/* object                                                            */
/*********************************************************************/
double kd_pdf(const kernel_density *kd, const double *x,
	      const double reltol, const double abstol
#ifdef DIAGNOSTIC
	      , unsigned int *nodecheck, unsigned int *leafcheck,
	      unsigned int *termcheck
#endif
	      ) {
  unsigned int i, nnode, nalloc, ptr, lchild, rchild;
  unsigned int *nodelist;
  double *nodepdf;
  double maxerr, relerr, abserr;
  double pdf, leftpdf, rightpdf;

  /* Initialize for diagnostic mode */
#ifdef DIAGNOSTIC
  *nodecheck = *leafcheck = *termcheck = 0;
#endif

  /* Analyze root node */
  pdf = kd_pdf_node(kd, x, ROOT);

  /* Record diagnostic data */
#ifdef DIAGNOSTIC
  (*nodecheck)++;
  if (kd->tree->tree[ROOT].splitdim == -1) (*leafcheck)++;
  (*termcheck)++;
#endif

  /* Initialize from root node */
  if (kd->tree->tree[ROOT].splitdim == -1) {

    /* Special case: if root node is a leaf, we just got the exact
       value, so return it and exit. */
    return(pdf);

  } else {

    /* The usual case: root node is not a leaf */

    /* Allocate memory for node list */
    if (!(nodelist = (unsigned int *) 
	  calloc(NODEBLOCKSIZE, sizeof(unsigned int)))) {
      fprintf(stderr, "bayesphot: error: unable to allocate memory in kd_pdf\n");
      exit(1);
    }
    if (!(nodepdf = (double *) 
	  calloc(NODEBLOCKSIZE, sizeof(double)))) {
      fprintf(stderr, "bayesphot: error: unable to allocate memory in kd_pdf\n");
      exit(1);
    }
    nalloc = NODEBLOCKSIZE;  

    /* Initialize error estimates */
    abserr = pdf;
    relerr = abserr / (pdf + DBL_MIN);

    /* Push root node onto the node list */
    nodelist[0] = ROOT;
    nodepdf[0] = pdf;
    nnode = 1;
  }

  /* Now work recursively through the tree, identifying the node in
     the tree that contributes most to the error and opening it until
     the error estimate is within our tolerance */
  while (1) {

    /* Find the node that is contributing the most error. */
    maxerr = nodepdf[0];
    ptr = 0;
    for (i=1; i<nnode; i++) {
      if (nodepdf[i] > maxerr) {
	ptr = i;
	maxerr = nodepdf[i];
      }
    }

    /* Subtract this node's contribution to the global
       estimates. Enforce positivity to avoid spurious negative
       results coming from roundoff error. */
    pdf -= nodepdf[ptr];
    abserr -= nodepdf[ptr];
    if (pdf < 0) pdf = 0.0;
    if (abserr < 0) abserr = 0.0;

    /* Compute estimates for that node's children */
    lchild = LEFT(nodelist[ptr]);
    rchild = RIGHT(nodelist[ptr]);
    leftpdf = kd_pdf_node(kd, x, lchild);
    rightpdf = kd_pdf_node(kd, x, rchild);

    /* Get new PDF and error estimates */
    pdf += leftpdf + rightpdf;
    if (kd->tree->tree[lchild].splitdim != -1) abserr += leftpdf;
    if (kd->tree->tree[rchild].splitdim != -1) abserr += rightpdf;
    relerr = abserr / (pdf + DBL_MIN);

    /* Check the termination condition */
    if ((abserr <= abstol) || (relerr <= reltol)) break;

    /* Remove the node we just analyzed from the list */
    for (i=ptr; i<nnode-1; i++) {
      nodelist[i] = nodelist[i+1];
      nodepdf[i] = nodepdf[i+1];
    }
    nnode--;

    /* Allocate more memory to hold child nodes if necessary */
    if (nnode+2 >= nalloc) {
      if (!(nodelist = (unsigned int *) 
	    realloc(nodelist, 2*nalloc*sizeof(unsigned int)))) {
	fprintf(stderr, "bayesphot: error: unable to allocate memory in kd_pdf\n");
	exit(1);
      }
      if (!(nodepdf = (double *) 
	    realloc(nodepdf, 2*nalloc*sizeof(double)))) {
	fprintf(stderr, "bayesphot: error: unable to allocate memory in kd_pdf\n");
	exit(1);
      }
      nalloc *= 2;
    }

    /* If children are not leaves, push them onto the node list */
    if (kd->tree->tree[lchild].splitdim != -1) {
      nodelist[nnode] = lchild;
      nodepdf[nnode] = leftpdf;
      nnode++;
    }
    if (kd->tree->tree[rchild].splitdim != -1) {
      nodelist[nnode] = rchild;
      nodepdf[nnode] = rightpdf;
      nnode++;
    }

    /* Bail out if there are no non-leaf nodes left to be
       analyzed. This should also give abserr = relerr = 0, and thus
       we should have exited the loop when we checked the termination
       condition. However, this can sometimes fail due to roundoff
       error if abstol and reltol are both very small, so this is a
       backstop. */
    if (nnode == 0) break;

    /* Record diagnostic data */
#ifdef DIAGNOSTIC
    *nodecheck += 2;
    if (kd->tree->tree[lchild].splitdim != -1) (*leafcheck)++;
    if (kd->tree->tree[rchild].splitdim != -1) (*leafcheck)++;
    (*termcheck)++;
#endif

  }

  /* Free memory */
  free(nodelist);
  free(nodepdf);

  /* Return */
  return(pdf);
}


/*********************************************************************/
/* Function to evaluate the PDF integrated over certain dimensions   */
/* approximately using a kernel_density object                       */
/*********************************************************************/
double kd_pdf_int(const kernel_density *kd, const double *x,
		  const unsigned int *dims, const unsigned int ndim,
		  const double reltol, const double abstol
#ifdef DIAGNOSTIC
		  , unsigned int *nodecheck, unsigned int *leafcheck,
		  unsigned int *termcheck
#endif
		  ) {
  unsigned int i, ndim_int, nnode, nalloc, ptr, lchild, rchild;
  unsigned int *nodelist;
  double hprod, ds_n, fac;
  double *nodepdf;
  double maxerr, relerr, abserr;
  double pdf, leftpdf, rightpdf;

  /* Pre-compute constant factor in the integrals we're evaluating;
     this is the part that depends only on h and the number of
     dimensions. */
  ndim_int = kd->tree->ndim - ndim;
  ds_n = ds(ndim_int);
  hprod = 1.0;
  for (i=0; i<kd->tree->ndim; i++) hprod *= kd->h[i];
  for (i=0; i<ndim; i++) hprod /= kd->h[dims[i]];
  fac = hprod * ds_n;
  switch (kd->ktype) {
  case epanechnikov: {
    fac *= 2.0 / (ndim_int * (ndim_int + 2));
    break;
  }
  case tophat: {
    fac /= ndim_int;
    break;
  }
  case gaussian: {
    fac *= pow(2.0, ndim_int/2.0 - 1) * gsl_sf_gamma(0.5*ndim_int);
    break;
  }
  }

  /* Initialize for diagnostic mode */
#ifdef DIAGNOSTIC
  *nodecheck = *leafcheck = *termcheck = 0;
#endif

  /* Estimate integral from root node */
  pdf = kd_pdf_node_int(kd, x, dims, ndim, ndim_int, fac, ROOT);

  /* Initialize from root node */
  if (kd->tree->tree[ROOT].splitdim == -1) {

    /* Special case: if root node is a leaf, we just got the exact
       value, so return it and exit. */
    return(pdf);

  } else {

    /* The usual case: root node is not a leaf */

    /* Allocate memory for node list */
    if (!(nodelist = (unsigned int *) 
	  calloc(NODEBLOCKSIZE, sizeof(unsigned int)))) {
      fprintf(stderr, "bayesphot: error: unable to allocate memory in kd_pdf_int\n");
      exit(1);
    }
    if (!(nodepdf = (double *) 
	  calloc(NODEBLOCKSIZE, sizeof(double)))) {
      fprintf(stderr, "bayesphot: error: unable to allocate memory in kd_pdf_int\n");
      exit(1);
    }
    nalloc = NODEBLOCKSIZE;  

    /* Initialize error estimates */
    abserr = pdf;
    relerr = abserr / (pdf + DBL_MIN);

    /* Push root node onto the node list */
    nodelist[0] = ROOT;
    nodepdf[0] = pdf;
    nnode = 1;
  }

  /* Record diagnostic data */
#ifdef DIAGNOSTIC
  (*nodecheck)++;
  if (kd->tree->tree[ROOT].splitdim == -1) (*leafcheck)++;
  (*termcheck)++;
#endif

  /* Now work recursively through the tree, identifying the node in
     the tree that contributes most to the error and opening it until
     the error estimate is within our tolerance */
  while (1) {

    /* Find the node that is contributing the most error. */
    maxerr = nodepdf[0];
    ptr = 0;
    for (i=1; i<nnode; i++) {
      if (nodepdf[i] > maxerr) {
	ptr = i;
	maxerr = nodepdf[i];
      }
    }

    /* Subtract this node's contribution to the global
       estimates. Enforce positivity to avoid spurious negative
       results coming from roundoff error. */
    pdf -= nodepdf[ptr];
    abserr -= nodepdf[ptr];
    if (pdf < 0) pdf = 0.0;
    if (abserr < 0) abserr = 0.0;

    /* Compute estimates for that node's children */
    lchild = LEFT(nodelist[ptr]);
    rchild = RIGHT(nodelist[ptr]);
    leftpdf = kd_pdf_node_int(kd, x, dims, ndim, ndim_int, fac, lchild);
    rightpdf = kd_pdf_node_int(kd, x, dims, ndim, ndim_int, fac, rchild);

    /* Get new PDF and error estimates */
    pdf += leftpdf + rightpdf;
    if (kd->tree->tree[lchild].splitdim != -1) abserr += leftpdf;
    if (kd->tree->tree[rchild].splitdim != -1) abserr += rightpdf;
    relerr = abserr / (pdf + DBL_MIN);

    /* Check the termination condition */
    if ((abserr <= abstol) || (relerr <= reltol)) break;

    /* Remove the node we just analyzed from the list */
    for (i=ptr; i<nnode-1; i++) {
      nodelist[i] = nodelist[i+1];
      nodepdf[i] = nodepdf[i+1];
    }
    nnode--;

    /* Allocate more memory to hold child nodes if necessary */
    if (nnode+2 >= nalloc) {
      if (!(nodelist = (unsigned int *) 
	    realloc(nodelist, 2*nalloc*sizeof(unsigned int)))) {
	fprintf(stderr, "bayesphot: error: unable to allocate memory in kd_pdf_int\n");
	exit(1);
      }
      if (!(nodepdf = (double *) 
	    realloc(nodepdf, 2*nalloc*sizeof(double)))) {
	fprintf(stderr, "bayesphot: error: unable to allocate memory in kd_pdf_int\n");
	exit(1);
      }
      nalloc *= 2;
    }

    /* If children are not leaves, push them onto the node list */
    if (kd->tree->tree[lchild].splitdim != -1) {
      nodelist[nnode] = lchild;
      nodepdf[nnode] = leftpdf;
      nnode++;
    }
    if (kd->tree->tree[rchild].splitdim != -1) {
      nodelist[nnode] = rchild;
      nodepdf[nnode] = rightpdf;
      nnode++;
    }

    /* Bail out if there are no non-leaf nodes left to be
       analyzed. This should also give abserr = relerr = 0, and thus
       we should have exited the loop when we checked the termination
       condition. However, this can sometimes fail due to roundoff
       error if abstol and reltol are both very small, so this is a
       backstop. */
    if (nnode == 0) break;

    /* Record diagnostic data */
#ifdef DIAGNOSTIC
    *nodecheck += 2;
    if (kd->tree->tree[lchild].splitdim != -1) (*leafcheck)++;
    if (kd->tree->tree[rchild].splitdim != -1) (*leafcheck)++;
    (*termcheck)++;
#endif

  }

  /* Free memory */
  free(nodelist);
  free(nodepdf);

  /* Return */
  return(pdf);
}


/*********************************************************************/
/* Function to estimate the contribution to a PDF from a node; if    */
/* specified node is a leaf, the result returned is its exact        */
/* contribution to the PDF. If the node is not a leaf, the value     */
/* returned is half the upper limit on its contribution. If this is  */
/* then used as the actual estimate, the value returned serves is    */
/* both the estimated value and an absolute upper limit on the       */
/* error in that estimate. More accurate estimates are possible, but */
/* testing shows that, in high dimensions, they're not worth the     */
/* extra time that is required to compute them.                      */
/*********************************************************************/
inline
double kd_pdf_node(const kernel_density *kd, const double *x,
		   const unsigned int curnode) {
  unsigned int i;
  unsigned int ndim = kd->tree->ndim;
  double d2, pdf;

  /* Is this node a leaf? If so, just sum over it and return that */
  if (kd->tree->tree[curnode].splitdim == -1) {

    /* Loop over points, summing their contribution */
    pdf = 0.0;
    for (i=0; i<kd->tree->tree[curnode].npt; i++) {

      /* Get distance in units of the kernel size */
      d2 = dist2(&(kd->tree->tree[curnode].x[ndim*i]), x,
		 ndim, ndim, NULL, NULL, kd->h, ndim);

      /* Compute contribution based on kernel type and presence or
	 absence of weights */
      if (kd->tree->tree[curnode].dptr != NULL) {
	switch (kd->ktype) {
	case epanechnikov: {
	  if (d2 < 1)
	    pdf += ((double *) kd->tree->tree[curnode].dptr)[i] * 
	      (1.0 - d2);
	  break;
	}
	case tophat: {
	  if (d2 < 1)
	    pdf += ((double *) kd->tree->tree[curnode].dptr)[i];
	  break;
	}
	case gaussian: {
	  pdf += ((double *) kd->tree->tree[curnode].dptr)[i] *
	    exp(-d2/2.0);
	  break;
	}
	}
      } else {
	switch (kd->ktype) {
	case epanechnikov: {
	  if (d2 < 1) pdf += 1.0 - d2;
	  break;
	}
	case tophat: {
	  if (d2 < 1) pdf += 1.0;
	  break;
	}
	case gaussian: {
	  pdf += exp(-d2/2.0);
	  break;
	}
	}
      }

    }

    /* Normalize and return */
    return(pdf * kd->norm_tot);

  } else {

    /* This node is not a leaf, so return half the maximum possible
       contribution to the PDF */

    /* Get minimum distance in units of the kernel size */
    d2 = box_min_dist2(x, 
		       (const double **) kd->tree->tree[curnode].xbnd,
		       ndim, ndim, NULL, NULL, kd->h, ndim);

    /* Return PDF at minimum distance */
    switch (kd->ktype) {
    case epanechnikov: {
      return(d2 < 1.0 ? 0.5 * kd->nodewgt[curnode] * kd->norm_tot 
	     * (1.0 - d2) : 0.0);
    }
    case tophat: {
      return(d2 < 1.0 ? 0.5 * kd->nodewgt[curnode] * kd->norm_tot : 0.0);
    }
    case gaussian: {
      return(0.5 * kd->nodewgt[curnode] * kd->norm_tot * exp(-d2/2.0));
    }
    }
  }
}


/*********************************************************************/
/* Function to estimate the contribution to an integrated PDF from a */
/* node. Functionality is identical to kdf_pdf_node, except that     */
/* certain dimensions are being integrated out.                      */
/*********************************************************************/
inline
double kd_pdf_node_int(const kernel_density *kd, const double *x,
		       const unsigned int *dims, const unsigned int ndim,
		       const unsigned int ndim_int, const double fac,
		       const unsigned int curnode) {
  unsigned int i;
  unsigned int ndim_tot = kd->tree->ndim;
  double d2, pdf;

  /* Is this node a leaf? If so, just sum over it and return that */
  if (kd->tree->tree[curnode].splitdim == -1) {

    /* Loop over points, summing their contribution */
    pdf = 0.0;
    for (i=0; i<kd->tree->tree[curnode].npt; i++) {

      /* Get distance in units of the kernel size */
      d2 = dist2(&(kd->tree->tree[curnode].x[ndim_tot*i]), x,
		 ndim_tot, ndim, NULL, dims, kd->h, ndim_tot);

      /* Now compute contribution from each point */
      if (kd->tree->tree[curnode].dptr != NULL) {
	switch (kd->ktype) {
	case epanechnikov: {
	  if (d2 < 1)
	    pdf += ((double *) kd->tree->tree[curnode].dptr)[i] * 
	      pow(1.0-d2, 1.0+0.5*ndim_int);
	  break;
	}
	case tophat: {
	  if (d2 < 1)
	    pdf += ((double *) kd->tree->tree[curnode].dptr)[i] *
	      pow(1.0-d2, 0.5*ndim_int);
	  break;
	}
	case gaussian: {
	  pdf += ((double *) kd->tree->tree[curnode].dptr)[i] *
	    exp(-d2/2.0);
	  break;
	}
	}
      } else {
	switch (kd->ktype) {
	case epanechnikov: {
	  if (d2 < 1) pdf += pow(1.0-d2, 1.0+0.5*ndim_int);
	  break;
	}
	case tophat: {
	  if (d2 < 1) pdf += pow(1.0-d2, 0.5*ndim_int);
	  break;
	}
	case gaussian: {
	  pdf += exp(-d2/2.0);
	  break;
	}
	}
      }

    }

    /* Normalize and return */
    return(pdf * fac * kd->norm_tot);

  } else {

    /* This node is not a leaf, so return half the maximum possible
       contribution to the PDF */

    /* Get minimum distance in units of the kernel size */
    d2 = box_min_dist2(x, 
		       (const double **) kd->tree->tree[curnode].xbnd,
		       ndim, ndim_tot, dims, NULL, kd->h, ndim_tot);

    /* Return half the integral assuming all points are at the minimum
       distance. */
    switch (kd->ktype) {
    case epanechnikov: {
      return(d2 < 1.0 ? 0.5 * kd->nodewgt[curnode] * kd->norm_tot * fac * 
	      pow(1.0-d2, 1.0+0.5*ndim_int) : 0.0);
    }
    case tophat: {
      return(d2 < 1.0 ? 0.5 * kd->nodewgt[curnode] * kd->norm_tot * fac * 
	     pow(1.0-d2, 0.5*ndim_int) : 0.0);
    }
    case gaussian: {
      return(0.5 * kd->nodewgt[curnode] * kd->norm_tot * fac * exp(-d2/2.0));
    }
    }
  }
}


/*********************************************************************/
/* Vectorized version of kd_pdf_int                                  */
/*********************************************************************/
void kd_pdf_int_vec(const kernel_density *kd, const double *x, 
		    const unsigned int *dims, const unsigned int ndim,
		    const unsigned int npt, const double reltol, 
		    const double abstol, double *pdf
#ifdef DIAGNOSTIC
		    , unsigned int *nodecheck, unsigned int *leafcheck,
		    unsigned int *termcheck
#endif
		    ) {

  unsigned int i;
#pragma omp parallel private(i)
  {
    #pragma omp for
    for (i=0; i<npt; i++) {
      pdf[i] = kd_pdf_int(kd, x+i*ndim, dims, ndim, reltol, abstol
#ifdef DIAGNOSTIC
			  , nodecheck+i, leafcheck+i, termcheck+i
#endif
			  );
    }
  }
}


/*********************************************************************/
/* Function to evaluate the PDF using a kernel_density object for a  */
/* vector of inputs positions                                        */
/*********************************************************************/
void kd_pdf_vec(const kernel_density *kd, const double *x, 
		const unsigned int npt, const double reltol, 
		const double abstol, double *pdf
#ifdef DIAGNOSTIC
		, unsigned int *nodecheck, unsigned int *leafcheck,
		unsigned int *termcheck
#endif
		) {

  unsigned int i;
#pragma omp parallel private(i)
  {
  #pragma omp for
    for (i=0; i<npt; i++) {
      pdf[i] = kd_pdf(kd, x+i*kd->tree->ndim, reltol, abstol
#ifdef DIAGNOSTIC
		      , nodecheck+i, leafcheck+i, termcheck+i
#endif
		      );
    }
  }
}

