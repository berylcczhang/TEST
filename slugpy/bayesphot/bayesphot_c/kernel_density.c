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
#include "kernel_density_util.h"

/*********************************************************************/
/* Useful macros                                                     */
/*********************************************************************/

#define SQR(x) ((x)*(x))
#define NODEBLOCKSIZE 16384

/*********************************************************************/
/* Static functions                                                  */
/*********************************************************************/

/* Functions to compute PDFs and integrals thereof on single tree nodes */
static inline
void kd_pdf_node(const kernel_density *kd, const double *x, 
		 const unsigned long curnode, double *pdf,
		 double *pdferr);

static inline
double kd_pdf_node_grid(const kernel_density *kd, const double *xfixed, 
			const unsigned long *dimfixed,
			const unsigned long ndimfixed,
			const double *xgrid,
			const unsigned long *dimgrid,
			const unsigned long ndimgrid,
			const unsigned long ngrid,
			const unsigned long curnode,
			double *pdf);

static inline
void kd_pdf_node_int(const kernel_density *kd, const double *x, 
		     const unsigned long *dims, const unsigned long ndim,
		     const unsigned long ndim_int, const double fac,
		     const unsigned long curnode,
		     double *pdf, double *pdferr);

static inline
double kd_pdf_node_int_grid(const kernel_density *kd, const double *xfixed, 
			    const unsigned long *dimfixed,
			    const unsigned long ndimfixed,
			    const double *xgrid,
			    const unsigned long *dimgrid,
			    const unsigned long ndimgrid,
			    const unsigned long ngrid,
			    const unsigned long ndim_int,
			    const double fac,
			    const unsigned long curnode,
			    double *pdf);

static inline
double kd_pdf_node_int_reggrid(const kernel_density *kd, 
			       const double *xfixed, 
			       const unsigned long *dimfixed,
			       const unsigned long ndimfixed,
			       const double *xgridlo,
			       const double *dxgrid,
			       const unsigned long *ngrid,
			       const unsigned long *nstencil,
			       const unsigned long ngridtot,
			       const unsigned long nstenciltot,
			       const unsigned long *dimgrid,
			       const unsigned long ndimgrid,
			       const unsigned long ndim_int,
			       const double fac,
			       const unsigned long curnode,
			       double *pdf);

static inline
double kd_pdf_node_reggrid(const kernel_density *kd, 
			   const double *xfixed, 
			   const unsigned long *dimfixed,
			   const unsigned long ndimfixed,
			   const double *xgridlo,
			   const double *dxgrid,
			   const unsigned long *ngrid,
			   const unsigned long *nstencil,
			   const unsigned long ngridtot,
			   const unsigned long nstenciltot,
			   const unsigned long *dimgrid,
			   const unsigned long ndimgrid,
			   const unsigned long curnode,
			   double *pdf);

/**********************************************************************/
/* Function to evaluate the PDF approximately using a kernel_density */
/* object                                                            */
/*********************************************************************/
double kd_pdf(const kernel_density *kd, const double *x,
	      const double reltol, const double abstol
#ifdef DIAGNOSTIC
	      , unsigned long *nodecheck, unsigned long *leafcheck,
	      unsigned long *termcheck
#endif
	      ) {
  unsigned long i, nnode, nalloc, ptr, lchild, rchild;
  unsigned long *nodelist;
  double *nodepdf, *nodeerr;
  double pdf, maxerr, relerr, abserr;
  double leafpdf, lpdf, lerr, rpdf, rerr;

  /* Initialize for diagnostic mode */
#ifdef DIAGNOSTIC
  *nodecheck = *leafcheck = *termcheck = 0;
#endif

  /* Analyze root node */
  kd_pdf_node(kd, x, ROOT, &pdf, &abserr);

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
    if (!(nodelist = (unsigned long *) 
	  calloc(NODEBLOCKSIZE, sizeof(unsigned long)))) {
      fprintf(stderr, "bayesphot: error: unable to allocate memory in kd_pdf\n");
      exit(1);
    }
    if (!(nodepdf = (double *) 
	  calloc(NODEBLOCKSIZE, sizeof(double)))) {
      fprintf(stderr, "bayesphot: error: unable to allocate memory in kd_pdf\n");
      exit(1);
    }
    if (!(nodeerr = (double *) 
	  calloc(NODEBLOCKSIZE, sizeof(double)))) {
      fprintf(stderr, "bayesphot: error: unable to allocate memory in kd_pdf\n");
      exit(1);
    }
    nalloc = NODEBLOCKSIZE;  

    /* Push root node onto the node list */
    nodelist[0] = ROOT;
    nodepdf[0] = pdf;
    nodeerr[0] = abserr;
    nnode = 1;
    leafpdf = 0.0;
  }

  /* Now work recursively through the tree, identifying the node in
     the tree that contributes most to the error and opening it until
     the error estimate is within our tolerance */
  while (1) {

    /* Compute the current estimate and total error, and find the node
       that is contributing the most to the error budget. Note that in
       principle we could do fewer arithmetic operations by keeping a
       running tally of abserr and pdf, and updating these each time
       we go through the loop, rather than recomputing them from
       scratch each time. However, in practice this involves a lot of
       cancelling additions and subtractions that can create problems
       with roundoff when the error tolerance is tight and a lot of
       node evaluations are needed to achieve it. Doing the
       calculation fresh each time through the loop avoids this
       problem, which is worth the small amount of extra work required. */
    abserr = maxerr = nodeerr[0];
    pdf = leafpdf + nodepdf[0];
    ptr = 0;
    for (i=1; i<nnode; i++) {
      pdf += nodepdf[i];
      abserr += nodeerr[i];
      if (nodeerr[i] > maxerr) {
	ptr = i;
	maxerr = nodeerr[i];
      }
    }

    /* Check for convergence */
    relerr = abserr / pdf;
    if ((abserr <= abstol) || (relerr <= reltol)) break;
    
    /* Compute estimates for that node's children */
    lchild = LEFT(nodelist[ptr]);
    rchild = RIGHT(nodelist[ptr]);
    kd_pdf_node(kd, x, lchild, &lpdf, &lerr);
    kd_pdf_node(kd, x, rchild, &rpdf, &rerr);

    /* If children are leaves, add their contribution to the PDF from
       leaves */
    if (kd->tree->tree[lchild].splitdim == -1) leafpdf += lpdf;
    if (kd->tree->tree[rchild].splitdim == -1) leafpdf += rpdf;

    /* Remove the node we just analyzed from the list */
    for (i=ptr; i<nnode-1; i++) {
      nodelist[i] = nodelist[i+1];
      nodepdf[i] = nodepdf[i+1];
      nodeerr[i] = nodeerr[i+1];
    }
    nnode--;

    /* Allocate more memory to hold child nodes if necessary */
    if (nnode+2 >= nalloc) {
      if (!(nodelist = (unsigned long *) 
	    realloc(nodelist, 2*nalloc*sizeof(unsigned long)))) {
	fprintf(stderr, "bayesphot: error: unable to allocate memory in kd_pdf\n");
	exit(1);
      }
      if (!(nodepdf = (double *) 
	    realloc(nodepdf, 2*nalloc*sizeof(double)))) {
	fprintf(stderr, "bayesphot: error: unable to allocate memory in kd_pdf\n");
	exit(1);
      }
      if (!(nodeerr = (double *) 
	    realloc(nodeerr, 2*nalloc*sizeof(double)))) {
	fprintf(stderr, "bayesphot: error: unable to allocate memory in kd_pdf\n");
	exit(1);
      }
      nalloc *= 2;
    }

    /* If children are not leaves, push them onto the node list */
    if (kd->tree->tree[lchild].splitdim != -1) {
      nodelist[nnode] = lchild;
      nodepdf[nnode] = lpdf;
      nodeerr[nnode] = lerr;
      nnode++;
    }
    if (kd->tree->tree[rchild].splitdim != -1) {
      nodelist[nnode] = rchild;
      nodepdf[nnode] = rpdf;
      nodeerr[nnode] = rerr;
      nnode++;
    }

    /* Record diagnostic data */
#ifdef DIAGNOSTIC
    *nodecheck += 2;
    if (kd->tree->tree[lchild].splitdim != -1) (*leafcheck)++;
    if (kd->tree->tree[rchild].splitdim != -1) (*leafcheck)++;
    (*termcheck)++;
#endif

    /* Safety check: bail out if no nodes left */
    if (nnode == 0) break;
  }

  /* Free memory */
  free(nodelist);
  free(nodepdf);
  free(nodeerr);

  /* Return */
  return(pdf);
}


/*********************************************************************/
/* Function to evaluate the PDF on a grid where some dimensions are  */
/* fixed and others are varying                                      */
/*********************************************************************/

void kd_pdf_grid(const kernel_density *kd, const double *xfixed,
		 const unsigned long *dimfixed, 
		 const unsigned long ndimfixed,
		 const unsigned long nfixed,
		 const double *xgrid,
		 const unsigned long *dimgrid,
		 const unsigned long ndimgrid,
		 const unsigned long ngrid,
		 const double reltol, const double abstol,
		 double *pdf) {

  unsigned long i, j;
  unsigned long nalloc=0, nnode, ptr, lchild, rchild;
  unsigned long fixedptr, outptr;
  unsigned long *nodelist = NULL;
  double *nodepdf = NULL;
  double relerr = 0.0, maxerr;
  double pdfmax, leftpdf, rightpdf;

  /* Loop over the input fixed points */
  for (i=0; i<nfixed; i++) {

    /* Move the pointers in the input and output arrays */
    fixedptr = i*ndimfixed;
    outptr = i*ngrid;

    /* Initialize output array to zero */
    for (j=0; j<ngrid; j++) pdf[i*ngrid+j] = 0.0;

    /* Analyze root node */
    pdfmax = kd_pdf_node_grid(kd, xfixed+fixedptr, dimfixed, ndimfixed,
			      xgrid, dimgrid, ndimgrid, ngrid, ROOT,
			      pdf+outptr);

    /* If root node is a leaf, we're done. Go to next fixed point. */
    if (kd->tree->tree[ROOT].splitdim == -1) continue;

    /* More common case where root node is not a leaf */

    /* Allocate memory for node list if necessary */
    if (nalloc == 0) {
      if (!(nodelist = (unsigned long *) 
	    calloc(NODEBLOCKSIZE, sizeof(unsigned long)))) {
	fprintf(stderr, "bayesphot: error: unable to allocate memory in kd_pdf_grid\n");
	exit(1);
      }
      if (!(nodepdf = (double *) 
	    calloc(NODEBLOCKSIZE, sizeof(double)))) {
	fprintf(stderr, "bayesphot: error: unable to allocate memory in kd_pdf_grid\n");
	exit(1);
      }
      nalloc = NODEBLOCKSIZE;  
    }

    /* Push root node onto the node list */
    nodelist[0] = ROOT;
    nodepdf[0] = pdfmax;
    nnode = 1;

    /* Now proceed through the tree, identifying the nodes that are
       contributing the most error and recursively opening them. */
    while (1) {

      /* Find the node that is contributing the most error. */
      maxerr = nodepdf[0];
      ptr = 0;
      for (j=1; j<nnode; j++) {
	if (nodepdf[j] > maxerr) {
	  ptr = j;
	  maxerr = nodepdf[j];
	}
      }

      /* Subtract this node's contribution to the maximum possible PDF
	 contribution from unopened nodes. Enforce positivity to avoid
	 spurious negative results coming from roundoff error. */
      pdfmax -= nodepdf[ptr];
      if (pdfmax < 0) pdfmax = 0.0;

      /* Compute estimates for this node's children */
      lchild = LEFT(nodelist[ptr]);
      rchild = RIGHT(nodelist[ptr]);
      leftpdf = kd_pdf_node_grid(kd, xfixed+fixedptr, dimfixed, ndimfixed,
				 xgrid, dimgrid, ndimgrid, ngrid, lchild,
				 pdf+outptr);
      rightpdf = kd_pdf_node_grid(kd, xfixed+fixedptr, dimfixed, ndimfixed,
				  xgrid, dimgrid, ndimgrid, ngrid, rchild,
				  pdf+outptr);

      /* Get new maximum contribution to PDF from unopened nodes */
      pdfmax += leftpdf + rightpdf;

      /* Check for termination on absolute error. Note that we can
	 save a factor of 2 by adding half the upper limit on the
	 contribution to every point. */
      if (0.5*pdfmax < abstol) {
	for (j=0; j<ngrid; j++) pdf[outptr+j] += 0.5*pdfmax;
	break;
      }

      /* Check for termination on relative error. Again, we save a
	 factor of 2 by adding half the upper limit on to each point. */
      for (j=0; j<ngrid; j++) {
	relerr = 0.5*pdfmax / (pdf[outptr+j] + 0.5*pdfmax + DBL_MIN);
	if (relerr >= reltol) break;
      }
      if (relerr < reltol) {
	for (j=0; j<ngrid; j++) pdf[outptr+j] += 0.5*pdfmax;
	break;
      }

      /* If we're here, we haven't converged yet. Remove the node we
	 just analyzed from the list. */
      for (j=ptr; j<nnode-1; j++) {
	nodelist[j] = nodelist[j+1];
	nodepdf[j] = nodepdf[j+1];
      }
      nnode--;

      /* Allocate more memory to hold child nodes if necessary */
      if (nnode+2 >= nalloc) {
	if (!(nodelist = (unsigned long *) 
	      realloc(nodelist, 2*nalloc*sizeof(unsigned long)))) {
	  fprintf(stderr, "bayesphot: error: unable to allocate memory in kd_pdf_grid\n");
	  exit(1);
	}
	if (!(nodepdf = (double *) 
	      realloc(nodepdf, 2*nalloc*sizeof(double)))) {
	  fprintf(stderr, "bayesphot: error: unable to allocate memory in kd_pdf_grid\n");
	  exit(1);
	}
	nalloc *= 2;
      }

      /* If children are not leaves, and their maximum contribution to
	 the PDF is not 0 (possible due to roundoff, or if the kernel
	 is compact), push them onto the node list */
      if (leftpdf > 0.0) {
	nodelist[nnode] = lchild;
	nodepdf[nnode] = leftpdf;
	nnode++;
      }
      if (rightpdf > 0.0) {
	nodelist[nnode] = rchild;
	nodepdf[nnode] = rightpdf;
	nnode++;
      }

      /* Bail out if there are no non-leaf nodes left to be
	 analyzed. This should also give abserr = 0, and thus
	 we should have exited the loop when we checked the termination
	 condition. However, this can sometimes fail due to roundoff
	 error if abstol and reltol are both very small, so this is a
	 backstop. */
      if (nnode == 0) break;
    }
  }

  /* Free memory */
  if (nodelist != NULL) free(nodelist);
  if (nodepdf != NULL) free(nodepdf);
}


/*********************************************************************/
/* Function to evaluate the PDF integrated over certain dimensions   */
/* approximately using a kernel_density object                       */
/*********************************************************************/
double kd_pdf_int(const kernel_density *kd, const double *x,
		  const unsigned long *dims, const unsigned long ndim,
		  const double reltol, const double abstol
#ifdef DIAGNOSTIC
		  , unsigned long *nodecheck, unsigned long *leafcheck,
		  unsigned long *termcheck
#endif
		  ) {
  unsigned long i, ndim_int, nnode, nalloc, ptr, lchild, rchild;
  unsigned long *nodelist = NULL;
  double hprod, ds_n, fac;
  double *nodepdf = NULL, *nodeerr = NULL;
  double maxerr, relerr, abserr;
  double pdf, leafpdf, lpdf, rpdf, lerr, rerr;
  
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

  /* Process the root node */
  kd_pdf_node_int(kd, x, dims, ndim, ndim_int, fac, ROOT,
		  &pdf, &abserr);

  /* Initialize node list */
  if (kd->tree->tree[ROOT].splitdim == -1) {

    /* Special case: root node is a leaf, so just return the exact value */
    return(pdf);

  } else {

    /* The usual case: root node is not a leaf */

    /* Allocate memory for node list */
    if (!(nodelist = (unsigned long *) 
	  calloc(NODEBLOCKSIZE, sizeof(unsigned long)))) {
      fprintf(stderr, "bayesphot: error: unable to allocate memory in kd_pdf_int\n");
      exit(1);
    }
    if (!(nodepdf = (double *) 
	  calloc(NODEBLOCKSIZE, sizeof(double)))) {
      fprintf(stderr, "bayesphot: error: unable to allocate memory in kd_pdf_int\n");
      exit(1);
    }
    if (!(nodeerr = (double *) 
	  calloc(NODEBLOCKSIZE, sizeof(double)))) {
      fprintf(stderr, "bayesphot: error: unable to allocate memory in kd_pdf_int\n");
      exit(1);
    }
    nalloc = NODEBLOCKSIZE;  

    /* Push root node onto the node list */
    nodelist[0] = ROOT;
    nodepdf[0] = pdf;
    nodeerr[0] = abserr;
    nnode = 1;
    leafpdf = 0.0;
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

    /* Compute the current estimate and total error, and find the node
       that makes the largest contribution to the total error */
    abserr = maxerr = nodeerr[0];
    pdf = leafpdf + nodepdf[0];
    ptr = 0;
    for (i=1; i<nnode; i++) {
      pdf += nodepdf[i];
      abserr += nodeerr[i];
      if (nodeerr[i] > maxerr) {
	ptr = i;
	maxerr = nodeerr[i];
      }
    }

    /* Check termination condition */
    relerr = abserr / pdf;
    if ((abserr <= abstol) || (relerr <= reltol)) break;

    /* Analyze the children of the node contributing the most error */
    lchild = LEFT(nodelist[ptr]);
    rchild = RIGHT(nodelist[ptr]);
    kd_pdf_node_int(kd, x, dims, ndim, ndim_int, fac, lchild, &lpdf, &lerr);
    kd_pdf_node_int(kd, x, dims, ndim, ndim_int, fac, rchild, &rpdf, &rerr);

    /* If children are leaves, add their contribution to the PDF from
       leaves */
    if (kd->tree->tree[lchild].splitdim == -1) leafpdf += lpdf;
    if (kd->tree->tree[rchild].splitdim == -1) leafpdf += rpdf;
    
    /* Remove the node we just analyzed from the list */
    for (i=ptr; i<nnode-1; i++) {
      nodelist[i] = nodelist[i+1];
      nodepdf[i] = nodepdf[i+1];
      nodeerr[i] = nodeerr[i+1];
    }
    nnode--;

    /* Allocate more memory to hold child nodes if necessary */
    if (nnode+2 >= nalloc) {
      if (!(nodelist = (unsigned long *) 
	    realloc(nodelist, 2*nalloc*sizeof(unsigned long)))) {
	fprintf(stderr, "bayesphot: error: unable to allocate memory in kd_pdf_int\n");
	exit(1);
      }
      if (!(nodepdf = (double *) 
	    realloc(nodepdf, 2*nalloc*sizeof(double)))) {
	fprintf(stderr, "bayesphot: error: unable to allocate memory in kd_pdf_int\n");
	exit(1);
      }
      if (!(nodeerr = (double *) 
	    realloc(nodeerr, 2*nalloc*sizeof(double)))) {
	fprintf(stderr, "bayesphot: error: unable to allocate memory in kd_pdf_int\n");
	exit(1);
      }
      nalloc *= 2;
    }

    /* If children are not leaves, push them onto the node list */
    if (kd->tree->tree[lchild].splitdim != -1) {
      nodelist[nnode] = lchild;
      nodepdf[nnode] = lpdf;
      nodeerr[nnode] = lerr;
      nnode++;
    }
    if (kd->tree->tree[rchild].splitdim != -1) {
      nodelist[nnode] = rchild;
      nodepdf[nnode] = rpdf;
      nodeerr[nnode] = rerr;
      nnode++;
    }

    /* Record diagnostic data */
#ifdef DIAGNOSTIC
    *nodecheck += 2;
    if (kd->tree->tree[lchild].splitdim != -1) (*leafcheck)++;
    if (kd->tree->tree[rchild].splitdim != -1) (*leafcheck)++;
    (*termcheck)++;
#endif
    
    /* Safety check: bail out if no nodes left */
    if (nnode == 0) break;
  }

  /* Free memory */
  free(nodelist);
  free(nodepdf);
  free(nodeerr);

  /* Return */
  return(pdf);
}


/*********************************************************************/
/* Function to evaluate the PDF on a grid where some dimensions are  */
/* fixed, others are varying, and still others are being integrated  */
/* out (i.e., marginalized over)                                     */
/*********************************************************************/

void kd_pdf_int_grid(const kernel_density *kd, const double *xfixed,
		     const unsigned long *dimfixed, 
		     const unsigned long ndimfixed,
		     const unsigned long nfixed,
		     const double *xgrid,
		     const unsigned long *dimgrid,
		     const unsigned long ndimgrid,
		     const unsigned long ngrid,
		     const double reltol, const double abstol,
		     double *pdf) {

  unsigned long i, j;
  unsigned long nalloc=0, nnode, ptr, lchild, rchild, ndim_int;
  unsigned long fixedptr, outptr;
  unsigned long *nodelist = NULL;
  double *nodepdf = NULL;
  double hprod, ds_n, fac;
  double maxerr, relerr = 0.0, pdfmax, leftpdf, rightpdf;

  /* Pre-compute constant factor in the integrals we're evaluating;
     this is the part that depends only on h and the number of
     dimensions. */
  ndim_int = kd->tree->ndim - ndimfixed - ndimgrid;
  ds_n = ds(ndim_int);
  hprod = 1.0;
  for (i=0; i<kd->tree->ndim; i++) hprod *= kd->h[i];
  for (i=0; i<ndimfixed; i++) hprod /= kd->h[dimfixed[i]];
  for (i=0; i<ndimgrid; i++) hprod /= kd->h[dimgrid[i]];
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

  /* Loop over the input fixed points */
  for (i=0; i<nfixed; i++) {

    /* Move the pointers in the input and output arrays */
    fixedptr = i*ndimfixed;
    outptr = i*ngrid;

    /* Initialize output array to zero */
    for (j=0; j<ngrid; j++) pdf[i*ngrid+j] = 0.0;

    /* Analyze root node */
    pdfmax = kd_pdf_node_int_grid(kd, xfixed+fixedptr, dimfixed, 
				  ndimfixed, xgrid, dimgrid, ndimgrid, 
				  ngrid, ndim_int, fac, ROOT, pdf+outptr);

    /* If root node is a leaf, we're done. Go to next fixed point. */
    if (kd->tree->tree[ROOT].splitdim == -1) continue;

    /* More common case where root node is not a leaf */

    /* Allocate memory for node list if necessary */
    if (nalloc == 0) {
      if (!(nodelist = (unsigned long *) 
	    calloc(NODEBLOCKSIZE, sizeof(unsigned long)))) {
	fprintf(stderr, "bayesphot: error: unable to allocate memory in kd_pdf_int_grid\n");
	exit(1);
      }
      if (!(nodepdf = (double *) 
	    malloc(NODEBLOCKSIZE*sizeof(double)))) {
	fprintf(stderr, "bayesphot: error: unable to allocate memory in kd_pdf_int_grid\n");
	exit(1);
      }
      nalloc = NODEBLOCKSIZE;  
    }

    /* Push root node onto the node list */
    nodelist[0] = ROOT;
    nodepdf[0] = pdfmax;
    nnode = 1;

    /* Now proceed through the tree, identifying the nodes that are
       contributing the most error and recursively opening them. */
    while (1) {

      /* Find the node that is contributing the most error. */
      maxerr = nodepdf[0];
      ptr = 0;
      for (j=1; j<nnode; j++) {
	if (nodepdf[j] > maxerr) {
	  ptr = j;
	  maxerr = nodepdf[j];
	}
      }

      /* Subtract this node's contribution to the maximum possible PDF
	 contribution from unopened nodes. Enforce positivity to avoid
	 spurious negative results coming from roundoff error. */
      pdfmax -= nodepdf[ptr];
      if (pdfmax < 0) pdfmax = 0.0;

      /* Compute estimates for this node's children */
      lchild = LEFT(nodelist[ptr]);
      rchild = RIGHT(nodelist[ptr]);
      leftpdf = kd_pdf_node_int_grid(kd, xfixed+fixedptr, dimfixed, 
				     ndimfixed, xgrid, dimgrid, ndimgrid, 
				     ngrid, ndim_int, fac, lchild, 
				     pdf+outptr);
      rightpdf = kd_pdf_node_int_grid(kd, xfixed+fixedptr, dimfixed, 
				      ndimfixed, xgrid, dimgrid, ndimgrid, 
				      ngrid, ndim_int, fac, rchild, 
				      pdf+outptr);

      /* Get new maximum contribution to PDF from unopened nodes */
      pdfmax += leftpdf + rightpdf;

      /* Check for termination on absolute error. Note that we can
	 save a factor of 2 by adding half the upper limit on the
	 contribution to every point. */
      if (0.5*pdfmax < abstol) {
	for (j=0; j<ngrid; j++) pdf[outptr+j] += 0.5*pdfmax;
	break;
      }

      /* Check for termination on relative error. Again, we save a
	 factor of 2 by adding half the upper limit on to each point. */
      for (j=0; j<ngrid; j++) {
	relerr = 0.5*pdfmax / (pdf[outptr+j] + 0.5*pdfmax + DBL_MIN);
	if (relerr >= reltol) break;
      }
      if (relerr < reltol) {
	for (j=0; j<ngrid; j++) pdf[outptr+j] += 0.5*pdfmax;
	break;
      }

      /* If we're here, we haven't converged yet. Remove the node we
	 just analyzed from the list. */
      for (j=ptr; j<nnode-1; j++) {
	nodelist[j] = nodelist[j+1];
	nodepdf[j] = nodepdf[j+1];
      }
      nnode--;

      /* Allocate more memory to hold child nodes if necessary */
      if (nnode+2 >= nalloc) {
	if (!(nodelist = (unsigned long *) 
	      realloc(nodelist, 2*nalloc*sizeof(unsigned long)))) {
	  fprintf(stderr, "bayesphot: error: unable to allocate memory in kd_pdf_grid\n");
	  exit(1);
	}
	if (!(nodepdf = (double *) 
	      realloc(nodepdf, 2*nalloc*sizeof(double)))) {
	  fprintf(stderr, "bayesphot: error: unable to allocate memory in kd_pdf_grid\n");
	  exit(1);
	}
	nalloc *= 2;
      }

      /* If children are not leaves, and their maximum contribution to
	 the PDF is not 0 (possible due to roundoff, or if the kernel
	 is compact), push them onto the node list */
      if (leftpdf > 0.0) {
	nodelist[nnode] = lchild;
	nodepdf[nnode] = leftpdf;
	nnode++;
      }
      if (rightpdf > 0.0) {
	nodelist[nnode] = rchild;
	nodepdf[nnode] = rightpdf;
	nnode++;
      }

      /* Bail out if there are no non-leaf nodes left to be
	 analyzed. This should also give abserr = 0, and thus
	 we should have exited the loop when we checked the termination
	 condition. However, this can sometimes fail due to roundoff
	 error if abstol and reltol are both very small, so this is a
	 backstop. */
      if (nnode == 0) break;
    }
  }

  /* Free memory */
  if (nodelist != NULL) free(nodelist);
  if (nodepdf != NULL) free(nodepdf);
}


/*********************************************************************/
/* Same as kd_pdf_int_grid, but for regular grids                    */
/*********************************************************************/
void kd_pdf_int_reggrid(const kernel_density *kd, const double *xfixed,
			const unsigned long *dimfixed, 
			const unsigned long ndimfixed,
			const unsigned long nfixed,
			const double *xgridlo,
			const double *xgridhi,
			const unsigned long *ngrid,
			const unsigned long *dimgrid,
			const unsigned long ndimgrid,
			const double reltol, const double abstol,
			double *pdf)
{
  unsigned long i, j;
  unsigned long nalloc=0, nnode, ptr, lchild, rchild, ndim_int;
  unsigned long ngridtot, nstenciltot, fixedptr, outptr;
  unsigned long *nodelist = NULL, *nstencil;
  double *dxgrid;
  double *nodepdf = NULL;
  double hprod, ds_n, fac;
  double maxerr, relerr = 0.0, pdfmax, leftpdf, rightpdf;

  /* Pre-compute constant factor in the integrals we're evaluating;
     this is the part that depends only on h and the number of
     dimensions. */
  ndim_int = kd->tree->ndim - ndimfixed - ndimgrid;
  ds_n = ds(ndim_int);
  hprod = 1.0;
  for (i=0; i<kd->tree->ndim; i++) hprod *= kd->h[i];
  for (i=0; i<ndimfixed; i++) hprod /= kd->h[dimfixed[i]];
  for (i=0; i<ndimgrid; i++) hprod /= kd->h[dimgrid[i]];
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

  /* Figure out the total number of points in the grid */
  ngridtot = ngrid[0];
  for (i=1; i<ndimgrid; i++) ngridtot *= ngrid[i];

  /* Figure out the grid spacing */
  if (!(dxgrid = (double *) 
	calloc(ndimgrid, sizeof(double)))) {
    fprintf(stderr, 
	    "bayesphot: error: unable to allocate memory in kd_pdf_int_reggrid\n");
    exit(1);
  }
  for (i=0; i<ndimgrid; i++) {
    if (ngrid[i] > 1) 
      dxgrid[i] = (xgridhi[i]-xgridlo[i]) / (ngrid[i]-1);
    else
      dxgrid[i] = 0.0;
  }

  /* Figure out how big a stencil we need */
  if (!(nstencil = (unsigned long *) 
	calloc(ndimgrid, sizeof(unsigned long)))) {
    fprintf(stderr, 
	    "bayesphot: error: unable to allocate memory in kd_pdf_int_reggrid\n");
    exit(1);
  }
  nstenciltot = 1;
  for (i=0; i<ndimgrid; i++) {
    if (ngrid[i] > 1) {
      switch (kd->ktype) {
      case epanechnikov: {
	nstencil[i] = (int) 
	  floor((kd->h[dimgrid[i]]-dxgrid[i]/2.0)/dxgrid[i]);
	break;
      }
      case tophat: {
	nstencil[i] = (int) 
	  floor((kd->h[dimgrid[i]]-dxgrid[i]/2.0)/dxgrid[i]);
	break;
      }
      case gaussian: {
	nstencil[i] = (int) 
	  floor(0.5 * (1.0+sqrt(1.0 - 8.0*log(reltol) * 
			       (kd->h[dimgrid[i]]/dxgrid[i]) *
			       (kd->h[dimgrid[i]]/dxgrid[i]))));
	break;
      }
      }
      if (nstencil[i] >= ngrid[i]) nstencil[i] = ngrid[i]-1;
    } else {
      nstencil[i] = 0;
    }
    nstenciltot *= (2*nstencil[i]+1);
  }

  /* Loop over the input fixed points */
  for (i=0; i<nfixed; i++) {

    /* Move the pointers in the input and output arrays */
    fixedptr = i*ndimfixed;
    outptr = i*ngridtot;

    /* Initialize output array to zero */
    for (j=0; j<ngridtot; j++) pdf[i*ngridtot+j] = 0.0;

    /* Analyze root node */
    pdfmax = 
      kd_pdf_node_int_reggrid(kd, xfixed+fixedptr, dimfixed, ndimfixed, 
			      xgridlo, dxgrid, ngrid, nstencil,
			      ngridtot, nstenciltot, dimgrid,
			      ndimgrid, ndim_int, fac, ROOT, 
			      pdf+outptr);

    /* If root node is a leaf, we're done. Go to next fixed point. */
    if (kd->tree->tree[ROOT].splitdim == -1) continue;

    /* More common case where root node is not a leaf */

    /* Allocate memory for node list if necessary */
    if (nalloc == 0) {
      if (!(nodelist = (unsigned long *) 
	    calloc(NODEBLOCKSIZE, sizeof(unsigned long)))) {
	fprintf(stderr, "bayesphot: error: unable to allocate memory in kd_pdf_int_grid\n");
	exit(1);
      }
      if (!(nodepdf = (double *) 
	    malloc(NODEBLOCKSIZE*sizeof(double)))) {
	fprintf(stderr, "bayesphot: error: unable to allocate memory in kd_pdf_int_grid\n");
	exit(1);
      }
      nalloc = NODEBLOCKSIZE;  
    }

    /* Push root node onto the node list */
    nodelist[0] = ROOT;
    nodepdf[0] = pdfmax;
    nnode = 1;

    /* Now proceed through the tree, identifying the nodes that are
       contributing the most error and recursively opening them. */
    while (1) {

      /* Find the node that is contributing the most error. */
      maxerr = nodepdf[0];
      ptr = 0;
      for (j=1; j<nnode; j++) {
	if (nodepdf[j] > maxerr) {
	  ptr = j;
	  maxerr = nodepdf[j];
	}
      }

      /* Subtract this node's contribution to the maximum possible PDF
	 contribution from unopened nodes. Enforce positivity to avoid
	 spurious negative results coming from roundoff error. */
      pdfmax -= nodepdf[ptr];
      if (pdfmax < 0) pdfmax = 0.0;

      /* Compute estimates for this node's children */
      lchild = LEFT(nodelist[ptr]);
      rchild = RIGHT(nodelist[ptr]);
      leftpdf = 
	kd_pdf_node_int_reggrid(kd, xfixed+fixedptr, dimfixed, ndimfixed, 
				xgridlo, dxgrid, ngrid, nstencil,
				ngridtot, nstenciltot, dimgrid,
				ndimgrid, ndim_int, fac, lchild, 
				pdf+outptr);
      rightpdf =
	kd_pdf_node_int_reggrid(kd, xfixed+fixedptr, dimfixed, ndimfixed, 
				xgridlo, dxgrid, ngrid, nstencil,
				ngridtot, nstenciltot, dimgrid,
				ndimgrid, ndim_int, fac, rchild, 
				pdf+outptr);

      /* Get new maximum contribution to PDF from unopened nodes */
      pdfmax += leftpdf + rightpdf;

      /* Check for termination on absolute error. Note that we can
	 save a factor of 2 by adding half the upper limit on the
	 contribution to every point. */
      if (0.5*pdfmax < abstol) {
	for (j=0; j<ngridtot; j++) pdf[outptr+j] += 0.5*pdfmax;
	break;
      }

      /* Check for termination on relative error. Again, we save a
	 factor of 2 by adding half the upper limit on to each point. */
      for (j=0; j<ngridtot; j++) {
	relerr = 0.5*pdfmax / (pdf[outptr+j] + 0.5*pdfmax + DBL_MIN);
	if (relerr >= reltol) break;
      }
      if (relerr < reltol) {
	for (j=0; j<ngridtot; j++) pdf[outptr+j] += 0.5*pdfmax;
	break;
      }

      /* If we're here, we haven't converged yet. Remove the node we
	 just analyzed from the list. */
      for (j=ptr; j<nnode-1; j++) {
	nodelist[j] = nodelist[j+1];
	nodepdf[j] = nodepdf[j+1];
      }
      nnode--;

      /* Allocate more memory to hold child nodes if necessary */
      if (nnode+2 >= nalloc) {
	if (!(nodelist = (unsigned long *) 
	      realloc(nodelist, 2*nalloc*sizeof(unsigned long)))) {
	  fprintf(stderr, "bayesphot: error: unable to allocate memory in kd_pdf_grid\n");
	  exit(1);
	}
	if (!(nodepdf = (double *) 
	      realloc(nodepdf, 2*nalloc*sizeof(double)))) {
	  fprintf(stderr, "bayesphot: error: unable to allocate memory in kd_pdf_grid\n");
	  exit(1);
	}
	nalloc *= 2;
      }

      /* If children are not leaves, and their maximum contribution to
	 the PDF is not 0 (possible due to roundoff, or if the kernel
	 is compact), push them onto the node list */
      if (leftpdf > 0.0) {
	nodelist[nnode] = lchild;
	nodepdf[nnode] = leftpdf;
	nnode++;
      }
      if (rightpdf > 0.0) {
	nodelist[nnode] = rchild;
	nodepdf[nnode] = rightpdf;
	nnode++;
      }

      /* Bail out if there are no non-leaf nodes left to be
	 analyzed. This should also give abserr = 0, and thus
	 we should have exited the loop when we checked the termination
	 condition. However, this can sometimes fail due to roundoff
	 error if abstol and reltol are both very small, so this is a
	 backstop. */
      if (nnode == 0) break;
    }
  }

  /* Free memory */
  if (nodelist != NULL) free(nodelist);
  if (nodepdf != NULL) free(nodepdf);
  free(dxgrid);
  free(nstencil);
}


/*********************************************************************/
/* Routine to compute the PDF in the case where some dimensions are  */
/* held fixed and the others are varied along a regular grid.        */
/*********************************************************************/
void kd_pdf_reggrid(const kernel_density *kd, const double *xfixed,
		    const unsigned long *dimfixed, 
		    const unsigned long ndimfixed,
		    const unsigned long nfixed,
		    const double *xgridlo,
		    const double *xgridhi,
		    const unsigned long *ngrid,
		    const unsigned long *dimgrid,
		    const unsigned long ndimgrid,
		    const double reltol, const double abstol,
		    double *pdf)
{
  unsigned long i, j;
  unsigned long nalloc=0, nnode, ptr, lchild, rchild;
  unsigned long ngridtot, nstenciltot, fixedptr, outptr;
  unsigned long *nodelist = NULL, *nstencil;
  double *dxgrid;
  double *nodepdf = NULL;
  double relerr = 0.0, maxerr;
  double pdfmax, leftpdf, rightpdf;

  /* Figure out the total number of points in the grid */
  ngridtot = ngrid[0];
  for (i=1; i<ndimgrid; i++) ngridtot *= ngrid[i];

  /* Figure out the grid spacing */
  if (!(dxgrid = (double *) 
	calloc(ndimgrid, sizeof(double)))) {
    fprintf(stderr, 
	    "bayesphot: error: unable to allocate memory in kd_pdf_reggrid\n");
    exit(1);
  }
  for (i=0; i<ndimgrid; i++) {
    if (ngrid[i] > 1) 
      dxgrid[i] = (xgridhi[i]-xgridlo[i]) / (ngrid[i]-1);
    else
      dxgrid[i] = 0.0;
  }

  /* Figure out how big a stencil we need */
  if (!(nstencil = (unsigned long *) 
	calloc(ndimgrid, sizeof(unsigned long)))) {
    fprintf(stderr, 
	    "bayesphot: error: unable to allocate memory in kd_pdf_reggrid\n");
    exit(1);
  }
  nstenciltot = 1;
  for (i=0; i<ndimgrid; i++) {
    if (ngrid[i] > 1) {
      switch (kd->ktype) {
      case epanechnikov: {
	nstencil[i] = (int) 
	  floor((kd->h[dimgrid[i]]-dxgrid[i]/2.0)/dxgrid[i]);
	break;
      }
      case tophat: {
	nstencil[i] = (int) 
	  floor((kd->h[dimgrid[i]]-dxgrid[i]/2.0)/dxgrid[i]);
	break;
      }
      case gaussian: {
	nstencil[i] = (int) 
	  floor(0.5 * (1.0+sqrt(1.0 - 8.0*log(reltol) * 
			       (kd->h[dimgrid[i]]/dxgrid[i]) *
			       (kd->h[dimgrid[i]]/dxgrid[i]))));
	break;
      }
      }
      if (nstencil[i] >= ngrid[i]) nstencil[i] = ngrid[i]-1;
    } else {
      nstencil[i] = 0;
    }
    nstenciltot *= (2*nstencil[i]+1);
  }

  /* Loop over the input fixed points */
  for (i=0; i<nfixed; i++) {

    /* Move the pointers in the input and output arrays */
    fixedptr = i*ndimfixed;
    outptr = i*ngridtot;

    /* Initialize output array to zero */
    for (j=0; j<ngridtot; j++) pdf[i*ngridtot+j] = 0.0;

    /* Analyze root node */
    pdfmax = kd_pdf_node_reggrid(kd, xfixed+fixedptr, dimfixed, 
				 ndimfixed, xgridlo, dxgrid, ngrid,
				 nstencil, ngridtot, nstenciltot, dimgrid, 
				 ndimgrid, ROOT, pdf+outptr);

    /* If root node is a leaf, we're done. Go to next fixed point. */
    if (kd->tree->tree[ROOT].splitdim == -1) continue;

    /* More common case where root node is not a leaf */

    /* Allocate memory for node list if necessary */
    if (nalloc == 0) {
      if (!(nodelist = (unsigned long *) 
	    calloc(NODEBLOCKSIZE, sizeof(unsigned long)))) {
	fprintf(stderr, "bayesphot: error: unable to allocate memory in kd_pdf_grid\n");
	exit(1);
      }
      if (!(nodepdf = (double *) 
	    calloc(NODEBLOCKSIZE, sizeof(double)))) {
	fprintf(stderr, "bayesphot: error: unable to allocate memory in kd_pdf_grid\n");
	exit(1);
      }
      nalloc = NODEBLOCKSIZE;  
    }

    /* Push root node onto the node list */
    nodelist[0] = ROOT;
    nodepdf[0] = pdfmax;
    nnode = 1;

    /* Now proceed through the tree, identifying the nodes that are
       contributing the most error and recursively opening them. */
    while (1) {

      /* Find the node that is contributing the most error. */
      maxerr = nodepdf[0];
      ptr = 0;
      for (j=1; j<nnode; j++) {
	if (nodepdf[j] > maxerr) {
	  ptr = j;
	  maxerr = nodepdf[j];
	}
      }

      /* Subtract this node's contribution to the maximum possible PDF
	 contribution from unopened nodes. Enforce positivity to avoid
	 spurious negative results coming from roundoff error. */
      pdfmax -= nodepdf[ptr];
      if (pdfmax < 0) pdfmax = 0.0;

      /* Compute estimates for this node's children */
      lchild = LEFT(nodelist[ptr]);
      rchild = RIGHT(nodelist[ptr]);
      leftpdf = kd_pdf_node_reggrid(kd, xfixed+fixedptr, dimfixed, 
				    ndimfixed, xgridlo, dxgrid, ngrid,
				    nstencil, ngridtot, nstenciltot, 
				    dimgrid, ndimgrid, lchild, pdf+outptr);
      rightpdf = kd_pdf_node_reggrid(kd, xfixed+fixedptr, dimfixed, 
				     ndimfixed, xgridlo, dxgrid, ngrid,
				     nstencil, ngridtot, nstenciltot, 
				     dimgrid, ndimgrid, rchild, pdf+outptr);

      /* Get new maximum contribution to PDF from unopened nodes */
      pdfmax += leftpdf + rightpdf;

      /* Check for termination on absolute error. Note that we can
	 save a factor of 2 by adding half the upper limit on the
	 contribution to every point. */
      if (0.5*pdfmax < abstol) {
	for (j=0; j<ngridtot; j++) pdf[outptr+j] += 0.5*pdfmax;
	break;
      }

      /* Check for termination on relative error. Again, we save a
	 factor of 2 by adding half the upper limit on to each point. */
      for (j=0; j<ngridtot; j++) {
	relerr = 0.5*pdfmax / (pdf[outptr+j] + 0.5*pdfmax + DBL_MIN);
	if (relerr >= reltol) break;
      }
      if (relerr < reltol) {
	for (j=0; j<ngridtot; j++) pdf[outptr+j] += 0.5*pdfmax;
	break;
      }

      /* If we're here, we haven't converged yet. Remove the node we
	 just analyzed from the list. */
      for (j=ptr; j<nnode-1; j++) {
	nodelist[j] = nodelist[j+1];
	nodepdf[j] = nodepdf[j+1];
      }
      nnode--;

      /* Allocate more memory to hold child nodes if necessary */
      if (nnode+2 >= nalloc) {
	if (!(nodelist = (unsigned long *) 
	      realloc(nodelist, 2*nalloc*sizeof(unsigned long)))) {
	  fprintf(stderr, "bayesphot: error: unable to allocate memory in kd_pdf_grid\n");
	  exit(1);
	}
	if (!(nodepdf = (double *) 
	      realloc(nodepdf, 2*nalloc*sizeof(double)))) {
	  fprintf(stderr, "bayesphot: error: unable to allocate memory in kd_pdf_grid\n");
	  exit(1);
	}
	nalloc *= 2;
      }

      /* If children are not leaves, and their maximum contribution to
	 the PDF is not 0 (possible due to roundoff, or if the kernel
	 is compact), push them onto the node list */
      if (leftpdf > 0.0) {
	nodelist[nnode] = lchild;
	nodepdf[nnode] = leftpdf;
	nnode++;
      }
      if (rightpdf > 0.0) {
	nodelist[nnode] = rchild;
	nodepdf[nnode] = rightpdf;
	nnode++;
      }

      /* Bail out if there are no non-leaf nodes left to be
	 analyzed. This should also give abserr = 0, and thus
	 we should have exited the loop when we checked the termination
	 condition. However, this can sometimes fail due to roundoff
	 error if abstol and reltol are both very small, so this is a
	 backstop. */
      if (nnode == 0) break;
    }
  }

  /* Free memory */
  if (nodelist != NULL) free(nodelist);
  if (nodepdf != NULL) free(nodepdf);
  free(dxgrid);
  free(nstencil);
}


/*********************************************************************/
/* Function to estimate the contribution to a PDF from a node, and   */
/* the error on it; if the input node is a leaf, the routine         */
/* computes the exact PDF contribution, and sets the error to zero.  */
/* If it is not a leaf, the routine computes the minimum and maximum */
/* possible contributions from the node, returns their average as    */
/* the central estimate, and half their difference as the maximum    */
/* possible error.                                                   */
/*********************************************************************/
inline
void kd_pdf_node(const kernel_density *kd, const double *x,
		 const unsigned long curnode, double *pdf,
		 double *pdferr) {
  unsigned long i;
  unsigned long ndim = kd->tree->ndim;
  double d2, pdfmin = 0.0, pdfmax = 0.0;

  /* Is this node a leaf? If so, just sum over it and return that */
  if (kd->tree->tree[curnode].splitdim == -1) {

    /* Loop over points, summing their contribution */
    *pdf = 0.0;
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
	    *pdf += ((double *) kd->tree->tree[curnode].dptr)[i] * 
	      (1.0 - d2);
	  break;
	}
	case tophat: {
	  if (d2 < 1)
	    *pdf += ((double *) kd->tree->tree[curnode].dptr)[i];
	  break;
	}
	case gaussian: {
	  *pdf += ((double *) kd->tree->tree[curnode].dptr)[i] *
	    exp(-d2/2.0);
	  break;
	}
	}
      } else {
	switch (kd->ktype) {
	case epanechnikov: {
	  if (d2 < 1) *pdf += 1.0 - d2;
	  break;
	}
	case tophat: {
	  if (d2 < 1) *pdf += 1.0;
	  break;
	}
	case gaussian: {
	  *pdf += exp(-d2/2.0);
	  break;
	}
	}
      }

    }

    /* Normalize pdf, set error to zero */
    *pdf *= kd->norm_tot;
    *pdferr = 0.0;

  } else {

    /* This node is not a leaf, so set the PDF to the average of the
       minimum and maximum possible contributions, and set the error */

    /* Get minimum distance in units of the kernel size */
    d2 = box_min_dist2(x, 
		       (const double **) kd->tree->tree[curnode].xbnd,
		       ndim, ndim, NULL, NULL, kd->h, ndim);

    /* Get PDF evaluated at minimum distance */
    switch (kd->ktype) {
    case epanechnikov: {
      pdfmax = d2 < 1.0 ? kd->nodewgt[curnode] * kd->norm_tot 
	* (1.0 - d2) : 0.0;
      break;
    }
    case tophat: {
      pdfmax = d2 < 1.0 ? kd->nodewgt[curnode] * kd->norm_tot : 0.0;
      break;
    }
    case gaussian: {
      pdfmax = kd->nodewgt[curnode] * kd->norm_tot * exp(-d2/2.0);
      break;
    }
    }

    /* Get maximum distance in units of the kernel size, and PDF
       evaluated at that distance */
    d2 = box_max_dist2(x, 
		       (const double **) kd->tree->tree[curnode].xbnd,
		       ndim, ndim, NULL, NULL, kd->h, ndim);
    switch (kd->ktype) {
    case epanechnikov: {
      pdfmin = d2 < 1.0 ? kd->nodewgt[curnode] * kd->norm_tot 
	* (1.0 - d2) : 0.0;
      break;
    }
    case tophat: {
      pdfmin = d2 < 1.0 ? kd->nodewgt[curnode] * kd->norm_tot : 0.0;
      break;
    }
    case gaussian: {
      pdfmin = kd->nodewgt[curnode] * kd->norm_tot * exp(-d2/2.0);
      break;
    }
    }

    /* Set estimate to average of min and max, and error to half the
       difference between them. */
    *pdf = (pdfmin + pdfmax) / 2.0;
    *pdferr = (pdfmax - pdfmin) / 2.0;

  }
}


/*********************************************************************/
/* Function to compute the contribution of a node to the PDF on a    */
/* grid of input points, with a single fixed point. If the node is   */
/* not a leaf, we return the maximum possible contribution to any    */
/* point on the grid, without examining the grid points individualy. */
/* If it is a leaf, we add the contribution to each point on the     */
/* grid, and return 0.                                               */
/*********************************************************************/
inline
double kd_pdf_node_grid(const kernel_density *kd, const double *xfixed, 
			const unsigned long *dimfixed,
			const unsigned long ndimfixed,
			const double *xgrid,
			const unsigned long *dimgrid,
			const unsigned long ndimgrid,
			const unsigned long ngrid,
			const unsigned long curnode,
			double *pdf) {

  unsigned long i, j;
  unsigned long ndim_tot = kd->tree->ndim;
  double d2, d2fixed, d2grid;

  /* Is this node a leaf? If so, just sum over the points in it */
  if (kd->tree->tree[curnode].splitdim == -1) {

    /* Loop over points in this leaf */
    for (i=0; i<kd->tree->tree[curnode].npt; i++) {

      /* Get distance in the fixed dimensions */
      d2fixed = dist2(&(kd->tree->tree[curnode].x[ndim_tot*i]),
		      xfixed, ndim_tot, ndimfixed,
		      NULL, dimfixed, kd->h, ndim_tot);

      /* Loop over the points in the grid */
      for (j=0; j<ngrid; j++) {

	/* Get distance in the grid dimensions */
	d2grid = dist2(&(kd->tree->tree[curnode].x[ndim_tot*i]),
		       xgrid+j*ndimgrid, ndim_tot, ndimgrid,
		       NULL, dimgrid, kd->h, ndim_tot);

	/* Sum distances and get contribution to the PDF of this
	   point */
	d2 = d2fixed + d2grid;
	if (kd->tree->tree[curnode].dptr != NULL) {
	  switch (kd->ktype) {
	  case epanechnikov: {
	    if (d2 < 1)
	      pdf[j] += ((double *) kd->tree->tree[curnode].dptr)[i] * 
		(1.0 - d2);
	    break;
	  }
	  case tophat: {
	    if (d2 < 1)
	      pdf[j] += ((double *) kd->tree->tree[curnode].dptr)[i];
	    break;
	  }
	  case gaussian: {
	    pdf[j] += ((double *) kd->tree->tree[curnode].dptr)[i] *
	      exp(-d2/2.0);
	    break;
	  }
	  }
	} else {
	  switch (kd->ktype) {
	  case epanechnikov: {
	    if (d2 < 1) pdf[j] += 1.0 - d2;
	    break;
	  }
	  case tophat: {
	    if (d2 < 1) pdf[j] += 1.0;
	    break;
	  }
	  case gaussian: {
	    pdf[j] += exp(-d2/2.0);
	    break;
	  }
	  }
	}
      }
    }

    /* Return */
    return(0.0);

  } else {

    /* This node is not a leaf */

    /* Get minimum distance in units of the kernel size */
    d2 = box_min_dist2(xfixed, 
		       (const double **) kd->tree->tree[curnode].xbnd,
		       ndimfixed, ndim_tot, dimfixed, NULL, 
		       kd->h, ndim_tot);

    /* Return 0.5 * PDF at minimum distance */
    switch (kd->ktype) {
    case epanechnikov: {
      return(d2 < 1.0 ? kd->nodewgt[curnode] * kd->norm_tot 
	     * (1.0 - d2) : 0.0);
    }
    case tophat: {
      return(d2 < 1.0 ? kd->nodewgt[curnode] * kd->norm_tot : 0.0);
    }
    case gaussian: {
      return(kd->nodewgt[curnode] * kd->norm_tot * exp(-d2/2.0));
    }
    }
    
    /*  Never get here -- exists to suppress compiler warning */
    return(0.0);
  }
}


/*********************************************************************/
/* Function to estimate the contribution to an integrated PDF from a */
/* node. Functionality is identical to kdf_pdf_node, except that     */
/* certain dimensions are being integrated out.                      */
/*********************************************************************/
inline
void kd_pdf_node_int(const kernel_density *kd, const double *x,
		     const unsigned long *dims, const unsigned long ndim,
		     const unsigned long ndim_int, const double fac,
		     const unsigned long curnode, double *pdf,
		     double *pdferr) {
  unsigned long i;
  unsigned long ndim_tot = kd->tree->ndim;
  double d2, pdfmin = 0.0, pdfmax = 0.0;

  /* Is this node a leaf? If so, just sum over it and return that */
  if (kd->tree->tree[curnode].splitdim == -1) {

    /* Loop over points, summing their contribution */
    *pdf = 0.0;
    for (i=0; i<kd->tree->tree[curnode].npt; i++) {

      /* Get distance in units of the kernel size */
      d2 = dist2(&(kd->tree->tree[curnode].x[ndim_tot*i]), x,
		 ndim_tot, ndim, NULL, dims, kd->h, ndim_tot);

      /* Now compute contribution from each point */
      if (kd->tree->tree[curnode].dptr != NULL) {
	switch (kd->ktype) {
	case epanechnikov: {
	  if (d2 < 1)
	    *pdf += ((double *) kd->tree->tree[curnode].dptr)[i] * 
	      pow(1.0-d2, 1.0+0.5*ndim_int);
	  break;
	}
	case tophat: {
	  if (d2 < 1)
	    *pdf += ((double *) kd->tree->tree[curnode].dptr)[i] *
	      pow(1.0-d2, 0.5*ndim_int);
	  break;
	}
	case gaussian: {
	  *pdf += ((double *) kd->tree->tree[curnode].dptr)[i] *
	    exp(-d2/2.0);
	  break;
	}
	}
      } else {
	switch (kd->ktype) {
	case epanechnikov: {
	  if (d2 < 1) *pdf += pow(1.0-d2, 1.0+0.5*ndim_int);
	  break;
	}
	case tophat: {
	  if (d2 < 1) *pdf += pow(1.0-d2, 0.5*ndim_int);
	  break;
	}
	case gaussian: {
	  *pdf += exp(-d2/2.0);
	  break;
	}
	}
      }

    }

    /* Normalize and return */
    *pdf *= fac * kd->norm_tot;
    *pdferr = 0.0;

  } else {

    /* This node is not a leaf, so return half the maximum possible
       contribution to the PDF */

    /* Get minimum distance in units of the kernel size */
    d2 = box_min_dist2(x, 
		       (const double **) kd->tree->tree[curnode].xbnd,
		       ndim, ndim_tot, dims, NULL, kd->h, ndim_tot);

    /* Get PDF evaluated at minimum distance */
    switch (kd->ktype) {
    case epanechnikov: {
      pdfmax = d2 < 1.0 ? kd->nodewgt[curnode] * kd->norm_tot * fac * 
	pow(1.0-d2, 1.0+0.5*ndim_int) : 0.0;
      break;
    }
    case tophat: {
      pdfmax = d2 < 1.0 ? kd->nodewgt[curnode] * kd->norm_tot * fac * 
	pow(1.0-d2, 0.5*ndim_int) : 0.0;
      break;
    }
    case gaussian: {
      pdfmax = kd->nodewgt[curnode] * kd->norm_tot * fac * exp(-d2/2.0);
      break;
    }
    }

      /* Get maximum distance in units of the kernel size */
    d2 = box_max_dist2(x, 
		       (const double **) kd->tree->tree[curnode].xbnd,
		       ndim, ndim_tot, dims, NULL, kd->h, ndim_tot);

    /* Get PDF evaluated at maximum distance */
    switch (kd->ktype) {
    case epanechnikov: {
      pdfmin = d2 < 1.0 ? kd->nodewgt[curnode] * kd->norm_tot * fac * 
	pow(1.0-d2, 1.0+0.5*ndim_int) : 0.0;
      break;
    }
    case tophat: {
      pdfmin = d2 < 1.0 ? kd->nodewgt[curnode] * kd->norm_tot * fac * 
	pow(1.0-d2, 0.5*ndim_int) : 0.0;
      break;
    }
    case gaussian: {
      pdfmin = kd->nodewgt[curnode] * kd->norm_tot * fac * exp(-d2/2.0);
      break;
    }
    }

    /* Set estimate to average of min and max, and error to half the
       difference between them. */
    *pdf = (pdfmin + pdfmax) / 2.0;
    *pdferr = (pdfmax - pdfmin) / 2.0;
  }
}



/*********************************************************************/
/* Function to compute the contribution of a node to the PDF on a    */
/* grid of input points, with a single fixed point, integrating out  */
/* some of the dimensions. If the node is not a leaf, we return the  */
/* maximum possible contribution to any point on the grid, without   */
/* examining the grid points individualy. If it is a leaf, we add    */
/* the contribution to each point on the grid, and return 0.         */
/*********************************************************************/
inline
double kd_pdf_node_int_grid(const kernel_density *kd, 
			    const double *xfixed, 
			    const unsigned long *dimfixed,
			    const unsigned long ndimfixed,
			    const double *xgrid,
			    const unsigned long *dimgrid,
			    const unsigned long ndimgrid,
			    const unsigned long ngrid,
			    const unsigned long ndim_int,
			    const double fac,
			    const unsigned long curnode,
			    double *pdf) {

  unsigned long i, j;
  unsigned long ndim_tot = kd->tree->ndim;
  double d2, d2fixed, d2grid;

  /* Is this node a leaf? If so, just sum over the points in it */
  if (kd->tree->tree[curnode].splitdim == -1) {

    /* Loop over points in this leaf */
    for (i=0; i<kd->tree->tree[curnode].npt; i++) {

      /* Get distance in the fixed dimensions */
      d2fixed = dist2(&(kd->tree->tree[curnode].x[ndim_tot*i]),
		      xfixed, ndim_tot, ndimfixed,
		      NULL, dimfixed, kd->h, ndim_tot);

      /* Loop over the points in the grid */
      for (j=0; j<ngrid; j++) {

	/* Get distance in the grid dimensions */
	d2grid = dist2(&(kd->tree->tree[curnode].x[ndim_tot*i]),
		       xgrid+j*ndimgrid, ndim_tot, ndimgrid,
		       NULL, dimgrid, kd->h, ndim_tot);

	/* Sum distances and get contribution to the PDF of this
	   point */
	d2 = d2fixed + d2grid;
	if (kd->tree->tree[curnode].dptr != NULL) {
	  switch (kd->ktype) {
	  case epanechnikov: {
	    if (d2 < 1)
	      pdf[j] += ((double *) kd->tree->tree[curnode].dptr)[i] * 
		pow(1.0-d2, 1.0+0.5*ndim_int) * fac * kd->norm_tot;
	    break;
	  }
	  case tophat: {
	    if (d2 < 1)
	      pdf[j] += ((double *) kd->tree->tree[curnode].dptr)[i] *
		pow(1.0-d2, 0.5*ndim_int) * fac * kd->norm_tot;
	    break;
	  }
	  case gaussian: {
	    pdf[j] += ((double *) kd->tree->tree[curnode].dptr)[i] *
	      exp(-d2/2.0) * fac * kd->norm_tot;
	    break;
	  }
	  }
	} else {
	  switch (kd->ktype) {
	  case epanechnikov: {
	    if (d2 < 1) 
	      pdf[j] += pow(1.0-d2, 1.0+0.5*ndim_int) * fac * kd->norm_tot;
	    break;
	  }
	  case tophat: {
	    if (d2 < 1) 
	      pdf[j] += pow(1.0-d2, 0.5*ndim_int) * fac * kd->norm_tot;
	    break;
	  }
	  case gaussian: {
	    pdf[j] += exp(-d2/2.0) * fac * kd->norm_tot;
	    break;
	  }
	  }
	}
      }
    }

    /* Return */
    return(0.0);

  } else {

    /* This node is not a leaf */

    /* Get minimum distance in units of the kernel size */
    d2 = box_min_dist2(xfixed, 
		       (const double **) kd->tree->tree[curnode].xbnd,
		       ndimfixed, ndim_tot, dimfixed, NULL, 
		       kd->h, ndim_tot);

    /* Return PDF at minimum distance */
    switch (kd->ktype) {
    case epanechnikov: {
      return(d2 < 1.0 ? kd->nodewgt[curnode] * kd->norm_tot * fac * 
	      pow(1.0-d2, 1.0+0.5*ndim_int) : 0.0);
    }
    case tophat: {
      return(d2 < 1.0 ? kd->nodewgt[curnode] * kd->norm_tot * fac * 
	     pow(1.0-d2, 0.5*ndim_int) : 0.0);
    }
    case gaussian: {
      return(kd->nodewgt[curnode] * kd->norm_tot * fac * exp(-d2/2.0));
    }
    }
    
    /*  Never get here -- exists to suppress compiler warning */
    return(0.0);
  }
}

/*********************************************************************/
/* Same as kd_pdf_node_int_grid, but for regular grids               */
/*********************************************************************/
inline
double kd_pdf_node_int_reggrid(const kernel_density *kd, 
			       const double *xfixed, 
			       const unsigned long *dimfixed,
			       const unsigned long ndimfixed,
			       const double *xgridlo,
			       const double *dxgrid,
			       const unsigned long *ngrid,
			       const unsigned long *nstencil,
			       const unsigned long ngridtot,
			       const unsigned long nstenciltot,
			       const unsigned long *dimgrid,
			       const unsigned long ndimgrid,
			       const unsigned long ndim_int,
			       const double fac,
			       const unsigned long curnode,
			       double *pdf) {

  unsigned long i, j;
  long k;
  long *ctr, *offset;
  long pdf_offset;
  long on_grid;
  unsigned long ndim_tot = kd->tree->ndim;
  double d2, d2fixed;
  double *d2grid;

  /* Is this node a leaf? If so, just sum over the points in it */
  if (kd->tree->tree[curnode].splitdim == -1) {

    /* Allocate memory we'll need */
    if (!(ctr = (long *) calloc(ndimgrid, sizeof(long)))) {
      fprintf(stderr, 
	      "bayesphot: error: unable to allocate memory in kd_pdf_int_reggrid\n");
      exit(1);
    }
    if (!(offset = (long *) calloc(ndimgrid, sizeof(long)))) {
      fprintf(stderr, 
	      "bayesphot: error: unable to allocate memory in kd_pdf_int_reggrid\n");
      exit(1);
    }
    if (!(d2grid = (double *) calloc(ndimgrid, sizeof(double)))) {
      fprintf(stderr, 
	      "bayesphot: error: unable to allocate memory in kd_pdf_int_reggrid\n");
      exit(1);
    }

    /* Loop over points in this leaf */
    for (i=0; i<kd->tree->tree[curnode].npt; i++) {

      /* Get distance in the fixed dimensions */
      d2fixed = dist2(&(kd->tree->tree[curnode].x[ndim_tot*i]),
		      xfixed, ndim_tot, ndimfixed,
		      NULL, dimfixed, kd->h, ndim_tot);

      /* Initialize the offset array */
      for (j=0; j<ndimgrid; j++) offset[j] = -nstencil[j];

      /* Figure out which grid cell this point lands in */
      for (j=0; j<ndimgrid; j++)
	ctr[j] = (int) 
	  round((kd->tree->tree[curnode].x[ndim_tot*i+dimgrid[j]] - 
		 xgridlo[j]) / dxgrid[j]);

      /* Initialize the d2grid array at the lower left corner of the
	 stencil grid */
      for (j=0; j<ndimgrid; j++)
	d2grid[j] = 
	  ((kd->tree->tree[curnode].x[ndim_tot*i+dimgrid[j]] -
	    (xgridlo[j] + (ctr[j]+offset[j])*dxgrid[j]))
	   / kd->h[dimgrid[j]]) *
	  ((kd->tree->tree[curnode].x[ndim_tot*i+dimgrid[j]] -
	    (xgridlo[j] + (ctr[j]+offset[j])*dxgrid[j]))
	   / kd->h[dimgrid[j]]);

      /* Loop over the stencil */
      for (j=0; j<nstenciltot; j++) {

	/* Make sure we're on the grid; if not, skip to updating the
	   offsets */
	on_grid = 1;
	for (k=0; k<ndimgrid; k++)
	  on_grid = on_grid & (ctr[k] + offset[k] >= 0)
	    & (ctr[k] + offset[k] < ngrid[k]);
	if (on_grid) {

	  /* Figure out where in the PDF array this point in the stencil
	     corresponds to */
	  pdf_offset = ctr[0] + offset[0];
	  for (k=1; k<ndimgrid; k++) {
	    pdf_offset *= ngrid[k];
	    pdf_offset += ctr[k] + offset[k];
	  }

	  /* Get total distance to this point */
	  d2 = d2fixed;
	  for (k=0; k<ndimgrid; k++) d2 += d2grid[k];

	  /* Add the appropriate contribution to the PDF */
	  if (kd->tree->tree[curnode].dptr != NULL) {
	    switch (kd->ktype) {
	    case epanechnikov: {
	      if (d2 < 1)
		pdf[pdf_offset] += ((double *) 
				    kd->tree->tree[curnode].dptr)[i] * 
		  pow(1.0-d2, 1.0+0.5*ndim_int) * fac * kd->norm_tot;
	      break;
	    }
	    case tophat: {
	      if (d2 < 1)
		pdf[pdf_offset] += 
		  ((double *) kd->tree->tree[curnode].dptr)[i] *
		  pow(1.0-d2, 0.5*ndim_int) * fac * kd->norm_tot;
	      break;
	    }
	    case gaussian: {
	      pdf[pdf_offset] += 
		((double *) kd->tree->tree[curnode].dptr)[i] *
		exp(-d2/2.0) * fac * kd->norm_tot;
	      break;
	    }
	    }
	  } else {
	    switch (kd->ktype) {
	    case epanechnikov: {
	      if (d2 < 1) 
		pdf[pdf_offset] += 
		  pow(1.0-d2, 1.0+0.5*ndim_int) * fac * kd->norm_tot;
	      break;
	    }
	    case tophat: {
	      if (d2 < 1) 
		pdf[pdf_offset] += 
		  pow(1.0-d2, 0.5*ndim_int) * fac * kd->norm_tot;
	      break;
	    }
	    case gaussian: {
	      pdf[pdf_offset] += exp(-d2/2.0) * fac * kd->norm_tot;
	      break;
	    }
	    }
	  }
	}

	/* Update the offsets and distances */
	for (k=ndimgrid-1; k>=0; k--) {
	  offset[k]++;
	  if (offset[k] > (int) nstencil[k]) offset[k] = -(int) nstencil[k];
	  d2grid[k] = 
	    ((kd->tree->tree[curnode].x[ndim_tot*i+dimgrid[k]] -
	      (xgridlo[k] + (ctr[k]+offset[k])*dxgrid[k]))
	     / kd->h[dimgrid[k]]) *
	    ((kd->tree->tree[curnode].x[ndim_tot*i+dimgrid[k]] -
	      (xgridlo[k] + (ctr[k]+offset[k])*dxgrid[k]))
	     / kd->h[dimgrid[k]]);
	  if (offset[k] != -(int) nstencil[k]) break;
	}
      }
    }

    /* Free memory */
    free(offset);
    free(ctr);
    free(d2grid);

    /* Return */
    return(0.0);

  } else {

    /* This node is not a leaf */

    /* Get minimum distance in units of the kernel size */
    d2 = box_min_dist2(xfixed, 
		       (const double **) kd->tree->tree[curnode].xbnd,
		       ndimfixed, ndim_tot, dimfixed, NULL, 
		       kd->h, ndim_tot);

    /* Return PDF at minimum distance */
    switch (kd->ktype) {
    case epanechnikov: {
      return(d2 < 1.0 ? kd->nodewgt[curnode] * kd->norm_tot * fac * 
	      pow(1.0-d2, 1.0+0.5*ndim_int) : 0.0);
    }
    case tophat: {
      return(d2 < 1.0 ? kd->nodewgt[curnode] * kd->norm_tot * fac * 
	     pow(1.0-d2, 0.5*ndim_int) : 0.0);
    }
    case gaussian: {
      return(kd->nodewgt[curnode] * kd->norm_tot * fac * exp(-d2/2.0));
    }
    }
    
    /*  Never get here -- exists to suppress compiler warning */
    return(0.0);
  }
}


/*********************************************************************/
/* Function to compute the contribution of a node to the PDF on a    */
/* regular grid of input points, with a single fixed point. Same as  */
/* kd_pdf_node_grid, except that the grid is required to be regular. */
/*********************************************************************/
inline
double kd_pdf_node_reggrid(const kernel_density *kd, 
			   const double *xfixed, 
			   const unsigned long *dimfixed,
			   const unsigned long ndimfixed,
			   const double *xgridlo,
			   const double *dxgrid,
			   const unsigned long *ngrid,
			   const unsigned long *nstencil,
			   const unsigned long ngridtot,
			   const unsigned long nstenciltot,
			   const unsigned long *dimgrid,
			   const unsigned long ndimgrid,
			   const unsigned long curnode,
			   double *pdf) {

  unsigned long i, j;
  long k;
  long *ctr, *offset;
  int on_grid;
  long pdf_offset;
  unsigned long ndim_tot = kd->tree->ndim;
  double d2, d2fixed;
  double *d2grid;

  /* Is this node a leaf? If so, just sum over the points in it */
  if (kd->tree->tree[curnode].splitdim == -1) {

    /* Allocate memory we'll need */
    if (!(ctr = (long *) calloc(ndimgrid, sizeof(long)))) {
      fprintf(stderr, 
	      "bayesphot: error: unable to allocate memory in kd_pdf_reggrid\n");
      exit(1);
    }
    if (!(offset = (long *) calloc(ndimgrid, sizeof(long)))) {
      fprintf(stderr, 
	      "bayesphot: error: unable to allocate memory in kd_pdf_reggrid\n");
      exit(1);
    }
    if (!(d2grid = (double *) calloc(ndimgrid, sizeof(double)))) {
      fprintf(stderr, 
	      "bayesphot: error: unable to allocate memory in kd_pdf_reggrid\n");
      exit(1);
    }

    /* Loop over points in this leaf */
    for (i=0; i<kd->tree->tree[curnode].npt; i++) {

      /* Get distance in the fixed dimensions */
      d2fixed = dist2(&(kd->tree->tree[curnode].x[ndim_tot*i]),
		      xfixed, ndim_tot, ndimfixed,
		      NULL, dimfixed, kd->h, ndim_tot);

      /* Initialize the offset array */
      for (j=0; j<ndimgrid; j++) offset[j] = -nstencil[j];

      /* Figure out which grid cell this point lands in */
      for (j=0; j<ndimgrid; j++)
	ctr[j] = (int) 
	  round((kd->tree->tree[curnode].x[ndim_tot*i+dimgrid[j]] - 
		 xgridlo[j]) / dxgrid[j]);

      /* Initialize the d2grid array at the lower left corner of the
	 stencil grid */
      for (j=0; j<ndimgrid; j++)
	d2grid[j] = 
	  ((kd->tree->tree[curnode].x[ndim_tot*i+dimgrid[j]] -
	    (xgridlo[j] + (ctr[j]+offset[j])*dxgrid[j]))
	   / kd->h[dimgrid[j]]) *
	  ((kd->tree->tree[curnode].x[ndim_tot*i+dimgrid[j]] -
	    (xgridlo[j] + (ctr[j]+offset[j])*dxgrid[j]))
	   / kd->h[dimgrid[j]]);

      /* Loop over the stencil */
      for (j=0; j<nstenciltot; j++) {

	/* Make sure we're on the grid; if not, skip to updating the
	   offsets */
	on_grid = 1;
	for (k=0; k<ndimgrid; k++)
	  on_grid = on_grid & (ctr[k] + offset[k] >= 0)
	    & (ctr[k] + offset[k] < ngrid[k]);
	if (on_grid) {

	  /* Figure out where in the PDF array this point in the stencil
	     corresponds to */
	  pdf_offset = ctr[0] + offset[0];
	  for (k=1; k<ndimgrid; k++) {
	    pdf_offset *= ngrid[k];
	    pdf_offset += ctr[k] + offset[k];
	  }

	  /* Get total distance to this point */
	  d2 = d2fixed;
	  for (k=0; k<ndimgrid; k++) d2 += d2grid[k];

	  /* Add the appropriate contribution to the PDF */
	  if (kd->tree->tree[curnode].dptr != NULL) {
	    switch (kd->ktype) {
	    case epanechnikov: {
	      if (d2 < 1)
		pdf[pdf_offset] += 
		  ((double *) kd->tree->tree[curnode].dptr)[i] * 
		  (1.0 - d2);
	    break;
	    }
	    case tophat: {
	      if (d2 < 1)
		pdf[pdf_offset] += 
		  ((double *) kd->tree->tree[curnode].dptr)[i];
	      break;
	    }
	    case gaussian: {
	      pdf[pdf_offset] += 
		((double *) kd->tree->tree[curnode].dptr)[i] *
		exp(-d2/2.0);
	      break;
	    }
	    }
	  } else {
	    switch (kd->ktype) {
	    case epanechnikov: {
	      if (d2 < 1) pdf[pdf_offset] += 1.0 - d2;
	      break;
	    }
	    case tophat: {
	      if (d2 < 1) pdf[pdf_offset] += 1.0;
	      break;
	    }
	    case gaussian: {
	      pdf[pdf_offset] += exp(-d2/2.0);
	      break;
	    }
	    }
	  }
	}

	/* Update the offsets and distances */
	for (k=ndimgrid-1; k>=0; k--) {
	  offset[k]++;
	  if (offset[k] > (int) nstencil[k]) offset[k] = -(int) nstencil[k];
	  d2grid[k] = 
	    ((kd->tree->tree[curnode].x[ndim_tot*i+dimgrid[k]] -
	      (xgridlo[k] + (ctr[k]+offset[k])*dxgrid[k]))
	     / kd->h[dimgrid[k]]) *
	    ((kd->tree->tree[curnode].x[ndim_tot*i+dimgrid[k]] -
	      (xgridlo[k] + (ctr[k]+offset[k])*dxgrid[k]))
	     / kd->h[dimgrid[k]]);
	  if (offset[k] != -(int) nstencil[k]) break;
	}
      }
    }

    /* Free memory */
    free(offset);
    free(ctr);
    free(d2grid);

    /* Return */
    return(0.0);

  } else {

    /* This node is not a leaf */

    /* Get minimum distance in units of the kernel size */
    d2 = box_min_dist2(xfixed, 
		       (const double **) kd->tree->tree[curnode].xbnd,
		       ndimfixed, ndim_tot, dimfixed, NULL, 
		       kd->h, ndim_tot);
    
    /* Return 0.5 * PDF at minimum distance */
    switch (kd->ktype) {
    case epanechnikov: {
      return(d2 < 1.0 ? kd->nodewgt[curnode] * kd->norm_tot 
	     * (1.0 - d2) : 0.0);
    }
    case tophat: {
      return(d2 < 1.0 ? kd->nodewgt[curnode] * kd->norm_tot : 0.0);
    }
    case gaussian: {
      return(kd->nodewgt[curnode] * kd->norm_tot * exp(-d2/2.0));
    }
    }
    
    /*  Never get here -- exists to suppress compiler warning */
    return(0.0);
  }
}

/*********************************************************************/
/* Vectorized version of kd_pdf_int                                  */
/*********************************************************************/
void kd_pdf_int_vec(const kernel_density *kd,
		    const double *x,
		    const double *bandwidth,
		    const unsigned long *dims,
		    const unsigned long ndim,
		    const unsigned long npt,
		    const double reltol, 
		    const double abstol,
		    double *pdf
#ifdef DIAGNOSTIC
		    , unsigned long *nodecheck,
		    unsigned long *leafcheck,
		    unsigned long *termcheck
#endif
		    ) {

  unsigned long i, j;
  kernel_density kd_tmp;
  double *bw_tmp;
  
  if (!bandwidth) {

    /* Case with constant bandwidth */
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

  } else {
    
    /* Case with varying bandwidth. The only part that is tricky here
       is that, in parallel, each thread needs to have a version of
       the kernel_density object that is identical to the original
       except for the values of kd->h, kd->norm, and kd->norm_tot. To
       accomplish this, we first copy kd to a temporary variable
       kd_tmp outside the loop; we declare that this variable is
       private inside the loop. Then each time we go through the loop
       we point kd_tmp.h to a newly-allocated array, whose value we
       are free to change without interfering with the other threads. */
#pragma omp parallel for private(i, j, kd_tmp, bw_tmp)
    for (i=0; i<npt; i++) {
	
      /* Make our private copy of kd, and allocate memory for its
	 bandwidth array; also allocate memory for a temporary
	 bandwidth array, and initialize it. */
      kd_tmp = *kd;
      kd_tmp.h = calloc(kd->tree->ndim, sizeof(double));
      bw_tmp = calloc(kd->tree->ndim, sizeof(double));
      for (j=0; j<ndim; j++) bw_tmp[j] = kd->h[j];

      /* Change bandwidth in our private copy of kd; note that we
	 have to properly map the input bandwidth, which is missing
	 certain dimensions, to the dimensions that remain */
      for (j=0; j<ndim; j++) bw_tmp[dims[j]] = bandwidth[i*ndim+j];
      kd_change_bandwidth(bw_tmp, &kd_tmp);

      /* Evauluate the PDF with the modified bandwidth */
      pdf[i] = kd_pdf_int(&kd_tmp, x+i*ndim, dims, ndim, reltol, abstol
#ifdef DIAGNOSTIC
			  , nodecheck+i, leafcheck+i, termcheck+i
#endif
			  );

      /* Free the local storage we allocated */
      free(kd_tmp.h);
      free(bw_tmp);
    }
  }
}


/*********************************************************************/
/* Function to evaluate the PDF using a kernel_density object for a  */
/* vector of inputs positions                                        */
/*********************************************************************/
void kd_pdf_vec(const kernel_density *kd,
		const double *x,
		const double *bandwidth,
		const unsigned long npt,
		const double reltol,
		const double abstol,
		double *pdf
#ifdef DIAGNOSTIC
		, unsigned long *nodecheck,
		unsigned long *leafcheck,
		unsigned long *termcheck
#endif
		) {

  kernel_density kd_tmp;

  if (!bandwidth) {

    /* Case with constant bandwidth */
    /*#pragma omp parallel private(i)*/
    {
#pragma omp parallel for
      for (unsigned long i=0; i<npt; i++) {
	pdf[i] = kd_pdf(kd, x+i*kd->tree->ndim, reltol, abstol
#ifdef DIAGNOSTIC
			, nodecheck+i, leafcheck+i, termcheck+i
#endif
			);
      }
    }

  } else {

    /* Case with varying bandwidth. The only part that is tricky here
       is that, in parallel, each thread needs to have a version of
       the kernel_density object that is identical to the original
       except for the values of kd->h, kd->norm, and kd->norm_tot. To
       accomplish this, we maintain a privde kd_tmp variable, which we
       populate by copying from the global one, allocating memory for
       a private bandwidth copy, and then changing that private copy
       to have whatever bandwidth we need. */
#pragma omp parallel for private(kd_tmp)
    for (unsigned long i=0; i<npt; i++) {
	
      /* Make our private copy of kd, and allocate memory for its
	 bandwidth array */
      kd_tmp = *kd;
      kd_tmp.h = calloc(kd->tree->ndim, sizeof(double));

      /* Change bandwidth in our private copy of kd  */
      kd_change_bandwidth(bandwidth+i*kd->tree->ndim, &kd_tmp);

      /* Evauluate the PDF with the modified bandwidth */
      pdf[i] = kd_pdf(&kd_tmp, x+i*kd->tree->ndim, reltol, abstol
#ifdef DIAGNOSTIC
		      , nodecheck+i, leafcheck+i, termcheck+i
#endif
		      );
      /*
      printf("x =");
      for (unsigned long j=0; j<kd->tree->ndim; j++)
	printf(" %f", x[i*kd->tree->ndim+j]);
      printf(", pdf = %f\n", pdf[i]);
      */      
      /* Free the local storage we allocated */
      free(kd_tmp.h);
    }
  }
}

