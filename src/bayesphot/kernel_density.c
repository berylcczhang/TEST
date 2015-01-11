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
#include "kernel_density.h"

/*********************************************************************/
/* Useful macros                                                     */
/*********************************************************************/

#define SQR(x) ((x)*(x))
#define NODEBLOCKSIZE 1024

/*********************************************************************/
/* Static functions                                                  */
/*********************************************************************/

/* Function to compute the surface area element for an N-sphere */
static inline
double ds(unsigned int n) {
  unsigned int m;
  if (n % 2) {
    /* Case for n odd */
    m = n/2;
    return 2 * n * gsl_pow_uint(2.0*M_PI, m) / gsl_sf_doublefact(2*m+1);
  } else {
    /* Case for n even */
    return n * gsl_pow_uint(M_PI, n/2) / gsl_sf_fact(n/2);
  }
}

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
/* Function to build a kernel_density object                         */
/*********************************************************************/
kernel_density* build_kd(double *x, unsigned int ndim, 
			 unsigned int npt, double *wgt,
			 unsigned int leafsize, double *bandwidth, 
			 kernel_type ktype) {
  unsigned int i, j, curnode;
  double ds_n, hprod, wgttot=0.0;
  kernel_density *kd;

  /* Allocate memory */
  if (!(kd = malloc(sizeof(kernel_density)))) {
    fprintf(stderr, "bayesphot: error: unable to allocate memory in build_kernel_density\n");
    exit(1);
  }
  if (!(kd->h = calloc(ndim, sizeof(double)))) {
    fprintf(stderr, "bayesphot: error: unable to allocate memory in build_kernel_density\n");
    exit(1);
  }

  /* Record data describing kernel */
  for (i=0; i<ndim; i++) kd->h[i] = bandwidth[i];
  kd->ktype = ktype;

  /* Surface element factor for an n-sphere */
  ds_n = ds(ndim);

  /* Compute the normalization factor for the kernel around each
     point; this is defined as norm = 1/ \int K(z,h) dV */
  hprod = 1.0;
  for (i=0; i<ndim; i++) hprod *= bandwidth[i];
  switch (ktype) {
  case epanechnikov: {
    /* K(z, h) = 1 - z^2/h^2, z < h */
    kd->norm = ndim*(ndim+2) / (2.0*ds_n*hprod);
    break;
  }
  case tophat: {
    /* K(z, h) = 1, z < h */
    kd->norm = ndim / (hprod * ds_n);
    break;
  }
  case gaussian: {
    /* K(z, h) = exp[ -z^2/(2h^2) ], all z */
    kd->norm = 1.0 / (hprod * pow(sqrt(2.0*M_PI), ndim));
    break;
  }
  }

  /* Compute the normalization factor for the entire PDF */
  if (wgt == NULL) {
    kd->norm_tot = kd->norm / npt;
  } else {
    for (i=0; i<npt; i++) wgttot += wgt[i];
    kd->norm_tot = kd->norm / wgttot;
  }

  /* Build the KD tree around the data */
  if (wgt == NULL) {
    kd->tree = build_tree(x, ndim, npt, leafsize, NULL, 0);
  } else {
    kd->tree = build_tree(x, ndim, npt, leafsize, wgt, sizeof(double));
  }

  /* Allocate memory to hold summed weights of the nodes of the tree */
  if (!(kd->nodewgt = calloc(kd->tree->nodes, sizeof(double)))) {
    fprintf(stderr, "bayesphot: error: unable to allocate memory in build_kernel_density\n");
    exit(1);
  }
  kd->nodewgt--; /* Change to 1 offset instead of 0 offset */

  /* Breadth-first traversal of tree */

  /* Go to leftmost leaf */
  curnode = ROOT;
  while (kd->tree->tree[curnode].splitdim != -1) 
    curnode = LEFT(curnode);

  /* Traverse the leaves of the tree, adding up the weight in each */
  for (i=curnode; i<2*curnode; i++) {
    if (wgt != NULL) {
      kd->nodewgt[i] = 0.0;  
      for (j=0; j<kd->tree->tree[i].npt; j++)
	kd->nodewgt[i] += 
	  ((double *) kd->tree->tree[i].dptr)[j];
    } else {
      kd->nodewgt[i] = ((double) kd->tree->tree[i].npt);
    }
  }

  /* Now compute the weights in the rest of the tree by summing up the
     children */
  curnode = PARENT(curnode);
  while (curnode != 0) {
    for (i=curnode; i<2*curnode; i++) {
      if (kd->tree->tree[i].splitdim != -1) {
	kd->nodewgt[i] = kd->nodewgt[LEFT(i)] + kd->nodewgt[RIGHT(i)];
      } else {
	if (wgt != NULL)
	  for (j=0; j<kd->tree->tree[i].npt; j++)
	    kd->nodewgt[i] += 
	      ((double *) kd->tree->tree[i].dptr)[j];
	else
	  kd->nodewgt[i] = ((double) kd->tree->tree[i].npt);
      }
    }
    curnode = PARENT(curnode);
  }

  /* Return */
  return(kd);
}

/*********************************************************************/
/* De-allocates a kernel_density object                              */
/*********************************************************************/
void free_kd(kernel_density *kd) {
  free(kd->h);
  kd->nodewgt++;
  free(kd->nodewgt);
  free_tree(kd->tree);
  free(kd);
}


/*********************************************************************/
/* Function to change the weights in a kernel_density object         */
/*********************************************************************/
void kd_change_wgt(const double *wgt, kernel_density *kd) {
  unsigned int i, j, curnode;
  double wgttot;

  /* Change the weights */
  kd->tree->tree[1].dptr = (void *) wgt;

  /* Compute the normalization factor for the entire PDF */
  if (wgt == NULL) {
    kd->norm_tot = kd->norm / kd->tree->tree[1].npt;
  } else {
    wgttot = 0;
    for (i=0; i<kd->tree->tree[1].npt; i++) wgttot += wgt[i];
    kd->norm_tot = kd->norm / wgttot;
  }

  /* Reassign the individual node weight pointers and weight sums in a
     breadth-first traversal of the tree */

  /* Go to leftmost leaf */
  curnode = ROOT;
  while (kd->tree->tree[curnode].splitdim != -1) 
    curnode = LEFT(curnode);

  /* Traverse the leaves of the tree; at each leaf, set the offset in
     the weight array to match that in the data array, and adding up
     the weight of each point to get a node weight */
  for (i=curnode; i<2*curnode; i++) {
    if (wgt != NULL) {

      /* Set pointer to weight for this node */
      if (kd->tree->tree[i].npt != 0)
	kd->tree->tree[i].dptr = (void *) 
	  (wgt + (kd->tree->tree[i].x - kd->tree->tree[ROOT].x) / 
	   kd->tree->ndim);

      /* Get weight */
      kd->nodewgt[i] = 0.0;  
      for (j=0; j<kd->tree->tree[i].npt; j++)
	kd->nodewgt[i] += 
	  ((double *) kd->tree->tree[i].dptr)[j];

    } else {

      /* Set pointer to weight for this node */
      if (kd->tree->tree[i].npt != 0)
	kd->tree->tree[i].dptr = NULL;      

      /* Get weight */
      kd->nodewgt[i] = ((double) kd->tree->tree[i].npt);
    }
  }

  /* Now repeat for the rest of the tree summing up the contribution
     from the children to get the summed weight */
  curnode = PARENT(curnode);
  while (curnode != 0) {
    for (i=curnode; i<2*curnode; i++) {

      /* Set pointer to weight for this node */
      if (wgt != NULL) {
	if (kd->tree->tree[i].npt != 0)
	  kd->tree->tree[i].dptr = (void *) 
	    (wgt + (kd->tree->tree[i].x - kd->tree->tree[ROOT].x) / 
	     kd->tree->ndim);
      } else {
	kd->tree->tree[i].dptr = NULL;
      }

      /* Get weight for this node */
      if (kd->tree->tree[i].splitdim != -1) {
	kd->nodewgt[i] = kd->nodewgt[LEFT(i)] + kd->nodewgt[RIGHT(i)];
      } else {
	if (wgt != NULL)
	  for (j=0; j<kd->tree->tree[i].npt; j++)
	    kd->nodewgt[i] += 
	      ((double *) kd->tree->tree[i].dptr)[j];
	else
	  kd->nodewgt[i] = ((double) kd->tree->tree[i].npt);
      }
    }
    curnode = PARENT(curnode);
  }

}

/*********************************************************************/
/* Function to change the bandwidth in a kernel_density object       */
/*********************************************************************/
void kd_change_bandwidth(const double *bandwidth, kernel_density *kd) {
  unsigned int i;
  double ds_n, hprod, old_norm;

  /* Change bandwidth */
  for (i=0; i<kd->tree->ndim; i++) kd->h[i] = bandwidth[i];

  /* Recompute the normalization of each kernel */
  old_norm = kd->norm;
  ds_n = ds(kd->tree->ndim);
  hprod = 1.0;
  for (i=0; i<kd->tree->ndim; i++) hprod *= bandwidth[i];
  switch (kd->ktype) {
  case epanechnikov: {
    /* K(z, h) = 1 - z^2/h^2, z < h */
    kd->norm = kd->tree->ndim*(kd->tree->ndim+2) / (2.0*ds_n*hprod);
    break;
  }
  case tophat: {
    /* K(z, h) = 1, z < h */
    kd->norm = kd->tree->ndim / (hprod * ds_n);
    break;
  }
  case gaussian: {
    /* K(z, h) = exp[ -z^2/(2h^2) ], all z */
    kd->norm = 1.0 / (hprod * pow(sqrt(2.0*M_PI), kd->tree->ndim));
    break;
  }
  }
  kd->norm_tot *= kd->norm / old_norm;
}


/*********************************************************************/
/* Finds nearest neighbors to a specified input point                */
/*********************************************************************/
void kd_neighbors(const kernel_density *kd, const double *xpt, 
		  const unsigned int *dims, const unsigned int ndim, 
		  const unsigned int nneighbor,
		  const bool bandwidth_units, double *pos,
		  void *dptr, double *d2) {
  /* Call KDtree neighbor finding routine */
  if (bandwidth_units) {
    neighbors(kd->tree, xpt, dims, ndim, nneighbor, kd->h, pos, dptr, d2);
  } else {
    neighbors(kd->tree, xpt, dims, ndim, nneighbor, NULL, pos, dptr, d2);
  }
}

/*********************************************************************/
/* Find the N nearest neighbors to all points in the KDtree          */
/*********************************************************************/
void kd_neighbors_all(const kernel_density *kd, 
		      const unsigned int nneighbor, 
		      const bool bandwidth_units, unsigned int *idx, 
		      double *d2) {
  /* Call the KDtree neighbor finding routine */
  if (bandwidth_units) {
    neighbors_all(kd->tree, nneighbor, kd->h, idx, d2);
  } else {
    neighbors_all(kd->tree, nneighbor, NULL, idx, d2);
  }
}

/*********************************************************************/
/* Find N nearest neighbors to a single point in the KD tree         */
/*********************************************************************/
void kd_neighbors_point(const kernel_density *kd, 
			const unsigned int idxpt, 
			const unsigned int nneighbor,
			const bool bandwidth_units,
			unsigned int *idx, double *d2) {
  /* Just call the KDtree routine */
  if (bandwidth_units) {
    neighbors_point(kd->tree, idxpt, nneighbor, kd->h, idx, d2);
  } else {
    neighbors_point(kd->tree, idxpt, nneighbor, NULL, idx, d2);
  }
}

/*********************************************************************/
/* Vectorized version of kd_neighbors_point                          */
/*********************************************************************/
void kd_neighbors_point_vec(const kernel_density *kd, 
			    const unsigned int *idxpt, 
			    const unsigned int npt,
			    const unsigned int nneighbor,
			    const bool bandwidth_units,
			    unsigned int *idx, double *d2) {
  unsigned int i;
  for (i=0; i<npt; i++) {
    if (bandwidth_units) {
      neighbors_point(kd->tree, idxpt[i], nneighbor, kd->h, 
		      idx+i*nneighbor, d2+i*nneighbor);
    } else {
      neighbors_point(kd->tree, idxpt[i], nneighbor, NULL, 
		      idx+i*nneighbor, d2+i*nneighbor);
    }
  }
}

/*********************************************************************/
/* Same as kd_neighbors, but for a vector of input points            */
/*********************************************************************/
void kd_neighbors_vec(const kernel_density *kd, const double *xpt, 
		      const unsigned int *dims, const unsigned int ndim, 
		      const unsigned int npt, const unsigned int nneighbor,
		      const bool bandwidth_units, double *pos,
		      void *dptr, double *d2) {
  unsigned int i;

  /* Loop over input points, calling KDtree neighbor find routine on each */
  for (i=0; i<npt; i++) {
    if (bandwidth_units) {
      neighbors(kd->tree, 
		xpt + i*ndim, 
		dims, ndim, nneighbor, kd->h, 
		pos + i*ndim*nneighbor, 
		(void *) (((double *) dptr) + i*nneighbor),
		d2 + i*nneighbor);
    } else {
      neighbors(kd->tree, 
		xpt + i*ndim, 
		dims, ndim, nneighbor, NULL, 
		pos + i*ndim*nneighbor, 
		(void *) (((double *) dptr) + i*nneighbor),
		d2 + i*nneighbor);
    }
  }
}


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
  for (i=0; i<npt; i++)
    pdf[i] = kd_pdf_int(kd, x+i*ndim, dims, ndim, reltol, abstol
#ifdef DIAGNOSTIC
			, nodecheck+i, leafcheck+i, termcheck+i
#endif
			);
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
  for (i=0; i<npt; i++)
    pdf[i] = kd_pdf(kd, x+i*kd->tree->ndim, reltol, abstol
#ifdef DIAGNOSTIC
		    , nodecheck+i, leafcheck+i, termcheck+i
#endif
		    );
}

/*********************************************************************/
/* Function to report if we were compiled in diagnotic mode          */
/*********************************************************************/
bool diagnostic_mode() {
#ifdef DIAGNOSTIC
  return true;
#else
  return false;
#endif
}
