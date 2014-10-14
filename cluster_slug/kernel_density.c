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

static inline
void kd_pdf_tol_node(const double *x, const unsigned int curnode, 
		     const kernel_density *kd, double *pdf,
		     double *pdferr);

/*********************************************************************/
/* Function to build a kernel_density object                         */
/*********************************************************************/
kernel_density* build_kd(double *x, unsigned int ndim, 
			 unsigned int npt, double *wgt, 
			 unsigned int leafsize, double bandwidth, 
			 kernel_type ktype) {

  unsigned int i, j, k, curnode;
  double ds_n, wgttot=0.0;
  kernel_density *kd;

  /* Allocate memory */
  if (!(kd = malloc(sizeof(kernel_density)))) {
    fprintf(stderr, "cluster_slug: error: unable to allocate memory in build_kernel_density\n");
    exit(1);
  }

  /* Record miscellaneous data */
  kd->h = bandwidth;
  kd->scalefac = NULL;
  kd->ktype = ktype;
  kd->copy_pos = false;

  /* Surface element factor for an n-sphere */
  ds_n = ds(ndim);

  /* Compute the normalization factor for the kernel around each
     point; this is defined as norm = 1/ \int K(z,h) dV */
  switch (ktype) {
  case epanechnikov: {
    /* K(z, h) = 1 - z^2/h^2, z < h */
    kd->norm = ndim*(ndim+2) / (2.0*ds_n*pow(kd->h, ndim));
    break;
  }
    /* K(z, h) = 1, z < h */
  case tophat: {
    kd->norm = ndim / (pow(kd->h, ndim) * ds_n);
    break;
  }
    /* K(z, h) = exp[ -z^2/(2h^2) ], all z */
  case gaussian: {
    kd->norm = 1.0 / pow(sqrt(2.0*M_PI)*kd->h, ndim);
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

  /* Allocate memory to hold summed weights and weighted positions of
     the nodes of the tree */
  if (!(kd->nodecm = calloc(kd->tree->nodes*ndim, sizeof(double)))) {
    fprintf(stderr, "cluster_slug: error: unable to allocate memory in build_kernel_density\n");
    exit(1);
  }
  if (!(kd->nodewgt = calloc(kd->tree->nodes, sizeof(double)))) {
    fprintf(stderr, "cluster_slug: error: unable to allocate memory in build_kernel_density\n");
    exit(1);
  }
  kd->nodewgt--; /* Change to 1 offset instead of 0 offset */
  kd->nodecm -= ndim;  /* Change to 1 offset instead of 0 offset */

  /* Breadth-first traversal of tree */

  /* Go to leftmost leaf */
  curnode = ROOT;
  while (kd->tree->tree[curnode].splitdim != -1) 
    curnode = LEFT(curnode);

  /* Traverse the leaves of the tree, adding up the weight in each */
  for (i=curnode; i<2*curnode; i++) {
    for (k=0; k<ndim; k++) kd->nodecm[i*ndim+k] = 0.0;
    if (wgt != NULL) {
      kd->nodewgt[i] = 0.0;  
      for (j=0; j<kd->tree->tree[i].npt; j++) {
	kd->nodewgt[i] += 
	  ((double *) kd->tree->tree[i].dptr)[j];
	for (k=0; k<ndim; k++)
	  kd->nodecm[i*ndim+k] += kd->tree->tree[i].x[j*ndim+k] *
	    ((double *) kd->tree->tree[i].dptr)[j];
      }
      kd->nodewgt[i] *= kd->norm / wgttot;
      for (k=0; k<ndim; k++) kd->nodecm[i*ndim+k] /= wgttot;
    } else {
      kd->nodewgt[i] = 
	((double) kd->tree->tree[i].npt) / 
	((double) npt) * kd->norm;
      for (j=0; j<kd->tree->tree[i].npt; j++) {
	for (k=0; k<ndim; k++) 
	  kd->nodecm[i*ndim+k] += kd->tree->tree[i].x[j*ndim+k];
      }
      for (k=0; k<ndim; k++) 
	kd->nodecm[i*ndim+k] /= kd->tree->tree[i].npt;
    }
  }

  /* Now compute the centers of mass and weights in the rest of the
     tree by summing up the children */
  curnode = PARENT(curnode);
  while (curnode != 0) {
    for (i=curnode; i<2*curnode; i++) {
      for (j=0; j<ndim; j++)
	kd->nodecm[i*ndim+j] = 
	  (kd->nodewgt[LEFT(i)] * kd->nodecm[LEFT(i)*ndim+j] +
	   kd->nodewgt[RIGHT(i)] * kd->nodecm[RIGHT(i)*ndim+j]) /
	  (kd->nodewgt[LEFT(i)] + kd->nodewgt[RIGHT(i)] );
      kd->nodewgt[i] = kd->nodewgt[LEFT(i)] + kd->nodewgt[RIGHT(i)];
    }
    curnode = PARENT(curnode);
  }

  /* Return */
  return(kd);
}


/*********************************************************************/
/* Function to build a kernel_density object where the bandwidth is  */
/* not the same in all dimensions                                    */
/*********************************************************************/
kernel_density* build_kd_bvec(double *x, unsigned int ndim, 
			      unsigned int npt, double *wgt, 
			      unsigned int leafsize, 
			      double *bandwidth, 
			      kernel_type ktype) {

  unsigned int i, j;
  unsigned long ndata;
  double *xcopy, *scalefac;
  kernel_density *kd;

  /* Make a copy of the positions */
  ndata = ndim*npt;
  if (!(xcopy = calloc(ndata, sizeof(double)))) {
    fprintf(stderr, "cluster_slug: error: unable to allocate memory in build_kernel_density_bvec\n");
    exit(1);
  }
  memcpy(xcopy, x, ndata*sizeof(double));

  /* Compute scale factors that will be used to render the bandwidth
     equal size in all dimensions */
  if (!(scalefac = calloc(ndim, sizeof(double)))) {
    fprintf(stderr, "cluster_slug: error: unable to allocate memory in build_kernel_density_bvec\n");
    exit(1);
  }
  scalefac[0] = 1.0;
  for (i=1; i<ndim; i++) scalefac[i] = bandwidth[0] / bandwidth[i];

  /* Rescale positions */
  for (i=0; i<npt; i++)
    for (j=1; j<ndim; j++) xcopy[i*ndim+j] *= scalefac[j];

  /* Now call the scalar build routine on the rescaled data */
  kd = build_kd(xcopy, ndim, npt, wgt, leafsize, bandwidth[0], ktype);

  /* Store the scale factors, and that we have copied the positions */
  kd->scalefac = scalefac;
  kd->copy_pos = true;

  /* Return */
  return(kd);
}


/*********************************************************************/
/* Finds nearest neighbors to a specified input point                */
/*********************************************************************/
void find_neighbors(const double *xpt, const unsigned int *dims,
		    const unsigned int ndim, 
		    const unsigned int nneighbor, 
		    kernel_density *kd, double *pos,
		    double *wgt, double *dist2) {

  unsigned int i, j;
  double dist2_tmp;
  double *xpt_copy;

  /* Rescale coordinates if necessary */
  if (kd->scalefac != NULL) {
    if (!(xpt_copy = calloc(ndim, sizeof(double)))) {
      fprintf(stderr, "cluster_slug: error: unable to allocate memory in find_neighbors\n");
      exit(1);
    }
    for (i=0; i<ndim; i++) xpt_copy[i] = xpt[i] * kd->scalefac[dims[i]];
  } else {
    xpt_copy = (double *) xpt;
  }

  /* Call KDtree neighbor finding routine */
  neighbors(kd->tree, xpt_copy, dims, ndim, nneighbor, pos, wgt, dist2);

  /* Rescale coordinates if necessary */
  if (kd->scalefac != NULL) {
    for (i=0; i<nneighbor; i++) {
      for (j=0; j<kd->tree->ndim; j++)
	pos[i*kd->tree->ndim+j] /= kd->scalefac[j];
      for (dist2_tmp=0.0, j=0; j<ndim; j++) 
	dist2_tmp += SQR(xpt[j] - pos[i*kd->tree->ndim+dims[j]]);
      dist2[i] = dist2_tmp;
    }
    free(xpt_copy);
  }
}


/*********************************************************************/
/* De-allocates a kernel_density object                              */
/*********************************************************************/
void free_kd(kernel_density **kd) {
  if ((*kd)->copy_pos) free((*kd)->tree->tree[1].x);
  if ((*kd)->scalefac != NULL) free((*kd)->scalefac);
  (*kd)->nodewgt++;
  (*kd)->nodecm += (*kd)->tree->ndim;
  free((*kd)->nodewgt);
  free((*kd)->nodecm);
  free_tree(&((*kd)->tree));
  free(*kd);
  *kd = NULL;
}

/*********************************************************************/
/* Function to evaluate the PDF using a kernel_density object        */
/*********************************************************************/
double kd_pdf(const double *x, const kernel_density* kd) {

  unsigned int i, npt;
  double *xsample, *wsample, *dist2, *x_copy;
  double h2, sum;

  /* If we're using a non-compact kernel, call the approximate PDF
     routine instead, using a default tolerance */
  if (kd->ktype == gaussian)
    return kd_pdf_tol(x, kd, 1e-6, 0.0);

  /* Rescale coordinates if necessary */
  if (kd->scalefac != NULL) {
    if (!(x_copy = calloc(kd->tree->ndim, sizeof(double)))) {
      fprintf(stderr, "cluster_slug: error: unable to allocate memory in kd_pdf\n");
      exit(1);
    }
    for (i=0; i<kd->tree->ndim; i++) x_copy[i] = x[i] * kd->scalefac[i];
  } else {
    x_copy = (double *) x;
  }

  /* Find all the points within a distance h of the requested point */
  npt = query_sphere(kd->tree, x_copy, kd->h, &xsample, 
		     (void *) &wsample, &dist2);

  /* Evaluate the kernel at the returned points */
  sum = 0.0;
  h2 = SQR(kd->h);
  if (kd->tree->tree[1].dptr != NULL) {
    /* Case with weights */
    if (kd->ktype == epanechnikov) {
      for (i=0; i<npt; i++) sum += wsample[i] * (1.0 - dist2[i]/h2);
    } else if (kd->ktype == tophat) {
      for (i=0; i<npt; i++) sum += wsample[i];
    }
  } else {
    /* Case without weights */
    if (kd->ktype == epanechnikov) {
      for (i=0; i<npt; i++) sum += 1.0 - dist2[i]/h2;
    } else if (kd->ktype == tophat) {
      for (i=0; i<npt; i++) sum += 1.0;
    }
  }
  sum = sum * kd->norm_tot;

  /* Renormalize if using scaled coordinates */
  if (kd->scalefac != NULL)
    for (i=1; i<kd->tree->ndim; i++) sum *= kd->scalefac[i];

  /* Free */
  if (npt > 0) {
    free(xsample);
    if (kd->tree->tree[1].dptr != NULL) free(wsample);
    free(dist2);
  }
  if (kd->scalefac != NULL) free(x_copy);

  /* Return */
  return(sum);
}

/*********************************************************************/
/* Function to evaluate the PDF approximately using a kernel_density */
/* object                                                            */
/*********************************************************************/
#define NODEBLOCKSIZE 1024
double kd_pdf_tol(const double *x, const kernel_density *kd,
		  const double reltol, const double abstol) {
  unsigned int i, nnode, nalloc, ptr;
  unsigned int *nodelist;
  double *x_copy, *nodepdf, *nodeerr;
  double maxerr, relerr, abserr;
  double pdf, pdferr;

  /* Rescale coordinates if necessary */
  if (kd->scalefac != NULL) {
    if (!(x_copy = (double *) calloc(kd->tree->ndim, sizeof(double)))) {
      fprintf(stderr, "cluster_slug: error: unable to allocate memory in kd_pdf_tol\n");
      exit(1);
    }
    for (i=0; i<kd->tree->ndim; i++) x_copy[i] = x[i] * kd->scalefac[i];
  } else {
    x_copy = (double *) x;
  }

  /* Allocate memory to start */
  if (!(nodelist = (unsigned int *) 
	calloc(NODEBLOCKSIZE, sizeof(unsigned int)))) {
    fprintf(stderr, "cluster_slug: error: unable to allocate memory in kd_pdf_tol\n");
    exit(1);
  }
  if (!(nodepdf = (double *) 
	calloc(NODEBLOCKSIZE, sizeof(double)))) {
    fprintf(stderr, "cluster_slug: error: unable to allocate memory in kd_pdf_tol\n");
    exit(1);
  }
  if (!(nodeerr = (double *) 
	calloc(NODEBLOCKSIZE, sizeof(double)))) {
    fprintf(stderr, "cluster_slug: error: unable to allocate memory in kd_pdf_tol\n");
    exit(1);
  }
  nalloc = NODEBLOCKSIZE;

  /* Analyze root node */
  kd_pdf_tol_node(x_copy, ROOT, kd, nodepdf, nodeerr);
  nodelist[0] = ROOT;

  /* Initialize the global estimates */
  pdf = nodepdf[0];
  pdferr = nodeerr[0];
  nnode = 1;
  abserr = pdferr;
  relerr = pdferr / (pdf+DBL_MIN);

  /* Now work recursively through the tree, identifying the node in
     the tree that contributes most to the error and opening it until
     the error estimate is within our tolerance */
  while (abserr > abstol && relerr > reltol) {

    /* Find the node that is contributing the most error */
    maxerr = 0.0;
    for (i=0; i<nnode; i++) {
      if (nodeerr[i] > maxerr) {
	ptr = i;
	maxerr = nodeerr[i];
      }
    }

    /* Allocate more memory if necessary */
    if (nnode == nalloc) {
      if (!(nodelist = (unsigned int *) 
	    realloc(nodelist, 2*nalloc*sizeof(unsigned int)))) {
	fprintf(stderr, "cluster_slug: error: unable to allocate memory in kd_pdf_tol\n");
	exit(1);
      }
      if (!(nodepdf = (double *) 
	    realloc(nodepdf, 2*nalloc*sizeof(double)))) {
	fprintf(stderr, "cluster_slug: error: unable to allocate memory in kd_pdf_tol\n");
	exit(1);
      }
      if (!(nodeerr = (double *) 
	    realloc(nodeerr, 2*nalloc*sizeof(double)))) {
	fprintf(stderr, "cluster_slug: error: unable to allocate memory in kd_pdf_tol\n");
	exit(1);
      }
      nalloc *= 2;
    }

    /* Subtract this node's contribution to the global estimates */
    pdf -= nodepdf[ptr];
    pdferr -= nodeerr[ptr];

    /* Compute estimates for that node's children, and add them to
       the list */
    kd_pdf_tol_node(x_copy, LEFT(nodelist[ptr]), kd, 
		    nodepdf+ptr, nodeerr+ptr);
    kd_pdf_tol_node(x_copy, RIGHT(nodelist[ptr]), kd,
		    nodepdf+nnode, nodeerr+nnode);
    nodelist[nnode] = RIGHT(nodelist[ptr]);
    nodelist[ptr] = LEFT(nodelist[ptr]);

    /* Recompute totals and error estimates */
    pdf += nodepdf[ptr] + nodepdf[nnode];
    pdferr += nodeerr[ptr] + nodeerr[nnode];
    abserr = pdferr;
    relerr = pdferr / pdf;

    /* Increment the node counter */
    nnode++;
  }

  /* Renormalize if using scaled coordinates */
  if (kd->scalefac != NULL)
    for (i=1; i<kd->tree->ndim; i++) pdf *= kd->scalefac[i];

  /* Free memory */
  free(nodelist);
  free(nodepdf);
  free(nodeerr);
  if (kd->scalefac != NULL) free(x_copy);

  /* Return */
  return(pdf);
}

/*********************************************************************/
/* Function to estimate the contribution to a PDF from a single node */
/*********************************************************************/
void kd_pdf_tol_node(const double *x, const unsigned int curnode, 
		     const kernel_density *kd, double *pdf,
		     double *pdferr) {
  unsigned int i;
  unsigned int ndim = kd->tree->ndim;
  double d2min, d2max, d2cm, pdfmin, pdfmax, pdfcm;
  double d2, h2 = SQR(kd->h);

  /* Is this node a leaf? If so, just sum over it and return that */
  if (kd->tree->tree[curnode].splitdim == -1) {

    *pdf = 0.0;
    *pdferr = 0.0;

    /* Loop over points */
    for (i=0; i<kd->tree->tree[curnode].npt; i++) {

      /* Get distance */
      d2 = euclidean_metric2(x, 
			     &(kd->tree->tree[curnode].x[ndim*i]),
			     ndim);

      /* Compute contribution based on kernel type and presence or
	 absence of weights */
      if (kd->tree->tree[curnode].dptr != NULL) {
	switch (kd->ktype) {
	case epanechnikov: {
	  if (d2 < h2)
	    *pdf += ((double *) kd->tree->tree[curnode].dptr)[i] * 
	      (1.0 - d2/h2);
	  break;
	}
	case tophat: {
	  if (d2 < h2)
	    *pdf += ((double *) kd->tree->tree[curnode].dptr)[i];
	  break;
	}
	case gaussian: {
	  *pdf += ((double *) kd->tree->tree[curnode].dptr)[i] *
	    exp(-d2/h2);
	  break;
	}
	}
      } else {
	switch (kd->ktype) {
	case epanechnikov: {
	  if (d2<h2)
	    *pdf += 1.0 - d2/h2;
	  break;
	}
	case tophat: {
	  if (d2<h2)
	    *pdf += 1.0;
	  break;
	}
	case gaussian: {
	  *pdf += exp(-d2/h2);
	  break;
	}
	}
      }

    }

    /* Normalize */
    *pdf *= kd->norm_tot;
    *pdferr *= kd->norm_tot;

    /* Return */
    return;
  }

  /* If we're here, this node is not a leaf, so return an estimate
     using the min and max distances, as well as the cm value */

  /* Get min, max, and cm distance */
  d2min = box_min_dist2(x, 
			(const double **) kd->tree->tree[curnode].xlim,
			ndim, ndim, NULL);
  d2max = box_max_dist2(x, 
			(const double **) kd->tree->tree[curnode].xlim,
			ndim, ndim, NULL, kd->tree->xswap);
  d2cm = euclidean_metric2(x, &(kd->nodecm[curnode*ndim]), ndim);

  /* Get min, max, and cm evaluation of PDF */
  switch (kd->ktype) {
  case epanechnikov: {
    pdfmin = d2min < h2 ? kd->nodewgt[ROOT] * (1.0 - d2min/h2) : 0.0;
    pdfmax = d2max < h2 ? kd->nodewgt[ROOT] * (1.0 - d2max/h2) : 0.0;
    pdfcm = d2cm < h2 ? kd->nodewgt[ROOT] * (1.0 - d2cm/h2) : 0.0;
    break;
  }
  case tophat: {
    pdfmin = d2min < h2 ? kd->nodewgt[ROOT] : 0.0;
    pdfmax = d2max < h2 ? kd->nodewgt[ROOT] : 0.0;
    pdfcm = d2cm < h2 ? kd->nodewgt[ROOT] : 0.0;
    break;
  }
  case gaussian: {
    pdfmin = kd->nodewgt[ROOT] * exp(-d2min/h2);
    pdfmax = kd->nodewgt[ROOT] * exp(-d2max/h2);
    pdfcm = kd->nodewgt[ROOT] * exp(-d2cm/h2);
    break;
  }
  }

  /* Store estimates */
  *pdf = pdfcm * kd->norm_tot;
  *pdferr = fmax(pdfmin - pdfcm, pdfcm - pdfmax) * kd->norm_tot;
}


/*********************************************************************/
/* Function to evaluate the PDF using a kernel_density object for a  */
/* vector of inputs positions                                        */
/*********************************************************************/
void kd_pdf_vec(const double *x, const unsigned int npt, 
		const kernel_density *kd, double *pdf) {
  unsigned int i;
  for (i=0; i<npt; i++)
    pdf[i] = kd_pdf(x+i*kd->tree->ndim, kd);
}

/*********************************************************************/
/* Function to evaluate the PDF at a specified point in some of the  */
/* dimensions, integrating over other dimensions                     */
/*********************************************************************/
double kd_pdf_int(const double *q, const unsigned int *qdim, 
		  unsigned int nqdim, const kernel_density *kd) {
  unsigned int i, npt, nint;
  double *xsample, *wsample, *dist2, *q_copy;
  double ds_n, h2, r2, sum;

  /* Safety assertion */
  assert(nqdim > 0);
  assert(nqdim < kd->tree->ndim);

  /* Rescale coordinates if necessary */
  if (kd->scalefac != NULL) {
    if (!(q_copy = calloc(nqdim, sizeof(double)))) {
      fprintf(stderr, "cluster_slug: error: unable to allocate memory in kd_pdf_int\n");
      exit(1);
    }
    for (i=0; i<nqdim; i++) q_copy[i] = q[i] * kd->scalefac[qdim[i]];
  } else {
    q_copy = (double *) q;
  }

  /* Find all the sample points within a distance h of the requested
     projection plane */
  npt = query_slab(kd->tree, q_copy, qdim, nqdim, kd->h, &xsample, 
		     (void *) &wsample, &dist2);

  /* Get the volume element in spherical coordinates for the
     dimensions we're integrating out */
  nint = kd->tree->ndim - nqdim;  /* Num. dim. in integral */
  ds_n = ds(nint); /* N-sphere surf. area*/

  /* Sum the kernel over the points found */
  sum = 0.0;
  h2 = SQR(kd->h);
  if (kd->tree->tree[1].dptr != NULL) {
    /* Case with weights */
    if (kd->ktype == epanechnikov) {
      for (i=0; i<npt; i++) {
	r2 = h2 - dist2[i];
	sum += wsample[i] * ds_n * pow(r2, 0.5*nint+1.0)/h2 *
	  (1.0/nint - 1.0/(nint+2.0));
      }
    } else if (kd->ktype == tophat) {
      for (i=0; i<npt; i++) {
	r2 = h2 - dist2[i];
	sum += wsample[i] * ds_n/nint * pow(r2, 0.5*nint);
      }
    }
  } else {
    /* Case without weights */
    if (kd->ktype == epanechnikov) {
      for (i=0; i<npt; i++) {
	r2 = h2 - dist2[i];
	sum += ds_n * pow(r2, 0.5*nint+1.0)/h2 *
	  (1.0/nint - 1.0/(nint+2.0));
      }
    } else if (kd->ktype == tophat) {
      for (i=0; i<npt; i++) {
	r2 = h2 - dist2[i];
	sum += ds_n/nint * pow(r2, 0.5*nint);
      }
    }
  }
  sum = sum * kd->norm_tot;

  /* Renormalize if using scaled coordinates */
  if (kd->scalefac != NULL)
    for (i=1; i<nqdim; i++) sum *= kd->scalefac[qdim[i]];

  /* Free */
  if (npt > 0) {
    free(xsample);
    if (kd->tree->tree[1].dptr != NULL) free(wsample);
    free(dist2);
  }
  if (kd->scalefac != NULL) free(q_copy);

  /* Return */
  return(sum);
}


/*********************************************************************/
/* Function to change the weights in a kernel_density object         */
/*********************************************************************/
void reweight_kd(const double *wgt, kernel_density *kd) {

  int i;
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
}
