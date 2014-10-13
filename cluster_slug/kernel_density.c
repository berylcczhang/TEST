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


/*********************************************************************/
/* Function to build a kernel_density object                         */
/*********************************************************************/
kernel_density* build_kd(double *x, unsigned int ndim, 
			 unsigned int npt, double *wgt, 
			 unsigned int leafsize, double bandwidth, 
			 kernel_type ktype) {

  unsigned int i;
  double ds_n, wgttot;
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
    wgttot = 0;
    for (i=0; i<npt; i++) wgttot += wgt[i];
    kd->norm_tot = kd->norm / wgttot;
  }

  /* Build the KD tree around the data */
  if (wgt == NULL) {
    kd->tree = build_tree(x, ndim, npt, leafsize, NULL, 0);
  } else {
    kd->tree = build_tree(x, ndim, npt, leafsize, wgt, sizeof(double));
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
