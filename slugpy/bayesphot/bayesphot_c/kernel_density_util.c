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
#include "geometry.h"
#include "kernel_density_util.h"

/*********************************************************************/
/* Function to build a kernel_density object                         */
/*********************************************************************/
kernel_density* build_kd(double *x, unsigned long ndim, 
			 unsigned long npt, double *wgt,
			 unsigned long leafsize, double *bandwidth, 
			 kernel_type ktype, unsigned long minsplit,
			 unsigned long *sortmap) {
  unsigned long curnode;
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
  for (unsigned long i=0; i<ndim; i++) kd->h[i] = bandwidth[i];
  kd->ktype = ktype;

  /* Surface element factor for an n-sphere */
  ds_n = ds(ndim);

  /* Compute the normalization factor for the kernel around each
     point; this is defined as norm = 1/ \int K(z,h) dV */
  hprod = 1.0;
  for (unsigned long i=0; i<ndim; i++) hprod *= bandwidth[i];
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
#pragma omp parallel for reduction(+:wgttot)
    for (unsigned long i=0; i<npt; i++) wgttot += wgt[i];
    kd->norm_tot = kd->norm / wgttot;
  }

  /* Build the KD tree around the data */
  if (wgt == NULL) {
    kd->tree = build_tree(x, ndim, npt, leafsize, NULL, 0, minsplit,
			  sortmap);
  } else {
    kd->tree = build_tree(x, ndim, npt, leafsize, wgt, 
			  sizeof(double), minsplit, sortmap);
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
#pragma omp parallel for
  for (unsigned long i=curnode; i<2*curnode; i++) {
    if (wgt != NULL) {
      kd->nodewgt[i] = 0.0;  
      for (unsigned long j=0; j<kd->tree->tree[i].npt; j++)
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
#pragma omp parallel for
    for (unsigned long i=curnode; i<2*curnode; i++) {
      if (kd->tree->tree[i].splitdim != -1) {
	kd->nodewgt[i] = kd->nodewgt[LEFT(i)] + kd->nodewgt[RIGHT(i)];
      } else {
	if (wgt != NULL)
	  for (unsigned long j=0; j<kd->tree->tree[i].npt; j++)
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
/* Function to build a kernel_density object, with a tree sorted     */
/* only on certain dimensions                                        */
/*********************************************************************/
kernel_density* build_kd_sortdims(double *x, unsigned long ndim, 
				  unsigned long npt, double *wgt,
				  unsigned long leafsize, 
				  double *bandwidth, 
				  kernel_type ktype, int *nosort,
				  unsigned long *sortmap) {
  unsigned long i, j, curnode;
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
    kd->tree = build_tree_sortdims(x, ndim, npt, leafsize, NULL, 0,
				   nosort, sortmap);
  } else {
    kd->tree = build_tree_sortdims(x, ndim, npt, leafsize, wgt, 
				   sizeof(double), nosort, sortmap);
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
/* Copies a kernel_density object                                    */
/*********************************************************************/
kernel_density *copy_kd(kernel_density *kd) {

  unsigned long i;
  kernel_density *kdcopy;

  /* Allocate memory */
  if (!(kdcopy = malloc(sizeof(kernel_density)))) {
    fprintf(stderr, "bayesphot: error: unable to allocate memory in build_kernel_density\n");
    exit(1);
  }
  if (!(kdcopy->h = calloc(kd->tree->ndim, sizeof(double)))) {
    fprintf(stderr, "bayesphot: error: unable to allocate memory in build_kernel_density\n");
    exit(1);
  }

  /* Copy data */
  kdcopy->tree = kd->tree;
  kdcopy->norm = kd->norm;
  kdcopy->norm_tot = kd->norm_tot;
  kdcopy->ktype = kd->ktype;
  kdcopy->nodewgt = kd->nodewgt;
  for (i=0; i<kd->tree->ndim; i++) kdcopy->h[i] = kd->h[i];

  /* Return */
  return kdcopy;
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
/* De-allocates a kernel_density object without deleting the KDtree  */
/*********************************************************************/
void free_kd_copy(kernel_density *kd) {
  free(kd->h);
  free(kd);
}


/*********************************************************************/
/* Function to change the weights in a kernel_density object         */
/*********************************************************************/
void kd_change_wgt(const double *wgt, kernel_density *kd) {
  unsigned long curnode;
  double wgttot;

  /* Change the weights */
  kd->tree->tree[1].dptr = (void *) wgt;

  /* Compute the normalization factor for the entire PDF */
  if (wgt == NULL) {
    kd->norm_tot = kd->norm / kd->tree->tree[1].npt;
  } else {
    wgttot = 0;
#pragma omp parallel for shared(wgt) reduction (+:wgttot)
    for (unsigned long i=0; i<kd->tree->tree[1].npt; i++)
      wgttot += wgt[i];
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
#pragma omp parallel for
  for (unsigned long i=curnode; i<2*curnode; i++) {
    if (wgt != NULL) {

      /* Set pointer to weight for this node */
      if (kd->tree->tree[i].npt != 0)
	kd->tree->tree[i].dptr = (void *) 
	  (wgt + (kd->tree->tree[i].x - kd->tree->tree[ROOT].x) / 
	   kd->tree->ndim);

      /* Get weight */
      kd->nodewgt[i] = 0.0;  
      for (unsigned long j=0; j<kd->tree->tree[i].npt; j++)
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
#pragma omp parallel for
    for (unsigned long i=curnode; i<2*curnode; i++) {

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
	  for (unsigned long j=0; j<kd->tree->tree[i].npt; j++)
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
  unsigned long i;
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
/* Function to report if we were compiled in diagnotic mode          */
/*********************************************************************/
bool diagnostic_mode() {
#ifdef DIAGNOSTIC
  return true;
#else
  return false;
#endif
}
