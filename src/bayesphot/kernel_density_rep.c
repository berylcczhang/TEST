/*********************************************************************
Copyright (C) 2015 Mark Krumholz
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
#include <gsl/gsl_sort.h>
#include "geometry.h"
#include "kernel_density_rep.h"

/*********************************************************************/
/* Functions to build the representation                             */
/*********************************************************************/

#define NODEBLOCKSIZE 512
#define PTSIZE 8192

unsigned long kd_rep(const kernel_density *kd, const double *x,
		     const unsigned long *dims, 
		     const unsigned long ndim,
		     const double reltol,
		     double **xpt, double **wgts) {

  unsigned long i, j, k, ptr, curnode, nlist, nlist_alloc, npt_alloc, npt;
  unsigned long dimptr, npt_final, ndim_ret = kd->tree->ndim - ndim;
  unsigned long *nodelist;
  size_t *idx;
  double wgt_in, wgt_out, maxwgt, d2;
  double *nodewgt;
  double *xtmp, *wgttmp;

  /* Make sure our kernel is Gaussian */
  assert(kd->ktype == gaussian);

  /* Initialize the outputs */
  npt = 0;
  if (!(xtmp = (double *) 
	calloc(ndim_ret*PTSIZE, sizeof(double)))) {
    fprintf(stderr, "bayesphot: error: unable to allocate memory in kd_rep_nodes\n");
    exit(1);
  }
  if (!(wgttmp = (double *) calloc(PTSIZE, sizeof(double)))) {
    fprintf(stderr, "bayesphot: error: unable to allocate memory in kd_rep_nodes\n");
    exit(1);
  }
  npt_alloc = PTSIZE;

  /* Initialize the internal lists */
  if (!(nodelist = (unsigned long *) 
	calloc(NODEBLOCKSIZE, sizeof(unsigned long)))) {
    fprintf(stderr, "bayesphot: error: unable to allocate memory in kd_pdf\n");
    exit(1);
  }
  if (!(nodewgt = (double *) 
	calloc(NODEBLOCKSIZE, sizeof(double)))) {
    fprintf(stderr, "bayesphot: error: unable to allocate memory in kd_pdf\n");
    exit(1);
  }
  nlist_alloc = NODEBLOCKSIZE;  

  /* Put the root node into the internal list */
  nodelist[0] = ROOT;
  nodewgt[0] = kd->nodewgt[ROOT];
  nlist = 1;

  /* Initialize the sums that keep track of the weight of nodes that
     we're going to return, and the maximum weight of those we've left
     out */
  wgt_in = 0.0;
  wgt_out = nodewgt[0];

  /* Main loop */
  while (1) {

    /* Find the node that is contributing the most error */
    maxwgt = nodewgt[0];
    ptr = 0;
    for (i=1; i<nlist; i++) {
      if (nodewgt[i] > maxwgt) {
	ptr = i;
	maxwgt = nodewgt[i];
      }
    }
    curnode = nodelist[ptr];

    /* Subtract this node's contribution from the excluded weight;
       avoid roundoff error */
    wgt_out -= maxwgt;
    if (wgt_out < 0) wgt_out = 0.0;

    /* If this node is a leaf, add its points to the return list, and
       add its tally to the weight we've accounted for */
    if (kd->tree->tree[curnode].splitdim == -1) {

      /* Add to return list, allocating more memory if needed */
      if (npt + kd->tree->tree[curnode].npt >= npt_alloc) {
	npt_alloc *= 2;
	if (!(xtmp = (double *) 
	      realloc(xtmp, ndim_ret*npt_alloc*sizeof(double)))) {
	  fprintf(stderr, "bayesphot: error: unable to allocate memory in kd_rep_nodes\n");
	  exit(1);
	}
	if (!(wgttmp = (double *)
	      realloc(wgttmp, npt_alloc*sizeof(double)))) {
	  fprintf(stderr, "bayesphot: error: unable to allocate memory in kd_rep_nodes\n");
	  exit(1);
	}
      }

      /* Add the contributions of points in this leaf, and store its
	 points in the output arrays */
      for (i=0; i<kd->tree->tree[curnode].npt; i++) {
      
	/* Get distance */
	d2 = dist2(&(kd->tree->tree[curnode].x[kd->tree->ndim*i]), x,
		   kd->tree->ndim, ndim, NULL, NULL, kd->h, 
		   kd->tree->ndim);
      
	/* Compute weight of point */
	wgttmp[npt] = exp(-d2/2.0);
	if (kd->tree->tree[ptr].dptr != NULL) 
	  wgttmp[npt] *= ((double *) kd->tree->tree[curnode].dptr)[i];

	/* Store point in output array */
	dimptr = 0;
	for (j=0; j<kd->tree->ndim; j++) {
	  for (k=0; k<ndim; k++) if (j == dims[k]) break;
	  if (k == ndim) {
	    xtmp[npt*ndim_ret+dimptr] = 
	      kd->tree->tree[curnode].x[kd->tree->ndim*i+j];
	    dimptr++;
	  }
	}
 
        /* Add to sum and increment counter */
	wgt_in += wgttmp[npt];
	npt++;
      }

      /* Remove this node from list of nodes to check */
      for (i=ptr; i<nlist-1; i++) {
	nodelist[i] = nodelist[i+1];
	nodewgt[i] = nodewgt[i+1];
      }
      nlist--;

    } else {

      /* This node is not a leaf */

      /* Allocate more memory in the node list if needed */
      if (nlist+1 == nlist_alloc) {
	nlist_alloc *= 2;
	if (!(nodelist = (unsigned long *) 
	      realloc(nodelist, nlist_alloc*sizeof(unsigned long)))) {
	  fprintf(stderr, "bayesphot: error: unable to allocate memory in kd_pdf\n");
	  exit(1);
	}
	if (!(nodewgt = (double *) 
	      realloc(nodewgt, nlist_alloc*sizeof(double)))) {
	  fprintf(stderr, "bayesphot: error: unable to allocate memory in kd_pdf\n");
	  exit(1);
	}
      }

      /* Replace this node in the list with its children */
      nodelist[ptr] = LEFT(curnode);
      nodelist[nlist] = RIGHT(curnode);

      /* Get weight for the child nodes */
      d2 = box_min_dist2(x, (const double **) 
			 kd->tree->tree[nodelist[ptr]].xbnd,
			 ndim, kd->tree->ndim, dims, NULL, kd->h, ndim);
      nodewgt[ptr] = kd->nodewgt[nodelist[ptr]] * exp(-d2/2.0);
      d2 = box_min_dist2(x, (const double **) 
			 kd->tree->tree[nodelist[nlist]].xbnd,
			 ndim, kd->tree->ndim, dims, NULL, kd->h, ndim);
      nodewgt[nlist] = kd->nodewgt[nodelist[nlist]] * exp(-d2/2.0);

      /* Add weight of child nodes to sum, and increment list pointer */
      wgt_out += nodewgt[ptr] + nodewgt[nlist];
      nlist++;
    }

    /* Check if we're converged */
    if (wgt_in / (wgt_in + wgt_out) > 1.0-reltol) break;
  }

  /* Free some memory */
  free(nodelist);
  free(nodewgt);

  /* Allocate memory to hold the final result; we will reduce the
     size of this later as needed */
  if (!(*xpt = (double *) calloc(ndim_ret*npt, sizeof(double)))) {
    fprintf(stderr, "bayesphot: error: unable to allocate memory in kd_pdf\n");
    exit(1);
  }
  if (!(*wgts = (double *) calloc(npt, sizeof(double)))) {
    fprintf(stderr, "bayesphot: error: unable to allocate memory in kd_pdf\n");
    exit(1);
  }

  /* Indirectly sort points by weight */
  if (!(idx = (size_t *) calloc(npt, sizeof(size_t)))) {
    fprintf(stderr, "bayesphot: error: unable to allocate memory in kd_pdf\n");
    exit(1);
  }
  gsl_sort_index(idx, wgttmp, 1, npt);

  /* Traverse the list of points, from highest weight to lowest; only
     keep the points we need to achieve the required level of accuracy */
  wgt_out += wgt_in;
  wgt_in = 0.0;
  for (i = npt-1, npt_final = 0; 
       wgt_in / (wgt_in + wgt_out) > 1.0-reltol;
       i--, npt_final++) {
    /* Copy point */
    (*wgts)[npt_final] = wgttmp[idx[i]];
    for (j=0; j<ndim_ret; j++) 
      (*xpt)[npt_final*ndim_ret+j] = xtmp[idx[i]*ndim_ret+j];
    /* Recompute weights */
    wgt_in += (*wgts)[npt_final];
    wgt_out -= (*wgts)[npt_final];
    /* Safety assertion: shouldn't be tripped, because if we add in
       all the points we should get back to an acceptable value of
       wgt_in */
    assert(i>0);
  };

  /* Free memory */
  free(xtmp);
  free(wgttmp);
  *xpt = (double *) realloc(*xpt, ndim_ret*npt_final*sizeof(double));
  *wgts = (double *) realloc(*wgts, npt_final*sizeof(double));

  /* Normalize weights to have sum of unity */
  for (i=0; i<npt_final; i++) (*wgts)[i] /= wgt_in;

  /* Return */
  return npt_final;
}

#undef NODEBLOCKSIZE
#undef PTSIZE

/*********************************************************************/
/* Memory freeing convenience function                               */
/*********************************************************************/

void free_kd_rep(double **xpt, double **wgts) {
  free(*xpt);
  free(*wgts);
  *xpt = NULL;
  *wgts = NULL;
}

