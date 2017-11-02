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
#include "kernel_density_util.h"
#include "kernel_density_rep.h"
#include "kernel_density_rep_data.h"



/*********************************************************************/
/* Declarations of local functions                                   */
/*********************************************************************/

static inline
double node_merge_error(double w, double logdx);

/*********************************************************************/
/* Functions to build the representation                             */
/*********************************************************************/

#define NODEBLOCKSIZE 512
#define PTSIZE 8192

unsigned long kd_rep(const kernel_density *kd, const double *x,
		     const unsigned long *dims, 
		     const unsigned long ndim, const double reltol,
		     const unsigned long *dim_return, 
		     const unsigned long ndim_return,
		     double **xpt, double **wgts) {

  unsigned long i, j, k, ptr, curnode, nlist, nlist_alloc, npt_alloc, npt;
  unsigned long dimptr, npt_final, ndim_ret;
  unsigned long *nodelist;
  size_t *idx;
  double wgt_in, wgt_out, maxwgt, d2;
  double *nodewgt;
  double *xtmp, *wgttmp;

  /* Make sure our kernel is Gaussian */
  assert(kd->ktype == gaussian);

  /* Set the number of dimensions we're returning to the default if
     dim_return is NULL */
  if (dim_return) ndim_ret = kd->tree->ndim - ndim;
  else ndim_ret = ndim_return;

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

  /* Main loop */
  while (1) {

    /* Find the node that is contributing the most error */
    maxwgt = nodewgt[0];
    ptr = 0;
    wgt_out = nodewgt[0];
    for (i=1; i<nlist; i++) {
      if (nodewgt[i] > maxwgt) {
	ptr = i;
	maxwgt = nodewgt[i];
      }
      wgt_out += nodewgt[i];
    }
    curnode = nodelist[ptr];
    wgt_out -= maxwgt;

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
		   kd->tree->ndim, ndim, NULL, dims, kd->h, 
		   kd->tree->ndim);
      
	/* Compute weight of point */
	wgttmp[npt] = exp(-d2/2.0);
	if (kd->tree->tree[curnode].dptr != NULL) 
	  wgttmp[npt] *= ((double *) kd->tree->tree[curnode].dptr)[i];

	/* Store point in output array */
	dimptr = 0;
	for (j=0; j<kd->tree->ndim; j++) {
	  /* If this dimension is in dims, skip it */
	  for (k=0; k<ndim; k++) if (j == dims[k]) break;
	  if (k != ndim) continue;
	  /* If this dimension is not in dim_return, and dim_return is
	     set, skip it */
	  if (dim_return) {
	    for (k=0; i<ndim_return; k++) if (j == dim_return[k]) break;
	    if (k == ndim) continue;
	  }
	  /* If we're here, copy the dimension */ 
	  xtmp[npt*ndim_ret+dimptr] = 
	    kd->tree->tree[curnode].x[kd->tree->ndim*i+j];
	  dimptr++;
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
      if (nlist+1 >= nlist_alloc) {
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

      /* Increment list pointer */
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
  for (npt_final = 0;  wgt_in / (wgt_in + wgt_out) < 1.0-reltol;
       npt_final++) {
    /* Copy point */
    i = npt - npt_final - 1;
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
/* Function to squeeze a kernel density representation by merging    */
/* Gaussians that lie close to one antoher                           */
/*********************************************************************/
unsigned long squeeze_rep(const unsigned long npt, 
			  const unsigned int ndim, double *h, 
			  const double tol, double **x, 
			  double **wgts) {
  kernel_density *kd;
  unsigned long i, j, levptr, npt_final;
  double d2, w, log_hmin, logdx, *node_err;

  /* Get minimum bandwidth */
  log_hmin = log10(h[0]);
  for (i=1; i<ndim; i++) 
    log_hmin = log_hmin < log10(h[i]) ? log_hmin : log10(h[i]);

  /* Build a kernel density representation of the data, with leaves
     consisting of two points */
  kd = build_kd(*x, ndim, npt, *wgts, 2, h, gaussian, 0, NULL);

  /* Allocate temporaries */
  if (!(node_err = (double *) calloc(kd->tree->nodes, sizeof(double)))) {
    fprintf(stderr, "bayesphot: error: unable to allocate memory in squeeze_rep\n");
    exit(1);
  }
  for (i=0; i<kd->tree->nodes; i++) node_err[i] = 0.0;

  /* Go to leftmost node of the tree */
  levptr = ROOT;
  while (kd->tree->tree[levptr].splitdim != -1) 
    levptr = LEFT(levptr);

  /* Breadth-first traversal of the tree, starting with the leaves
     and working up */
  for (; levptr != 0; levptr=PARENT(levptr)) {

    /* Traverse the list of nodes from least weight to greatest */
    for (i=levptr; i<2*levptr; i++) {

      /* Skip nodes with 0 points ; we don't need to do anything with
	 them */
      if (kd->tree->tree[i].npt <= 1) continue;

      /* Case where the node is a leaf */
      if (kd->tree->tree[i].splitdim == -1) {

	/* Get relative weights of the two points */
	w = ((double *) kd->tree->tree[i].dptr)[0] /
	  (((double *) kd->tree->tree[i].dptr)[0] +
	   ((double *) kd->tree->tree[i].dptr)[1]);

	/* Get log separation of node points in units of the minimum
	   bandwidth */ 
	d2 = dist2(&(kd->tree->tree[i].x[0]),
		   &(kd->tree->tree[i].x[ndim]),
		   ndim, ndim, NULL, NULL, NULL, ndim);
	logdx = 0.5 * log10(d2) - log_hmin;

	/* Compute the error associated with merging the points in this
	   node */
	node_err[i] = kd->nodewgt[i] * kd->norm * 
	  node_merge_error(w, logdx);

	/* If the error is too big, do nothing; go to next node */
	if (node_err[i] / (kd->nodewgt[i] * kd->norm) > tol)
	  continue;

	/* Compute the location of the merged point, and store in first
	   slot for this node */
	for (j=0; j<ndim; j++)
	  kd->tree->tree[i].x[j] = w * kd->tree->tree[i].x[j] +
	    (1.0-w) * kd->tree->tree[i].x[ndim+j];

	/* Set the weight of the merged point to the weight of the
	   entire node, and zero out the other point */
	((double *) kd->tree->tree[i].dptr)[0] = kd->nodewgt[i];
	((double *) kd->tree->tree[i].dptr)[1] = 0.0;

      } else {

	/* Case where node is not a leaf */

	/* Don't attempt to merge nodes if we've already rejected
	   merging their children */
	if (((double *) kd->tree->tree[LEFT(i)].dptr)[0] != 
	    kd->nodewgt[LEFT(i)])
	  continue;
	if (((double *) kd->tree->tree[RIGHT(i)].dptr)[0] != 
	    kd->nodewgt[RIGHT(i)])
	  continue;

	/* Get relative weights of the two points */
	w = kd->nodewgt[LEFT(i)] / 
	  (kd->nodewgt[LEFT(i)] + kd->nodewgt[RIGHT(i)]);

	/* Get log separation of node points in units of the minimum
	   bandwidth */ 
	d2 = dist2(kd->tree->tree[LEFT(i)].x,
		   kd->tree->tree[RIGHT(i)].x,
		   ndim, ndim, NULL, NULL, NULL, ndim);
	logdx = 0.5 * log10(d2) - log_hmin;

	/* Compute the error associated with merging the points in this
	   node */
	node_err[i] = kd->nodewgt[i] * kd->norm * 
	  node_merge_error(w, logdx);

	/* Add in errors from children */
	node_err[i] += node_err[LEFT(i)] + node_err[RIGHT(i)];

	/* If the error is too big, do nothing; go to next node */
	if (node_err[i] / (kd->nodewgt[i] * kd->norm) > tol)
	  continue;

	/* Compute the location of the merged point, and store in first
	   slot for this node */
	for (j=0; j<ndim; j++)
	  kd->tree->tree[i].x[j] = w * kd->tree->tree[LEFT(i)].x[j] +
	    (1.0-w) * kd->tree->tree[RIGHT(i)].x[j];

	/* Set the weight of the merged point to the weight of the
	   entire node, and zero out the other point */
	((double *) kd->tree->tree[i].dptr)[0] = kd->nodewgt[i];
	((double *) kd->tree->tree[RIGHT(i)].dptr)[0] = 0.0;
      }
    }
  }

  /* Last step: compress the final arrays, resize them in memory, then
     return the new size */
  for (i=0, npt_final=0; i<npt; i++) {
    if ((*wgts)[i] != 0.0) {
      /* Keep this point */
      if (npt_final < i) {
	for (j=0; j<ndim; j++) (*x)[npt_final*ndim+j] = (*x)[i*ndim+j];
	(*wgts)[npt_final] = (*wgts)[i];
      }
      npt_final++;
    }
  }
  (*x) = realloc(*x, npt_final*ndim*sizeof(double));
  (*wgts) = realloc(*wgts, npt_final*sizeof(double));
  return npt_final;
}


/*********************************************************************/
/* Memory freeing convenience function                               */
/*********************************************************************/

void free_kd_rep(double **xpt, double **wgts) {
  free(*xpt);
  free(*wgts);
  *xpt = NULL;
  *wgts = NULL;
}


/*********************************************************************/
/* Return the maximum error that results from merging the points in  */
/* a node; the quantity returns is normalized to the maximum of the  */
/* merged Gaussian                                                   */
/*********************************************************************/

inline
double node_merge_error(double w, double logdx) {
  double wwgt, dxwgt;
  int widx, dxidx;

  /* Get indices into the table */
  if (w>0.5) w=1.0-w;
  widx = (int) ((w-wmin) / dw);
  dxidx = (int) ((logdx-logdxmin) / dlogdx);

  /* Get weights in the table, being sure to catch cases off the table */
  if (widx < 0) {
    widx = 0;
    wwgt = 1.0;
  } else if (widx >= nw) {
    widx = nw-1;
    wwgt = -1.0;
  } else {
    wwgt = 1.0 + widx - (w-wmin)/dw; 
  }
  if (dxidx < 0) {
    dxidx = 0;
    dxwgt = 1.0;
  } else if (widx >= nw) {
    dxidx = ndx-1;
    dxwgt = -1.0;
  } else {
    dxwgt = 1.0 + dxidx - (logdx-logdxmin)/dlogdx; 
  }

  /* Do linear interpolation to get result */
  return pow(10, wwgt * dxwgt * gauss_err_tab[widx][dxidx] +
	     (1.0-wwgt) * dxwgt * gauss_err_tab[widx+1][dxidx] +
	     wwgt * (1.0-dxwgt) * gauss_err_tab[widx][dxidx+1] +
	     (1.0-wwgt) * (1.0-dxwgt) * gauss_err_tab[widx+1][dxidx+1]);
}
