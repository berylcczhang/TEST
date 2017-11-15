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

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include "kernel_density_draw.h"
#include "gsl/gsl_randist.h"
#if 0
#  include "gsl/gsl_roots.h"
#endif
#include "gsl/gsl_sort.h"

#define NODEBLOCKSIZE 16384

/*********************************************************************/
/* Static functions                                                  */
/*********************************************************************/

static inline
double get_node_wgt(const kernel_density *kd, const double *x, 
		    const unsigned long *dims, 
		    const unsigned long ndim, 
		    const unsigned long curnode,
		    double *pointwgt);

static inline
int assign_idx(const kernel_density *kd, 
	       const double *x, 
	       const unsigned long *dims, 
	       const unsigned long ndim, 
	       const unsigned long *leaflist,
	       const double *leafwgt, 
	       const double leafwgt_sum,
	       const double leafwgt_sum_lim,
	       const unsigned long nsample, 
	       const double *rnd, 
	       const size_t* rndsort, 
	       const double *pointwgt,
	       size_t *idxmap);

/*********************************************************************/
/* Helper functions for use with gsl                                 */
/*********************************************************************/

#if 0
/* These functions computes 
   f = (n+2)/2 r^n - (n/2) r^(n+2) - x
   dfdr = n(n+2)/2 [r^(n-1) - r^(n+1)]
   where n = number of dimensions, x = input random number. This is
   the function whose root must be found for the transformation method
   approach to getting the distance using an Epanechnikov
   kernel. Note that we make n a double instead of an integer type
   for convenience. */
struct ep_root_params { double x; double ndim; };

double ep_root(double r, void *p) {
  struct ep_root_params *params = (struct ep_root_params *) p;
  double x = params->x;
  double n = params->ndim;
  return 0.5*(n+2.0)*pow(r, n) - 0.5*n*pow(r, n+2.0) - x;
}

double ep_root_deriv(double r, void *p) {
  struct ep_root_params *params = (struct ep_root_params *) p;
  double x = params->x;
  double n = params->ndim;
  return 0.5*n*(n+2.0) * (pow(r, n-1.0) - pow(r, n+1.0));
}

void ep_root_fdf(double r, void *p, double *f, double *df) {
  *f = ep_root(r, p);
  *df = ep_root_deriv(r, p);
}
#endif

/*********************************************************************/
/* Functions to initialze and free the random number generator       */
/*********************************************************************/

gsl_rng* rng_init(unsigned long seed) {
  FILE *fp;
  gsl_rng *r = gsl_rng_alloc(gsl_rng_mt19937);
  if (seed > 0) {
    gsl_rng_set(r, seed);
  } else {
    if ((fp = fopen("/dev/urandom", "r"))) {
      fread(&seed, 1, sizeof(seed), fp);
      fclose(fp);
    }
    gsl_rng_set(r, seed);
  }
  return r;
}  

void rng_free(gsl_rng *r) {
  gsl_rng_free(r);
}

/*********************************************************************/
/* Function to assign indices to random draws using breadth-first    */
/* traversal of the tree                                             */
/*********************************************************************/
void kd_pdf_draw_breadth(const kernel_density *kd,
			 const double *x,
			 const unsigned long *dims,
			 const unsigned long ndim,
			 const unsigned long nsample,
			 const double *rnd,
			 const size_t *rndsort,
			 size_t *idxmap) {

  double lastwgt;
  double *nodewgt, *pointwgt;
  size_t offset;
  unsigned long i, j, curnode;

  /* Allocate memory for temporary arrays */
  if (!(nodewgt = (double *) calloc(kd->tree->nodes, sizeof(double)))) {
    fprintf(stderr, "bayesphot: error: unable to allocate memory in kd_pdf_draw_breadth\n");
    exit(1);
  }
  if (!(pointwgt = (double *) calloc(kd->tree->tree[ROOT].npt,
				     sizeof(double)))) {
    fprintf(stderr, "bayesphot: error: unable to allocate memory in kd_pdf_draw_breadth\n");
    exit(1);
  }

  /* Go to leftmost leaf */
  curnode = ROOT;
  while (kd->tree->tree[curnode].splitdim != -1) 
    curnode = LEFT(curnode);

  /* Traverse the tree from left to right, starting from leftmost
     leaf */
  while (curnode != 0) {

    for (i=curnode; i<2*curnode; i++) {

      /* Is this node a leaf */
      if (kd->tree->tree[i].splitdim != -1) {
	
	/* Non-leaf nodes just have a cumulative weight equal to that
	   of their right-most child */
	nodewgt[i] = nodewgt[RIGHT(i)];

      } else {

	/* For leaf nodes, we compute the weights of all points in the
	   leaf, and of the leaf as a whole, if there are any points */

	/* Get weight at end of previous leaf */
	if (i == curnode) lastwgt = 0.0;
	else lastwgt = nodewgt[i-1];

	/* Does this leaf have any points */
	if (kd->tree->tree[i].npt == 0) {

	  /* No points in this leaf, so just set its cumulative weight
	     equal to the last cumulative weight and move on */
	  nodewgt[i] = lastwgt;

	} else {

	  /* This leaf does have points, so get offset to them */
	  offset = (kd->tree->tree[i].x - kd->tree->tree[ROOT].x) /
	    kd->tree->ndim;

	  /* Get differential weights of points and of entire node */
	  nodewgt[i] = get_node_wgt(kd, x, dims, ndim, i,
				    pointwgt + offset);

	  
	  /* Change recorded weights in this leaf from differential to
	     cumulative */
	  nodewgt[i] += lastwgt;
	  if (kd->tree->tree[i].npt > 0) pointwgt[offset] += lastwgt;
	  for (j=1; j<kd->tree->tree[i].npt; j++) {
	    pointwgt[offset+j] += pointwgt[offset+j-1];
	  }
	}
      }
    }

    /* Go up one level in tree */
    curnode = PARENT(curnode);
  }

  /* Now traverse tree to assign indices to reach random number, using
     depth-first traversal */
  curnode = ROOT;
  for (i=0; i<nsample; i++) {

    /* Loop until we find the point we're looking for */
    while (1) {

      /* Does this node contain the point we are interested in? */
      if (rnd[rndsort[i]] < nodewgt[curnode]/nodewgt[ROOT]) {

	/* This node does contain point. See if it is a leaf */
	if (kd->tree->tree[curnode].splitdim != -1) {

	  /* It is not the leaf, so descend one level in the tree */
	  curnode = LEFT(curnode);
	  
	} else {

	  /* This is a leaf, so see which point in the leaf is the one
	     we want */
	  offset = (kd->tree->tree[curnode].x - kd->tree->tree[ROOT].x)
	    / kd->tree->ndim;
	  j = 0;
	  while (rnd[rndsort[i]] >= pointwgt[offset+j]/nodewgt[ROOT]) j++;

	  /* Record the point, then go to next random number */
	  idxmap[i] = offset+j;
	  break;
	}
	
      } else {

	/* This node does not contain the point we are interested in,
	   so go to next node */
	SETNEXT(curnode);

      }
    }
  }

  /* Free temporary storage */
  free(nodewgt);
  free(pointwgt);
}
  
/*********************************************************************/
/* Function to draw from the PDF using depth-first traversal         */
/*********************************************************************/
void kd_pdf_draw_depth(const kernel_density *kd,
		       const double *x,
		       const unsigned long *dims,
		       const unsigned long ndim,
		       const unsigned long nsample,
		       const double *rnd,
		       const size_t *rndsort,
		       size_t *idxmap) {

  unsigned long i, nnodealloc, nleafalloc, nnode, nleaf;
  unsigned long nleaf_old, nnode_old, ptr, lptr, rptr;
  unsigned long *leaflist, *nodelist;
  double rnd_max, maxwgt, leafwgt_sum, nodewgt_sum, lwgt, rwgt;
  double *leafwgt, *nodewgt, *pointwgt;

  /* Allocate temporary storage */
  if (!(nodelist = (unsigned long *) 
	calloc(NODEBLOCKSIZE, sizeof(unsigned long)))) {
    fprintf(stderr, "bayesphot: error: unable to allocate memory in kd_pdf_draw\n");
    exit(1);
  }
  if (!(nodewgt = (double *) 
	calloc(NODEBLOCKSIZE, sizeof(double)))) {
    fprintf(stderr, "bayesphot: error: unable to allocate memory in kd_pdf_draw\n");
    exit(1);
  }
  if (!(leaflist = (unsigned long *) 
	calloc(NODEBLOCKSIZE, sizeof(unsigned long)))) {
    fprintf(stderr, "bayesphot: error: unable to allocate memory in kd_pdf_draw\n");
    exit(1);
  }
  if (!(leafwgt = (double *) 
	calloc(NODEBLOCKSIZE, sizeof(double)))) {
    fprintf(stderr, "bayesphot: error: unable to allocate memory in kd_pdf_draw\n");
    exit(1);
  }
  if (!(pointwgt = (double *)
	calloc(kd->tree->leafsize*(NODEBLOCKSIZE+2), sizeof(double)))) {
    fprintf(stderr, "bayesphot: error: unable to allocate memory in kd_pdf_draw\n");
    exit(1);
  }
  nnodealloc = nleafalloc = NODEBLOCKSIZE;

  /* Shorthand to largest of random numbers */
  rnd_max = rnd[rndsort[nsample-1]];

  /* Analyze the root node, and push onto the appropriate list */
  if (kd->tree->tree[ROOT].splitdim == -1) {
    leaflist[0] = ROOT;
    leafwgt[0] = get_node_wgt(kd, x, dims, ndim, ROOT, pointwgt);
    nodewgt_sum = 0.0;
    leafwgt_sum = leafwgt[0];
    nnode = 0;
    nleaf = 1;
  } else {
    nodelist[0] = ROOT;
    nodewgt[0] = get_node_wgt(kd, x, dims, ndim, ROOT, pointwgt);
    nodewgt_sum = nodewgt[0];
    leafwgt_sum = 0.0;
    nnode = 1;
    nleaf = 0;
  }

  /* Now traverse the tree */
  while (1) {

    /* Are we done? There are two conditions for this: first, the
       weight checked must exceed the largest random number we
       have. Second, the assigments of random numbers to specific
       items in the tree must not depend on the unknown weight of the
       parts of the tree we have not yet checked. We also put in a
       safety to force the loop to exit if there are no nodes left
       to check; this is needed in case nodewgt_sum doesn't go to
       zero precisely enough to satisfy the weight condition as a
       result of finite precision arithmetic. */
    if (nnode == 0 || 
	leafwgt_sum / (leafwgt_sum + nodewgt_sum) > rnd_max) {
      if (assign_idx(kd, x, dims, ndim, leaflist, leafwgt, leafwgt_sum, 
		     leafwgt_sum + nodewgt_sum, nsample, rnd, 
		     rndsort, pointwgt, idxmap))
	break;
    }

    /* Find the node with the largest possible weight that is not a
       leaf */
    maxwgt = nodewgt[0];
    ptr = 0;
    for (i=1; i<nnode; i++) {
      if (nodewgt[i] > maxwgt) {
	maxwgt = nodewgt[i];
	ptr = i;
      }
    }

    /* Subtract this node's weight from the total */
    nodewgt_sum -= maxwgt;
    nnode--;

    /* Analyze this node's children */
    nleaf_old = nleaf;
    nnode_old = nnode;
    lptr = LEFT(nodelist[ptr]);
    rptr = RIGHT(nodelist[ptr]);
    lwgt = get_node_wgt(kd, x, dims, ndim, lptr, 
			pointwgt+nleaf*kd->tree->leafsize);
    if (kd->tree->tree[lptr].splitdim == -1) nleaf++; else nnode++;
    rwgt = get_node_wgt(kd, x, dims, ndim, rptr,
			pointwgt+nleaf*kd->tree->leafsize);
    if (kd->tree->tree[rptr].splitdim == -1) nleaf++; else nnode++;

    /* Allocate more memory to hold child nodes if necessary */
    if (nnode > nnodealloc) {
      nnodealloc *= 2;
      if (!(nodelist = (unsigned long *) 
	    realloc(nodelist, nnodealloc*sizeof(unsigned long)))) {
	fprintf(stderr, "bayesphot: error: unable to allocate memory in kd_pdf_draw\n");
	exit(1);
      }
      if (!(nodewgt = (double *) 
	    realloc(nodewgt, nnodealloc*sizeof(double)))) {
	fprintf(stderr, "bayesphot: error: unable to allocate memory in kd_pdf_draw\n");
	exit(1);
      }
    }
    if (nleaf > nleafalloc) {
      nleafalloc *= 2;
      if (!(leaflist = (unsigned long *) 
	    realloc(leaflist, nleafalloc*sizeof(unsigned long)))) {
	fprintf(stderr, "bayesphot: error: unable to allocate memory in kd_pdf_draw\n");
	exit(1);
      }
      if (!(leafwgt = (double *) 
	    realloc(leafwgt, nleafalloc*sizeof(double)))) {
	fprintf(stderr, "bayesphot: error: unable to allocate memory in kd_pdf_draw\n");
	exit(1);
      }
      if (!(pointwgt = (double *) 
	    realloc(pointwgt, 
		    (nleafalloc+2)*kd->tree->leafsize*sizeof(double)))) {
	fprintf(stderr, "bayesphot: error: unable to allocate memory in kd_pdf_draw\n");
	exit(1);
      }
    }

    /* Put child nodes on the appropriate list, and update running
       tallies */
    if (kd->tree->tree[lptr].splitdim == -1) {
      leaflist[nleaf_old] = lptr;
      leafwgt[nleaf_old] = lwgt;
      leafwgt_sum += lwgt;
      nleaf_old++;
    } else {
      nodelist[ptr] = lptr;
      nodewgt[ptr] = lwgt;
      nodewgt_sum += lwgt;
      nnode_old++;
    }
    if (kd->tree->tree[rptr].splitdim == -1) {
      leaflist[nleaf_old] = rptr;
      leafwgt[nleaf_old] = rwgt;
      leafwgt_sum += rwgt;
    } else {
      if (kd->tree->tree[lptr].splitdim == -1) {
	nodelist[ptr] = rptr;
	nodewgt[ptr] = rwgt;
      } else {
	nodelist[nnode_old] = rptr;
	nodewgt[nnode_old] = rwgt;
      }
      nodewgt_sum += rwgt;
    }

    /* If we haven't already overwritten this node in the node list
       (which is the case if both the children are leaves) do so now */
    if ((kd->tree->tree[lptr].splitdim == -1) &&
	(kd->tree->tree[rptr].splitdim == -1)) {
      for (i=ptr; i<nnode_old; i++) {
	nodelist[i] = nodelist[i+1];
	nodewgt[i] = nodewgt[i+1];
      }
    }
  }

  /* Free temporaries */
  free(nodelist);
  free(nodewgt);
  free(leaflist);
  free(leafwgt);
  free(pointwgt);
}


/*********************************************************************/
/* Function to draw from the PDF                                     */
/*********************************************************************/
void kd_pdf_draw(const kernel_density *kd, const double *x,
		 const unsigned long *dims, const unsigned long ndim,
		 const unsigned long nsample, const draw_method method,
		 const gsl_rng *r, double *out) {

  size_t *rndsort, *idxmap;
  double *rnd, *rndgauss = NULL;
  double norm, dist;
  unsigned long *dims_out;
  unsigned long ndim_out, i, j, k;
#if 0
  const gsl_root_fsolver_type *T;
  gsl_root_fsolver *s = NULL;
  struct ep_root_params ep_params;
  gsl_function_fdf FDF;
#endif

  /* Epanechnikov not done yet */
  assert(kd->ktype != epanechnikov);

  /* Allocate memory for temporary arrays */
  if (!(rndsort = (size_t *) calloc(nsample, sizeof(size_t)))) {
    fprintf(stderr, "bayesphot: error: unable to allocate memory in kd_pdf_draw\n");
    exit(1);
  }
  if (!(rnd = (double *) calloc(nsample, sizeof(double)))) {
    fprintf(stderr, "bayesphot: error: unable to allocate memory in kd_pdf_draw\n");
    exit(1);
  }
  ndim_out = kd->tree->ndim - ndim;
  if (!(dims_out = (unsigned long *) 
	calloc(ndim_out, sizeof(unsigned long)))) {
    fprintf(stderr, "bayesphot: error: unable to allocate memory in kd_pdf_draw\n");
    exit(1);
  }
  if (!(idxmap = (size_t *) calloc(nsample, sizeof(size_t)))) {
    fprintf(stderr, "bayesphot: error: unable to allocate memory in kd_pdf_draw\n");
    exit(1);
  }

  /* Figure out which dimensions we'll be outputting */
  for (i=0, j=0; i<kd->tree->ndim; i++) {
    for (k=0; k<ndim; k++) if (dims[k] == i) break;
    if (k==ndim) dims_out[j++] = i;
  }

  /* Generate a random number in [0,1) for each sample; make an
     indirectly sorted version of the random number list for later
     convenience */
  for (i=0; i<nsample; i++) rnd[i] = gsl_rng_uniform(r);
  gsl_sort_index(rndsort, rnd, 1, nsample);

  /* Traverse tree to map the random numbers we've just picked into
     sample indices */
  if (method == draw_depth_first) {
    kd_pdf_draw_depth(kd, x, dims, ndim, nsample, rnd, rndsort, idxmap);
  } else {
    kd_pdf_draw_breadth(kd, x, dims, ndim, nsample, rnd, rndsort, idxmap);
  }

  /* Now draw actual values of samples */
  for (i=0; i<nsample; i++) {
    switch (kd->ktype) {
    case gaussian: {
      /* Gaussian kernel: just pick from a Gaussian distribution in
	 each dimension, since the probabilities are nicely separable */
      for (j=0; j<ndim_out; j++)
	out[i*ndim_out+j] = 
	  kd->tree->tree[ROOT].x[idxmap[i]*kd->tree->ndim+dims_out[j]] +
	  gsl_ran_gaussian_ziggurat(r, kd->h[dims_out[j]]);
      break;
    }
    case tophat: {
      /* Tophat kernel: in this case we can compute the distance using a
	 transformation method, then get the angles using the method
	 of Marsiglia (1972) -- see
	 http://mathworld.wolfram.com/HyperspherePointPicking.html. */
      dist = pow(gsl_rng_uniform(r), 1.0/ndim_out);
      if (rndgauss == NULL) {
	if (!(rndgauss = (double *) calloc(ndim_out, sizeof(double)))) {
	  fprintf(stderr, "bayesphot: error: unable to allocate memory in kd_pdf_draw\n");
	  exit(1);
	}
      }
      norm = 0.0;
      for (j=0; j<ndim_out; j++) {
	rndgauss[j] = gsl_ran_gaussian_ziggurat(r, 1.0);
	norm += rndgauss[j] * rndgauss[j];
      }
      norm = sqrt(norm);
      for (j=0; j<ndim_out; j++)
	out[i*ndim_out+j] = 
	  kd->tree->tree[ROOT].x[idxmap[i]*kd->tree->ndim+dims_out[j]] +
	  dist * kd->h[dims_out[j]] * rndgauss[j] / norm;
      break;
    }
    case epanechnikov: {
      /* Epanechnikov kernel: in this case we get a distance using a
	 transformation method as for the tophat kernel, but in n
	 dimensions the inversion that is the final step of the
	 transformation method requires finding the root of a
	 higher-order polynomial equation numerically. This is then
	 followed up with the same Marsiglia method to get the
	 angles. */
      /* Note done yet */
      exit(1);
      break;
    }
    }
  }
  if ((kd->ktype == tophat) || (kd->ktype == epanechnikov)) free(rndgauss);

  /* Free temporaries */
  free(idxmap);
  free(rndsort);
  free(rnd);
  free(dims_out);
}


/*********************************************************************/
/* Static functions                                                  */
/*********************************************************************/

/* This function computes the maximum total weight of a node; if the
   node is a leaf it computes the actual weight of each point in the
   leaf, and of the entire node */
inline
double get_node_wgt(const kernel_density *kd, const double *x, 
		    const unsigned long *dims,
		    const unsigned long ndim, 
		    const unsigned long curnode,
		    double *pointwgt) {

  unsigned long i;
  unsigned long ndim_tree = kd->tree->ndim;
  double wgt, d2;

  /* Is this node a leaf? */
  if (kd->tree->tree[curnode].splitdim == -1) {

    /* Loop over points, summing their contribution */
    wgt = 0.0;
    for (i=0; i<kd->tree->tree[curnode].npt; i++) {

      /* Get distance in units of the kernel size */
      d2 = dist2(&(kd->tree->tree[curnode].x[ndim_tree*i]), x,
		 ndim_tree, ndim, NULL, dims, kd->h, ndim_tree);

      /* Compute the weight of this point */
      if (kd->tree->tree[curnode].dptr != NULL) {
	/* Case where points have individual weights */
	switch (kd->ktype) {
	case epanechnikov: {
	  if (d2 < 1)
	    pointwgt[i] = ((double *) kd->tree->tree[curnode].dptr)[i] * 
	      (1.0 - d2) * kd->norm_tot;
	  break;
	}
	case tophat: {
	  if (d2 < 1)
	    pointwgt[i] = ((double *) kd->tree->tree[curnode].dptr)[i] 
	      * kd->norm_tot;
	  break;
	}
	case gaussian: {
	  pointwgt[i] = ((double *) kd->tree->tree[curnode].dptr)[i] *
	    exp(-d2/2.0) * kd->norm_tot;
	  break;
	}
	}
      } else {
	/* Case where all points have equal weights */
	switch (kd->ktype) {
	case epanechnikov: {
	  if (d2 < 1) pointwgt[i] = (1.0 - d2) * kd->norm_tot;
	  break;
	}
	case tophat: {
	  if (d2 < 1) pointwgt[i] = kd->norm_tot;
	  break;
	}
	case gaussian: {
	  pointwgt[i] = exp(-d2/2.0) * kd->norm_tot;
	  break;
	}
	}
      }

      /* Add contribution of this point to total for this leaf */
      wgt += pointwgt[i];

    }

    /* Normalize and return */
    return(wgt);

  } else {

    /* This node is not a leaf, so return an upper limit on its weight */
    /* Get minimum distance in units of the kernel size */
    d2 = box_min_dist2(x, 
		       (const double **) kd->tree->tree[curnode].xbnd,
		       ndim, ndim_tree, dims, NULL, kd->h, ndim_tree);

    /* Return value at minimum distance */
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

    /* Dummy return statement to suppress compiler warning; never executed */
    return 0.0;
  }
}


/* This function maps from a list of random numbers to a list of
   indices in the kernel density tree, using a set of nodes and their
   associated weights */
inline
int assign_idx(const kernel_density *kd, 
	       const double *x, 
	       const unsigned long *dims, 
	       const unsigned long ndim, 
	       const unsigned long *leaflist,
	       const double *leafwgt, 
	       const double leafwgt_sum,
	       const double leafwgt_sum_lim,
	       const unsigned long nsample, 
	       const double *rnd, 
	       const size_t *rndsort,
	       const double *pointwgt,
	       size_t *idxmap) {
  unsigned long i, j, jsave, k, ksave = 0;
  double wgt, cumwgt = 0, cumwgt1 = 0;

  /* Loop over the random numbers, in sorted order */
  for (i=0, jsave=-1, j=0; i<nsample; i++) {

    /* Move the pointer in the leaf list until the total weight is
       bigger than the random number we're working on */
    while ((cumwgt+leafwgt[j])/leafwgt_sum < rnd[rndsort[i]]) 
      cumwgt += leafwgt[j++];

    /* Now loop through the points in this node */
    if (j != jsave) {
      ksave = 0;
      cumwgt1 = cumwgt;
      jsave = j;
    }
    for (k=ksave; k<kd->tree->tree[leaflist[j]].npt; k++) {

      /* Record point weight */
      wgt = pointwgt[j*kd->tree->leafsize+k];
      cumwgt1 += wgt;

      /* Is this the point we want? */
      if (cumwgt1/leafwgt_sum > rnd[rndsort[i]]) {

	/* Yes; make sure we get the same point using the limit, and,
	   if not, bail */
	if (cumwgt1/leafwgt_sum_lim <= rnd[rndsort[i]]) {
	  return 0;
	}

	/* Record the index */
	idxmap[rndsort[i]] = k + (unsigned long)
	  (kd->tree->tree[leaflist[j]].x -
	   kd->tree->tree[ROOT].x)/kd->tree->ndim;
	ksave = k+1;

	/* Check if there are any more sample points that should also
	   be assigned this data point; if so, do the same operations
	   on them */
	for (; i<nsample-1; i++) {
	  if (cumwgt1/leafwgt_sum <= rnd[rndsort[i+1]]) break;
	  if (cumwgt1/leafwgt_sum_lim <= rnd[rndsort[i+1]])
	    return 0;
	  idxmap[rndsort[i+1]] = k + (unsigned long)
	    (kd->tree->tree[leaflist[j]].x -
	     kd->tree->tree[ROOT].x)/kd->tree->ndim;
	}

	/* We found our point, so stop looping through the leaf */
	break;
      }
    }
  }

  /* If we're here, we successfully assigned all points */
  return 1;
}

#undef NODEBLOCKSIZE

