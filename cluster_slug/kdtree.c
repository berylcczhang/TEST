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
#include <math.h>
#include <string.h>
#include <float.h>
#include "kdtree.h"

/*********************************************************************/
/* Block size for memory allocation                                  */
/*********************************************************************/
#define MEMBLOCKSIZE 256

/*********************************************************************/
/*********************************************************************/
/*********************************************************************/

/*********************************************************************/
/* Static function forward declarations                              */
/*********************************************************************/


static inline 
void exchange_points(double *pt1, double *pt2, unsigned int ndim,
		     void *dptr1, void *dptr2, void *dswap,
		     size_t dsize);
/* Trivial routine to exchange two points.

   Parameters:
      INPUT/OUTPUT pt1
         position of first point to be swapped
      INPUT/OUTPUT pt2
         position of second point to be swapped
      INPUT ndim
         number of dimensions for each point
      INPUT/OUTPUT dptr1
         pointer to extra data for first point to be swapped
      INPUT/OUTPUT dptr2
         pointer to extra data for second point to be swapped
      INPUT/OUTPUT dswap
         swap space for exchanging dptr data; must be of at least size
         dsize
      INPUT dsize
         size of dptr1 and dptr2

   Returns:
      Nothing
*/

static inline
void partition_node(KDnode *node, unsigned int ndim, void *dswap,
		    size_t dsize);
/* This routine partitions the points within a node into two equal
   sub-lists. Parition happens along node->splitdim.

   Parameters:
      INPUT/OUTPUT node
         pointer to the node to be partitioned
      INPUT ndim
         number of dimensions in the data set
      INPUT/OUTPUT dswap
         swap space for exchanging dptr data; must be of at least size
         dsize
      INPUT dsize
         size of extra associated to points

   Returns:
      Nothing
*/

static inline
void xdata_alloc(unsigned int nnew, unsigned int *ncur, double **x, 
		 void **dptr, double **dist2, unsigned int ndim,
		 size_t dsize, bool shrink);
/* This is a utility routine to allocate data for returning lists of
   data points found by tree searches.

   Parameters:
      INPUT nnew
         new number of allocated elements
      INPUT/OUTPUT ncur
         pointer to current number of allocated elements; on return,
         this is >= nnew
      INPUT/OUTPUT x
         pointer to current position data
      INPUT/OUTPUT dptr
         pointer to current extra data
      INPUT ndim
         number of dimensions in the data set
      INPUT dsize
         size of each element of extra data
      INPUT shrink
         if shrink is false, then if nnew < *ncur, the routine does
         nothing; if it is true, then if nnew < *ncur, memory is
         shunken to the requested size

   Returns
      Nothing
*/


/*********************************************************************/
/*********************************************************************/
/*********************************************************************/

/*********************************************************************/
/* Functions definitions                                             */
/*********************************************************************/

/*********************************************************************/
/* Routine to get maximum of the squared distance between a point,   */
/* line, plane, or analogous higher-dimensional object and the       */
/* nearest point in a box                                            */
/*********************************************************************/
double box_max_dist2(const double *x, const double *xbox[2],
		     const unsigned int ndim_x, 
		     const unsigned int ndim_xbox,
		     const unsigned int *dim_x,
		     double *swap) {
  /* Here we make use of the fact that the most distant point is
     always a corner, so we just need to check all the corners */
  unsigned int i, j, mask, ptr;
  unsigned int ncorner=1;
  double d2max = 0.0;
  for (i=0; i<ndim_xbox; i++) ncorner = ncorner<<1;
  /* Loop over all corners */
  for (i=0; i<ncorner; i++) {
    mask = 1;
    for (j=0; j<ndim_xbox; j++, mask = mask<<1) {
      ptr = ((i & mask) != 0);
      swap[j] = xbox[ptr][j];
    }
    if (dim_x != NULL) {
      d2max = fmax(d2max, euclidean_metric2_dim(swap, x, ndim_xbox,
						ndim_x, dim_x));
    } else {
      d2max = fmax(d2max, euclidean_metric2(swap, x, ndim_xbox));
    }
  }
  return(d2max);
}


/*********************************************************************/
/* Routine to get minimum of the squared distance between a point,   */
/* line, plane, or analogous higher-dimensional object and the       */
/* nearest point in a box                                            */
/*********************************************************************/
double box_min_dist2(const double *x, const double *xbox[2],
		     const unsigned int ndim_x, 
		     const unsigned int ndim_xbox,
		     const unsigned int *dim_x) {
  unsigned int i;
  double tmp, dist2 = 0.0;
  if (dim_x != NULL) {
    for (i=0; i<ndim_x; i++) {
      if (x[i] < xbox[0][dim_x[i]]) {
	tmp = x[i] - xbox[0][dim_x[i]];
	dist2 += tmp*tmp;
      } else if (x[i] > xbox[1][dim_x[i]]) {
	tmp = x[i] - xbox[1][dim_x[i]];
	dist2 += tmp*tmp;
      }
    }
  } else {
    for (i=0; i<ndim_x; i++) {
      if (x[i] < xbox[0][i]) {
	tmp = x[i] - xbox[0][i];
	dist2 += tmp*tmp;
      } else if (x[i] > xbox[1][i]) {
	tmp = x[i] - xbox[1][i];
	dist2 += tmp*tmp;
      }
    }
  }
  return(dist2);
}

/*********************************************************************/
/* Routine to build the tree                                         */
/*********************************************************************/
KDtree* build_tree(double *x, unsigned int ndim, unsigned int npt,
		   unsigned int leafsize, void *dptr, size_t dsize) {

  unsigned int i, n, idx, curnode;
  double width;
  KDtree *tree;

  /* Allocate memory for the tree */
  if (!(tree = (KDtree *) malloc(sizeof(KDtree)))) {
    fprintf(stderr, "cluster_slug: error: unable to allocate memory in build_tree\n");
    exit(1);
  }

  /* Figure out the dimensions of the tree */
  tree->ndim = ndim;
  for (n=npt, tree->nodes=1, tree->leaves=1, tree->levels=1; 
       n>=leafsize; n=n>>1) {
    tree->leaves = tree->leaves << 1;
    tree->nodes += tree->leaves;
    tree->levels++;
  }
  tree->dsize = dsize;

  /* Allocate memory for the tree and for swap */
  if (!(tree->tree = (KDnode *) calloc(tree->nodes, sizeof(KDnode)))) {
    fprintf(stderr, "cluster_slug: error: unable to allocate memory in build_tree\n");
    exit(1);
  }
  tree->tree--; /* Offset by 1 so tree->tree[1] is the root node */
  for (n=1; n<=tree->nodes; n++) {
    if (!(tree->tree[n].xlim[0] = 
	  (double *) calloc(tree->ndim, sizeof(double)))) {
      fprintf(stderr, "cluster_slug: error: unable to allocate memory in build_tree\n");
      exit(1);
    }
    if (!(tree->tree[n].xlim[1] = 
	  (double *) calloc(tree->ndim, sizeof(double)))) {
      fprintf(stderr, "cluster_slug: error: unable to allocate memory in build_tree\n");
      exit(1);
    }
  }
  if (!(tree->xswap = (double *) calloc(tree->ndim, sizeof(double)))) {
    fprintf(stderr, "cluster_slug: error: unable to allocate memory in build_tree\n");
    exit(1);
  }
  if (dsize > 0) {
    if (!(tree->dswap = malloc(dsize))) {
      fprintf(stderr, "cluster_slug: error: unable to allocate memory in build_tree\n");
      exit(1);
    }
  } else {
    tree->dswap = NULL;
  }

  /* As a safety measure, initialize all nodes to have no points;
     this is important for breadth-first traversals of the tree,
     where we may encounter nodes with no points */
  for (i=1; i<=tree->nodes; i++) tree->tree[i].npt = 0;

  /* Initialize data in the root node */
  curnode = ROOT;
  tree->tree[curnode].npt = npt;
  tree->tree[curnode].x = x;
  tree->tree[curnode].dptr = dptr;
  /* Set bounding box for root node */
  for (i=0; i<ndim; i++)
    tree->tree[curnode].xlim[0][i] = tree->tree[curnode].xlim[1][i] = x[i];
  for (n=1; n<npt; n++) {
    for (i=0; i<ndim; i++) {
      tree->tree[curnode].xlim[0][i] 
	= fmin(tree->tree[curnode].xlim[0][i], x[n*tree->ndim+i]);
      tree->tree[curnode].xlim[1][i] 
	= fmax(tree->tree[curnode].xlim[1][i], x[n*tree->ndim+i]);
    }
  }

  /* Now build the rest of the tree */
  while (1) {

    /* See if we should be a leaf or a parent */
    if (tree->tree[curnode].npt <= leafsize) {

      /* We are a leaf; mark us as having no split, then move to next node */
      tree->tree[curnode].splitdim = -1;
      SETNEXT(curnode);

    } else {

      /* We are a parent */

      /* Figure out which dimension to split along */
      width = tree->tree[curnode].xlim[1][0] - 
	tree->tree[curnode].xlim[0][0];
      tree->tree[curnode].splitdim = 0;
      for (i=1; i<=tree->ndim; i++) {
	if (tree->tree[curnode].xlim[1][i] - 
	    tree->tree[curnode].xlim[0][i] > width) {
	  width = tree->tree[curnode].xlim[1][i] - 
	    tree->tree[curnode].xlim[0][i];
	  tree->tree[curnode].splitdim = i;
	}
      }

      /* Partition the points along the split dimension */
      partition_node(&(tree->tree[curnode]), tree->ndim, 
		     tree->dswap, dsize);

      /* Set the bounding box on the child nodes */
      for (i=0; i<tree->ndim; i++) {
	if (i==tree->tree[curnode].splitdim) {
	  /* In splitting dimension, split at middle index */
	  idx = ndim*((tree->tree[curnode].npt-1)/2) + 
	    tree->tree[curnode].splitdim;
	  tree->tree[LEFT(curnode)].xlim[0][i] = 
	    tree->tree[curnode].xlim[0][i];
	  tree->tree[LEFT(curnode)].xlim[1][i] = 
	    tree->tree[curnode].x[idx];
	  tree->tree[RIGHT(curnode)].xlim[0][i] = 
	    tree->tree[curnode].x[idx];
	  tree->tree[RIGHT(curnode)].xlim[1][i] =
	    tree->tree[curnode].xlim[1][i];
	} else {
	  /* Every other dimension just copies values */
	  tree->tree[LEFT(curnode)].xlim[0][i] = 
	    tree->tree[curnode].xlim[0][i];
	  tree->tree[LEFT(curnode)].xlim[1][i] = 
	    tree->tree[curnode].xlim[1][i];
	  tree->tree[RIGHT(curnode)].xlim[0][i] = 
	    tree->tree[curnode].xlim[0][i];
	  tree->tree[RIGHT(curnode)].xlim[1][i] = 
	    tree->tree[curnode].xlim[1][i];
	}
      }

      /* Set counters and pointers for child nodes */
      tree->tree[LEFT(curnode)].npt = (tree->tree[curnode].npt+1) / 2;
      tree->tree[RIGHT(curnode)].npt = tree->tree[curnode].npt / 2;
      tree->tree[LEFT(curnode)].x = tree->tree[curnode].x;
      tree->tree[RIGHT(curnode)].x = tree->tree[curnode].x + 
	tree->tree[LEFT(curnode)].npt*tree->ndim;
      tree->tree[LEFT(curnode)].dptr = tree->tree[curnode].dptr;
      tree->tree[RIGHT(curnode)].dptr = 
	((char *) tree->tree[curnode].dptr) +
	dsize*tree->tree[LEFT(curnode)].npt;

      /* Continue to next node */
      curnode = LEFT(curnode);
    }

    /* Check if we're done */
    if (curnode == ROOT) break;
  }

  /* Return */
  return tree;
}


/*********************************************************************/
/* Routine to free a tree                                            */
/*********************************************************************/
void free_tree(KDtree **tree) {
  int i;
  for (i=1; i<=(*tree)->nodes; i++) {
    free((*tree)->tree[i].xlim[0]);
    free((*tree)->tree[i].xlim[1]);
  }
  (*tree)->tree++; /* Undo offset by 1 */
  free((*tree)->tree);
  free((*tree)->xswap);
  if ((*tree)->dsize > 0) free((*tree)->dswap);
  free(*tree);
  *tree = NULL;
}


/*********************************************************************/
/* Routine to do nearest neighbor search of tree                     */
/*********************************************************************/
void neighbors(const KDtree *tree, const double *xpt, 
	       const unsigned int *dims, const unsigned int ndim, 
	       const unsigned int nneighbor, double *pos,
	       void *dptr, double *dist2) {
  unsigned int i, j, k;
  unsigned int curnode, ptr;
  unsigned int *dims_tmp;
  double d2;

  /* If necessary, create a temporary array to hold dimensions */
  if (ndim < tree->ndim) {
    dims_tmp = (unsigned int *) dims;  /* Cast to avoid compiler
					  warning */
  } else {
    if (!(dims_tmp = (unsigned int *) 
	  calloc(ndim, sizeof(unsigned int)))) {
      fprintf(stderr, "cluster_slug: error: unable to allocate memory in kdtree:neighbors\n");
      exit(1);
    }
    for (i=0; i<ndim; i++) dims_tmp[i] = i;
  }

  /* Initialize all distances to something big */
  for (i=0; i<nneighbor; i++) dist2[i] = DBL_MAX;

  /* Now begin to search the tree */
  curnode = ROOT;
  while (1) {

    /* See if any part of this node is within a distance less than the
       current maximum distance from the search point; if not, got to
       next node; side CS subtlety: the typecast in the next line is
       needed to avoid compiler warnings; see 
       http://c-faq.com/ansi/constmismatch.html */
    d2 = box_min_dist2(xpt, (const double **) tree->tree[curnode].xlim, 
		       ndim, tree->ndim, dims_tmp);
    if (d2 > dist2[nneighbor-1]) {
      SETNEXT(curnode);
      if (curnode==ROOT) break;
      continue;
    }

    /* If we're here, we need to investigate this node. If this node
       is not a leaf, just go to its left child. */
    if (tree->tree[curnode].splitdim != -1) {
      curnode = LEFT(curnode);
      continue;
    }

    /* If we're here, we have found a leaf that we need to
       investigate, so check each of the points it contains */
    for (i=0; i<tree->tree[curnode].npt; i++) {

      /* Is this point more distance than the most distance of the
	 current neighbors? If so, continue to next point. */
      d2 = euclidean_metric2_dim(&(tree->tree[curnode].x[tree->ndim*i]),
				 xpt, tree->ndim, ndim, dims_tmp);
      if (d2  > dist2[nneighbor-1]) continue;

      /* If we're here, we've found a point in this leaf that is
	 closer than any of our current neighbors. Figure out where
	 in the current neighbor list to insert it. */
      ptr = 0;
      while (dist2[ptr] < d2) ptr++;

      /* Move the rest of the current neighbor list back one space */
      for (j=nneighbor-1; j>ptr; j--) {
	for (k=0; k<tree->ndim; k++)
	  pos[tree->ndim*j+k] = pos[tree->ndim*(j-1)+k];
	dist2[j] = dist2[j-1];
	if (tree->dsize > 0) {
	  memcpy(((char *) dptr) + tree->dsize*j,
		 ((char *) dptr) + tree->dsize*(j-1),
		 tree->dsize);
	}
      }

      /* Insert the point we just found into the neighbor list */
      for (k=0; k<tree->ndim; k++)
	pos[tree->ndim*ptr+k] = tree->tree[curnode].x[tree->ndim*i+k];
      dist2[ptr] = d2;
      if (tree->dsize > 0) {
	memcpy(((char *) dptr) + tree->dsize*ptr,
	       ((char *) tree->tree[curnode].dptr) + tree->dsize*i,
	       tree->dsize);
      }
    }

    /* We're done with this leaf, so go on to next one */
    SETNEXT(curnode);
    if (curnode==ROOT) break;
  }

  /* Free memory */
  if (ndim > tree->ndim) free(dims_tmp);
}


/*********************************************************************/
/*********************************************************************/
/*********************************************************************/


/*********************************************************************/
/* Definitions of static functions                                   */
/*********************************************************************/

/*********************************************************************/
/* Routine to exchange a pair of points                              */
/*********************************************************************/
static inline
void exchange_points(double *pt1, double *pt2, unsigned int ndim,
		     void *dptr1, void *dptr2, void *dswap,
		     size_t dsize) {
  double tmp;
  unsigned int i;
  for (i=0; i<ndim; i++) {
    tmp = pt1[i];
    pt1[i] = pt2[i];
    pt2[i] = tmp;
  }
  if (dsize>0) {
    memcpy(dswap, dptr1, dsize);
    memcpy(dptr1, dptr2, dsize);
    memcpy(dptr2, dswap, dsize);
  }
}

/*********************************************************************/
/* Euclidean metric                                                  */
/*********************************************************************/

double euclidean_metric2(const double *x1, const double *x2,
			 const unsigned int ndim) {
  unsigned int i;
  double tmp, dist = 0.0;
  for (i=0; i<ndim; i++) {
    tmp = x1[i]-x2[i];
    dist += tmp*tmp;
  }
  return dist;
}

double euclidean_metric2_dim(const double *x1, const double *x2,
			     const unsigned int ndim1,
			     const unsigned int ndim2,
			     const unsigned int *dim2) {
  unsigned int i;
  double tmp, dist = 0.0;
  for (i=0; i<ndim2; i++) {
    tmp = x1[dim2[i]] - x2[i];
    dist += tmp*tmp;
  }
  return dist;
}

/*********************************************************************/
/* Routine to partition a node                                       */
/*********************************************************************/
void partition_node(KDnode *node, unsigned int ndim, void *dswap, 
		    size_t dsize) {

  unsigned int l, r, m, lptr, rptr;
  double splitval;

  /* Initialize pointers */
  l = 0;
  r = node->npt-1;
  m = (node->npt-1) >> 1;

  /* Iterate until partition is complete */
  while (l < r) {

    /* Initialize pointers */
    lptr = l;
    rptr = r;

    /* Choose a partition value for this iteration */
    splitval = node->x[ndim*((lptr+rptr)/2) + node->splitdim];

    /* Put points smaller than the partition value to the left, and
       points larger than the partition value to the right */
    while (1) {
      while ((node->x[ndim*lptr + node->splitdim] < splitval) &&
	     (lptr < rptr))
	lptr++;
      while ((node->x[ndim*rptr + node->splitdim] > splitval) &&
	     (lptr <= rptr))
	rptr--;
      if (lptr >= rptr) break;
      exchange_points(&(node->x[ndim*lptr]), 
		      &(node->x[ndim*rptr]), 
		      ndim,
		      ((char *) node->dptr) + dsize*lptr,
		      ((char *) node->dptr) + dsize*rptr,
		      dswap, dsize);
      lptr++;
      rptr--;
    }

    /* Set rptr to point to the first element larger than the
       partition value */
    if ((node->x[ndim*rptr + node->splitdim] > splitval) &&
	(lptr <= rptr)) 
      rptr--;
    rptr++;

    /* Set l and r to point to enclose the part of the list containing
       the middle index. If the middle index is the only thing left,
       we're done. */
    if (rptr>m) r = rptr-1;
    else l = rptr;
  }
}

/*********************************************************************/
/* Routine to check if a slab entirely contains a given node         */
/*********************************************************************/
bool query_contained_by_slab(const KDnode *node, const double *xslab, 
			     const unsigned int *dim_slab, 
			     unsigned int ndim_slab,
			     double offset, unsigned int ndim) {
  unsigned int i;
  bool contain = true;
  for (i=0; i<ndim_slab; i++) {
    contain = contain && (xslab[i] + offset 
			  >= node->xlim[1][dim_slab[i]]);
    contain = contain && (xslab[i] - offset
			  <= node->xlim[0][dim_slab[i]]);
  }
  return contain;
}

/*********************************************************************/
/* Routine to check if a sphere entirely contains a given node       */
/*********************************************************************/
bool query_contained_by_sphere(const KDnode *node, const double *xpt,
			       double radius, unsigned int ndim, 
			       double *swap) {
  /* Note: we make use of the fact that a cube is entirely contained
     within a sphere if and only if all its corners are inside the 
     sphere. */
  unsigned int i, j, mask, ptr;
  unsigned int ncorner=1;
  double r2 = radius*radius;
  for (i=0; i<ndim; i++) ncorner = ncorner<<1;
  /* Loop over all corners */
  for (i=0; i<ncorner; i++) {
    mask = 1;
    for (j=0; j<ndim; j++, mask = mask<<1) {
      ptr = ((i & mask) != 0);
      swap[j] = node->xlim[ptr][j];
    }
    if (euclidean_metric2(swap, xpt, ndim) > r2) return(false);
  }
  return(true);
}

/*********************************************************************/
/* Routine to check if a node entirely contains a sphere             */
/*********************************************************************/
bool query_contains_sphere(const KDnode *node, const double *xpt,
			   double radius, unsigned int ndim, 
			   double *swap) {
  /* Note: we make use of the fact that a sphere is entirely contained
     within a cube if and only if the cube that curcumscribes the
     sphere is also contained within the given cube. */
  unsigned int i;
  bool contain = true;
  for (i=0; i<ndim; i++) {
    contain = contain && (xpt[i] + radius <= node->xlim[1][i]);
    contain = contain && (xpt[i] - radius >= node->xlim[0][i]);
  }
  return contain;
}


/*********************************************************************/
/* Routine to check for intersections between a hyperslab and the    */
/* hypercube associated with this node.                              */
/*********************************************************************/
bool query_intersect_slab(const KDnode *node, const double *xslab, 
			  const unsigned int *dim_slab,
			  unsigned int ndim_slab,
			  double offset, unsigned int ndim) {
  unsigned int i;
  bool intersect = true;
  for (i=0; i<ndim_slab; i++) {
    intersect = intersect && (xslab[i] + offset
			      >= node->xlim[0][dim_slab[i]]);
    intersect = intersect && (xslab[i] - offset
			      <= node->xlim[1][dim_slab[i]]);
  }
  return intersect;
}


/*********************************************************************/
/* Routine to check for intersections of a sphere with the cube      */
/* associated with a node                                            */
/*********************************************************************/
bool query_intersect_sphere(const KDnode *node, const double *xpt, 
			    double radius, unsigned int ndim) {
  unsigned int i;
  bool intersect = true;
  for (i=0; i<ndim; i++) {
    intersect = intersect && (xpt[i] + radius >= node->xlim[0][i]);
    intersect = intersect && (xpt[i] - radius <= node->xlim[1][i]);
  }
  return intersect;
}

/*********************************************************************/
/* Routine to find all the points that lie within a certain distance */
/* of a specified hyperslab.                                        */
/*********************************************************************/
unsigned int query_slab(const KDtree *tree, const double *xslab, 
			const unsigned int *dim_slab, 
			unsigned int ndim_slab, 
			double offset, double **x, void **dptr,
			double **dist2) {

  unsigned int i, j, ptr;
  unsigned int npt = 0, nalloc = 0;
  unsigned int curnode = ROOT;
  char *dptr1, *dptr2;
  double d2;
  double offset2 = offset*offset;

  /* Begin search */
  while (1) {

    /* Check if anything within this node is within the required
       distance of the hyperslab.  If no intersection is found, go on
       to next node. If next node is root, we're done searching the
       entire tree. */
    if (!query_intersect_slab(&(tree->tree[curnode]), xslab,
			      dim_slab, ndim_slab, offset, 
			      tree->ndim)) {
      SETNEXT(curnode);
      if (curnode==ROOT) break;
      else continue;
    }

    /* If we're here, we have a potential intersection between the
       search region and this node. Now see if the entire node is
       within the search region. If it is, add all the points in this
       node to the list. */
    if (query_contained_by_slab(&(tree->tree[curnode]), xslab,
				dim_slab, ndim_slab, offset, 
				tree->ndim)) {

      /* Increment number of points and allocate memory */
      ptr = npt;
      npt += tree->tree[curnode].npt;
      xdata_alloc(npt, &nalloc, x, dptr, dist2, tree->ndim,
		  tree->dsize, false);

      /* Store data */
      for (i=0; i<tree->tree[curnode].npt; i++) {
	for (j=0; j<tree->ndim; j++) {
	  (*x)[(ptr+i)*tree->ndim + j] = 
	    tree->tree[curnode].x[tree->ndim*i+j];
	}
	for (j=0; j<ndim_slab; j++)
	  tree->xswap[j] = 
	    tree->tree[curnode].x[tree->ndim*i + dim_slab[j]];
	(*dist2)[ptr+1] = 
	  euclidean_metric2(tree->xswap, xslab, ndim_slab);
	if (tree->dsize > 0) {
	  dptr1 = (char *) (*dptr) + tree->dsize*(ptr+i);
	  dptr2 = (char *) tree->tree[curnode].dptr + tree->dsize*i;
	  memcpy(dptr1, dptr2, tree->dsize);
	}
      }

      /* Go to next node. If next node is root, we're done. */
      SETNEXT(curnode);
      if (curnode==ROOT) break;
      else continue;
    }

    /* If we're here, we've found a node that overlaps the search
       region but is not entirely contained by it. If this node is not
       a leaf, we simply move on to its children. If it is a leaf, we
       search all the points it contains to check if they are in the
       search region. */
    if (tree->tree[curnode].splitdim != -1) {

      /* This is not a leaf */
      curnode = LEFT(curnode);

    } else {

      /* This is a leaf, so loop over its points */
      for (i=0; i<tree->tree[curnode].npt; i++) {

	/* Is this point in the search region? */
	for (j=0; j<ndim_slab; j++)
	  tree->xswap[j] = 
	    tree->tree[curnode].x[tree->ndim*i + dim_slab[j]];
	d2 = euclidean_metric2(tree->xswap, xslab, ndim_slab);
	if (d2 < offset2) {

	  /* This point is in the search region, so we need to add it */
	  ptr = npt;
	  npt++;
	  xdata_alloc(npt, &nalloc, x, dptr, dist2, tree->ndim, 
		      tree->dsize, false);
	  for (j=0; j<tree->ndim; j++) {
	    (*x)[ptr*tree->ndim + j] = 
	      tree->tree[curnode].x[tree->ndim*i+j];
	  }
	  (*dist2)[ptr] = d2;
	  if (tree->dsize > 0) {
	    dptr1 = (char *) (*dptr) + tree->dsize*ptr;
	    dptr2 = (char *) tree->tree[curnode].dptr + tree->dsize*i;
	    memcpy(dptr1, dptr2, tree->dsize);
	  }
	}
      }

      /* Point to next node. If next node is ROOT, we're done */
      SETNEXT(curnode);
      if (curnode==ROOT) break;
      else continue;
    }
  }

  /* Shrink the memory allocated for the return values to the correct
     size */
  xdata_alloc(npt, &nalloc, x, dptr, dist2, tree->ndim, tree->dsize, 
	      true);

  /* Return */
  return npt;
}

/*********************************************************************/
/* Routine to find all points in the tree within a sphere of a       */
/* specified radius                                                  */
/*********************************************************************/
unsigned int query_sphere(const KDtree *tree, const double *xpt,
			  double radius, double **x, void **dptr, 
			  double **dist2) {

  unsigned int i, j, ptr;
  unsigned int npt = 0, nalloc = 0;
  unsigned int curnode=ROOT;
  char *dptr1, *dptr2;
  double r2;
  double rad2 = radius*radius;

  /* Begin search */
  while (1) {

    /* Check if anything within this node is within the required
       distance. If no intersection is found, go on to next node. If
       next node is root, we're done searching the entire tree. */
    if (!query_intersect_sphere(&(tree->tree[curnode]), xpt, radius, 
				tree->ndim)) {
      SETNEXT(curnode);
      if (curnode==ROOT) break;
      else continue;
    }

    /* If we're here, we have a potential intersection between the
       search region and this node. Now see if the entire node is
       within the search region. If it is, add all the points in this
       node to the list. */
    if (query_contained_by_sphere(&(tree->tree[curnode]), xpt, radius,
				  tree->ndim, tree->xswap)) {

      /* Increment number of points and allocate memory */
      ptr = npt;
      npt += tree->tree[curnode].npt;
      xdata_alloc(npt, &nalloc, x, dptr, dist2, tree->ndim, 
		  tree->dsize, false);

      /* Store data */
      for (i=0; i<tree->tree[curnode].npt; i++) {
	for (j=0; j<tree->ndim; j++) {
	  (*x)[(ptr+i)*tree->ndim + j] = 
	    tree->tree[curnode].x[tree->ndim*i+j];
	}
	(*dist2)[ptr+i] = 
	  euclidean_metric2(&(tree->tree[curnode].x[tree->ndim*i]),
			    xpt, tree->ndim);
	if (tree->dsize > 0) {
	  dptr1 = (char *) (*dptr) + tree->dsize*(ptr+i);
	  dptr2 = (char *) tree->tree[curnode].dptr + tree->dsize*i;
	  memcpy(dptr1, dptr2, tree->dsize);
	}
      }

      /* Go to next node. If next node is root, we're done. */
      SETNEXT(curnode);
      if (curnode==ROOT) break;
      else continue;
    }

    /* If we're here, we've found a node that overlaps the search
       region but is not entirely contained by it. If this node is not
       a leaf, we simply move on to its children. If it is a leaf, we
       search all the points it contains to check if they are in the
       search region. */
    if (tree->tree[curnode].splitdim != -1) {

      /* This is not a leaf */
      curnode = LEFT(curnode);

    } else {

      /* This is a leaf, so loop over its points */
      for (i=0; i<tree->tree[curnode].npt; i++) {

	/* Is this point in the search region? */
	r2 = euclidean_metric2(&(tree->tree[curnode].x[tree->ndim*i]),
			       xpt, tree->ndim);
	if (r2 < rad2) {

	  /* This point is in the search region, so we need to add it */
	  ptr = npt;
	  npt++;
	  xdata_alloc(npt, &nalloc, x, dptr, dist2, tree->ndim,
		      tree->dsize, false);
	  for (j=0; j<tree->ndim; j++) {
	    (*x)[ptr*tree->ndim + j] = 
	      tree->tree[curnode].x[tree->ndim*i+j];
	  }
	  (*dist2)[ptr] = r2;
	  if (tree->dsize > 0) {
	    dptr1 = (char *) (*dptr) + tree->dsize*ptr;
	    dptr2 = (char *) tree->tree[curnode].dptr + tree->dsize*i;
	    memcpy(dptr1, dptr2, tree->dsize);
	  }
	}
      }

      /* Point to next node. If next node is ROOT, we're done */
      SETNEXT(curnode);
      if (curnode==ROOT) break;
      else continue;
    }
  }

  /* Shrink the memory allocated for the return values to the correct
     size */
  xdata_alloc(npt, &nalloc, x, dptr, dist2, tree->ndim, tree->dsize,
	      true);

  /* Return */
  return npt;
}

/*********************************************************************/
/* Utility memory allocation routine                                 */
/*********************************************************************/
static inline
void xdata_alloc(unsigned int nnew, unsigned int *ncur, double **x, 
		 void **dptr, double **dist2, unsigned int ndim,
		 size_t dsize, bool shrink) {

  unsigned int nblock;
  size_t memsize;

  /* Do nothing if enough memory is already available and we aren't
     being requested to shrink */
  if ((nnew < *ncur) && (!shrink)) return;

  /* Figure out how many blocks we need */
  if (nnew > 0) {
    nblock = (nnew-1) / MEMBLOCKSIZE + 1;
  } else {
    nblock = 0;
  }

  /* Allocate */
  if (*ncur == 0) {
    /* Call malloc if no data is allocated yet */
    if (nblock > 0) {
      memsize = nblock*MEMBLOCKSIZE*ndim*sizeof(double);
      if (!(*x = malloc(memsize))) {
	fprintf(stderr, "cluster_slug: error: unable to allocate memory for KDtree return value\n");
	exit(1);
      }
      memsize = nblock*MEMBLOCKSIZE*sizeof(double);
      if (!(*dist2 = malloc(memsize))) {
	fprintf(stderr, "cluster_slug: error: unable to allocate memory for KDtree return value\n");
	exit(1);
      }
      if (dsize > 0) {
	memsize = nblock*MEMBLOCKSIZE*dsize;
	if (!(*dptr = malloc(memsize))) {
	  fprintf(stderr, "cluster_slug: error: unable to allocate memory for KDtree return value\n");
	  exit(1);
	}
      }
    }
  } else {
    /* Call realloc to resize existing memory */
    memsize = nblock*MEMBLOCKSIZE*ndim*sizeof(double);
    *x = realloc(*x, memsize);
    if ((x==NULL) && (memsize>0)) {
      fprintf(stderr, "cluster_slug: error: unable to allocate memory for KDtree search\n");
      exit(1);
    }
    memsize = nblock*MEMBLOCKSIZE*sizeof(double);
    *dist2 = realloc(*dist2, memsize);
    if ((dist2==NULL) && (memsize>0)) {
      fprintf(stderr, "cluster_slug: error: unable to allocate memory for KDtree search\n");
      exit(1);
    }
    if (dsize > 0) {
      memsize = nblock*MEMBLOCKSIZE*dsize;
      *dptr = realloc(*dptr, memsize);
      if ((dptr==NULL) && (memsize>0)) {
	fprintf(stderr, "cluster_slug: error: unable to allocate memory for KDtree search\n");
	exit(1);
      }
    }
  }

  /* Record amount now allocated */
  *ncur = nblock*MEMBLOCKSIZE;
}
