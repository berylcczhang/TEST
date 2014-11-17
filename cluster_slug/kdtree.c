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
#include "geometry.h"
#include "kdtree.h"


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
/* Routine to partition a node                                       */
/*********************************************************************/
static inline
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
/* Utility memory allocation routine                                 */
/*********************************************************************/
#define MEMBLOCKSIZE 256
static inline
void xdata_alloc(unsigned int nnew, unsigned int *ncur, double **x, 
		 void **dptr, double **d2, unsigned int ndim,
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
      if (!(*x = (double *) malloc(memsize))) {
	fprintf(stderr, "cluster_slug: error: unable to allocate memory for KDtree return value\n");
	exit(1);
      }
      memsize = nblock*MEMBLOCKSIZE*sizeof(double);
      if (!(*d2 = (double *) malloc(memsize))) {
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
    *x = (double *) realloc(*x, memsize);
    if ((x==NULL) && (memsize>0)) {
      fprintf(stderr, "cluster_slug: error: unable to allocate memory for KDtree search\n");
      exit(1);
    }
    memsize = nblock*MEMBLOCKSIZE*sizeof(double);
    *d2 = (double *) realloc(*d2, memsize);
    if ((d2==NULL) && (memsize>0)) {
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
#undef MEMBLOCKSIZE



/*********************************************************************/
/* Routine to build the tree                                         */
/*********************************************************************/
KDtree* build_tree(double *x, unsigned int ndim, unsigned int npt,
		   unsigned int leafsize, void *dptr, size_t dsize) {

  unsigned int i, n, idx, curnode;
  KDtree *tree;
  void *dswap = NULL;

  /* Allocate memory for the tree */
  if (!(tree = (KDtree *) malloc(sizeof(KDtree)))) {
    fprintf(stderr, "cluster_slug: error: unable to allocate memory in build_tree\n");
    exit(1);
  }

  /* Allocate swap space */
  if (dsize > 0) {
    if (!(dswap = malloc(dsize))) {
      fprintf(stderr, "cluster_slug: error: unable to allocate memory in build_tree\n");
      exit(1);
    }
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

  /* Allocate memory for the tree  */
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
    if (!(tree->tree[n].xbnd[0] = 
	  (double *) calloc(tree->ndim, sizeof(double)))) {
      fprintf(stderr, "cluster_slug: error: unable to allocate memory in build_tree\n");
      exit(1);
    }
    if (!(tree->tree[n].xbnd[1] = 
	  (double *) calloc(tree->ndim, sizeof(double)))) {
      fprintf(stderr, "cluster_slug: error: unable to allocate memory in build_tree\n");
      exit(1);
    }
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
  for (i=0; i<ndim; i++) {
    tree->tree[curnode].xbnd[0][i] = tree->tree[curnode].xlim[0][i];
    tree->tree[curnode].xbnd[1][i] = tree->tree[curnode].xlim[1][i];
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
      if (curnode != ROOT)
	tree->tree[curnode].splitdim = 
	  (tree->tree[PARENT(curnode)].splitdim + 1) % tree->ndim;
      else
	tree->tree[curnode].splitdim = 0;

      /* Partition the points along the split dimension */
      partition_node(&(tree->tree[curnode]), tree->ndim, dswap, dsize);

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

      /* Get actual data limits on child nodes */
      if (tree->tree[LEFT(curnode)].npt > 0) {
	for (i=0; i<tree->ndim; i++) 
	  tree->tree[LEFT(curnode)].xbnd[0][i] =
	    tree->tree[LEFT(curnode)].xbnd[1][i] = 
	    tree->tree[LEFT(curnode)].x[i];
	for (n=1; n<tree->tree[LEFT(curnode)].npt; n++) {
	  for (i=0; i<ndim; i++) {
	    tree->tree[LEFT(curnode)].xbnd[0][i] 
	      = fmin(tree->tree[LEFT(curnode)].xbnd[0][i], 
		     tree->tree[LEFT(curnode)].x[n*tree->ndim+i]);
	    tree->tree[LEFT(curnode)].xbnd[1][i] 
	      = fmax(tree->tree[LEFT(curnode)].xbnd[1][i], 
		     tree->tree[LEFT(curnode)].x[n*tree->ndim+i]);
	  }
	}
      }
      if (tree->tree[RIGHT(curnode)].npt > 0) {
	for (i=0; i<tree->ndim; i++) 
	  tree->tree[RIGHT(curnode)].xbnd[0][i] =
	    tree->tree[RIGHT(curnode)].xbnd[1][i] = 
	    tree->tree[RIGHT(curnode)].x[i];
	for (n=1; n<tree->tree[RIGHT(curnode)].npt; n++) {
	  for (i=0; i<ndim; i++) {
	    tree->tree[RIGHT(curnode)].xbnd[0][i] 
	      = fmin(tree->tree[RIGHT(curnode)].xbnd[0][i], 
		     tree->tree[RIGHT(curnode)].x[n*tree->ndim+i]);
	    tree->tree[RIGHT(curnode)].xbnd[1][i] 
	      = fmax(tree->tree[RIGHT(curnode)].xbnd[1][i], 
		     tree->tree[RIGHT(curnode)].x[n*tree->ndim+i]);
	  }
	}
      }

      /* Continue to next node */
      curnode = LEFT(curnode);
    }

    /* Check if we're done */
    if (curnode == ROOT) break;
  }

  /* Free memory */
  if (dswap != NULL) free(dswap);

  /* Return */
  return tree;
}


/*********************************************************************/
/* Routine to free a tree                                            */
/*********************************************************************/
void free_tree(KDtree *tree) {
  unsigned int i;
  for (i=1; i<=tree->nodes; i++) {
    free(tree->tree[i].xlim[0]);
    free(tree->tree[i].xlim[1]);
    free(tree->tree[i].xbnd[0]);
    free(tree->tree[i].xbnd[1]);
  }
  tree->tree++; /* Undo offset by 1 */
  free(tree->tree);
  free(tree);
}


/*********************************************************************/
/* Neighbor search routine                                           */
/*********************************************************************/
void neighbors(const KDtree *tree, const double *xpt, 
	       const unsigned int *dims, const unsigned int ndim, 
	       const unsigned int nneighbor,
	       const double *scale, double *pos,
	       void *dptr, double *d2) {
  unsigned int i, j, k;
  unsigned int curnode, ptr;
  double r2;

  /* Initialize all distances to something big */
  for (i=0; i<nneighbor; i++) d2[i] = DBL_MAX;

  /* Now begin to search the tree */
  curnode = ROOT;
  while (1) {

    /* See if any part of this node is within a distance less than the
       current maximum distance from the search point; if not, got to
       next node; side CS subtlety: the typecast in the next line is
       needed to avoid compiler warnings; see 
       http://c-faq.com/ansi/constmismatch.html */
    r2 = box_min_dist2(xpt, (const double **) tree->tree[curnode].xbnd, 
		       ndim, tree->ndim, dims, NULL, scale, ndim);

    if (r2 > d2[nneighbor-1]) {
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

      /* Is this point more distant than the most distance of the
	 current neighbors? If so, continue to next point. */
      r2 = dist2(xpt, &(tree->tree[curnode].x[tree->ndim*i]), ndim,
		 tree->ndim, dims, NULL, scale, tree->ndim);
      if (r2  > d2[nneighbor-1]) continue;

      /* If we're here, we've found a point in this leaf that is
	 closer than any of our current neighbors. Figure out where
	 in the current neighbor list to insert it. */
      ptr = 0;
      while (d2[ptr] < r2) ptr++;

      /* Move the rest of the current neighbor list back one space */
      for (j=nneighbor-1; j>ptr; j--) {
	for (k=0; k<tree->ndim; k++)
	  pos[tree->ndim*j+k] = pos[tree->ndim*(j-1)+k];
	d2[j] = d2[j-1];
	if (tree->dsize > 0) {
	  memcpy(((char *) dptr) + tree->dsize*j,
		 ((char *) dptr) + tree->dsize*(j-1),
		 tree->dsize);
	}
      }

      /* Insert the point we just found into the neighbor list */
      for (k=0; k<tree->ndim; k++)
	pos[tree->ndim*ptr+k] = tree->tree[curnode].x[tree->ndim*i+k];
      d2[ptr] = r2;
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
}


/*********************************************************************/
/* Routine to search a box                                           */
/*********************************************************************/
unsigned int query_box(const KDtree *tree, const double *xbox[2], 
		       unsigned int ndim, const unsigned int *dim, 
		       const double *scale, double **x, void **dptr,
		       double **d2) {
  unsigned int i, j, ptr;
  unsigned int npt = 0, nalloc = 0;
  unsigned int curnode = ROOT;
  char *dptr1, *dptr2;
  double *xcen;
  bool contained;

  /* Allocate temporary space and store box center in it */
  if (!(xcen = (double *) calloc(tree->ndim, sizeof(double)))) {
    fprintf(stderr, "cluster_slug: error: unable to allocate memory in query_box\n");
    exit(1);
  }
  for (j=0; j<ndim; j++) xcen[j] = 0.5*(xbox[0][j] + xbox[1][j]);

  /* Begin search */
  while (1) {

    /* Check if anything within this node is within the box. If no
       intersection is found, go on to next node. If next node is
       root, we're done searching the entire tree. */
    if (!box_intersects_box((const double **) tree->tree[curnode].xbnd,
			    xbox, tree->ndim, ndim, NULL, dim)) {
      SETNEXT(curnode);
      if (curnode==ROOT) break;
      else continue;
    }

    /* If we're here, we have a potential intersection between the
       search region and this node. Now see if the entire node is
       within the search region. If it is, add all the points in this
       node to the list. */
    if (box_in_box((const double **) tree->tree[curnode].xbnd,
		   xbox, tree->ndim, ndim, NULL, dim)) {

      /* Increment number of points and allocate memory */
      ptr = npt;
      npt += tree->tree[curnode].npt;
      xdata_alloc(npt, &nalloc, x, dptr, d2, tree->ndim,
		  tree->dsize, false);

      /* Store data */
      for (i=0; i<tree->tree[curnode].npt; i++) {
	for (j=0; j<tree->ndim; j++) {
	  (*x)[(ptr+i)*tree->ndim + j] = 
	    tree->tree[curnode].x[tree->ndim*i+j];
	}
	(*d2)[ptr+1] = dist2(&(tree->tree[curnode].x[tree->ndim*i+j]),
			     xcen, tree->ndim, ndim, NULL, dim,
			     scale, tree->ndim);
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
	contained = true;
	for (j=0; j<ndim; j++) {
	  contained = contained &&
	    (tree->tree[curnode].x[tree->ndim*i + dim[j]] > xbox[0][j]);
	  contained = contained &&
	    (tree->tree[curnode].x[tree->ndim*i + dim[j]] < xbox[1][j]);
	}
	if (!contained) continue;

	/* This point is in the search region, so we need to add it */
	ptr = npt;
	npt++;
	xdata_alloc(npt, &nalloc, x, dptr, d2, tree->ndim, 
		    tree->dsize, false);
	for (j=0; j<tree->ndim; j++) {
	  (*x)[ptr*tree->ndim + j] = 
	    tree->tree[curnode].x[tree->ndim*i+j];
	}
	(*d2)[ptr] = dist2(&(tree->tree[curnode].x[tree->ndim*i+j]),
			   xcen, tree->ndim, ndim, NULL, dim,
			   scale, tree->ndim);
	if (tree->dsize > 0) {
	  dptr1 = (char *) (*dptr) + tree->dsize*ptr;
	  dptr2 = (char *) tree->tree[curnode].dptr + tree->dsize*i;
	  memcpy(dptr1, dptr2, tree->dsize);
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
  xdata_alloc(npt, &nalloc, x, dptr, d2, tree->ndim, tree->dsize, 
	      true);

  /* Free temporary */
  free(xcen);

  /* Return */
  return npt;		   
}

/*********************************************************************/
/* Routine to search a sphere                                        */
/*********************************************************************/
unsigned int query_sphere(const KDtree *tree, const double *xcen,
			  unsigned int ndim, const unsigned int *dim, 
			  const double radius, const double *scale,
			  double **x, void **dptr, double **d2) {
  unsigned int i, j, ptr;
  unsigned int npt = 0, nalloc = 0;
  unsigned int curnode=ROOT;
  char *dptr1, *dptr2;
  double r2;
  double *xbox[2];
  double rad2 = radius*radius;

  /* For convenience, define a box that encloses the search search
     sphere. This will be used to discard nodes that don't need to be
     searched. */
  if (!(xbox[0] = calloc(ndim, sizeof(double)))) {
    fprintf(stderr, "cluster_slug: error: unable to allocate memory in query_sphere\n");
    exit(1);
  }
  if (!(xbox[1] = calloc(ndim, sizeof(double)))) {
    fprintf(stderr, "cluster_slug: error: unable to allocate memory in query_sphere\n");
    exit(1);
  }
  if (scale == NULL) {
    /* Case with no dimensional scaling */
    for (i=0; i<ndim; i++) {
      xbox[0][i] = xcen[i] - radius;
      xbox[1][i] = xcen[i] + radius;
    }
  } else {
    /* Case with dimensional scaling */
    if (dim == NULL) {
      /* No dimensional mapping specified */
      for (i=0; i<ndim; i++) {
	xbox[0][i] = xcen[i] - radius/scale[i];
	xbox[1][i] = xcen[i] + radius/scale[i];
      }
    } else {
      /* Dimensional mapping specified */
      for (i=0; i<ndim; i++) {
	xbox[0][i] = xcen[i] - radius/scale[dim[i]];
	xbox[1][i] = xcen[i] + radius/scale[dim[i]];
      }
    }
  }

  /* Begin search */
  while (1) {

    /* Check if the box for this node intersects the search region. If
       no intersection is found, go on to next node. If next node is
       root, we're done searching the entire tree. */
    if (!(box_intersects_box((const double **) tree->tree[curnode].xbnd, 
			     (const double **) xbox, 
			     tree->ndim, ndim, NULL, dim))) {
      SETNEXT(curnode);
      if (curnode==ROOT) break;
      else continue;
    }

    /* If we're here, we have a potential intersection between the
       search region and this node. Now see if the entire node is
       within the search region. If it is, add all the points in this
       node to the list. */
    if (box_in_sphere((const double **) &(tree->tree[curnode].xbnd),
		      xcen, tree->ndim, ndim, NULL, dim, radius,
		      scale, tree->ndim)) {

      /* Increment number of points and allocate memory */
      ptr = npt;
      npt += tree->tree[curnode].npt;
      xdata_alloc(npt, &nalloc, x, dptr, d2, tree->ndim, 
		  tree->dsize, false);

      /* Store data */
      for (i=0; i<tree->tree[curnode].npt; i++) {
	for (j=0; j<tree->ndim; j++) {
	  (*x)[(ptr+i)*tree->ndim + j] = 
	    tree->tree[curnode].x[tree->ndim*i+j];
	}
	(*d2)[ptr+i] = 
	  dist2(&(tree->tree[curnode].x[tree->ndim*i]), xcen,
		tree->ndim, ndim, NULL, dim, scale, tree->ndim);
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
	r2 = dist2(&(tree->tree[curnode].x[tree->ndim*i]), xcen,
		   tree->ndim, ndim, NULL, dim, scale, tree->ndim);
	if (r2 < rad2) {

	  /* This point is in the search region, so we need to add it */
	  ptr = npt;
	  npt++;
	  xdata_alloc(npt, &nalloc, x, dptr, d2, tree->ndim,
		      tree->dsize, false);
	  for (j=0; j<tree->ndim; j++) {
	    (*x)[ptr*tree->ndim + j] = 
	      tree->tree[curnode].x[tree->ndim*i+j];
	  }
	  (*d2)[ptr] = r2;
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
  xdata_alloc(npt, &nalloc, x, dptr, d2, tree->ndim, tree->dsize,
	      true);

  /* Free the inscribing box */
  free(xbox[0]);
  free(xbox[1]);

  /* Return */
  return npt;
}
