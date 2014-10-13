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

/*********************************************************************/
/* This file contains routines to build and search KD trees          */
/*********************************************************************/

#ifndef _KDTREE_H_
#define _KDTREE_H_

#include <stdbool.h>
#include <stdlib.h>

/*********************************************************************/
/* Macros for traversing KD trees                                    */
/*********************************************************************/
#define ROOT 1
#define LEFT(i) (i<<1)
#define RIGHT(i) ((i<<1)+1)
#define PARENT(i) (i>>1)
#define SIBLING(i) ((i&1)?i-1:i+1)
#define SETNEXT(i) { while (i&1) i=i>>1; ++i; }

/*********************************************************************/
/* Structures for a node and a tree                                  */
/*********************************************************************/
typedef struct {
  /* Number of points in this node */
  unsigned int npt;
  /* Corners of box bounding region occupied by this node and its
     children. */
  double *xlim[2];
  /* Dimension along which this node splits. -1 for a leaf. */
  int splitdim;
  /* Pointer to the positions for this node */
  double *x;
  /* Extra data associated to the positions */
  void *dptr;
} KDnode;

typedef struct {
  /* Number of dimensions in the tree */
  unsigned int ndim;
  /* Number of levels in the tree */
  unsigned int levels;
  /* Number of leaves in the tree */
  unsigned int leaves;
  /* Number of nodes in the tree */
  unsigned int nodes;
  /* Size of extra data elements associated with positions */
  size_t dsize;
  /* Swap space */
  double *xswap;
  void *dswap;
  /* Pointer to tree root */
  KDnode *tree;
} KDtree;

/*********************************************************************/
/* Function definitions                                              */
/*********************************************************************/

KDtree* build_tree(double *x, unsigned int ndim, unsigned int npt, 
		   unsigned int leafsize, void *dptr, size_t dsize);
/* This routine builds a KD tree from the input data.

   Parameters:
      INPUT/OUTPUT x
         array of npt * ndim elements containing the positions,
         ordered so that element x[j + i*ndim] is the jth coordinate
         of point i; on return, this will be sorted into a tree
      INPUT ndim
         number of dimensions in the data set
      INPUT npt
         number of data points
      INPUT leafsize
         number of points to place in a leaf of the tree
      INPUT/OUTPUT dptr
         pointer to extra data associated to each position
      INPUT dsize
         size of each element of dptr

   Returns:
      tree : OUTPUT
         a KD tree decomposition of the data
*/


void free_tree(KDtree **tree);
/* Frees the memory associated with a KD tree.

   Parameters
      INPUT/OUTPUT tree
         The KDtree to be de-allocated

   Returns
      Nothing
*/


void neighbors(const KDtree *tree, const double *xpt, 
	       const unsigned int *dims, const unsigned int ndim, 
	       const unsigned int nneighbor, double *pos,
	       void *dptr, double *dist2);
/* Routine to find the N nearest neighbors to an input point; input
   points can have fewer dimensions than the search space, in which
   case the routine searches for the nearest neighbors to a line,
   plane, or higher-dimensional object.

   Parameters:
      INPUT tree
         the KD tree to be searched
      INPUT xpt
         array of ndim elements giving the position of the search
         point
      INPUT dims
         array specifying which dimensions in the kd tree are given in
         xpt; e.g., if the KD tree has tree dimensions, and dims = [0,
         2], then the two elements of x will be interpreted as
         specifying x and z coordinates, and the routine will search
         for the closest neighbors to a line as the specified x and
         z. If ndim is equal to the number of dimensions in the KD
         tree, this argument may be left as NULL.
      INPUT ndim
         number of elements in x and dims
      INPUT nneighbor
         the number of neighbors to find
      OUTPUT x
         positions of the nearest neighbor points; on entry, this
         pointer must point to a block of at least
         tree->ndim*nneighbor elements, and on return element
	 x[i*tree->ndim+j] contains the jth coordinate for the ith
         neighbor found; points are sorted by distance from xpt
      OUTPUT dptr
         extra data associated with each of the nearest neighbor
         points in x; element dptr[i] is the extra data for the ith
         point; must point to a block of valid memory at least i
         elements long
      OUTPUT dist2
         squared distances of all particles found from xpt; on entry,
         this pointer mut point to a block of nneighbor elements, and
         on return dist2[i] gives the distance from the ith point
         found to xpt
*/


inline
bool query_contained_by_slab(const KDnode *node, const double *xslab, 
			     const unsigned int *dim_slab, 
			     unsigned int ndim_slab,
			     double offset, unsigned int ndim);
/* Returns true if this node is entirely contained within a slab of
   the specified position and thickness.

   Parameters:
      INPUT node
         pointer to the node being queried
      INPUT xslab
         An array giving the coordinates of the center of the
         hyperslab. The dimensions whose coordinates are being
         specified is given by the dim_slab parameter.
      INPUT dim_slab
         An array specifying the dimensions in which the hyperslab
         lies. These refer to the dimensions of the data set indexed
         by the KD tree, and start at 0. Each element of dim_slab
         must be in the range 0 ... ndim-1, where ndim is the number
         of dimensions in the KD tree data set.
      INPUT ndim_slab
         The number of dimensions of the hyperslab being searched,
         i.e. the number of elements in xslab and dim_slab.
      INPUT offset
         half-thickness of the slab around the slab
      INPUT ndim
         number of dimensions in the data

   Returns:
      OUTPUT contain
         true if the node is entirely contained inside the slab, false
         otherwise
*/


inline
bool query_contained_by_sphere(const KDnode *node, const double *xpt, 
			       double radius, unsigned int ndim, 
			       double *swap);
/* Returns true if this node is entirely contained within a sphere of
   the specified radius centered on the input point.

   Parameters:
      INPUT node
         pointer to the node being queried
      INPUT xpt
         center of the search sphere
      INPUT radius
         radius of the search sphere
      INPUT ndim
         number of dimensions in the data
      INPUT/OUTPUT swap
         an array of ndim elements to use as swap space

   Returns:
      OUTPUT contain
         true if the node is entirely contained in the sphere, false
         otherwise
*/


inline
bool query_contains_sphere(const KDnode *node, const double *xpt, 
			   double radius, unsigned int ndim, 
			   double *swap);
/* Returns true if this node is entirely contains a sphere of the
   specified radius centered on the input point.

   Parameters:
      INPUT node
         pointer to the node being queried
      INPUT xpt
         center of the search sphere
      INPUT radius
         radius of the search sphere
      INPUT ndim
         number of dimensions in the data
      INPUT/OUTPUT swap
         an array of ndim elements to use as swap space

   Returns:
      OUTPUT contain
         true if the node entirely contains the sphere, false
         otherwise
*/


inline
bool query_intersect_slab(const KDnode *node, const double *xslab, 
			  const unsigned int *dim_slab,
			  unsigned int ndim_slab,
			  double offset, unsigned int ndim);
/* Returns true if the cube associatd with this node intersects a slab
   of finite thickness centered on a specified slab, false otherwise.

   Parameters:
      INPUT node
         pointer to the node being queried
      INPUT xslab
         An array giving the coordinates of the center of the
         hyperslab. The dimensions whose coordinates are being
         specified is given by the dim_slab parameter.
      INPUT dim_slab
         An array specifying the dimensions in which the hyperslab
         lies. These refer to the dimensions of the data set indexed
         by the KD tree, and start at 0. Each element of dim_slab
         must be in the range 0 ... ndim-1, where ndim is the number
         of dimensions in the KD tree data set.
      INPUT ndim_slab
         The number of dimensions of the hyperslab being searched,
         i.e. the number of elements in xslab and dim_slab.
      INPUT offset
         half-thickness of the slab around the slab
      INPUT ndim
         number of dimensions in the data

   Returns:
      OUTPUT contain
         true if the node intersects the slab, false otherwise
*/

inline
bool query_intersect_sphere(const KDnode *node, const double *xpt, 
			    double radius, unsigned int ndim);
/* Returns true if the cube associatd with this node intersects a sphere of
   the specified radius centered on the input point.

   Parameters:
      INPUT node
         pointer to the node being queried
      INPUT xpt
         center of the search sphere
      INPUT radius
         radius of the search sphere
      INPUT ndim
         number of dimensions in the data

   Returns:
      OUTPUT contain
         true if the node is intersects the sphere, false otherwise
*/

unsigned int query_slab(const KDtree *tree, const double *xslab, 
			const unsigned int *dim_slab, 
			unsigned int ndim_slab, 
			double offset, double **x, void **dptr,
			double **dist2);
/* This routine searches the KD tree and finds all the points that lie
   within a certain distance of a specified hyperslab, aligned with
   the coordinate axes of the data in the KD tree.

   Parameters:
      INPUT tree
         The KD tree to be searched
      INPUT xslab
         An array giving the coordinates of the center of the
         hyperslab. The dimensions whose coordinates are being specified
         is given by the dim_slab parameter.
      INPUT dim_slab
         An array specifying the dimensions in which the hyperslab
         lies. These refer to the dimensions of the data set indexed
         by the KD tree, and start at 0. Each element of dim_slab
         must be in the range 0 ... ndim-1, where ndim is the number
         of dimensions in the KD tree data set.
      INPUT ndim_slab
         The number of dimensions of the hyperslab being searched,
         i.e. the number of elements in xslab and dim_slab.
      INPUT offset
         The extent of the region around the hyperslab to be
         searched. Points are found if their perpendicular distance to
         the hyperslab is <= offset.
      OUTPUT x
         pointer to an array containing the locations of all the
         points found; on output (*x)[tree->ndim*i + j] is the jth
         coordinate of point i; (*x) should not point to valid memory
         on input, as the required memory is allocated within the
         routine.
      OUTPUT dptr
         pointer to an array of extra data associated with the points
         found; on output ((char *) (*dptr)) + i*tree->dsize is a
         pointer to the start of the data for point i; (*dptr) should
         not point to valid memory on input, as the required memory is
         allocated within the routine
      OUTPUT dist2
         pointer to an array containing the square of the distance
         from each point returned by the query to the center of the
         hyperslab; (*dist2) should not point to valid memory on
         input, as the required memory is allocated within the routine

   Returns:
      OUTPUT npt
         the number of points returned by the search
*/


unsigned int query_sphere(const KDtree *tree, const double *xpt, 
			  double radius, double **x, void **dptr,
			  double **dist2);
/* This routine searchers a KD tree to find all points within a
   specified distance of an input data point. It returns the positions
   all the points found and copies of their associated data.

   Parameters:
      INPUT tree
         pointer to the KDtree to be searched
      INPUT xpt
         center of the search region; must be an array of tree->ndim
         elements
      INPUT radius
         search radius
      OUTPUT x
         pointer to an array containing the locations of all the
         points found; on output (*x)[tree->ndim*i + j] is the jth
         coordinate of point i; (*x) should not point to valid memory
         on input, as the required memory is allocated within the
         routine.
      OUTPUT dptr
         pointer to an array of extra data associated with the points
         found; on output ((char *) (*dptr)) + i*tree->dsize is a
         pointer to the start of the data for point i; (*dptr) should
         not point to valid memory on input, as the required memory is
         allocated within the routine
      OUTPUT dist2
         pointer to an array containing the square of the distance
         from each point returned by the query to the center of the
         sphere; (*dist2) should not point to valid memory on input,
         as the required memory is allocated within the routine
   Returns:
      OUTPUT npt
         the number of points returned by the search
*/

#endif
/* _KDTREE_H_ */
