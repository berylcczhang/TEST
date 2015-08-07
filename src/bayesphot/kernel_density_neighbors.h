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
/* This module contains routines to do neighbor searching on kernel  */
/* density estimation structures.                                    */
/*********************************************************************/

#ifndef _KERNEL_DENSITY_NEIGHBORS_H_
#define _KERNEL_DENSITY_NEIGHBORS_H_

#include <stdbool.h>
#include "kernel_density.h"

void kd_neighbors(const kernel_density *kd, const double *xpt, 
		  const unsigned long *dims, const unsigned long ndim, 
		  const unsigned long nneighbor,
		  const bool bandwidth_units, double *pos,
		  void *dptr, double *d2);
/* Routine to find the N nearest neighbors to an input point; input
   points can have fewer dimensions than the search space, in which
   case the routine searches for the nearest neighbors to a line,
   plane, or higher-dimensional object.

   Parameters:
      INPUT kd
         the kernel density object to be searched
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
      INPUT bandwidth_units
         if true, the metric used to determine relative distance is
         normalized to the dimension-dependent bandwidth; if false,
         the metric is a simple Euclidean one
      OUTPUT pos
         positions of the nearest neighbor points; on entry, this
         pointer must point to a block of at least
         tree->ndim*nneighbor elements, and on return element
	 x[i*tree->ndim+j] contains the jth coordinate for the ith
         neighbor found; points are sorted by distance from xpt
      OUTPUT dptr
         extra data associated with each of the nearest neighbor
         points in x; element dptr[i] is the extra data for the ith
         point; must point to a block of valid memory at least
         nneighbor elements long
      OUTPUT d2
         squared distances of all particles found from xpt; on entry,
         this pointer must point to a block of nneighbor elements, and
         on return dist2[i] gives the distance from the ith point
         found to xpt

   Returns
      Nothing
*/

void kd_neighbors_all(const kernel_density *kd, 
		      const unsigned long nneighbor, 
		      const bool bandwidth_units, unsigned long *idx, 
		      double *d2);
/* This routine returns the indices of the N nearest neighbors for
   every point in the data set indexed by the kernel density
   object. Points are not considered their own neighbors.

   Parameters:
      INPUT kd
         the kernel density object to be searched
      INPUT nneighbor
         the number of neighbors to find
      INPUT bandwidth_units
         if true, the metric used to determine relative distance is
         normalized to the dimension-dependent bandwidth; if false,
         the metric is a simple Euclidean one
      OUTPUT idx
         indices of the neighbors found element idx[i*nneighbor+j]
         gives the index of the jth nearest neighbor of the ith point,
	 where point numbers and indices refer to the array x that is
         indexed by the kernel density object; must point to
         nneighbor * kd->tree->npt valid memory elements on entry
      OUTPUT d2
         squared distances of all neighbors found, indexed in the same
         way as idx; as with idx, this array must point to nneighbor *
         kd->tree->npt valid memory elements on entry

   Returns
      Nothing
*/

void kd_neighbors_point(const kernel_density *kd, 
			const unsigned long idxpt, 
			const unsigned long nneighbor,
			const bool bandwidth_units,
			unsigned long *idx, double *d2);
/* This routine finds the N nearest neighbors of a data point in the
   data set. Points are not considered their own neighbors.

   Parameters:
      INPUT kd
         the kernel density object to be searched
      INPUT idxpt
         the index of the point whose neighbors are to be found, in
         the array x that holds the kernel density data set
      INPUT nneighbor
         the number of neighbors to find
      INPUT bandwidth_units
         if true, the metric used to determine relative distance is
         normalized to the dimension-dependent bandwidth; if false,
         the metric is a simple Euclidean one
      OUTPUT idx
         indices of the neighbors found, where indices refer to the
         array x that is indexed by the kernel density object; must
         point to nneighbor valid memory elements on entry
      OUTPUT d2
         squared distances of all neighbors found; must point to
         nneighbor valid memory elements on entry

   Returns
      Nothing
*/

void kd_neighbors_point_vec(const kernel_density *kd, 
			    const unsigned long *idxpt, 
			    const unsigned long npt,
			    const unsigned long nneighbor,
			    const bool bandwidth_units,
			    unsigned long *idx, double *d2);
/* This routine is identical to kd_neighbors_point, except that it
   operates on a vector of input points instead of a single point.

   Parameters:
      INPUT kd
         the kernel density object to be searched
      INPUT idxpt
         the indices of the points whose neighbors are to be found, in
         the array x that holds the kernel density data set
      INPUT npt
         number of points whose neighbors are to be found
      INPUT nneighbor
         the number of neighbors to find
      INPUT bandwidth_units
         if true, the metric used to determine relative distance is
         normalized to the dimension-dependent bandwidth; if false,
         the metric is a simple Euclidean one
      OUTPUT idx
         indices of the neighbors found, where indices refer to the
         array x that is indexed by the kernel density object; element
         idx[i*nneighbor+j] is the jth neighbor of the ith input point;
	 must point to nneighbor * npt valid memory elements on entry
      OUTPUT d2
         squared distances of all neighbors found, indexed in the same
         way as idx; must point to nneighbor * npt valid memory
         elements on entry

   Returns
      Nothing
*/


void kd_neighbors_vec(const kernel_density *kd, const double *xpt, 
		      const unsigned long *dims, const unsigned long ndim, 
		      const unsigned long npt, const unsigned long nneighbor,
		      const bool bandwidth_units, double *pos,
		      void *dptr, double *d2);
/* This routine is identical to kd_neighbors, except that it operates
   on a vector of input points instead of a single input point.

   Parameters:
      INPUT kd
         the kernel density object to be searched
      INPUT xpt
         array of ndim * npt elements giving the position of the search
         points; xpt[ndim*i+j] is the jth coordinate of the ith point
      INPUT dims
         array specifying which dimensions in the kd tree are given in
         xpt; e.g., if the KD tree has tree dimensions, and dims = [0,
         2], then the two elements of x will be interpreted as
         specifying x and z coordinates, and the routine will search
         for the closest neighbors to a line at the specified x and
         z. If ndim is equal to the number of dimensions in the KD
         tree, this argument may be left as NULL.
      INPUT ndim
         number of elements in dims, and in each set of coordinates in
         x
      INPUT npt
         number of points whose neighbors are to be found
      INPUT nneighbor
         the number of neighbors to find
      INPUT bandwidth_units
         if true, the metric used to determine relative distance is
         normalized to the dimension-dependent bandwidth; if false,
         the metric is a simple Euclidean one
      OUTPUT pos
         positions of the nearest neighbor points; on entry, this
         pointer must point to a block of at least
         tree->ndim*nneighbor*npt elements, and on return element
	 x[i*tree->ndim*npt+j*tree->ndim+k] contains the kth
         coordinate of the jth neighbor point for the ith input point;
	 points are sorted by distance from xpt
      OUTPUT dptr
         extra data associated with each of the nearest neighbor
         points in x; element dptr[i*npt+j] is the extra data for the
         jth neighbor of the ith input point; must point to a block of
         valid memory at least nneighbor * npt elements long
      OUTPUT d2
         squared distances of all particles found from xpt; on entry,
         this pointer mut point to a block of nneighbor*npt elements,
         and on return dist2[i*nneighbor+j] gives the distance from
         the jth nearest neighbor to the ith point found to xpt

   Returns
      Nothing
*/

#endif
/* _KERNEL_DENSITY_NEIGHBORS_H_ */
