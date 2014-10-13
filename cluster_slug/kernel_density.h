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
/* This module contains routines to compute kernel density estimates */
/* from KD trees                                                     */
/*********************************************************************/

#ifndef _KERNEL_DENSITY_H_
#define _KERNEL_DENSITY_H_

#include <stdbool.h>
#include "kdtree.h"

/*********************************************************************/
/* Types of kernels                                                  */
/*********************************************************************/
typedef enum kernelType { epanechnikov, tophat, gaussian } 
  kernel_type;

/*********************************************************************/
/* Kernel density structure                                          */
/*********************************************************************/
typedef struct {
  /* The KD tree describing the data */
  KDtree *tree;
  /* The bandwidth of the kernel density estimate */
  double h;
  /* Scale factors for relative sizes of different dimensions */
  double *scalefac;
  /* The normalization factor for the kernel around each point */
  double norm;
  /* The normalization factor for the entire PDF */
  double norm_tot;
  /* The type of kernel */
  kernel_type ktype;
  /* Do we have our own copy of the data positions? */
  bool copy_pos;
} kernel_density;


/*********************************************************************/
/* Function definitions                                              */
/*********************************************************************/

kernel_density* build_kd(double *x, unsigned int ndim, 
			 unsigned int npt, double *wgt, 
			 unsigned int leafsize, double bandwidth, 
			 kernel_type ktype);
/* This routine builds a kernel density object from a set of input
   samples and weights.

   Parameters:
      INPUT/OUTPUT x
         array of npt * ndim elements containing the positions,
         ordered so that element x[j + i*ndim] is the jth coordinate
         of point i; on return, this will be sorted into a tree; data
         are not copied, so altering x will result in bogus results
      INPUT ndim
         number of dimensions in the data set
      INPUT npt
         number of data points
      INPUT/OUTPUT wgt
         array of npt elements giving weights of all points; set to
         NULL for equal weighting; data are not copied, so altering
         wgt will lead to bogus results
      INPUT leafsize
         number of points to place in a leaf of the tree
      INPUT bandwidth
         bandwidth of the kernel
      INPUT ktype
         functional form of the kernel

   Returns:
      OUTPUT kd
         a kernel_density object
*/

kernel_density* build_kd_bvec(double *x, unsigned int ndim, 
			      unsigned int npt, double *wgt,
			      unsigned int leafsize, double *bandwidth, 
			      kernel_type ktype);
/* This routine builds a kernel density object from a set of input
   samples and weights; each dimension has its own bandwidth.

   Parameters:
      INPUT x
         array of npt * ndim elements containing the positions,
         ordered so that element x[j + i*ndim] is the jth coordinate
         of point i; unlike with build_kd, these data ARE copied, so
         the array x may be altered or disposed of once this function
         returns
      INPUT ndim
         number of dimensions in the data set
      INPUT npt
         number of data points
      INPUT/OUTPUT wgt
         array of npt elements giving weights of all points; set to
         NULL for equal weighting; data are not copied, so altering
         wgt will lead to bogus results
      INPUT leafsize
         number of points to place in a leaf of the tree
      INPUT bandwidth
         array of ndim elements giving the bandwidth of the kernel in
         each dimension
      INPUT ktype
         functional form of the kernel

   Returns:
      OUTPUT kd
         a kernel_density object
*/

void find_neighbors(const double *xpt, const unsigned int *dims,
		    const unsigned int ndim, 
		    const unsigned int nneighbor, 
		    kernel_density *kd, double *pos,
		    double *wgt, double *dist2);
/* Routine to find the N nearest neighbors to an input point; input
   points can have fewer dimensions than the search space, in which
   case the routine searches for the nearest neighbors to a line,
   plane, or higher-dimensional object.

   Parameters:
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
      INPUT kd
         the kernel_density object to be searched
      OUTPUT x
         positions of the nearest neighbor points; on entry, this
         pointer must point to a block of at least
         tree->ndim*nneighbor elements, and on return element
	 x[i*tree->ndim+j] contains the jth coordinate for the ith
         neighbor found; points are sorted by distance from xpt
      OUTPUT wgt
         weights of the points found, only set if the kernel_density
         object uses weights; on exit, wgt[i] is the weight of the ith
         point point; if weights are set, this must point to a block
         of memory at least nneighbor elements long
      OUTPUT dist2
         squared distances of all particles found from xpt; on entry,
         this pointer mut point to a block of nneighbor elements, and
         on return dist2[i] gives the distance from the ith point
         found to xpt
*/

void free_kd(kernel_density **kd);
/* Frees the memory associated with a kernel_density object.

   Parameters
      INPUT/OUTPUT kd
         The kernel_density object to be de-allocated

   Returns
      Nothing
*/

double kd_pdf(const double *x, const kernel_density *kd);
/* This routine returns the value of the probability distribution
   function for a kernel_density object evaluated at a specified
   position.

   Parameters:
      INPUT x
         an ndim element array giving the position at which the PDF is
         to be evaluated
      INPUT kd
         the kernel_density object to be used to evaluate the PDF

   Returns:
      OUT pdf
         the PDF evaluated at x
*/


void kd_pdf_vec(const double *x, const unsigned int npt, 
		const kernel_density *kd, double *pdf);
/* This routine returns the value of the probability distribution
   function for a kernel_density object evaluated at a serires of
   specified positions. The computation is identical to that in
   kd_pdf, just computed on a vector of points instead of a single
   one.

   Parameters:
      INPUT x
         an ndim*npt element array giving the positions at which the
         PDF is to be evaluated; element x[i*ndim+j] is the jth
         coordinate of the ith input data point
      INPUT npt
         number of input positions
      INPUT kd
         the kernel_density object to be used to evaluate the PDF
      OUTPUT pdf
         the computed values of the PDF; array must point to npt
         elements of allocated, writeable memory on input

   Returns:
      Nothing
*/


double kd_pdf_int(const double *q, const unsigned int *qdim, 
		  unsigned int nqdim, const kernel_density *kd);
/* This function returns value of the probability distribution
   function for a kernel_density object evaluated at a specified
   position over certain dimensions, and integrated over the other
   dimensions. That is, if p(x_0, x_1, ... x_{N-1}) is the PDF, this
   function returns the quantity:
   f(q_{a_0}, q_{a_1}, ... q_{x_{M-1}}) =
   \int p(x_0, x_1, ... x_{N-1}) * delta(x_{a_0}-q_{a_0}) *
       delta(x_{a_1}-q_{a_1}) * ... * delta(x_{a_{M-1}}-q_{a_{M-1}}) dV,
   where a_0 ... a_M are any set of indices in the range 0 ... N-1.

   Parameters
      INPUT q
         An array giving the coordinates q_{a_0}, q_{a,1} ... q_{a_M-1}
      INPUT qdim
         An array giving the dimensions a_0, a_1, ... a_{M-1} for the
         coordinates.
      INPUT nqdim
         Number of elements in q and qdim (in the notation above, M);
         must be in the range 1 ... N-1.
      INPUT kd
         The kernel density object over which the projection is to be
         carried out.

   Returns
      OUTPUT f
        The value of the function f(q_{a_0}, q_{a_1}, ... q_{x_{M-1}})
	as defined above.
*/

void reweight_kd(const double *wgt, kernel_density *kd);
/* This routine changes the weights in a kernel_density object, while
   leaving the point positions unchanged.

   Parameters
      INPUT wgts
         Any array giving the new weights
      INPUT kd
         The kernel density object to be re-weighted

   Returns
      Nothing
*/

#endif
/* _KERNEL_DENSITY_H_ */
