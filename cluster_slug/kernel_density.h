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
  /* The bandwidth of the kernel density estimate; can be different in
     every dimension */
  double *h;
  /* The normalization factor for the kernel around each point */
  double norm;
  /* The normalization factor for the entire PDF */
  double norm_tot;
  /* The type of kernel */
  kernel_type ktype;
  /* Sums of weights in nodes of the tree */
  double *nodewgt;
} kernel_density;

/*********************************************************************/
/* Function definitions                                              */
/*********************************************************************/

kernel_density* build_kd(double *x, unsigned int ndim, 
			 unsigned int npt, double *wgt,
			 unsigned int leafsize, double *bandwidth, 
			 kernel_type ktype);
/* This routine builds a kernel density object from a set of input
   samples and weights

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

void free_kd(kernel_density *kd);
/* Frees the memory associated with a kernel_density object.

   Parameters
      INPUT/OUTPUT kd
         The kernel_density object to be de-allocated

   Returns
      Nothing
*/

void kd_change_wgt(const double *wgt, kernel_density *kd);
/* This routine changes the weights in a kernel_density object, while
   leaving the point positions unchanged.

   Parameters
      INPUT wgts
         Array giving the new weights
      INPUT/OUTPUT kd
         The kernel density object to be re-weighted

   Returns
      Nothing
*/

void kd_change_bandwidth(const double *bandwidth, kernel_density *kd);
/* This routine changes the bandwidths in a kernel_density object,
   leaving the point positions and weights unchanged.

   Parameters
      INPUT bandwidth
         Array giving the new bandwidths; must have tree->ndim
         elements
      INPUT/OUTPUT kd
         The kernel density object whose bandwidths are to be changed

   Returns
      Nothing
*/

void kd_neighbors(const kernel_density *kd, const double *xpt, 
		  const unsigned int *dims, const unsigned int ndim, 
		  const unsigned int nneighbor,
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
         point; must point to a block of valid memory at least i
         elements long
      OUTPUT d2
         squared distances of all particles found from xpt; on entry,
         this pointer mut point to a block of nneighbor elements, and
         on return dist2[i] gives the distance from the ith point
         found to xpt
*/

double kd_pdf(const kernel_density *kd, const double *x,
	      const double reltol, const double abstol
#ifdef DIAGNOSTIC
	      , unsigned int *nodecheck, unsigned int *leafcheck,
	      unsigned int *termcheck
#endif
	      );
/* This routine returns the value of the probability distribution
   function for a kernel_density object evaluated at a specified
   position, evaluated with some specified relative and absolute
   tolerances.

   Parameters:
      INPUT kd
         the kernel_density object to be used to evaluate the PDF
      INPUT x
         an ndim element array giving the position at which the PDF is
         to be evaluated
      INPUT reltol
         Relative error tolerance in the computation. An approximate
         value pdf_approx will be returned once the estimated error
	 | pdf_approx - pdf_true | / pdf_true < reltol.
      INPUT abstol
         Absolute error tolerance in the computation. An approximate
         value pdf_approx will be returned once the estimated error
	 | pdf_approx - pdf_true | < abstol.
      OUTPUT nodecheck
         Number of individual nodes examined during the evaluation;
         only if compiled with DIAGNOSTIC set
      OUTPUT leafcheck
         Number of leaves examined during the evaluation; only if
         compiled with DIAGNOSTIC set
      OUTPUT termcheck
         Number of nodes examined during the evaluation which were not
         futher sub-divided; only if compiled with DIAGNOSTIC set

   Returns:
      OUT pdf_approx
         an approximation to the PDF evaluated at x, satisfying the
         input error tolerances
*/

double kd_pdf_int(const kernel_density *kd, const double *x,
		  const unsigned int *dims, const unsigned int ndim,
		  const double reltol, const double abstol
#ifdef DIAGNOSTIC
		  , unsigned int *nodecheck, unsigned int *leafcheck,
		  unsigned int *termcheck
#endif
		  );
/* This routine returns the value of the input probability distrbution
   function evaluated at a particular point x in certain dimensions,
   with all other dimensions integrated out. For example, if the PDF
   depends on n variables, and is written out as

   p(x(0), x(1), x(2), ... x(n-1)),

   the input data point is

   x = [0.1, 0.3, 0.5],

   and the input list of dimensions is

   dims = [0, 1, 3],

   then the value returned will be

   \int p(0.1, 0.3, x(2), 0.5, x(4), x(5) ... x(n-1)) 
       dx(2) dx(4) dx(5) ... dx(n-1)

   Parameters:
      INPUT kd
         the kernel_density object to be used to evaluate the PDF
      INPUT x
         an ndim element array giving the position at which the PDF is
         to be evaluated
      INPUT dims
         an ndim element array specifying the dimensions included in x
      INPUT ndim
         number of dimensions in x; must be less than the number in
         the kd tree
      INPUT reltol
         Relative error tolerance in the computation. An approximate
         value pdf_approx will be returned once the estimated error
	 | pdf_approx - pdf_true | / pdf_true < reltol.
      INPUT abstol
         Absolute error tolerance in the computation. An approximate
         value pdf_approx will be returned once the estimated error
	 | pdf_approx - pdf_true | < abstol.
      OUTPUT nodecheck
         Number of individual nodes examined during the evaluation;
         only if compiled with DIAGNOSTIC set
      OUTPUT leafcheck
         Number of leaves examined during the evaluation; only if
         compiled with DIAGNOSTIC set
      OUTPUT termcheck
         Number of nodes examined during the evaluation which were not
         futher sub-divided; only if compiled with DIAGNOSTIC set

   Returns:
      OUT pdf_approx
         an approximation to the output integral, accurate within the
         specified error tolerances
*/

void kd_pdf_vec(const kernel_density *kd, const double *x, 
		const unsigned int npt, const double reltol, 
		const double abstol, double *pdf
#ifdef DIAGNOSTIC
		, unsigned int *nodecheck, unsigned int *leafcheck,
		unsigned int *termcheck
#endif
		);
/* This routine returns the value of the probability distribution
   function for a kernel_density object evaluated at a serires of
   specified positions, with a specified relative tolerance. The
   computation is identical to that in kd_pdf, just computed on a
   vector of points instead of a single one.

   Parameters:
      INPUT kd
         the kernel_density object to be used to evaluate the PDF
      INPUT x
         an ndim*npt element array giving the positions at which the
         PDF is to be evaluated; element x[i*ndim+j] is the jth
         coordinate of the ith input data point
      INPUT npt
         number of input positions
      INPUT reltol
         the relative tolerance for the computation; see kd_pdf_tol
      INPUT abstol
         the absolute tolerance for the computation; see kd_pdf_tol
      OUTPUT pdf
         the computed values of the PDF; array must point to npt
         elements of allocated, writeable memory on input
      OUTPUT nodecheck
         Number of individual nodes examined during each evaluation;
         must point to npt elements of allocatd, writeable memory;
         only if compiled with DIAGNOSTIC set
      OUTPUT leafcheck
         Number of leaves examined during each evaluation; must point
         to npt elements of allocatd, writeable memory; only if
         compiled with DIAGNOSTIC set
      OUTPUT termcheck
         Number of nodes examined during each evaluation which were not
         futher sub-divided; must point to npt elements of allocatd,
         writeable memory; only if compiled with DIAGNOSTIC set

   Returns:
      Nothing
*/

bool diagnostic_mode(void);
/* This routine just returns true if the code was compiled in
   diagnostic mode, false if it was not. */

#endif
/* _KERNEL_DENSITY_H_ */
