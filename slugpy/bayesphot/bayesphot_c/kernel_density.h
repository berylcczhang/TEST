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

double kd_pdf(const kernel_density *kd, const double *x,
	      const double reltol, const double abstol
#ifdef DIAGNOSTIC
	      , unsigned long *nodecheck, unsigned long *leafcheck,
	      unsigned long *termcheck
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

void kd_pdf_grid(const kernel_density *kd, const double *xfixed,
		 const unsigned long *dimfixed, 
		 const unsigned long ndimfixed,
		 const unsigned long nfixed,
		 const double *xgrid,
		 const unsigned long *dimgrid,
		 const unsigned long ndimgrid,
		 const unsigned long ngrid,
		 const double reltol, const double abstol,
		 double *pdf);
/* This routine returns the value of the probability distribution
   function for a kernel_density object evaluated on a grid where some
   of the dimensions held at fixed points and others are varying. This
   can be used, for example, to specify a fixed set of photometric
   values, and then evaluate the PDF on a grid of physical values at
   that fixed value of the photometry.

   Parameters:
      INPUT kd
         the kernel_density object to be used to evaluate the PDF
      INPUT xfixed
         an ndimfixed * nfixed element array; each block of ndimfixed
         elements gives a position in the dimensions specified by
         dimfixed, and there are nfixed such blocks.
      INPUT dimfixed
         an ndimfixed element array specifying the dimensions in each
         of the nfixed blocks in xfixed; dimensions must not be
         repeated or appear in both dimfixed and dimgrid, and every
         dimension in the kernel density object must appear in either
         dimfixed or dimgrid
      INPUT nfixed
         number of fixed points
      INPUT xgrid
         an ndimgrid * ngrid element array; each block of ndimgrid
         elements gives the position of a grid point, and there are
         ngrid such blocks.
      INPUT dimgrid
         an ndimgrid element array specifying the dimensions in each
         of the ngrid blocks in xgrid; dimensions must not be
         repeated or appear in both dimfixed and dimgrid, and every
         dimension in the kernel density object must appear in either
         dimfixed or dimgrid
      INPUT ngrid
         number of grid points
      INPUT reltol
         Relative error tolerance in the computation. An approximate
         value pdf_approx will be returned once the estimated error
	 | pdf_approx - pdf_true | / pdf_true < reltol.
      INPUT abstol
         Absolute error tolerance in the computation. An approximate
         value pdf_approx will be returned once the estimated error
	 | pdf_approx - pdf_true | < abstol.
      OUT pdf
         an ngrid * nfixed element array giving an approximation to
         the PDF, satisfying the input error tolerances; element
         pdf[i*ngrid + j] gives the PDF for the jth grid point
         evaluated for the ith fixed point; this array must point to
         valid, allocated memory when it is passed

   Returns:
      Nothing
*/


double kd_pdf_int(const kernel_density *kd, const double *x,
		  const unsigned long *dims, const unsigned long ndim,
		  const double reltol, const double abstol
#ifdef DIAGNOSTIC
		  , unsigned long *nodecheck, unsigned long *leafcheck,
		  unsigned long *termcheck
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


void kd_pdf_int_grid(const kernel_density *kd, const double *xfixed,
		     const unsigned long *dimfixed, 
		     const unsigned long ndimfixed,
		     const unsigned long nfixed,
		     const double *xgrid,
		     const unsigned long *dimgrid,
		     const unsigned long ndimgrid,
		     const unsigned long ngrid,
		     const double reltol, const double abstol,
		     double *pdf);
/* This routine returns the value of the probability distribution
   function for a kernel_density object evaluated on a grid where some
   of the dimensions held at fixed points, others are varying, and
   still others are integrated out. This can be used, for example, to
   specify a fixed set of photometric values, and then evaluate the
   PDF on a grid of physical values, marginalizing over other physical
   values.

   Parameters:
      INPUT kd
         the kernel_density object to be used to evaluate the PDF
      INPUT xfixed
         an ndimfixed * nfixed element array; each block of ndimfixed
         elements gives a position in the dimensions specified by
         dimfixed, and there are nfixed such blocks.
      INPUT dimfixed
         an ndimfixed element array specifying the dimensions in each
         of the nfixed blocks in xfixed; dimensions must not be
         repeated or appear in both dimfixed and dimgrid; dimensions
         that are not included in either dimfixed or dimgrid will be
         integrated over (i.e., marginalized out)
      INPUT nfixed
         number of fixed points
      INPUT xgrid
         an ndimgrid * ngrid element array; each block of ndimgrid
         elements gives the position of a grid point, and there are
         ngrid such blocks.
      INPUT dimgrid
         an ndimgrid element array specifying the dimensions in each
         of the ngrid blocks in xgrid; dimensions must not be
         repeated or appear in both dimfixed and dimgrid;  dimensions
         that are not included in either dimfixed or dimgrid will be
         integrated over (i.e., marginalized out)
      INPUT ngrid
         number of grid points
      INPUT reltol
         Relative error tolerance in the computation. An approximate
         value pdf_approx will be returned once the estimated error
	 | pdf_approx - pdf_true | / pdf_true < reltol.
      INPUT abstol
         Absolute error tolerance in the computation. An approximate
         value pdf_approx will be returned once the estimated error
	 | pdf_approx - pdf_true | < abstol.
      OUT pdf
         an ngrid * nfixed element array giving an approximation to
         the PDF evaluated at x, satisfying the input error
         tolerances; element pdf[i*ngrid + j] gives the PDF for the
         jth grid point evaluated for the ith fixed point; this array
         must point to valid, allocated memory when it is passed

   Returns:
      Nothing
*/

void kd_pdf_int_reggrid(const kernel_density *kd, const double *xfixed,
			const unsigned long *dimfixed, 
			const unsigned long ndimfixed,
			const unsigned long nfixed,
			const double *xgridlo,
			const double *xgridhi,
			const unsigned long *ngrid,
			const unsigned long *dimgrid,
			const unsigned long ndimgrid,
			const double reltol, const double abstol,
			double *pdf);
/* This routine returns the value of the probability distribution
   function for a kernel_density object evaluated on a grid where some
   of the dimensions held at fixed points, others are varying, and
   still others are integrated out. This can be used, for example, to
   specify a fixed set of photometric values, and then evaluate the
   PDF on a grid of physical values, marginalizing over other physical
   values. This routine is identical to kd_pdf_int_grid, except that
   the grid is required to be regular, which allows for a significant
   computational speedup.

   Parameters:
      INPUT kd
         the kernel_density object to be used to evaluate the PDF
      INPUT xfixed
         an ndimfixed * nfixed element array; each block of ndimfixed
         elements gives a position in the dimensions specified by
         dimfixed, and there are nfixed such blocks.
      INPUT dimfixed
         an ndimfixed element array specifying the dimensions in each
         of the nfixed blocks in xfixed; dimensions must not be
         repeated or appear in both dimfixed and dimgrid; dimensions
         that are not included in either dimfixed or dimgrid will be
         integrated over (i.e., marginalized out)
      INPUT ndimfixed
         number of dimensions that are fixed; must be > 0
      INPUT nfixed
         number of fixed points; must be > 0
      INPUT xgridlo
         an ndimgrid element array, giving the coordinates of the
         lower left corner of the grid on which the PDF is to be
         evaluated
      INPUT xgridhi
         an ndimgrid element array, giving the coordinates of the
         upper right corner of the grid on which the PDF is to be
         evaluated
      INPUT ngrid
         an ndimgrid element array specifying the number of grid
         points in each dimension of the grid; must be > 0
      INPUT dimgrid
         an ndimgrid element array specifying the dimensions in each
         of the ngrid blocks in xgrid; dimensions must not be
         repeated or appear in both dimfixed and dimgrid, and every
         dimension in the kernel density object must appear in either
         dimfixed or dimgrid
      INPUT reltol
         Relative error tolerance in the computation. An approximate
         value pdf_approx will be returned once the estimated error
	 | pdf_approx - pdf_true | / max(pdf_true) < reltol.
      INPUT abstol
         Absolute error tolerance in the computation. An approximate
         value pdf_approx will be returned once the estimated error
	 | pdf_approx - pdf_true | < abstol.
      OUT pdf
         an ngrid[0] * ngrid[1] * ... ngrid[ndimgrid-1] * nfixed
         element array giving an approximation to the PDF, satisfying
         the input error tolerances; this array must point to valid,
         allocated memory when it is passed; data are packed in
         standard c order, with nfixed as the slowest varying
         dimension. For example, if ngrid = [5, 4, 2] and nfixed
         = 3, then the PDF for element (4, 3, 1) of the grid,
         evaluated for fixed point 2, is located at
	 pdf[1 + ngrid[2]*(3+ngrid[1]*(4+ngrid[0]*nfixed))]

   Returns:
      Nothing

   Notes:
      Grids points for PDF evaluations are placed so that the jth grid
      point along dimension dimgrid[i] will be placed at position
      xgridlo[i] + (xgridhi[i] - xgridlo[i]) / (ngrid[i]-1) * j
      for j = 0 ... ngrid[i]-1. Special case: if ngrid[i] = 1,
      then a single grid point will be placed at xgridlo[i], and
      xgridhi[i] is ignored.
*/


void kd_pdf_int_vec(const kernel_density *kd,
		    const double *x, 
		    const double *bandwidth,
		    const unsigned long *dims,
		    const unsigned long ndim,
		    const unsigned long npt,
		    const double reltol, 
		    const double abstol,
		    double *pdf
#ifdef DIAGNOSTIC
		    , unsigned long *nodecheck,
		    unsigned long *leafcheck,
		    unsigned long *termcheck
#endif
		    );
/* This routine is identical to kd_pdf_int, except that it operates
   on a vector of input points x instead of a single input point.

   Parameters:
      INPUT kd
         the kernel_density object to be used to evaluate the PDF
      INPUT x
         an ndim*npt element array giving the position at which the PDF is
         to be evaluated; element x[i*ndim+j] is the jth element of
         the ith point
      INPUT bandwidth
         if this is not NULL, it must be an ndim*npt element array
         giving the bandwidths to be used in evaluating the PDF for
         each input point; bandwidth[i*ndim+j] is the bandwidth
         that will be used for the jth dimension and the ith point;
         if this is NULL, then the evaluation will use the bandwidth
         stored in the kernel_density object for all points
      INPUT dims
         an ndim element array specifying the dimensions included in x
      INPUT ndim
         number of dimensions in x; must be less than the number in
         the kd tree
      INPUT npt
         number of input positions
      INPUT reltol
         Relative error tolerance in the computation. An approximate
         value pdf_approx will be returned once the estimated error
	 | pdf_approx - pdf_true | / pdf_true < reltol.
      INPUT abstol
         Absolute error tolerance in the computation. An approximate
         value pdf_approx will be returned once the estimated error
	 | pdf_approx - pdf_true | < abstol.
      OUT pdf
         an approximation to the output integral, accurate within the
         specified error tolerances
      OUTPUT nodecheck
         Number of individual nodes examined during the evaluation;
         only if compiled with DIAGNOSTIC set; must point to npt
         elements of valid, writeable memory
      OUTPUT leafcheck
         Number of leaves examined during the evaluation; only if
         compiled with DIAGNOSTIC set; must point to npt elements of
         valid, writeable memory
      OUTPUT termcheck
         Number of nodes examined during the evaluation which were not
         futher sub-divided; only if compiled with DIAGNOSTIC set;
         must point to npt elements of valid, writeable memory

   Returns:
      Nothing

   Notes:
      Setting the argument bandwidth to a value other than NULL has
      the same effect as calling kd_change_bandwidth followed by kd_pdf_int
      npt times, then calling change bandwidth to set the bandwidth
      back to its original value. However, this routine executes this
      procedure in a manner that is thread-safe.
*/


void kd_pdf_reggrid(const kernel_density *kd, const double *xfixed,
		    const unsigned long *dimfixed, 
		    const unsigned long ndimfixed,
		    const unsigned long nfixed,
		    const double *xgridlo,
		    const double *xgridhi,
		    const unsigned long *ngrid,
		    const unsigned long *dimgrid,
		    const unsigned long ndimgrid,
		    const double reltol, const double abstol,
		    double *pdf);
/* This routine returns the value of the probability distribution
   function for a kernel_density object evaluated on a regular grid
   where some of the dimensions held at fixed points and others are
   varying. This can be used, for example, to specify a fixed set of
   photometric values, and then evaluate the PDF on a grid of physical
   values at that fixed value of the photometry. This routine differs
   from kd_pdf_grid in that the grid for this routine is required to
   be regularly spaced in the varying dimensions. This requirement
   allows a significant computational speedup.

   Parameters:
      INPUT kd
         the kernel_density object to be used to evaluate the PDF
      INPUT xfixed
         an ndimfixed * nfixed element array; each block of ndimfixed
         elements gives a position in the dimensions specified by
         dimfixed, and there are nfixed such blocks.
      INPUT dimfixed
         an ndimfixed element array specifying the dimensions in each
         of the nfixed blocks in xfixed; dimensions must not be
         repeated or appear in both dimfixed and dimgrid, and every
         dimension in the kernel density object must appear in either
         dimfixed or dimgrid
      INPUT nfixed
         number of fixed points
      INPUT xgridlo
         an ndimgrid element array, giving the coordinates of the
         lower left corner of the grid on which the PDF is to be
         evaluated
      INPUT xgridlo
         an ndimgrid element array, giving the coordinates of the
         upper right corner of the grid on which the PDF is to be
         evaluated
      INPUT ngrid
         an ndimgrid element array specifying the number of grid
         points in each dimension of the grid; must be > 0
      INPUT dimgrid
         an ndimgrid element array specifying the dimensions in each
         of the ngrid blocks in xgrid; dimensions must not be
         repeated or appear in both dimfixed and dimgrid, and every
         dimension in the kernel density object must appear in either
         dimfixed or dimgrid
      INPUT reltol
         Relative error tolerance in the computation. The returned
         value pdf[i] for the PDF at gridpoint x[i] is considered
         converged if | pdf[i] - pdf_true(x[i]) | / max(pdf[i]) <
         reltol, where pdf_true(x[i]) is the true PDF evaluated at
         x[i]. Note that this means that relative tolerances are
         expressed as error relative to the maximum of the PDF, not
         with respect to each point.
      INPUT abstol
         Absolute error tolerance in the computation. The return value
         pdf[i] for the PDF at gridpoint x[i] is considered converged
         if | pdf[i] - pdf_true(x[i]) | < abstol, where pdf_true(x[i])
         is the true PDF evaluated at x[i].
      OUT pdf
         an ngrid[0] * ngrid[1] * ... ngrid[ndimgrid-1] * nfixed
         element array giving an approximation to the PDF, satisfying
         the input error tolerances; this array must point to valid,
         allocated memory when it is passed; data are packed in
         standard c order, with nfixed as the slowest varying
         dimension. For example, if ngrid = [5, 4, 2] and nfixed
         = 3, then the PDF for element (4, 3, 1) of the grid,
         evaluated for fixed point 2, is located at
	 pdf[1 + ngrid[2]*(3+ngrid[1]*(4+ngrid[0]*nfixed))]

   Returns:
      Nothing

   Notes:
      Grids points for PDF evaluationa are placed so that the jth grid
      point along dimension dimgrid[i] will be placed at position
      xgridlo[i] + (xgridhi[i] - xgridlo[i]) / (ngrid[i]-1) * j
      for j = 0 ... ngrid[i]-1. Special case: if ngrid[i] = 1,
      then a single grid point will be placed at xgridlo[i], and
      xgridhi[i] is ignored.
*/

void kd_pdf_vec(const kernel_density *kd,
		const double *x,
		const double *bandwidth,
		const unsigned long npt,
		const double reltol,
		const double abstol,
		double *pdf
#ifdef DIAGNOSTIC
		, unsigned long *nodecheck,
		unsigned long *leafcheck,
		unsigned long *termcheck
#endif
		);
/* This routine returns the value of the probability distribution
   function for a kernel_density object evaluated at a series of
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
      INPUT bandwidth
         if this is not NULL, it must be an ndim*npt element array
         giving the bandwidths to be used in evaluating the PDF for
         each input point; bandwidth[i*ndim+j] is the bandwidth
         that will be used for the jth dimension and the ith point;
         if this is NULL, then the evaluation will use the bandwidth
         stored in the kernel_density object for all points
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

   Notes:
      Setting the argument bandwidth to a value other than NULL has
      the same effect as calling kd_change_bandwidth followed by kd_pdf
      npt times, then calling change bandwidth to set the bandwidth
      back to its original value. However, this routine executes this
      procedure in a manner that is thread-safe.
*/


#endif
/* _KERNEL_DENSITY_H_ */
