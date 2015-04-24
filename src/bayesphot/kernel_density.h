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

void kd_pdf_grid(const kernel_density *kd, const double *xfixed,
		 const unsigned int *dimfixed, 
		 const unsigned int ndimfixed,
		 const unsigned int nfixed,
		 const double *xgrid,
		 const unsigned int *dimgrid,
		 const unsigned int ndimgrid,
		 const unsigned int ngrid,
		 const double reltol, const double abstol,
		 double *pdf);
/* This routine returns the value of the probability distribution
   function for a kernel_density object evaluated on a grid where some
   of the dimensions held fixed points and others are varying. This
   can be used, for example, to specify a fixed set of photometric
   values, and then evaluate the PDF on a grid of physical values to
   generate the marginal PDF that goes with the photometry.

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
         the PDF evaluated at x, satisfying the input error
         tolerances; element pdf[i*ngrid + j] gives the PDF for the
         jth grid point evaluated for the ith fixed point; this array
         must point to valid, allocated memory when it is passed

   Returns:
      Nothing
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

void kd_pdf_int_vec(const kernel_density *kd, const double *x, 
		    const unsigned int *dims, const unsigned int ndim,
		    const unsigned int npt, const double reltol, 
		    const double abstol, double *pdf
#ifdef DIAGNOSTIC
		    , unsigned int *nodecheck, unsigned int *leafcheck,
		    unsigned int *termcheck
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


#endif
/* _KERNEL_DENSITY_H_ */
