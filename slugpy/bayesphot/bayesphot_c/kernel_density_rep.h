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

/*********************************************************************/
/* This module contains routines that are used in construct          */
/* approximate representations of KD tree objects.                   */
/*********************************************************************/

#ifndef _KERNEL_DENSITY_REP_H_
#define _KERNEL_DENSITY_REP_H_

#include "kernel_density.h"

/*********************************************************************/
/* Function definitions                                              */
/*********************************************************************/

unsigned long kd_rep(const kernel_density *kd, const double *x,
		     const unsigned long *dims, 
		     const unsigned long ndim, const double reltol,
		     const unsigned long *dim_return,
		     const unsigned long ndim_return,
		     double **xpt, double **wgts);
/* This routine returns a list of points that can be used to do fast
   kernel density estimation for a given input point. Performing the
   estimate using the list of points returned is guaranteed to
   recover, when integrated over all space, a fraction 1.0-reltol of
   the probability that one would recover by using all the data
   available. This only works for Gaussian kernels, because only
   Gaussian kernels have the property that the weights resulting from
   distances in different dimensions can simply be multiplied by each
   other to get a total weight.

   Parameters:
      INPUT kd
         the kernel_density object to be used to evaluate the PDF
      INPUT x
         an ndim element array giving the position for which the
         approximate representation is to be constructed
      INPUT dims
         an ndim element array specifying which dimensions the elements
	 of x correspond to
      INPUT ndims
         the number of dimensions in x and dims
      INPUT reltol
         Tolerance on the integral of the returned approximation; the
         sample points returned are guaranteed to account for >
	 1-inttol of the integrated probability, i.e., if reltol =
         0.01 then if one integrates the PDF that results from
         computing a kernel density estimate just on the returned
         points, that will recover 99% of the total probability
      INPUT dim_return
         array of ndim_return elements specifying the dimensions to
         return in xpt; no element can appear in both this and dims;
         if left as NULL, all dimensions except those in dims will be
         returned
      INPUT ndim_return
         number of dimensions in dim_return; if dim_return is NULL,
         this argument is ignored, and the number of dimensions
         returned will automatically be kd->tree->ndim - ndim
      OUTPUT xpt
         Array of positions for the approximate representation, in the
         dimensions in ndim_return, and the array is ordered such that
         xpt[i*ndim_return + j] is the ith coordinate of the jth
         point; the array has npt * ndim_return, elements, where npt
         is the value returned by the function
      OUTPUT wgts
         Array of npt elements containing the weights of the returned
         points, adjusted based on their distance to the input point;
	 weights are normalized to have a sum of unity

   Returns:
      npt, the number of points in xpt and wgts. A return value of 0
      indicates that the routine terminated prematurely due to
      inadequate numerical precision to achieve the desired
      tolerance. This usually indicates that the input point x was too
      far from the library for the exp(-dist) kernel functions to be
      evaluated with adequate precision.

   Notes:
      The approximate PDF can be computed from the representation as
      p(x[j]) = sum_i wgts[i] * prod_j exp(-(x[j] - xpt[i,j]) / 
      (2*bandwidth[j]**2)),
      where the sum over i runs over all the returned points, the
      product over j runs over all the returned dimensions, and x[j]
      is the point at which the PDF is to be computed.
*/

void free_kd_rep(double **xpt, double **wgts);
/* This routine frees the xpt and wgts arrays returned by kd_rep. It
   is provided as a convenience for python, so that it can free memory
   that is returned by kd_rep.

   Parameters:
      INPUT xpt
         The xpt array to be freed
      INPUT wgts
         The wgts array to be freed

   Returns:
      Nothing
*/

unsigned long squeeze_rep(const unsigned long npts, 
			  const unsigned int ndim, double *h, 
			  const double tol, double **x, 
			  double **wgts);
/* This routine shrinks the size of a kernel density representation by
   approximating Gaussians that are close to one another with a single
   Gaussian.

   Parameters:
      INPUT npt
         Number of points in the representation to be squeezed
      INPUT ndim
         Number of dimensions in the representation
      INPUT h
         Array of ndim elements giving the bandwidth in each dimension
      INPUT tol
         Error tolerance of the squeezing procedure; the fractional
         difference between the original and squeezed representations
         is guaranteed to be < tol
      INPUT/OUTPUT x
         Array of npt * ndim elements giving the positions of the
         points to be squeezed; on exit, the array will have been
         resized to contain npt_out * ndim elements
      INPUT/OUTPUT wgts
         Array of npt elements giving the weights of the points to be
         squeezed; on exit, the array will have been resized to
         contain npt_out elements

   Returns:
      npt_out, the number of points in the squeezed
      representation. The squeezed representation may be used to
      reconstruct the PDF in exactly the same manner as the original
      one (see kd_rep).
*/

#endif
/* _KERNEL_DENSITY_REP_H_ */
