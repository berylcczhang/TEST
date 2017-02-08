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
/* This module contains routines to do basic operations on kernel    */
/* density estimation structures.                                    */
/*********************************************************************/

#ifndef _KERNEL_DENSITY_UTIL_H_
#define _KERNEL_DENSITY_UTIL_H_

#include <stdbool.h>
#include "kernel_density.h"

/*********************************************************************/
/* Function definitions                                              */
/*********************************************************************/

kernel_density* build_kd(double *x, unsigned long ndim, 
			 unsigned long npt, double *wgt,
			 unsigned long leafsize, double *bandwidth, 
			 kernel_type ktype, unsigned long minsplit,
			 unsigned long *sortmap);
/* This routine builds a kernel density object from a set of input
   samples and weights

   Parameters:
      INPUT/OUTPUT x
         array of npt * ndim elements containing the positions,
         ordered so that element x[j + i*ndim] is the jth coordinate
         of point i; the data are re-ordered but not copied, so
         altering x after calling build_kd will lead to bogus results
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
      INPUT minsplit
         when building the KD tree, no splittings along dimensions <
	 minsplit will be performed
      OUTPUT sortmap
         array of npt elements that, on return, will contain the
         mapping from the initial position of every input to its final
         position in the list; if this is not needed, set to NULL when
         calling

   Returns:
      OUTPUT kd
         a kernel_density object
*/

kernel_density* build_kd_sortdims(double *x, unsigned long ndim, 
				  unsigned long npt, double *wgt,
				  unsigned long leafsize, double *bandwidth, 
				  kernel_type ktype, int *nosort,
				  unsigned long *sortmap);
/* This routine builds a kernel density object from a set of input
   samples and weights. It differs from build_kd in that it allows
   users to specify that the tree is only to be sorted in certain
   dimensions. This can be useful in cases where it is known in
   advance that tree will be searched exclusively in certain
   dimensions.

   Parameters:
      INPUT/OUTPUT x
         array of npt * ndim elements containing the positions,
         ordered so that element x[j + i*ndim] is the jth coordinate
         of point i; the data are re-ordered but not copied, so
         altering x after calling build_kd will lead to bogus results
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
      INPUT nosort
         array of ndim values; for any dimension for which nosort is
         non-zero, the tree will not be partitioned along that
         dimension
      OUTPUT sortmap
         array of npt elements that, on return, will contain the
         mapping from the initial position of every input to its final
         position in the list; if this is not needed, set to NULL when
         calling

   Returns:
      OUTPUT kd
         a kernel_density object
*/

kernel_density *copy_kd(kernel_density *kd);
/* Makes a copy of a kernel density object, without copying the
   underlying KD tree. This method is useful for generating multiple
   kernel_density objects with different bandwidths but that index the
   same data.

   Parameters
      INPUT kd
         The kernel_density object to be copied

   Returns
      OUTPUT kd_copy
         a kernel_density object that is a copy of the input kd, but
         with the bandwidth, node weights, and normalizations stored
         in different memory locations
*/

void free_kd(kernel_density *kd);
/* Frees the memory associated with a kernel_density object. Note that
   this does not free any of the memory associated with the data
   indexed by the KD tree, just the metadata specific to the tree.

   Parameters
      INPUT/OUTPUT kd
         The kernel_density object to be de-allocated

   Returns
      Nothing
*/

void free_kd_copy(kernel_density *kd);
/* Frees the memory associated with a kernel_density object, but not
   the underlying data or the underlying KD tree. This differs from
   free_kd in that the KD tree is left unaffected by this routine,
   while it is deleted by free_kd. This routine is most useful in
   conjunction with copy_kd.

   Parameters
      INPUT/OUTPUT kd
         The kernel_density object to be de-allocated

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

void kd_change_wgt(const double *wgt, kernel_density *kd);
/* This routine changes the weights in a kernel_density object, while
   leaving the point positions unchanged.

   Parameters
      INPUT wgts
         Array giving the new weights; can be NULL, in which case all
         points are given equal weight
      INPUT/OUTPUT kd
         The kernel density object to be re-weighted

   Returns
      Nothing
*/

bool diagnostic_mode(void);
/* This routine just returns true if the code was compiled in
   diagnostic mode, false if it was not. */

#endif
/* _KERNEL_DENSITY_UTIL_H_ */


