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
/* This module contains routines that provide the ability to draw    */
/* samples from a kernel density PDF                                 */
/*********************************************************************/

#ifndef _KERNEL_DENSITY_DRAW_H_
#define _KERNEL_DENSITY_DRAW_H_

#include "kernel_density.h"
#include "gsl/gsl_rng.h"

typedef enum draw_method_enum
  { draw_depth_first, draw_breadth_first } draw_method;

/*********************************************************************/
/* Function definitions                                              */
/*********************************************************************/

void kd_pdf_draw(const kernel_density *kd, const double *x,
		 const unsigned long *dims, const unsigned long ndim,
		 const unsigned long nsample, const draw_method method,
		 const gsl_rng *r, double *out);
/* This routine draws from a kernel density PDF, subject to a
   constaint that certain dimensions have fixed values. For example,
   for a kernel density PDF representing the joint distribution of two
   variables, this can be used to draw random samples from the second
   variable at a particular value of the first variable.

   Parameters:
      INPUT kd
         The kernel density object from which to draw
      INPUT x
         An ndim element array giving the positions in the dimensions
         that are fixed; if ndim is 0, this is ignored
      INPUT dims
         An ndim element array specifying which dimensions are fixed;
         dimensions must not be repeated; if ndim is 0, this is
         ignored
      INPUT ndim
         The number of fixed dimensions; must be < kd->tree->ndim, and
         can be 0
      INPUT nsample
         Number of samples to draw; must be > 0
      INPUT method
         Method of drawing -- depth-first or breadth first
      INPUT/OUTPUT r
         The random number generator to use
      OUTPUT out
         An array of (kd->tree->ndim - ndim) * nfixed elements
         containing the returned samples. The data are packed so that
         the first (kd->tree->ndim - ndim) elements correspond to the
         first sample, the next (kd->tree->ndim - ndim) to the next
         sample, and so forth. Each sample contains values only for
         the dimensions not contained in dims, and the remaining
         dimensions are listed in increasing order. This array must be
         allocated to the correct size before calling this routine.

   Returns:
      Nothing
*/

gsl_rng* rng_init(unsigned long seed);
/* This routine initializes a random number generator for use by the
   drawing routines.

   Parameters:
      INPUT seed
         the seed used to initialize the random number generator; if
         set to 0, one will be generated automatically by reading
         /dev/urandom if available, or using the default seed if not

   Returns:
      OUT rng
         a pointer to a random number generator
*/

void rng_free(gsl_rng *r);
/* This routine frees a random number generator initialized by
   rng_init.

   Parameters:
      INPUT r
         The random number generator to be freed.

   Returns:
      Nothing
*/

#endif
/* _KERNEL_DENSITY_DRAW_H_ */
