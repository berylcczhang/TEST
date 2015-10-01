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

#include "kernel_density_neighbors.h"

/*********************************************************************/
/* Finds nearest neighbors to a specified input point                */
/*********************************************************************/
void kd_neighbors(const kernel_density *kd, const double *xpt, 
		  const unsigned long *dims, const unsigned long ndim, 
		  const unsigned long nneighbor,
		  const bool bandwidth_units, double *pos,
		  void *dptr, double *d2) {
  /* Call KDtree neighbor finding routine */
  if (bandwidth_units) {
    neighbors(kd->tree, xpt, dims, ndim, nneighbor, kd->h, pos, dptr, d2);
  } else {
    neighbors(kd->tree, xpt, dims, ndim, nneighbor, NULL, pos, dptr, d2);
  }
}

/*********************************************************************/
/* Find the N nearest neighbors to all points in the KDtree          */
/*********************************************************************/
void kd_neighbors_all(const kernel_density *kd, 
		      const unsigned long nneighbor, 
		      const bool bandwidth_units, unsigned long *idx, 
		      double *d2) {
  /* Call the KDtree neighbor finding routine */
  if (bandwidth_units) {
    neighbors_all(kd->tree, nneighbor, kd->h, idx, d2);
  } else {
    neighbors_all(kd->tree, nneighbor, NULL, idx, d2);
  }
}

/*********************************************************************/
/* Find N nearest neighbors to a single point in the KD tree         */
/*********************************************************************/
void kd_neighbors_point(const kernel_density *kd, 
			const unsigned long idxpt, 
			const unsigned long nneighbor,
			const bool bandwidth_units,
			unsigned long *idx, double *d2) {
  /* Just call the KDtree routine */
  if (bandwidth_units) {
    neighbors_point(kd->tree, idxpt, nneighbor, kd->h, idx, d2);
  } else {
    neighbors_point(kd->tree, idxpt, nneighbor, NULL, idx, d2);
  }
}

/*********************************************************************/
/* Vectorized version of kd_neighbors_point                          */
/*********************************************************************/
void kd_neighbors_point_vec(const kernel_density *kd, 
			    const unsigned long *idxpt, 
			    const unsigned long npt,
			    const unsigned long nneighbor,
			    const bool bandwidth_units,
			    unsigned long *idx, double *d2) {
  unsigned long i;
  for (i=0; i<npt; i++) {
    if (bandwidth_units) {
      neighbors_point(kd->tree, idxpt[i], nneighbor, kd->h, 
		      idx+i*nneighbor, d2+i*nneighbor);
    } else {
      neighbors_point(kd->tree, idxpt[i], nneighbor, NULL, 
		      idx+i*nneighbor, d2+i*nneighbor);
    }
  }
}

/*********************************************************************/
/* Same as kd_neighbors, but for a vector of input points            */
/*********************************************************************/
void kd_neighbors_vec(const kernel_density *kd, const double *xpt, 
		      const unsigned long *dims, const unsigned long ndim, 
		      const unsigned long npt, const unsigned long nneighbor,
		      const bool bandwidth_units, double *pos,
		      void *dptr, double *d2) {
  unsigned long i;

  /* Loop over input points, calling KDtree neighbor find routine on each */
  for (i=0; i<npt; i++) {
    if (bandwidth_units) {
      neighbors(kd->tree, 
		xpt + i*ndim, 
		dims, ndim, nneighbor, kd->h, 
		pos + i*kd->tree->ndim*nneighbor, 
		(void *) (((double *) dptr) + i*nneighbor),
		d2 + i*nneighbor);
    } else {
      neighbors(kd->tree, 
		xpt + i*ndim, 
		dims, ndim, nneighbor, NULL, 
		pos + i*kd->tree->ndim*nneighbor, 
		(void *) (((double *) dptr) + i*nneighbor),
		d2 + i*nneighbor);
    }
  }
}
