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

#include <assert.h>
#include <math.h>
#include <stdlib.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_sf_gamma.h>
#include "geometry.h"

#define MAX(X,Y) ((X) > (Y) ? (X) : (Y))
#define MIN(X,Y) ((X) < (Y) ? (X) : (Y))


bool box_in_box(const double *xbox1[2], const double *xbox2[2], 
		const unsigned long ndim1, const unsigned long ndim2,
		const unsigned long *dim1, const unsigned long *dim2) {
  unsigned long i, j;

  /* Handle cases of dimensional specifications being present or absent */
  if (dim1 == NULL) {

    /* No specification for dim1 */

    if (dim2 == NULL) {

      /* No specification for dim2 */
      for (i=0; i < MIN(ndim1, ndim2); i++) {
        if (xbox1[0][i] < xbox2[0][i]) return(false);
	if (xbox1[1][i] > xbox2[1][i]) return(false);
      }
      return(true);

    } else {

      /* Dimension specification for dim2 but not dim1 */
      for (i=0; i < MIN(ndim1, ndim2); i++) {
	if (xbox1[0][dim2[i]] < xbox2[0][i]) return(false);
	if (xbox1[1][dim2[i]] > xbox2[1][i]) return(false);
      }
      return(true);

    }
  } else {

    /* Dimension specification for dim1 */

    if (dim2 == NULL) {

      /* No specification for dim2, specification for dim1 */
      for (i=0; i < MIN(ndim1, ndim2); i++) {
        if (xbox1[0][i] < xbox2[0][dim1[i]]) return(false);
	if (xbox1[1][i] > xbox2[1][dim1[i]]) return(false);
      }
      return(true);

    } else {

      /* Both dimensions specified */

      /* Check if dim1 < dim2, meaning that box 1 is a
      higher-dimensional object than box 2. If that is the case, then
      containment is automatically false. */
      if (ndim1 < ndim2) return(false);

      /* Loop over dimensions of box 2 */
      for (j=0; j<ndim2; j++) {

	/* Find corresponding dimension in box 1, and check that it
	   fits inside the box 2 in that dimension. */
	for (i=0; i<ndim1; i++) {
	  if (dim1[i] == dim2[j]) {
	    if (xbox1[0][dim2[j]] < xbox2[0][dim1[i]]) return(false);
	    if (xbox2[1][dim2[j]] > xbox2[1][dim1[i]]) return(false);
	    break;
	  }
	}

	/* Make sure that we found an element of dim1 that matches
	   this element of dim2. If not, then box 1 is infinite in a
	   dimension where box 2 is finite, and containment fails. */
	if (i==ndim1) return(false);
      }

      /* If we made it to here, return true */
      return(true);

    }
  }

  /* Never get here */
  assert(false);
  return(false);
}



bool box_in_sphere(const double *xbox1[2], const double *xcen2, 
		   const unsigned long ndim1, const unsigned long ndim2,
		   const unsigned long *dim1, const unsigned long *dim2,
		   const double radius, const double *scale,
		   const unsigned long ndim) {
  /* Note: we make use of the fact that a cube is entirely contained
     within a sphere if and only if all its corners are inside the 
     sphere. Thus we only need to check if box_max_dist2 is less than
     our specified distance. */
  return (radius*radius >=
	  box_max_dist2(xcen2, xbox1, ndim2, ndim1, dim2, dim1, 
			scale, ndim));
}



bool box_intersects_box(const double *xbox1[2], const double *xbox2[2], 
			const unsigned long ndim1, const unsigned long ndim2,
			const unsigned long *dim1, const unsigned long *dim2) {
  unsigned long i, j;

  /* Handle cases of dimensional specifications being present or absent */
  if (dim1 == NULL) {

    /* No specification for dim1 */

    if (dim2 == NULL) {

      /* No specification for dim2 */
      for (i=0; i < MIN(ndim1, ndim2); i++) {
        if (xbox1[1][i] < xbox2[0][i]) return(false);
	if (xbox1[0][i] > xbox2[1][i]) return(false);
      }
      return(true);

    } else {

      /* Dimension specification for dim2 but not dim1 */
      for (i=0; i < MIN(ndim1, ndim2); i++) {
	if (xbox1[1][dim2[i]] < xbox2[0][i]) return(false);
	if (xbox1[0][dim2[i]] > xbox2[1][i]) return(false);
      }
      return(true);

    }
  } else {

    /* Dimension specification for dim1 */

    if (dim2 == NULL) {

      /* No specification for dim2, specification for dim1 */
      for (i=0; i < MIN(ndim1, ndim2); i++) {
        if (xbox1[1][i] < xbox2[0][dim1[i]]) return(false);
	if (xbox1[0][i] > xbox2[1][dim1[i]]) return(false);
      }
      return(true);

    } else {

      /* Both dimensions specified */

      /* Loop over dimensions of box 1 */
      for (i=0; i<ndim1; j++) {

	/* Find corresponding dimension in box 2, and check it */
	for (j=0; j<ndim2; i++) {
	  if (dim1[i] == dim2[j]) {
	    if (xbox1[1][dim2[j]] < xbox2[0][dim1[i]]) return(false);
	    if (xbox2[0][dim2[j]] > xbox2[1][dim1[i]]) return(false);
	    break;
	  }
	}
      }

      /* If we made it to here, return true */
      return(true);

    }
  }

  /* Never get here */
  assert(false);
  return(false);
}



double box_max_dist2(const double *x1, const double *xbox2[2],
		     const unsigned long ndim1, const unsigned long ndim2,
		     const unsigned long *dim1, const unsigned long *dim2,
		     const double *scale, const unsigned long ndim) {
  /* Here we make use of the fact that the most distant point is
     always a corner, so we just need to check all the corners */
  unsigned long i, j, dm1, dm2, mask, ptr;
  unsigned long ncorner=1;
  double *corner = NULL;
  double d2max = 0.0;

  /* First check if the box is semi-infinite in some dimension where
     the input object is finite; if so, the distance is automatically
     infinity. */
  if ((dim1 == NULL) && (dim2 == NULL)) {
    /* If no dimensions are specified, just compare number of
       dimensions */
    if (ndim2 < ndim1) return(INFINITY);
  } else {
    /* Dimensions are specified, so we need to check that there are no
       constraints on the point (specified by dim1) that are not
       matched by constraints on the box (specified by dim2). */
    for (i=0; i<ndim1; i++) {
      if (dim1 == NULL) dm1 = i;
      else dm1 = dim1[i];
      for (j=0; j<ndim2; j++) {
	if (dim2 == NULL) dm2 = j;
	else dm2 = dim2[j];
	if (dm1 == dm2) break;
      }
      if (j==ndim2) return(INFINITY);
    }
  }

  /* Allocate memory for swap space */
  corner = calloc(ndim, sizeof(double));

  /* Figure out how many corners our box has */
  for (i=0; i<ndim2; i++) ncorner = ncorner<<1;

  /* Loop over all corners */
  for (i=0; i<ncorner; i++) {
    mask = 1;
    for (j=0; j<ndim2; j++, mask = mask<<1) {
      ptr = ((i & mask) != 0);
      corner[j] = xbox2[ptr][j];
    }
    d2max = fmax(d2max, 
		 dist2(x1, corner, ndim1, ndim2, dim1, dim2, scale, ndim));
  }

  /* Free and return */
  free(corner);
  return(d2max);
}



double box_min_dist2(const double *x1, const double *xbox2[2],
		     const unsigned long ndim1, const unsigned long ndim2,
		     const unsigned long *dim1, const unsigned long *dim2,
		     const double *scale, const unsigned long ndim) {
  unsigned long i, j;
  double tmp, dist2 = 0.0;

  /* Handle various cases */
  if (scale == NULL) {

    /* No dimensional scaling */

    if (dim1 == NULL) {

      /* No dimensional specification for x1 */

      if (dim2 == NULL) {

	/* No dimensional specification for x2 */
	for (i=0; i < MIN(ndim1, ndim2); i++) {
	  if (x1[i] < xbox2[0][i]) {
	    tmp = x1[i] - xbox2[0][i];
	    dist2 += tmp*tmp;
	  } else if (x1[i] > xbox2[1][i]) {
	    tmp = x1[i] - xbox2[1][i];
	    dist2 += tmp*tmp;
	  }
	}

      } else {

	/* Dimensional specification for x2 but not for x1 */
	for (i=0; i < MIN(ndim1, ndim2); i++) {
	  if (x1[dim2[i]] < xbox2[0][i]) {
	    tmp = x1[dim2[i]] - xbox2[0][i];
	    dist2 += tmp*tmp;
	  } else if (x1[dim2[i]] > xbox2[1][i]) {
	    tmp = x1[dim2[i]] - xbox2[1][i];
	    dist2 += tmp*tmp;
	  }
	}
      }

    } else {

      /* Dimensional specification for x1 */

      if (dim2 == NULL) {

	/* No dimensional specification for x2 */
	for (i=0; i < MIN(ndim1, ndim2); i++) {
	  if (x1[i] < xbox2[0][dim1[i]]) {
	    tmp = x1[i] - xbox2[0][dim1[i]];
	    dist2 += tmp*tmp;
	  } else if (x1[i] > xbox2[1][dim1[i]]) {
	    tmp = x1[i] - xbox2[1][dim1[i]];
	    dist2 += tmp*tmp;
	  }
	}

      } else {

	/* Dimensional specification for both x1 and x2 */
	for (i=0; i<ndim1; i++) {
	  for (j=0; j<ndim2; j++) {
	    if (dim1[i] == dim2[j]) {
	      if (x1[dim2[j]] < xbox2[0][dim1[i]]) {
		tmp = x1[dim2[j]] - xbox2[0][dim1[i]];
		dist2 += tmp*tmp;
	      } else if (x1[dim2[j]] > xbox2[1][dim1[i]]) {
		tmp = x1[dim2[j]] - xbox2[1][dim1[i]];
		dist2 += tmp*tmp;
	      }
	      break;
	    }
	  }
	}
      }
    }

  } else {

    /* Dimensional scaling provided */

    if (dim1 == NULL) {

      /* No dimensional specification for x1 */

      if (dim2 == NULL) {

	/* No dimensional specification for x2 */
	for (i=0; i < MIN(ndim1, ndim2); i++) {
	  if (x1[i] < xbox2[0][i]) {
	    tmp = (x1[i] - xbox2[0][i]) / scale[i];
	    dist2 += tmp*tmp;
	  } else if (x1[i] > xbox2[1][i]) {
	    tmp = (x1[i] - xbox2[1][i]) / scale[i];
	    dist2 += tmp*tmp;
	  }
	}

      } else {

	/* Dimensional specification for x2 but not for x1 */
	for (i=0; i < MIN(ndim1, ndim2); i++) {
	  if (x1[dim2[i]] < xbox2[0][i]) {
	    tmp = (x1[dim2[i]] - xbox2[0][i]) / scale[dim2[i]];
	    dist2 += tmp*tmp;
	  } else if (x1[dim2[i]] > xbox2[1][i]) {
	    tmp = (x1[dim2[i]] - xbox2[1][i]) / scale[dim2[i]];
	    dist2 += tmp*tmp;
	  }
	}
      }

    } else {

      /* Dimensional specification for x1 */

      if (dim2 == NULL) {

	/* No dimensional specification for x2 */
	for (i=0; i < MIN(ndim1, ndim2); i++) {
	  if (x1[i] < xbox2[0][dim1[i]]) {
	    tmp = (x1[i] - xbox2[0][dim1[i]]) / scale[dim1[i]];
	    dist2 += tmp*tmp;
	  } else if (x1[i] > xbox2[1][dim1[i]]) {
	    tmp = (x1[i] - xbox2[1][dim1[i]]) / scale[dim1[i]];
	    dist2 += tmp*tmp;
	  }
	}

      } else {

	/* Dimensional specification for both x1 and x2 */
	for (i=0; i<ndim1; i++) {
	  for (j=0; j<ndim2; j++) {
	    if (dim1[i] == dim2[j]) {
	      if (x1[dim2[j]] < xbox2[0][dim1[i]]) {
		tmp = (x1[dim2[j]] - xbox2[0][dim1[i]]) / scale[dim1[i]];
		dist2 += tmp*tmp;
	      } else if (x1[dim2[j]] > xbox2[1][dim1[i]]) {
		tmp = (x1[dim2[j]] - xbox2[1][dim1[i]]) / scale[dim1[i]];
		dist2 += tmp*tmp;
	      }
	      break;
	    }
	  }
	}
      }
    }
  }
  return(dist2);
}



double dist2(const double *x1, const double *x2, 
	     const unsigned long ndim1, const unsigned long ndim2,
	     const unsigned long *dim1, const unsigned long *dim2, 
	     const double *scale, const unsigned long ndim) {
  unsigned long i, j;
  double tmp, dist = 0.0;

  /* Handle various cases */
  if (scale == NULL) {

    /* No dimensional scalings */

    if (dim1 == NULL) {

      /* No dimensional specification for x1 */

      if (dim2 == NULL) {

	/* No dimensional specification for x2 */
	for (i=0; i < MIN(ndim1, ndim2); i++) {
	  tmp = x1[i] - x2[i];
	  dist += tmp*tmp;
	}

      } else {

	/* Dimensional specification for x2 only */
	for (i=0; i < MIN(ndim1, ndim2); i++) {
	  tmp = x1[dim2[i]] - x2[i];
	  dist += tmp*tmp;
	}

      }

    } else {

      /* Dimensional specification for x1 */

      if (dim2 == NULL)  {

	/* No dimensional specification for x2 */
	for (i=0; i < MIN(ndim1, ndim2); i++) {
	  tmp = x1[i] - x2[dim1[i]];
	  dist += tmp*tmp;
	}

      } else {

	/* Dimensional specification for both x1 and x2 */
	for (i=0; i<ndim1; i++) {
	  for (j=0; j<ndim2; j++) {
	    if (dim1[i] == dim2[j]) {
	      tmp = x1[dim1[i]] - x2[dim2[j]];
	      dist += tmp*tmp;
	      break;
	    }
	  }
	}
      }
    }

  } else {

    /* Dimensional scalings provided */

    if (dim1 == NULL) {

      /* No dimensional specification for x1 */

      if (dim2 == NULL) {

	/* No dimensional specification for x2 */
	for (i=0; i < MIN(ndim1, ndim2); i++) {
	  tmp = (x1[i] - x2[i]) / scale[i];
	  dist += tmp*tmp;
	}

      } else {

	/* Dimensional specification for x2 only */
	for (i=0; i < MIN(ndim1, ndim2); i++) {
	  tmp = (x1[dim2[i]] - x2[i]) / scale[dim2[i]];
	  dist += tmp*tmp;
	}

      }

    } else {

      /* Dimensional specification for x1 */

      if (dim2 == NULL)  {

	/* No dimensional specification for x2 */
	for (i=0; i < MIN(ndim1, ndim2); i++) {
	  tmp = (x1[i] - x2[dim1[i]]) / scale[dim1[i]];
	  dist += tmp*tmp;
	}

      } else {

	/* Dimensional specification for both x1 and x2 */
	for (i=0; i<ndim1; i++) {
	  for (j=0; j<ndim2; j++) {
	    if (dim1[i] == dim2[j]) {
	      tmp = (x1[dim1[i]] - x2[dim2[j]]) / scale[dim1[i]];
	      dist += tmp*tmp;
	      break;
	    }
	  }
	}
      }
    }
  }

  return dist;
}    


double ds(unsigned int n) {
  /* Surface area element for an N-sphere */
  unsigned int m;
  if (n % 2) {
    /* Case for n odd */
    m = n/2;
    return 2 * n * gsl_pow_uint(2.0*M_PI, m) / gsl_sf_doublefact(2*m+1);
  } else {
    /* Case for n even */
    return n * gsl_pow_uint(M_PI, n/2) / gsl_sf_fact(n/2);
  }
}


bool sphere_in_box(const double *xcen1, const double *xbox2[2], 
		   const unsigned long ndim1, const unsigned long ndim2,
		   const unsigned long *dim1, const unsigned long *dim2,
		   const double radius, const double *scale,
		   const unsigned long ndim) {
  /* Note: we make use of the fact that a sphere x1 is entirely
     contained within cube x2 if and only if the cube in which x1 is
     inscribed is entirely contained in x2. */
  double *xbox1[2];
  unsigned long i;
  bool res;

  /* Allocate memory */
  xbox1[0] = calloc(ndim1, sizeof(double));
  xbox1[1] = calloc(ndim1, sizeof(double));

  /* Get coordinates of cube in which x1 is inscribed */
  for (i=0; i<ndim1; i++) {
    if (scale == NULL) {
      xbox1[0][i] = xcen1[i] - radius;
      xbox1[1][i] = xcen1[i] + radius;
    } else {
      if (dim1 == NULL) {
	xbox1[0][i] = xcen1[i] - radius/scale[i];
	xbox1[1][i] = xcen1[i] + radius/scale[i];
      } else {
	xbox1[0][i] = xcen1[i] - radius/scale[dim1[i]];
	xbox1[1][i] = xcen1[i] + radius/scale[dim1[i]];
      }
    }
  }

  /* Check if the inscribing cube is contained in x2 */
  res = box_in_box((const double **) xbox1, xbox2, ndim1, ndim2, dim1, dim2);

  /* Free memory */
  free(xbox1[0]);
  free(xbox1[1]);

  /* Return */
  return(res);
}
