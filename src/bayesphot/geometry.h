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

#include <stdbool.h>

/*********************************************************************/
/* This module defines some basic geometric functions                */
/*********************************************************************/

#ifndef _GEOMETRY_H_
#define _GEOMETRY_H_


bool box_in_box(const double *xbox1[2], const double *xbox2[2], 
		const unsigned long ndim1, const unsigned long ndim2,
		const unsigned long *dim1, const unsigned long *dim2);
/* Returns true if box 1 is entirely contained within box 2, false
   otherwise; note that the boxes need not match up in dimensions, in
   which case one or the other is infinite in extent in certain
   directions.

   Parameters:
      INPUT xbox1
         2D array of 2 x ndim1 elements giving positions of the
         corners of box 1
      INPUT xbox2
         2D array of 2 x ndim2 elements giving the positions of
         the corners of box 2
      INPUT ndim1
         Number of coordinates in xbox1
      INPUT ndim2
         Number of coordinates in xbox2
      INPUT dim1
         Array of ndim1 elements specifying which dimensions of the
         space are specified by the entries in xbox1; if left as NULL,
         dimensions are assumed to be [0, 1, 2, ... ndim - 1],
         corresponding to xbox1 being a closed box rather than a slab
      INPUT dim2
         Array of ndim2 elements specifying which dimensions of the
         space are specified by the entries in xbox2; if left as NULL,
         dimensions are assumed to be [0, 1, 2, ... ndim - 1],
         corresponding to xbox2 being a closed box rather than a slab

   Returns:
      OUTPUT contained
         true if xbox1 is contained in xbox2, false otherwise
*/



bool box_in_sphere(const double *xbox1[2], const double *xcen2, 
		   const unsigned long ndim1, const unsigned long ndim2,
		   const unsigned long *dim1, const unsigned long *dim2,
		   const double radius, const double *scale,
		   const unsigned long ndim);
/* Returns true if xbox1 is entirely contained within a sphere centered
   at xcen2 of specified radius. Also works for the more general cases
   where xbox1 has fewer dimensions than the space (so it is a slab,
   rectangular prism, or similar semi-infinite object), xcen has fewer
   dimensions than the space (so that the "sphere" is a cylinder
   centered on a line or higher-dimensional analog thereof), or the
   dimensions are not all scaled equally in the metric (so the sphere
   is an ellipsoid).

   Parameters:
      INPUT xbox1
         2D array of 2 x ndim1 elements giving positions of the
         corners of the box
      INPUT xcen2
         array of ndim2 elements specifying the center of the sphere
      INPUT ndim1
         Number of coordinates in xbox1
      INPUT ndim2
         Number of coordinates in xbox2
      INPUT dim1
         Array of ndim1 elements specifying which dimensions of the
         space are specified by the entries in xbox1; if left as NULL,
         dimensions are assumed to be [0, 1, 2, ... ndim - 1],
         corresponding to xbox1 being a closed box rather than a slab
      INPUT dim2
         Array of ndim2 elements specifying which dimensions of the
         space are specified by the entries in xcen2; if left as NULL,
         dimensions are assumed to be [0, 1, 2, ... ndim - 1],
         corresponding to xcen2 being a point
      INPUT radius
         radius of the sphere
      INPUT scale
         Array of ndim elements giving the scale factors for the
         Euclidean metric; if left as NULL, all scale factors are set to 1
      INPUT ndim
         number of dimensions in the space being used

   Returns:
      OUTPUT contain
         true if box is contained by sphere, false otherwise
*/



bool box_intersects_box(const double *xbox1[2], const double *xbox2[2], 
			const unsigned long ndim1, const unsigned long ndim2,
			const unsigned long *dim1, const unsigned long *dim2);
/* Returns true if box 1 and box2 have non-zero intersection, false
   otherwise; note that the boxes need not match up in dimensions, in
   which case one or the other is infinite in extent in certain
   directions.

   Parameters:
      INPUT xbox1
         2D array of 2 x ndim1 elements giving positions of the
         corners of box 1
      INPUT xbox2
         2D array of 2 x ndim2 elements giving the positions of
         the corners of box 2
      INPUT ndim1
         Number of coordinates in xbox1
      INPUT ndim2
         Number of coordinates in xbox2
      INPUT dim1
         Array of ndim1 elements specifying which dimensions of the
         space are specified by the entries in xbox1; if left as NULL,
         dimensions are assumed to be [0, 1, 2, ... ndim - 1],
         corresponding to xbox1 being a closed box rather than a slab
      INPUT dim2
         Array of ndim2 elements specifying which dimensions of the
         space are specified by the entries in xbox2; if left as NULL,
         dimensions are assumed to be [0, 1, 2, ... ndim - 1],
         corresponding to xbox2 being a closed box rather than a slab

   Returns:
      OUTPUT intersect
         true if xbox1 intersects xbox2, false otherwise
*/



double box_max_dist2(const double *x1, const double *xbox2[2],
		     const unsigned long ndim1, const unsigned long ndim2,
		     const unsigned long *dim1, const unsigned long *dim2,
		     const double *scale, const unsigned long ndim);
/* Returns the squared distance between an input object (point, line,
   plane, or higher-dimensional analog) and the most distant point in
   the box or semi-infinite slab, rectangular prism, or
   higher-dimensional analog whose corners are at positions xbox,
   using a scaled Euclidean metric.

   Parameters:
      INPUT x1
         Array of ndim1 elements giving the position of the point
      INPUT xbox2
         2D array of 2 x ndim2 elements giving positions of box
         corners
      INPUT ndim1
         number of dimensions in x1
      INPUT ndim2
         number of dimensions in xbox2
      INPUT dim1
         Array of ndim1 elements specifying which dimensions of the
         space are specified by the entries in x1; if left as NULL,
         dimensions are assumed to be [0, 1, 2, ... ndim - 1],
         corresponding to x1 being a point
      INPUT dim2
         Array of ndim2 elements specifying which dimensions of the
         space are specified by the entries in x1; if left as NULL,
         dimensions are assumed to be [0, 1, 2, ... ndim - 1],
         corresponding to x1 being a point
      INPUT scale
         Array of ndim elements giving scale factors to describe
         relative sizes of dimensions; if left as NULL, scale[i] = 1
         for all dimensions
      INPUT ndim
         Number of dimensions in the space

   Returns:
      OUTPUT dist2
         the squared distance between x and the farthest point in the
         box, where distances between points are defined by
	 dist2 = sum_i (x1[i] - x2[i])^2 / scale[i]^2

   Example 1:
      Suppose this routine is called with the following arguments:
         x1 = [0.0, 1.0]
         xbox2 = [[1.0, 3.0], [2.0, 4.0], [3.0, 5.0]]
         ndim1 = 2
         ndim2 = 3
         dim1 = [0, 2]
         dim2 = NULL
         ndim = 3
         scale = NULL
      In this case, the routine is being asked to find the distance
      between the line at (x = 0, z = 1) and the most distant point in
      a box that goes from x = 1 - 3, y = 2 - 4, and z = 3 - 5. The
      most distant corder of the box is at x = 3, y = 2 or 4 (both are
      equidistant) and z = 5, and the squared distance is (3 - 0)^2 +
      (5 - 1)^2 = 25. This is what the routine will return.

   Example 2:
      Suppose this routine is called with the following arguments:
         x1 = [0.0, 1.0, 2.0]
	 xbox2 = [[1.0, 3.0], [2.0, 4.0]]
	 ndim1 = 3
	 ndim2 = 2
	 dim1 = NULL
	 dim2 = [0, 2]
	 ndim = 3
	 scale = NULL
      In this case, the routine is being asked to find the distance
      between the point (0, 1, 2) and the most distant point in the
      rectangular prism that goes from x = 1 - 3 and z = 2 - 4. This
      distance is infinite, since the prism extends to +- infinity in
      the y direction, so the routine will return infinity.
*/


double box_min_dist2(const double *x1, const double *xbox2[2],
		     const unsigned long ndim1, const unsigned long ndim2,
		     const unsigned long *dim1, const unsigned long *dim2,
		     const double *scale, const unsigned long ndim);
/* Returns the squared distance between an input object (point, line,
   plane, or higher-dimensional analog) and the closest point in the
   box (or semi-infinite slab, rectangular prism, or
   higher-dimensional analog) whose corners are at positions xbox,
   using a scaled Euclidean metric.

   Parameters:
      INPUT x1
         Array of ndim1 elements giving the position of the point
      INPUT xbox2
         2D array of 2 x ndim2 elements giving positions of box
         corners
      INPUT ndim1
         number of dimensions in x1
      INPUT ndim2
         number of dimensions in xbox2
      INPUT dim1
         Array of ndim1 elements specifying which dimensions of the
         space are specified by the entries in x1; if left as NULL,
         dimensions are assumed to be [0, 1, 2, ... ndim - 1],
         corresponding to x1 being a point
      INPUT dim2
         Array of ndim2 elements specifying which dimensions of the
         space are specified by the entries in x2; if left as NULL,
         dimensions are assumed to be [0, 1, 2, ... ndim - 1],
         corresponding to x2 being a point
      INPUT scale
         Array of ndim elements giving scale factors to describe
         relative sizes of dimensions; if left as NULL, scale[i] = 1
         for all dimensions
      INPUT ndim
         Number of dimensions in the space

   Returns:
      OUTPUT dist2
         the squared distance between x and the farthest point in the
         box, where distances between points are defined by
	 dist2 = sum_i (x1[i] - x2[i])^2 / scale[i]^2

   Example 1:
      Suppose this routine is called with the following arguments:
         x1 = [0.0, 1.0]
         xbox2 = [[1.0, 3.0], [2.0, 4.0], [3.0, 5.0]]
         ndim1 = 2
         ndim2 = 3
         dim1 = [0, 2]
         dim2 = NULL
         ndim = 3
         scale = NULL
      In this case, the routine is being asked to find the distance
      between the line at (x = 0, z = 1) and the closest point in
      a box that goes from x = 1 - 3, y = 2 - 4, and z = 3 - 5. The
      closest point is on the edge defined by x = 1, z = 3, and the
      squared distance is (1 - 0)^2 + (3 - 1)^2 = 5. This is what the
      routine will return.

   Example 2:
      Suppose this routine is called with the following arguments:
         x1 = [0.0, 1.0, 2.0]
	 xbox2 = [[1.0, 3.0], [2.0, 4.0]]
	 ndim1 = 3
	 ndim2 = 2
	 dim1 = NULL
	 dim2 = [0, 2]
	 ndim = 3
	 scale = NULL
      In this case, the routine is being asked to find the distance
      between the point (0, 1, 2) and the nearest point in the
      rectangular prism that goes from x = 1 - 3 and z = 2 - 4. This
      distance is (1-0)^2 + (2-2)^2 = 1.
*/


double dist2(const double *x1, const double *x2, 
	     const unsigned long ndim1, const unsigned long ndim2,
	     const unsigned long *dim1, const unsigned long *dim2, 
	     const double *scale, const unsigned long ndim);
/* Returns the squared distance between a pair of points, lines,
   plane, or similar higher dimensional object using the N-dimensional
   Euclidean metric, with scale factors for different directions.

   Parameters:
      INPUT x1
         Array of ndim elements giving position of the first object
         (point, line, etc.)
      INPUT x2
         Array of ndim2 elements giving position of the second object
         (point, line, etc.)
      INPUT ndim1
         Number of coordinates in x1
      INPUT ndim2
         Number of coordinates in x2
      INPUT dim1
         Array of ndim1 elements specifying which dimensions of the
         space are specified by the entries in x1; if left as NULL,
         dimensions are assumed to be [0, 1, 2, ... ndim - 1],
         corresponding to x1 being a point
      INPUT dim2
         Array of ndim2 elements specifying which dimensions of the
         space are specified by the entries in x1; if left as NULL,
         dimensions are assumed to be [0, 1, 2, ... ndim - 1],
         corresponding to x1 being a point
      INPUT scale
         Array of ndim elements giving scale factors to describe
         relative sizes of dimensions; if left as NULL, scale[i] = 1
         for all dimensions
      INPUT ndim
         Number of dimensions in the space

   Returns:
      OUTPUT dist2
         The square of the scaled distance between the input objects,
         computing using a metric where the distance between points is
	 sum_i (x1[i] - x2[i])^2 / scale[i]^2

   Example:
      Suppose this routine is called with the following arguments:
         x1 = [1.0, 2.0, 3.0]
	 x2 = [0.0, -2.0]
	 ndim1 = 3
	 dim1 = NULL
	 ndim2 = 2
	 dim2 = [0, 2]
	 ndim = 3
	 scale = NULL
      In this case we are asking for the distance between the point
      (1.0, 2.0, 3.0) and the line (x = 0.0, z = -2.0), with no
      dimensional scaling. The value returned will be 25 = (3 - (-2))^2.
*/


double ds(unsigned int n);
/* Returns the surface area element for an N-sphere.

   Parameters:
      INPUT n
         number of dimensions

   Returns:
      OUTPUT ds
         the surface area element
*/


bool sphere_in_box(const double *xcen1, const double *xbox2[2], 
		   const unsigned long ndim1, const unsigned long ndim2,
		   const unsigned long *dim1, const unsigned long *dim2,
		   const double radius, const double *scale,
		   const unsigned long ndim);
/* Returns true if a sphere centered at xcen1 of of the specified
   radius is entirely contained in xbox2. Also works for the more
   general case where xbox2 has fewer dimensions than the space (so it
   is a slab, rectangular prism, or similar semi-infinite object)
   and/or xcen1 has fewer dimensions than the space (so that the
   "sphere" is a cylinder centered on a line or higher-dimensional
   analog thereof).

   Parameters:
      INPUT xcen1
         array of ndim1 elements specifying the center of the sphere
      INPUT xbox2
         2D array of 2 x ndim2 elements giving positions of the
         corners of the box
      INPUT ndim1
         Number of coordinates in xbox1
      INPUT ndim2
         Number of coordinates in xbox2
      INPUT dim1
         Array of ndim1 elements specifying which dimensions of the
         space are specified by the entries in xbox1; if left as NULL,
         dimensions are assumed to be [0, 1, 2, ... ndim - 1],
         corresponding to xbox1 being a closed box rather than a slab
      INPUT dim2
         Array of ndim2 elements specifying which dimensions of the
         space are specified by the entries in xcen2; if left as NULL,
         dimensions are assumed to be [0, 1, 2, ... ndim - 1],
         corresponding to xcen2 being a point
      INPUT radius
         radius of the sphere
      INPUT scale
         Array of ndim elements giving the scale factors for the
         Euclidean metric; if left as NULL, all scale factors are set to 1
      INPUT ndim
         number of dimensions in the space being used

   Returns:
      OUTPUT contain
         true if sphere is contained by box, false otherwise
*/

#endif
/* _GEOMETRY_H_ */
