/*********************************************************************
Copyright (C) 2017 Mark Krumholz
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

////////////////////////////////////////////////////////////////////////
// This header file just defines some aliases and defaults used by all
// the modules that deal with interpolation.
////////////////////////////////////////////////////////////////////////

#ifndef _slug_interpolation_H_
#define _slug_interpolation_H_

#include <boost/multi_array.hpp>
extern "C" {
#   include <gsl/gsl_interp.h>
#   include <gsl/gsl_spline.h>
}

// typedef's for boost arrays
typedef boost::multi_array<double, 1> array1d;
typedef boost::multi_array<double, 2> array2d;
typedef boost::multi_array<double, 3> array3d;
typedef boost::multi_array<double, 4> array4d;
typedef boost::multi_array<double, 5> array5d;
typedef boost::multi_array<gsl_spline *, 1> spl_arr_1d;
typedef boost::multi_array<gsl_spline *, 2> spl_arr_2d;
typedef boost::multi_array<gsl_interp_accel *, 1> acc_arr_1d;
typedef boost::multi_array<gsl_interp_accel *, 2> acc_arr_2d;
typedef array2d::array_view<1>::type view1d;
typedef array3d::array_view<2>::type view2d;
typedef array3d::array_view<3>::type view3d;
typedef spl_arr_2d::array_view<1>::type spl_arr_view_1d;
typedef acc_arr_2d::array_view<1>::type acc_arr_view_1d;
typedef boost::multi_array_types::index_range range_t;

// default interpolation types
#if GSLVERSION >= 2
const gsl_interp_type * const slug_default_interpolator = gsl_interp_steffen;
#else
const gsl_interp_type * const slug_default_interpolator = gsl_interp_akima;
#endif

#endif
// _slug_interpolation_H_
