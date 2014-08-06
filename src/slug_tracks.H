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

////////////////////////////////////////////////////////////////////////
// class slug_tracks
//
// This class holds a set of evolutionary tracks for stars, and
// provides methods to build isochrones from them usng the isocrhone
// synthesis method.
//
// The isochrone synthesis method currently implemented is as
// follows: the raw data are a set of points in the (log t, log m)
// plane, where log m is the initial mass and log t is the time. The
// points lie on lines of constant log m, but not constant log t; at
// each value of log m, there are an equal number of log t points,
// order in increasing time. Graphically, we may think of the data as
// defining a series of quadrilaterals in the (log t, log m)
// plane. The quadrilaterals form a regular array, and their edges in
// the log t direction are parallel to one another, but not in the log
// m direction. We wish to interpolate the data to predict quantities
// for a set of masses at constant t (an isochrone). To do so, we
//
// (1) Build Akima spline interpolations on the grid in the mass and
// time directions. In the mass direction, these are interpolations in
// time at constant mass. In the time direction, since the points are
// not all at the same time, these are interpolations along a sloped
// path in the (log m, log t) plane, connecting the corresponding
// points along parallel mass tracks. This step is done once, during
// initialization.
//
// (2) When an isochrone at a particular time is required, we walk a
// path that starts at the lowest mass track and continues to the
// highest mass track, or to the death line at that time. As this path
// intersects the edges of the quadrilaterials, we record the
// intersection points and use the stored spline interpolations to
// generate values of the desired quantities at them. The result is a
// set of (log m, Q) points, where Q is whatever quantity we're
// interested in.
//
// (3) Given the set of (log m, Q) points, we build a *new* Akima
// spline that interpolates in log m. This represents an isochrone.
//
// (4) We use use the spline representing the isochrone to evaluate
// any quantity Q of interest at whatever mass is desired. We save the
// interpolating function, so that subsequent calls to compute values
// at new masses but the same time do not require that steps (3) and
// (4) be reapeated.
////////////////////////////////////////////////////////////////////////

#ifndef _slug_tracks_H_
#define _slug_tracks_H_

#include <vector>
#include <string>
#include <cmath>
#include <boost/multi_array.hpp>
#include <gsl/gsl_interp.h>
#include <gsl/gsl_spline.h>
#include "slug.H"

typedef boost::multi_array<double, 2> array2d;
typedef std::vector<gsl_spline *> spl_arr;
typedef std::vector<gsl_interp_accel *> acc_arr;

class slug_tracks {

public:
  // Construct a set of tracks from a file containing them
  slug_tracks(const char *fname, double my_metallicity = -1.0,
	      double my_WR_mass = -1.0);

  // Destructor
  ~slug_tracks();

  // Maximum and minimum mass in the tracks
  double max_mass() const { return exp(logmass[0]); }
  double min_mass() const { return exp(logmass[ntrack-1]); }

  // Age of star dying at a particular time
  double death_mass(const double time) const;

  // Lifetime of a star of a specified mass
  double star_lifetime(const double mass) const;

  // Method to return the metallicity (normalized to Solar)
  double get_metallicity() const { return metallicity; }

  // Method to return the minimum mass for the star to have a WR phase
  double get_WR_mass() const { return WR_mass; }

  // Method to produce vector of stellar data from input vector of
  // mass and fixed time. Input mass m is in Msun, and the time t is
  // in yr.
  std::vector<slug_stardata>
  get_isochrone(const double t, const std::vector<double> &m) const;

private:

  // Method to set up the isochrone interpolator for a particular time
  void compute_isochrone(const double logt) const;

  // Private data
  std::string trackfileName;    // Name of the track file
  std::string trackDesc;        // Descriptor string on 1st line of file
  std::vector<std::string> tracktype; // Type of track (WR, RO, ML)
  std::vector<double> logmass;  // Log base e initial masses [Msun]
  array2d logtimes;             // Log base e of times [yr]
  array2d logcur_mass;          // Log base e of current stellar masses [Msun]
  array2d logL;                 // Log base 10 of luminosity [Lsun]
  array2d logTeff;              // Log base 10 of effective temperature [K]
  array2d h_surf;               // Surface abundance of H
  array2d he_surf;              // Surface abundance of He
  array2d c_surf;               // Surface abundance of C
  array2d n_surf;               // Surface abundance of N
  array2d o_surf;               // Surface abundance of O
  array2d logTstar;             // Log base 10 of T_* [K]
  array2d logmDot;              // Log base 10 of wind mass flux [Msun/yr]
  array2d slopes;               // Slopes connectiong points in (log
				// t, logm)
  double metallicity;           // Metallicity normalized to solar
  double WR_mass;               // Mass at which WR stars form
  unsigned int ntrack;          // Number of tracks
  unsigned int ntime;           // Number of times per track

  // Storage for the isochrone
  mutable double isochrone_logt;
  mutable gsl_spline *isochrone_logcur_mass, *isochrone_logL, 
    *isochrone_logTeff, *isochrone_h_surf, *isochrone_c_surf, 
    *isochrone_n_surf;
  mutable gsl_interp_accel *isochrone_logcur_mass_acc, 
    *isochrone_logL_acc, *isochrone_logTeff_acc, 
    *isochrone_h_surf_acc, *isochrone_c_surf_acc,
    *isochrone_n_surf_acc;

  // Interpolation functions for various quantities in the mass and
  // time directions; these correspond to the data arrays listed
  // above. Quantities ending in m_interp are interpolations in mass,
  // and quantities ending in t_interp are interpolations in time.
  spl_arr logcur_mass_m_interp; // Log current mass
  spl_arr logL_m_interp;        // Log luminosity
  spl_arr logTeff_m_interp;     // Log Teff
  spl_arr h_surf_m_interp;      // Surface abundance of H
  spl_arr c_surf_m_interp;      // Surface abundance of C
  spl_arr n_surf_m_interp;      // Surface abundance of N
  spl_arr logcur_mass_t_interp; // Log current mass
  spl_arr logL_t_interp;        // Log luminosity
  spl_arr logTeff_t_interp;     // Log Teff
  spl_arr h_surf_t_interp;      // Surface abundance of H
  spl_arr c_surf_t_interp;      // Surface abundance of C
  spl_arr n_surf_t_interp;      // Surface abundance of N

  // Distances along tracks in the time direction
  array2d trackdist;

  // Accelerators to go with the interpolation functions
  acc_arr logcur_mass_m_acc;    // Log current mass
  acc_arr logL_m_acc;           // Log luminosity
  acc_arr logTeff_m_acc;        // Log Teff
  acc_arr h_surf_m_acc;         // Surface abundance of H
  acc_arr c_surf_m_acc;         // Surface abundance of C
  acc_arr n_surf_m_acc;         // Surface abundance of N
  acc_arr logcur_mass_t_acc;    // Log current mass
  acc_arr logL_t_acc;           // Log luminosity
  acc_arr logTeff_t_acc;        // Log Teff
  acc_arr h_surf_t_acc;         // Surface abundance of H
  acc_arr c_surf_t_acc;         // Surface abundance of C
  acc_arr n_surf_t_acc;         // Surface abundance of N
};

#endif
// _slug_tracks_H_
