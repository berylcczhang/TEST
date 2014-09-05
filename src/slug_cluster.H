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
// class slug_cluster
//
// This class represents a single star cluster. It knows its age, and
// it knows the masses of all stars it contains. It knows how to draw
// stars from an IMF, how to evolve the star list in time, and how to
// write its proprties to a file.
////////////////////////////////////////////////////////////////////////
#ifndef _slug_cluster_H_
#define _slug_cluster_H_

#include "slug.H"
#include "slug_filter_set.H"
#include "slug_galaxy.H"
#include "slug_PDF.H"
#include "slug_specsyn.H"
#include "slug_tracks.H"
#include <iostream>
#include <fstream>
#ifdef ENABLE_FITS
#   include "fitsio.h"
#endif

class slug_cluster {

public:
  // Construct a cluster with specified target mass, birth time, IMF
  slug_cluster(const unsigned long my_id, const double my_mass, 
	       const double time, const slug_PDF *my_imf, 
	       const slug_tracks *my_tracks, 
	       const slug_specsyn *my_specsyn,
	       const slug_filter_set *my_filters,
	       const slug_PDF *my_clf=NULL);

  // Destructor
  ~slug_cluster() { }

  // Routine to return the id
  unsigned long get_id() const { return id; }

  // Set the cluster back to its initial state; this involves
  // re-drawing all the stars. If keep_id is false, the new ID will be
  // set to the old ID plus 1.
  void reset(bool keep_id = false);

  // Routine to advance to a specified time
  void advance(double time);

  // Routine to report most mass still-living star
  double get_most_massive_star() const { 
    if (stars.size() > 0) return stars[stars.size()-1];
    else return 0.0;
  }

  // Routine to report various masses and similar data
  double get_target_mass() const { return targetMass; }
  double get_birth_mass() const { return birthMass; }
  double get_alive_mass() const { return aliveMass; }
  double get_non_stoch_alive_mass() const { return nonStochAliveMass; }
  double get_age() const { return curTime - formationTime; }
  double get_stellar_death_mass() const { return stellarDeathMass; }

  // Routine to report if this cluster is disrupted
  bool disrupted() const { return is_disrupted; }

  // Routine to return lifetime of cluster against disruption
  double get_lifetime() const { return lifetime; }

  // Routine to return the number of stochastic stars
  std::vector<double>::size_type get_nstars() const 
  { return stars.size(); }

  // Routines to return the stellar masses, other stellar data, the
  // bolometric luminosity, the spectrum, the wavelength data to
  // go with the spectrum, and the photometry. Where something must be
  // computed, these functions just invoke the corresponding set
  // method and then return.
  const std::vector<double> &get_stars() const { return stars; }

  const std::vector<slug_stardata> &get_isochrone()
  { set_isochrone(); return stardata; }

  double get_Lbol() 
  { set_Lbol(); return Lbol; }

  const std::vector<double> &get_spectrum()
  { set_spectrum(); return L_lambda; }

  void get_spectrum(std::vector<double> &lambda_out, 
		    std::vector<double> &L_lambda_out)
  { set_spectrum(); 
    lambda_out = specsyn->lambda(); 
    L_lambda_out = L_lambda; }

  const std::vector<double> &get_photometry()
  { set_photometry(); return phot; }

  // Output routines. The optional arguemnt cluster_only is used to
  // flag that only this cluster is being written, rather than a
  // collection of them, and that this cluster should therefore write
  // out its time and the number 1, indicating that it is the only
  // cluster. This only makes a difference if out_mode == BINARY.
  void write_prop(std::ofstream& outfile, const outputMode out_mode, 
		  const bool cluster_only = false) const;
  void write_spectrum(std::ofstream& outfile, const outputMode out_mode,
		      const bool cluster_only = false);
  void write_photometry(std::ofstream& outfile, const outputMode out_mode,
			const bool cluster_only = false);

#ifdef ENABLE_FITS
  // These are identical to the previous three routines, but they
  // write to a FITS file instead of an ofstream
  void write_prop(fitsfile *out_fits, int trial);
  void write_spectrum(fitsfile *out_fits, int trial);
  void write_photometry(fitsfile *out_fits, int trial);
#endif

private:

  // Routines to compute stellar data, bolometric luminosity, spectrum
  void set_isochrone();
  void set_Lbol();
  void set_spectrum();
  void set_photometry();

  // Invariant data
  const double targetMass;            // Target mass
  const slug_PDF *imf;                // IMF
  const slug_PDF *clf;                // CLF
  const slug_tracks *tracks;          // Evolutionary track set
  const slug_specsyn *specsyn;        // Spectral synthesizer
  const slug_filter_set *filters;     // Photometric filters

  // Malleable data
  unsigned long id;                   // Unique ID
  double birthMass;                   // Actual mass at birth
  double aliveMass;                   // Mass of still-living stars
  double nonStochMass;                // Mass not being treated stochastically
  double nonStochAliveMass;           // Mass of still-living stars w/out stochasticity
  double formationTime;               // Cluster birth time
  double curTime;                     // Current time
  double lifetime;                    // Time at which this cluster disrupts
  double stellarDeathMass;            // Mass of star dying at current age
  double Lbol;                        // Bolometric luminosity
  std::vector<double> stars;          // List of stellar masses
  std::vector<slug_stardata> stardata; // Data on stellar properties
  std::vector<double> L_lambda;       // Spectrum
  std::vector<double> phot;           // Photometry
  bool is_disrupted;                  // Is this cluster disrupted?
  bool data_set, Lbol_set, spec_set, phot_set;  // Status indicators
};

#endif
// _slug_cluster_H_
