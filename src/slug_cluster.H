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
#include "slug_galaxy.H"
#include "slug_PDF.H"
#include "slug_specsyn.H"
#include "slug_tracks.H"
#include <iostream>
#include <fstream>

class slug_cluster {

public:
  // Construct a cluster with specified target mass, birth time, IMF
  slug_cluster(const long my_id, const double my_mass, 
	       const double time, slug_PDF *my_imf, 
	       slug_tracks *my_tracks, slug_specsyn *my_specsyn, 
	       slug_PDF *clf=NULL);

  // Destructor
  ~slug_cluster() { };

  // Routine to return the id
  unsigned long get_id() { return id; }

  // Routine to advance to a specified time
  void advance(double time);

  // Routine to report most mass still-living star
  double get_most_massive_star() { 
    if (stars.size() > 0) return stars[stars.size()-1];
    else return 0.0;
  }

  // Routine to report various masses and similar data
  double get_target_mass() { return targetMass; }
  double get_birth_mass() { return birthMass; }
  double get_alive_mass() { return aliveMass; }
  double get_non_stoch_alive_mass() { return nonStochAliveMass; }
  double get_age() { return curTime - formationTime; }
  double get_stellar_death_mass() { return stellarDeathMass; }

  // Routine to report if this cluster is disrupted
  bool disrupted() { return is_disrupted; }

  // Routine to return lifetime of cluster against disruption
  double get_lifetime() { return lifetime; }

  // Routine to return the number of stochastic stars
  unsigned int get_nstars() { return stars.size(); }

  // Routines to return the stellar masses, other stellar data, the
  // bolometric luminosity, the spectrum, and the wavelength data to
  // go with the spectrum. These functions just invoke the
  // corresponding set method and then return.
  const vector<double> &get_stars() { return stars; }
  const vector<slug_stardata> &get_isochrone()
  { set_isochrone(); return stardata; }
  double get_Lbol() 
  { set_Lbol(); return Lbol; }
  const vector<double> &get_spectrum()
  { set_spectrum(); return L_lambda; }
  void get_spectrum(vector<double> &lambda_out, 
		    vector<double> &L_lambda_out)
  { set_spectrum(); 
    lambda_out = specsyn->lambda(); 
    L_lambda_out = L_lambda; }

  // Output routines
  void write_prop(ofstream& outfile, const outputMode out_mode);
  void write_spectrum(ofstream& outfile, const outputMode out_mode);

private:

  // Routines to compute stellar data, bolometric luminosity, spectrum
  void set_isochrone();
  void set_Lbol();
  void set_spectrum();

  // Data
  const unsigned long id;             // Unique ID
  const double targetMass;            // Target mass
  double birthMass;                   // Actual mass at birth
  double aliveMass;                   // Mass of still-living stars
  double nonStochMass;                // Mass not being treated stochastically
  double nonStochAliveMass;           // Mass of still-living stars w/out stochasticity
  double formationTime;               // Cluster birth time
  double curTime;                     // Current time
  double lifetime;                    // Time at which this cluster disrupts
  double stellarDeathMass;            // Mass of star dying at current age
  vector<double> stars;               // List of stellar masses
  vector<slug_stardata> stardata;     // Data on stellar properties
  slug_PDF *imf;                      // IMF
  slug_tracks *tracks;                // Evolutionary track set
  slug_specsyn *specsyn;              // Spectral synthesizer
  bool is_disrupted;                  // Is this cluster disrupted?
  vector<double> L_lambda;            // Spectrum
  double Lbol;                        // Bolometric luminosity
  bool data_set, Lbol_set, spec_set;  // Status indicators
};

#endif
// _slug_cluster_H_
