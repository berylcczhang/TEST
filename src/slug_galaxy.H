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
// class slug_galaxy
//
// This class stores the list of clusters and field stars for a
// galaxy, and it knows how to evolve those lists in time. It also
// knows how to call the routines that take these lists as input and
// return spectra and photometry as outputs. It knows how to write its
// physical properties.
////////////////////////////////////////////////////////////////////////

#ifndef _slug_galaxy_H_
#define _slug_galaxy_H_

#include "slug.H"
#include "slug_specsyn.H"
#include <list>
#include <iostream>
#include <fstream>
#include <vector>

using namespace std;

class slug_galaxy {

public:
  // Creation, destruction
  slug_galaxy(slug_parmParser& pp, slug_PDF* my_imf, slug_PDF* my_cmf,
	      slug_PDF* my_clf, slug_PDF* my_sfh, slug_tracks* tracks, 
	      slug_specsyn* my_specsyn);
  ~slug_galaxy();

  // Set galaxy back to its initial state
  void reset(bool reset_cluster_id=false);  

  // Advance to specified time
  void advance(double time);

  // Routines to return the bolometric luminosity and the
  // spectrum. These just invoke the corresponding set method and
  // return.
  double get_Lbol() { set_Lbol(); return Lbol; }
  const vector<double> &get_spectrum()
  { set_spectrum(); return L_lambda; }
  void get_spectrum(vector<double>& lambda_out, 
		    vector<double>& L_lambda_out)
  { set_spectrum();
    lambda_out = specsyn->lambda();
    L_lambda_out = L_lambda; }

  // Functions to open output files and write headers
  void open_integrated_prop(slug_parmParser& pp);  
  void open_cluster_prop(slug_parmParser& pp);
  void open_integrated_spec(slug_parmParser& pp);
  void open_cluster_spec(slug_parmParser& pp);

  // Output functions
  void write_integrated_prop();       // Integrated galaxy properties
  void write_cluster_prop();          // Cluster properties
  void write_integrated_spec();       // Integrated galaxy spectrum
  void write_cluster_spec();          // Cluster spectra

private:

  // Routines to set stellar data, bolometric luminosity, spectrum
  void set_field_data();
  void set_Lbol();
  void set_spectrum();

  // Distributions and physical models used in simulation
  slug_PDF *imf;                      // IMF
  slug_PDF *cmf;                      // CMF
  slug_PDF *clf;                      // CLF
  slug_PDF *sfh;                      // SFH
  slug_tracks *tracks;                // Evolutionary track set
  slug_specsyn *specsyn;              // Spectral synthesizer

  // Simulation data
  double curTime;                     // Current time
  double mass;                        // Mass of stars created so far
  double targetMass;                  // Target mass at this time
  double aliveMass;                   // Mass of currently alive stars
  double clusterMass;                 // Mass in non-disrupted clusters
  double fc;                          // Fraction of stars in clusters
  vector<slug_star> field_stars;      // Field stars
  vector<slug_stardata> field_data;   // Field star data
  double nonStochFieldMass;        // Mass of field stars w/out stochasticity
  list<slug_cluster *> disrupted_clusters; // Disrupted clusters
  list<slug_cluster *> clusters;      // Cluster list
  long cluster_id;                    // Cluster ID counter
  double Lbol;                        // Bolometric luminosity of galaxy
  vector<double> L_lambda;            // Specific luminosity of galaxy

  // Status indicators
  bool Lbol_set;                      // Is bolometric luminosity current?
  bool spec_set;                      // Is the spectrum current?
  bool field_data_set;                // Is the field star data current?

  // Output information
  ofstream int_prop_file;             // Integrated properties file
  ofstream cluster_prop_file;         // Cluster properties file
  ofstream int_spec_file;             // Integrated spectra file
  ofstream cluster_spec_file;         // Cluster spectra file
  outputMode out_mode;                // Output mode
  
};

#endif
// _slug_galaxy_H_

