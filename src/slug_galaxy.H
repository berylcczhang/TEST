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
#include "slug_extinction.H"
#include "slug_IO.H"
#include "slug_nebular.H"
#include "filters/slug_filter_set.H"
#include "specsyn/slug_specsyn.H"
#include "yields/slug_yields.H"
#include <list>
#include <iostream>
#include <fstream>
#include <vector>
#ifdef ENABLE_FITS
extern "C" {
#   include "fitsio.h"
}
#endif

class slug_galaxy {

public:
  // Creation, destruction
  slug_galaxy(const slug_parmParser& pp, const slug_PDF* imf_, 
	      const slug_PDF* cmf_, const slug_PDF* clf_, 
	      const slug_PDF* sfh_, const slug_tracks* tracks_, 
	      const slug_specsyn* specsyn_,
	      const slug_filter_set* filters_,
	      const slug_extinction* extinct_,
	      const slug_nebular *nebular_,
	      const slug_yields *yields_,
	      slug_ostreams& ostreams_);
  ~slug_galaxy();

  // Set galaxy back to its initial state
  void reset(bool reset_cluster_id=false);  

  // Advance to specified time
  void advance(double time);

  // Routines to return the bolometric luminosity, the spectrum, and
  // the photometry. These just invoke the corresponding set method
  // and return.
  double get_Lbol() { set_Lbol(); return Lbol; }

  const std::vector<double> 
  &get_spectrum(const bool del_cluster = false)
  { set_spectrum(del_cluster); return L_lambda; }

  void get_spectrum(std::vector<double>& lambda_out, 
		    std::vector<double>& L_lambda_out,
		    const bool del_cluster = false)
  { set_spectrum(del_cluster);
    lambda_out = specsyn->lambda();
    L_lambda_out = L_lambda; }

  const std::vector<double>& 
  get_photometry(const bool del_cluster = false)
  { set_photometry(del_cluster); return phot; }

  const std::vector<double>&
  get_yield(const bool del_cluster = false)
  { set_yield(del_cluster); return all_yields; }
  
  // Output functions
  void write_integrated_prop(std::ofstream& int_prop_file, 
			     const outputMode out_mode, 
			     const unsigned long trial, const std::vector<double>& imfvp = {});
  void write_cluster_prop(std::ofstream& cluster_prop_file, 
			  const outputMode out_mode,
			  const unsigned long trial, const std::vector<double>& imfvp = {});
  void write_integrated_spec(std::ofstream& int_spec_file,
			     const outputMode out_mode,
			     const unsigned long trial,
			     const bool del_cluster = false);
  void write_cluster_spec(std::ofstream& cluster_spec_file,
			  const outputMode out_mode,
			  const unsigned long trial);
  void write_integrated_phot(std::ofstream& outfile,
			     const outputMode out_mode,
			     const unsigned long trial,
			     const bool del_cluster = false);
  void write_cluster_phot(std::ofstream& outfile,
			  const outputMode out_mode,
			  const unsigned long trial);
  void write_integrated_yield(std::ofstream& outfile,
			      const outputMode out_mode,
			      const unsigned long trial,
			      const bool del_cluster = false);
  void write_cluster_yield(std::ofstream& outfile,
			   const outputMode out_mode,
			   const unsigned long trial);

#ifdef ENABLE_FITS
  // FITS output functions
  void write_integrated_prop(fitsfile* int_prop_fits, unsigned long trial, const std::vector<double>& imfvp = {});
  void write_cluster_prop(fitsfile* cluster_prop_fits, unsigned long trial, const std::vector<double>& imfvp = {});
  void write_integrated_spec(fitsfile* int_spec_fits, unsigned long trial,
			     const bool del_cluster = false);
  void write_cluster_spec(fitsfile* cluster_spec_fits, unsigned long trial);
  void write_integrated_phot(fitsfile* int_phot_fits, unsigned long trial,
			     const bool del_cluster = false);
  void write_cluster_phot(fitsfile* cluster_phot_fits, unsigned long trial);
  void write_integrated_yield(fitsfile* int_yield_fits, unsigned long trial,
			      const bool del_cluster = false);
  void write_cluster_yield(fitsfile* cluster_yield_fits, unsigned long trial);
#endif

protected:
  // IO handler
  slug_ostreams &ostreams;

private:

  // Routines to set stellar data, bolometric luminosity, spectrum,
  // photometry
  void set_field_data();
  void set_Lbol();
  void set_spectrum(const bool del_cluster = false);
  void set_photometry(const bool del_cluster = false);
  void set_yield(const bool del_cluster = false);

  // Distributions and physical models used in simulation
  const slug_PDF *imf;                // IMF
  const slug_PDF *cmf;                // CMF
  const slug_PDF *clf;                // CLF
  const slug_PDF *sfh;                // SFH
  const slug_tracks *tracks;          // Evolutionary track set
  const slug_specsyn *specsyn;        // Spectral synthesizer
  const slug_filter_set *filters;     // Filter set
  const slug_extinction *extinct;     // Extinction object
  const slug_nebular *nebular;        // Nebular object
  const slug_yields *yields;          // Stellar yields

  // Integration helper classes, used to do integral over the IMF or
  // over time to handle non-stochastic stars
  const slug_imf_integrator<double> integ; // Scalar version
  slug_imf_integrator
  <std::vector<double> > v_integ;     // Vector version (used for yields)

  // Simulation data
  double curTime;                     // Current time
  double mass;                        // Mass of stars created so far
  double targetMass;                  // Target mass at this time
  double aliveMass;                   // Mass of currently alive stars
  double stellarMass;                 // Mass of all stars, including remnants
  double nonStochAliveMass;           // Non-stochastic component of aliveMass
  double nonStochRemnantMass;         // Non-stochastic remnant mass
  double clusterAliveMass;            // Cluster component of aliveMass
  double clusterStellarMass;          // Cluster component of stellarMass
  double fieldAliveMass;              // Stochastic field component of aliveMass
  double fieldRemnantMass;            // Mass of field star remnants
  double clusterMass;                 // Mass in non-disrupted clusters
  double fc;                          // Fraction of stars in clusters
  double nonStochFieldMass;           // Birth mass of field stars w/out stochasticity
  double Lbol;                        // Bolometric luminosity of galaxy
  double Lbol_ext;                    // Bolometric luminsoity w/extinction
  double last_yield_time;             // Last time yields were computed
  unsigned long cluster_id;           // Cluster ID counter
  std::vector<slug_star> field_stars; // Field stars
  std::vector<slug_star> dead_field_stars; // Field stars that died this time step
  std::vector<slug_stardata> field_data; // Field star data
  std::vector<double> field_star_AV;  // Field star extinctions
  std::vector<double> field_star_AV_neb;  // Field star nebular extinctions
  std::list<slug_cluster *> disrupted_clusters; // Disrupted clusters
  std::list<slug_cluster *> clusters; // Cluster list
  std::vector<double> L_lambda;       // Specific luminosity of galaxy
  std::vector<double> phot;           // Integrated photometry of galaxy
  std::vector<double> L_lambda_ext;   // Specific luminosity w/extinction
  std::vector<double> phot_ext;       // Integrated photometry w/extinction
  std::vector<double> L_lambda_neb;   // Specific star+nebular luminosity
  std::vector<double> phot_neb;       // Integrated star+nebular photometry
  std::vector<double> L_lambda_neb_ext; // Spec star+nebular lum w/extinction
  std::vector<double> phot_neb_ext;   // Integ star+nebular phot w/extinction
  std::vector<double> stoch_field_yields; // Yields from stochastic field stars
  std::vector<double> all_yields;     // Yields from all stars

  // Status indicators
  bool Lbol_set;                      // Is bolometric luminosity current?
  bool spec_set;                      // Is the spectrum current?
  bool field_data_set;                // Is the field star data current?
  bool phot_set;                      // Is photometry current
  bool yield_set;                     // Are the yields current?
  const slug_line_list *lines = nullptr;        // Line list
};

#endif
// _slug_galaxy_H_

