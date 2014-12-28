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

#ifdef __INTEL_COMPILER
// Need this to fix a bug in the intel compilers relating to c++11
namespace std
{
     typedef decltype(nullptr) nullptr_t;
}
#endif
#include "constants.H"
#include "int_tabulated.H"
#include "slug_cluster.H"
#include <algorithm>
#include <cassert>
#include <cmath>
#include <iomanip>

using namespace std;

////////////////////////////////////////////////////////////////////////
// The constructor
////////////////////////////////////////////////////////////////////////
slug_cluster::slug_cluster(const unsigned long my_id, 
			   const double my_mass, 
			   const double time, const slug_PDF *my_imf, 
			   const slug_tracks *my_tracks, 
			   const slug_specsyn *my_specsyn, 
			   const slug_filter_set *my_filters,
			   const slug_extinction *my_extinct,
			   const slug_nebular *my_nebular,
			   const slug_PDF *my_clf) :
  targetMass(my_mass), imf(my_imf), clf(my_clf), tracks(my_tracks), 
  specsyn(my_specsyn), filters(my_filters), extinct(my_extinct),
  nebular(my_nebular), id(my_id), formationTime(time), curTime(time)
{

  // Initialize to non-disrupted
  is_disrupted = false;

  // Populate with stars
  birthMass = imf->drawPopulation(targetMass, stars);

  // If the population only represents part of the mass range due to
  // restrictions on what range is being treated stochastically, be
  // sure to account for that
  nonStochMass = nonStochAliveMass = 
     targetMass * (1.0 - imf->mass_frac_restrict());
  birthMass += nonStochMass;

  // Initialize the living star mass
  aliveMass = birthMass;

  // Sort the stars
  sort(stars.begin(), stars.end());

  // If we were given a lifetime function, use it to draw a lifetime
  if (clf != NULL) {
    lifetime = clf->draw();
  } else {
    lifetime = constants::big;
  }

  // If we are using extinction, draw one
  if (extinct != NULL) {
    A_V = extinct->draw_AV();
  }

  // Initialize status flags for what data has been stored
  spec_set = Lbol_set = data_set = phot_set = false;
}

////////////////////////////////////////////////////////////////////////
// Routine to reset the cluster
////////////////////////////////////////////////////////////////////////
void
slug_cluster::reset(bool keep_id) {

  // Get new ID
  if (!keep_id) id++;

  // Reset the time, the disruption state, and all flags
  curTime = 0.0;
  is_disrupted = false;
  data_set = Lbol_set = spec_set = phot_set = false;

  // Delete current stellar masses and data
  stars.resize(0);
  stardata.resize(0);
  stars.shrink_to_fit();
  stardata.shrink_to_fit();

  // Re-populate with stars
  birthMass = imf->drawPopulation(targetMass, stars);

  // If the population only represents part of the mass range due to
  // restrictions on what range is being treated stochastically, be
  // sure to account for that
  nonStochMass = nonStochAliveMass = 
     targetMass * (1.0 - imf->mass_frac_restrict());
  birthMass += nonStochMass;

  // Initialize the living star mass
  aliveMass = birthMass;

  // Sort the stars
  sort(stars.begin(), stars.end());

  // If we were given a lifetime function, use it to draw a lifetime
  if (clf != NULL) {
    lifetime = clf->draw();
  } else {
    lifetime = constants::big;
  }

  // If we are using extinction, draw one
  if (extinct != NULL) {
    A_V = extinct->draw_AV();
  }
}


////////////////////////////////////////////////////////////////////////
// Advance routine. All we do is destroy all the stars whose lifetime
// is less than the current time, and set the flag to disrupted if
// we've exceeded the survival time.
////////////////////////////////////////////////////////////////////////
void
slug_cluster::advance(double time) {

  // Make sure we're not trying to go back into the past
  assert(time >= curTime);

  // Do nothing if the time hasn't changed
  if (time == curTime) return;

  // Get current age
  double clusterAge = time - formationTime;

  // Get stellar death mass corresponding to this age; save the
  // current one as a temporary in case we need it below
  stellarDeathMass = tracks->death_mass(clusterAge);

  // Traverse the list, popping off stars that have died, and
  // correcting the mass downward as we go
  while (stars.size() > 0) {
    if (stars.back() > stellarDeathMass) stars.pop_back();
    else break;
  }

  // If the maximum mass for non-stochastic treatment is smaller than
  // the stellar death mass, decrease the mass in the non-stochstic
  // bin
  if (stellarDeathMass < imf->get_xStochMin()) {
    nonStochAliveMass = nonStochMass * 
      imf->integral(imf->get_xStochMin(), stellarDeathMass) /
      imf->integral_restricted();
  }

  // Flag if we're disrupted
  if (clusterAge > lifetime) is_disrupted = true;

  // Set current time
  curTime = time;

  // Mark that data are not current
  data_set = spec_set = Lbol_set = phot_set = false;
}


////////////////////////////////////////////////////////////////////////
// Routine to get stellar data at this time
////////////////////////////////////////////////////////////////////////
void slug_cluster::set_isochrone() {

  // Do nothing if already set; if not set, refresh stellar data and
  // flag that it is now current
  if (data_set) return;
  stardata = tracks->get_isochrone(curTime-formationTime, stars);
  data_set = true;
}


////////////////////////////////////////////////////////////////////////
// Routine to get bolometric luminosity at this time
////////////////////////////////////////////////////////////////////////
void slug_cluster::set_Lbol() {

  // Do nothing if already set
  if (Lbol_set) return;

  // Initialize
  Lbol = 0.0;

  // Stochastic stars part
  if (stars.size() > 0) {

    // Refresh the stellar data
    set_isochrone();

    // Add bolometric luminosity from stochastic stars
    for (unsigned int i=0; i<stardata.size(); i++)
      Lbol += pow(10.0, stardata[i].logL);
  }

  // Non-stochastic part
  if (imf->has_stoch_lim())
    Lbol += specsyn->get_Lbol_cts(birthMass, curTime-formationTime);

  // If using extinction, to compute Lbol we need to first compute the
  // spectrum, then integrate it
  if (extinct != NULL) {
    set_spectrum();
    Lbol_ext = int_tabulated::
      integrate(extinct->lambda(), L_lambda_ext) / constants::Lsun;
  }

  // Flag that things are set
  Lbol_set = true;
}



////////////////////////////////////////////////////////////////////////
// Spectral synthesis routine. Note that this routine also sets Lbol
// in the process because the extra cost of computing it is negligible.
////////////////////////////////////////////////////////////////////////
void
slug_cluster::set_spectrum() {

  // Do nothing if already set
  if (spec_set) return;

  // Initialize
  vector<double>::size_type nl = specsyn->n_lambda();
  L_lambda.assign(nl, 0.0);
  Lbol = 0.0;

  // Stochastic stars part
  if (stars.size() > 0) {

    // Refresh the stellar data
    set_isochrone();

    // Get spectrum for stochastic stars
    L_lambda = specsyn->get_spectrum(stardata);

    // Add bolometric luminosity from stochastic stars
    for (unsigned int i=0; i<stardata.size(); i++)
      Lbol += pow(10.0, stardata[i].logL);
  }

  // Non-stochastic part
  if (imf->has_stoch_lim()) {
    double Lbol_tmp;
    vector<double> spec;
    specsyn->get_spectrum_cts(birthMass, curTime-formationTime, spec, 
			      Lbol_tmp);
    for (unsigned int i=0; i<nl; i++) L_lambda[i] += spec[i];
    Lbol += Lbol_tmp;
  }

  // If using extinction, compute the extincted spectrum and the
  // bolometric luminosity after extinction is applied
  if (extinct != NULL) {
    L_lambda_ext = extinct->spec_extinct(A_V, L_lambda);
    Lbol_ext = int_tabulated::
      integrate(extinct->lambda(), L_lambda_ext) 
      / constants::Lsun;
  }

  // Flag that things are set
  spec_set = Lbol_set = true;
}


////////////////////////////////////////////////////////////////////////
// Photometry calculation routine
////////////////////////////////////////////////////////////////////////
void
slug_cluster::set_photometry() {

  // Do nothing if already set
  if (phot_set) return;

  // Compute the spectrum
  set_spectrum();

  // Grab the wavelength table
  const vector<double>& lambda = specsyn->lambda();

  // Compute photometry
  phot = filters->compute_phot(lambda, L_lambda);

  // If any of the photometric values are -big, that indicates that we
  // want the bolometric luminosity, so insert that
  for (vector<double>::size_type i=0; i<phot.size(); i++)
    if (phot[i] == -constants::big) phot[i] = Lbol;

  // If using extinction, compute photometry on the extincted
  // spectrum; in the process, be careful to mask filters whose
  // response curve doesn't overlap with the extincted spectrum
  if (extinct != NULL) {
    phot_ext = filters->compute_phot(extinct->lambda(), 
				     L_lambda_ext);
    for (vector<double>::size_type i=0; i<phot_ext.size(); i++) {
      if (phot_ext[i] == -constants::big) {
	phot_ext[i] = Lbol_ext;
      } else if (filters->get_filter(i)->photon_filter() &&
		 (filters->get_filter(i)->get_wavelength_min() >
		  extinct->lambda_max())) {
	phot_ext[i] = nan("");
      } else if ((filters->get_filter(i)->get_wavelength_min() <
		  extinct->lambda_min()) ||
		 (filters->get_filter(i)->get_wavelength_max() >
		  extinct->lambda_max())) {
	phot_ext[i] = nan("");
      }
    }
  }
  
  // Flag that the photometry is set
  phot_set = true;
}  


////////////////////////////////////////////////////////////////////////
// Output physical properties
////////////////////////////////////////////////////////////////////////
void
slug_cluster::write_prop(ofstream& outfile, const outputMode out_mode,
			 const unsigned long trial, bool cluster_only) const {
  if (out_mode == ASCII) {
    outfile << setprecision(5) << scientific
	    << setw(11) << right << id << "   "
	    << setw(11) << right << curTime << "   "
	    << setw(11) << right << formationTime << "   "
	    << setw(11) << right << lifetime << "   "
	    << setw(11) << right << targetMass << "   "
	    << setw(11) << right << birthMass << "   "
	    << setw(11) << right << aliveMass << "   "
	    << setw(11) << right << stars.size() << "   ";
    if (stars.size() > 0)
      outfile << setw(11) << right << stars[stars.size()-1];
    else
      outfile << setw(11) << right << 0.0;
    if (extinct != NULL)
      outfile << "   " << setw(11) << right << A_V;
    outfile << endl;
  } else if (out_mode == BINARY) {
    if (cluster_only) {
      outfile.write((char *) &trial, sizeof trial);
      outfile.write((char *) &curTime, sizeof curTime);
      vector<double>::size_type n = 1;
      outfile.write((char *) &n, sizeof n);
    }
    outfile.write((char *) &id, sizeof id);
    outfile.write((char *) &formationTime, sizeof formationTime);
    outfile.write((char *) &lifetime, sizeof lifetime);
    outfile.write((char *) &targetMass, sizeof targetMass);
    outfile.write((char *) &birthMass, sizeof birthMass);
    outfile.write((char *) &aliveMass, sizeof aliveMass);
    vector<double>::size_type n = stars.size();
    outfile.write((char *) &n, sizeof n);
    if (stars.size() > 0)
      outfile.write((char *) &(stars[stars.size()-1]), 
		    sizeof stars[stars.size()-1]);
    else {
      double mstar = 0.0;
      outfile.write((char *) &mstar, sizeof mstar);
    }
    if (extinct != NULL)
      outfile.write((char *) &A_V, sizeof A_V);
  }
}


////////////////////////////////////////////////////////////////////////
// Output physical properties in FITS mode
////////////////////////////////////////////////////////////////////////
#ifdef ENABLE_FITS
void
slug_cluster::write_prop(fitsfile *out_fits, unsigned long trial) {

  // Get current number of entries
  int fits_status = 0;
  long nrows = 0;
  fits_get_num_rows(out_fits, &nrows, &fits_status);

  // Write a new entry
  fits_write_col(out_fits, TULONG, 1, nrows+1, 1, 1, &trial, 
		 &fits_status);
  fits_write_col(out_fits, TULONG, 2, nrows+1, 1, 1, &id,
		 &fits_status);
  fits_write_col(out_fits, TDOUBLE, 3, nrows+1, 1, 1, &curTime,
		 &fits_status);
  fits_write_col(out_fits, TDOUBLE, 4, nrows+1, 1, 1, &formationTime,
		 &fits_status);
  fits_write_col(out_fits, TDOUBLE, 5, nrows+1, 1, 1, &lifetime,
		 &fits_status);
  double tmass = targetMass; // Needed to avoid compiler complaint
  fits_write_col(out_fits, TDOUBLE, 6, nrows+1, 1, 1, &tmass,
		 &fits_status);
  fits_write_col(out_fits, TDOUBLE, 7, nrows+1, 1, 1, &birthMass,
		 &fits_status);
  fits_write_col(out_fits, TDOUBLE, 8, nrows+1, 1, 1, &aliveMass,
		 &fits_status);
  vector<double>::size_type n = stars.size();
  fits_write_col(out_fits, TULONG, 9, nrows+1, 1, 1, &n,
		 &fits_status);
  double mstar;
  if (n>0) mstar = stars.back();
  else mstar = 0.0;
  fits_write_col(out_fits, TDOUBLE, 10, nrows+1, 1, 1, &mstar,
		 &fits_status);
  if (extinct != NULL)
    fits_write_col(out_fits, TDOUBLE, 11, nrows+1, 1, 1, &A_V,
		   &fits_status);
}
#endif


////////////////////////////////////////////////////////////////////////
// Output spectrum
////////////////////////////////////////////////////////////////////////
void
slug_cluster::
write_spectrum(ofstream& outfile, const outputMode out_mode,
	       const unsigned long trial,
	       bool cluster_only) {

  // Make sure information is current
  if (!spec_set) set_spectrum();

  if (out_mode == ASCII) {
    const vector<double>& lambda = specsyn->lambda();
    for (unsigned int i=0; i<lambda.size(); i++) {
      outfile << setprecision(5) << scientific 
	      << setw(11) << right << id << "   "
	      << setw(11) << right << curTime << "   "
	      << setw(11) << right << lambda[i] << "   "
	      << setw(11) << right << L_lambda[i];
      if (extinct != NULL) {
	int j = i - extinct->off();
	if ((j >= 0) && (j < L_lambda_ext.size())) {
	  outfile << "   "
		  << setw(11) << right << L_lambda_ext[j];
	}
      }
      outfile << endl;
    }
  } else {
    if (cluster_only) {
      outfile.write((char *) &trial, sizeof trial);
      outfile.write((char *) &curTime, sizeof curTime);
      vector<double>::size_type n = 1;
      outfile.write((char *) &n, sizeof n);
    }
    outfile.write((char *) &id, sizeof id);
    outfile.write((char *) L_lambda.data(), 
		  L_lambda.size()*sizeof(double));
    if (extinct != NULL)
      outfile.write((char *) L_lambda_ext.data(), 
		    L_lambda_ext.size()*sizeof(double));
  }
}

////////////////////////////////////////////////////////////////////////
// Output spectrum in FITS mode
////////////////////////////////////////////////////////////////////////
#ifdef ENABLE_FITS
void
slug_cluster::
write_spectrum(fitsfile *out_fits, unsigned long trial) {

  // Make sure information is current
  if (!spec_set) set_spectrum();

  // Get current number of entries
  int fits_status = 0;
  long nrows = 0;
  fits_get_num_rows(out_fits, &nrows, &fits_status);

  // Write data
  fits_write_col(out_fits, TULONG, 1, nrows+1, 1, 1, &trial, 
		 &fits_status);
  fits_write_col(out_fits, TULONG, 2, nrows+1, 1, 1, &id, 
		 &fits_status);
  fits_write_col(out_fits, TDOUBLE, 3, nrows+1, 1, 1, &curTime, 
		 &fits_status);
  fits_write_col(out_fits, TDOUBLE, 4, nrows+1, 1, L_lambda.size(), 
		 L_lambda.data(), &fits_status);
  if (extinct != NULL)
    fits_write_col(out_fits, TDOUBLE, 5, nrows+1, 1, L_lambda_ext.size(), 
		   L_lambda_ext.data(), &fits_status);
}
#endif


////////////////////////////////////////////////////////////////////////
// Output photometry
////////////////////////////////////////////////////////////////////////
void
slug_cluster::
write_photometry(ofstream& outfile, const outputMode out_mode,
		 const unsigned long trial,
		 bool cluster_only) {

  // Make sure information is current
  if (!phot_set) set_photometry();

  if (out_mode == ASCII) {
    outfile << setprecision(5) << scientific 
	    << setw(18) << right << id << "   "
	    << setw(18) << right << curTime;
    for (vector<double>::size_type i=0; i<phot.size(); i++)
      outfile << "   " << setw(18) << right << phot[i];
    if (extinct != NULL) {
      for (vector<double>::size_type i=0; i<phot_ext.size(); i++) {
	if (!std::isnan(phot_ext[i]))
	  outfile << "   " << setw(18) << right << phot_ext[i];
	else
	  outfile << "   " << setw(18) << right << " ";
      }
    }
    outfile << endl;
  } else {
    if (cluster_only) {
      outfile.write((char *) &trial, sizeof trial);
      outfile.write((char *) &curTime, sizeof curTime);
      vector<double>::size_type n = 1;
      outfile.write((char *) &n, sizeof n);
    }
    outfile.write((char *) &id, sizeof id);
    outfile.write((char *) phot.data(), 
		  phot.size()*sizeof(double));
    if (extinct != NULL)
      outfile.write((char *) phot_ext.data(), 
		    phot_ext.size()*sizeof(double));
  }
}


////////////////////////////////////////////////////////////////////////
// Output photometry in FITS mode
////////////////////////////////////////////////////////////////////////
#ifdef ENABLE_FITS
void
slug_cluster::
write_photometry(fitsfile *out_fits, unsigned long trial) {

  // Make sure information is current
  if (!phot_set) set_photometry();

  // Get current number of entries
  int fits_status = 0;
  long nrows = 0;
  fits_get_num_rows(out_fits, &nrows, &fits_status);

  // Write data
  fits_write_col(out_fits, TULONG, 1, nrows+1, 1, 1, &trial, 
		 &fits_status);
  fits_write_col(out_fits, TULONG, 2, nrows+1, 1, 1, &id, 
		 &fits_status);
  fits_write_col(out_fits, TDOUBLE, 3, nrows+1, 1, 1, &curTime, 
		 &fits_status);
  for (int i=0; i<phot.size(); i++) {
    int colnum = i+4;
    fits_write_col(out_fits, TDOUBLE, colnum, nrows+1, 1, 1,
		   &(phot[i]), &fits_status);
  }
  if (extinct != NULL) {
    for (int i=0; i<phot_ext.size(); i++) {
      int colnum = i+4+phot_ext.size();
      fits_write_col(out_fits, TDOUBLE, colnum, nrows+1, 1, 1,
		     &(phot_ext[i]), &fits_status);
    }
  }
}
#endif
