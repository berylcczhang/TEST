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
#include <boost/bind.hpp>

using namespace std;

////////////////////////////////////////////////////////////////////////
// Trivial little helper function that just takes stardata and returns
// the current stellar mass. Used below.
////////////////////////////////////////////////////////////////////////
namespace cluster {
  double curMass(const slug_stardata &data) {
    return exp(data.logM/constants::loge);
  }
}

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
  nebular(my_nebular), integ(my_tracks, my_imf, nullptr), id(my_id),
  formationTime(time), curTime(time)
{

  // Initialize to non-disrupted
  is_disrupted = false;

  // Populate with stars
  stochBirthMass = stochAliveMass = stochStellarMass =
    imf->drawPopulation(targetMass, stars);

  // If the population only represents part of the mass range due to
  // restrictions on what range is being treated stochastically, be
  // sure to account for that
  nonStochBirthMass = nonStochAliveMass = nonStochStellarMass =
     targetMass * (1.0 - imf->mass_frac_restrict());

  // Set the various mass tallies
  birthMass = stochBirthMass + nonStochBirthMass;
  aliveMass = stochAliveMass + nonStochAliveMass;
  stellarMass = stochStellarMass + nonStochStellarMass;
  stochRemnantMass = nonStochRemnantMass = 0.0;

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
#ifndef __INTEL_COMPILER
  // At this time, intel does not support this c++ function
  stars.shrink_to_fit();
  stardata.shrink_to_fit();
#endif

  // Re-populate with stars
  stochBirthMass = stochAliveMass = stochStellarMass
    = imf->drawPopulation(targetMass, stars);

  // If the population only represents part of the mass range due to
  // restrictions on what range is being treated stochastically, be
  // sure to account for that
  nonStochBirthMass = nonStochAliveMass = nonStochStellarMass
     = targetMass * (1.0 - imf->mass_frac_restrict());

  // Reset the various mass tallies
  birthMass = stochBirthMass + nonStochBirthMass;
  aliveMass = stochAliveMass + nonStochAliveMass;
  stellarMass = stochStellarMass + nonStochStellarMass;
  stochRemnantMass = nonStochRemnantMass = 0.0;

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

  // Handle cases of monotonic and non-monotonic tracks differently
  if (tracks->check_monotonic()) {

    // Monotonic track case

    // Get stellar death mass corresponding to this age; save the
    // current one as a temporary in case we need it below
    stellarDeathMass = tracks->death_mass(clusterAge);

    // Traverse the list, popping off stars that have died, and adding
    // to the remannt mass tally as we go
    while (stars.size() > 0) {
      if (stars.back() > stellarDeathMass) {
	stochRemnantMass += tracks->remnant_mass(stars.back());
	stars.pop_back();
      }
      else break;
    }

  } else {

    // Non-monotonic track case
    vector<double> mass_cuts = tracks->live_mass_range(clusterAge);

    // Traverse the list to pop off stars that are more massive than
    // the upper end of the most massive "alive mass" interval; add to
    // the remnant mass total
    while (stars.size() > 0) {
      if (stars.back() > mass_cuts.back()) {
	stochRemnantMass += tracks->remnant_mass(stars.back());
	stars.pop_back();
      }
      else break;
    }

    // Now kill off stars outside every other alive mass interval
    int interval_ptr = mass_cuts.size()-2;
    int starptr = stars.size()-1;
    while (interval_ptr >= 0) {

      // Stop if we reach the bottom of the star list
      if (starptr < 0) break;

      // Find the lowest mass star that is still alive
      while (stars[starptr] > mass_cuts[interval_ptr]) {
	starptr--;
	if (starptr <= 0) break;
      }

      // If we've reached the bottom of the star list, we're done
      if (starptr <= 0) break;

      // If this star is below the minimum mass star we're killing
      // off, go to next iteration
      if (stars[starptr] < mass_cuts[interval_ptr-1]) {
	interval_ptr -= 2;
	continue;
      }

      // If we're here, starptr points to a star that should be dead
      // and this age. We now proceed through the star list to find
      // the minimum mass stars that is still alive.
      int starptr2 = starptr;
      for (; starptr2 > 0; starptr2--)
	if (stars[starptr2-1] < mass_cuts[interval_ptr-1]) break;

      // The interval from starptr2 to starptr represents stars whose
      // masses are such that they are dead. First add their
      // contribution to the remnant mass, then remove them from the
      // list of stars. Note that we need to add 1 to starptr because
      // the c++ vector erase method excludes the last element.
      for (vector<double>::size_type i=starptr2; i<=starptr; i++)
	stochRemnantMass += tracks->remnant_mass(stars[i]);
      stars.erase(stars.begin()+starptr2, stars.begin()+starptr+1);

      // Move the death mass point and star pointer
      interval_ptr -= 2;
      starptr = starptr2 - 1;

    }

    // Last step: if the minimum alive mass is non-zero, then we need
    // to kill all stars smaller than the smallest alive mass
    if (mass_cuts[0] > 0.0) {

      // Find all stars smaller than the minimum alive mass
      for (starptr=0; starptr<stars.size(); starptr++)
	if (stars[starptr] > mass_cuts[0]) break;

      // Kill those stars, adding to the remnant mass
      for (unsigned int i=0; i<starptr; i++) 
	stochRemnantMass += tracks->remnant_mass(stars[i]);
      stars.erase(stars.begin(), stars.begin()+starptr);

    }
  }

  // Flag if we're disrupted
  if (clusterAge > lifetime) is_disrupted = true;

  // Set current time
  curTime = time;

  // Mark that data are not current
  data_set = spec_set = Lbol_set = phot_set = false;

  // Update all the stellar data to the new isochrone
  set_isochrone();

  // Recompute the present-day masses of the surviving stochastic
  // stars using the new isochrone. We sum the masses of stars below
  // the minimum track mass assuming they have no mass loss, and for
  // all other stars we use the present-day mass from the tracks.
  stochAliveMass = 0.0;
  for (vector<double>::size_type i=0; i<stars.size(); i++) { 
    if (stars[i] >= tracks->min_mass()) break;
    stochAliveMass += stars[i];
  }
  for (vector<slug_stardata>::size_type i=0; i<stardata.size(); i++) 
    stochAliveMass += exp(stardata[i].logM / constants::loge);

  // Now do the same calculation for the non-stochastic stars
  nonStochAliveMass = 0.0;
  if (imf->get_xStochMin() > imf->get_xMin()) {
    if (imf->get_xMin() < tracks->min_mass()) {
      nonStochAliveMass += targetMass * 
	imf->mass_frac(imf->get_xMin(), 
		       min(tracks->min_mass(), imf->get_xStochMin()));
    }
    nonStochAliveMass += integ.integrate(targetMass, curTime-formationTime,
					 boost::bind(cluster::curMass, _1));
  }

  // Recompute the remnant mass from the non-stochastic stars; note a
  // bit of fancy footwork here: the slug_tracks::remnant_mass
  // function is overloaded, so we need to static_cast to the version
  // of it we want before passing to boost::bind
  nonStochRemnantMass = 
    integ.integrate_nt(targetMass, 
		       curTime-formationTime,
		       boost::bind(static_cast<double (slug_tracks::*)
				   (const double, const double) const> 
				   (&slug_tracks::remnant_mass), 
				   tracks, _1, _2));

  // Compute the new alive and total stellar masses
  aliveMass = nonStochAliveMass + stochAliveMass;
  stellarMass = aliveMass + stochRemnantMass + nonStochRemnantMass;
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

  // If using nebular emission, compute the stellar+nebular spectrum
  if (nebular != NULL) 
    L_lambda_neb = nebular->get_tot_spec(L_lambda, this->get_age());

  // If using extinction, compute the extincted spectrum and the
  // bolometric luminosity after extinction is applied
  if (extinct != NULL) {
    L_lambda_ext = extinct->spec_extinct(A_V, L_lambda);
    Lbol_ext = int_tabulated::
      integrate(extinct->lambda(), L_lambda_ext) 
      / constants::Lsun;
    if (nebular != NULL)
      L_lambda_neb_ext = 
	extinct->spec_extinct_neb(A_V, L_lambda_neb);
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

  // Compute photometry
  phot = filters->compute_phot(specsyn->lambda(), L_lambda);

  // If any of the photometric values are -big, that indicates that we
  // want the bolometric luminosity, so insert that
  for (vector<double>::size_type i=0; i<phot.size(); i++)
    if (phot[i] == -constants::big) phot[i] = Lbol;

  // Repeat for stellar+nebular spectrum
  if (nebular != NULL) {
    phot_neb = filters->compute_phot(nebular->lambda(), L_lambda_neb);
    // Special cases: force ionizing luminosity to be zero exactly for
    // this spectrum, and bolometric luminosity to be exactly the same
    // as for the non-nebular case
    for (vector<double>::size_type i=0; i<phot_neb.size(); i++) {
      if (phot_neb[i] == -constants::big) phot_neb[i] = Lbol;
      if (filters->get_filter(i)->photon_filter() &&
	  (filters->get_filter(i)->get_wavelength_max()
	   <= constants::lambdaHI))
	phot_neb[i] = 0.0;
    }
  }

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
		  extinct->lambda().back())) {
	phot_ext[i] = nan("");
      } else if ((filters->get_filter(i)->get_wavelength_min() <
		  extinct->lambda().front()) ||
		 (filters->get_filter(i)->get_wavelength_max() >
		  extinct->lambda().back())) {
	phot_ext[i] = nan("");
      }
    }

    // Repeat for stellar+nebular spectrum
    if (nebular != NULL) {
      phot_neb_ext = filters->compute_phot(extinct->lambda_neb(), 
					   L_lambda_neb_ext);
      for (vector<double>::size_type i=0; i<phot_neb_ext.size(); i++) {
	if (phot_neb_ext[i] == -constants::big) {
	  phot_neb_ext[i] = Lbol_ext;
	} else if (filters->get_filter(i)->photon_filter() &&
		   (filters->get_filter(i)->get_wavelength_min() >
		    extinct->lambda_neb().back())) {
	  phot_neb_ext[i] = nan("");
	} else if ((filters->get_filter(i)->get_wavelength_min() <
		    extinct->lambda_neb().front()) ||
		   (filters->get_filter(i)->get_wavelength_max() >
		    extinct->lambda_neb().back())) {
	  phot_neb_ext[i] = nan("");
	} else if (filters->get_filter(i)->photon_filter() &&
		   (filters->get_filter(i)->get_wavelength_max()
		    <= constants::lambdaHI)) {
	  phot_neb_ext[i] = 0.0;
	}
      }
    }
  }
  
  // Flag that the photometry is set
  phot_set = true;
}  


////////////////////////////////////////////////////////////////////////
// Routines to return the spectrum, without or without nebular
// contributions and extinction, and with our without the wavelength
// table
////////////////////////////////////////////////////////////////////////
const vector<double> &slug_cluster::get_spectrum() {
  set_spectrum();
  return L_lambda;
}

const vector<double> &slug_cluster::get_spectrum_neb() {
  set_spectrum();
  return L_lambda_neb;
}

const vector<double> &slug_cluster::get_spectrum_extinct() {
  set_spectrum();
  return L_lambda_ext;
}

const vector<double> &slug_cluster::get_spectrum_neb_extinct() { 
  set_spectrum(); 
  return L_lambda_neb_ext;
}

void slug_cluster::get_spectrum(vector<double> &lambda_out, 
				vector<double> &L_lambda_out,
				bool rest) { 
  set_spectrum(); 
  lambda_out = specsyn->lambda(rest);
  L_lambda_out = L_lambda;
}

void slug_cluster::get_spectrum_neb(vector<double> &lambda_out, 
				    vector<double> &L_lambda_out,
				    bool rest) {
  set_spectrum(); 
  lambda_out = nebular->lambda(rest); 
  L_lambda_out = L_lambda_neb;
}

void slug_cluster::get_spectrum_extinct(vector<double> &lambda_out, 
					vector<double> &L_lambda_out,
					bool rest) { 
  set_spectrum(); 
  lambda_out = extinct->lambda(rest); 
  L_lambda_out = L_lambda_ext;
}

void slug_cluster::get_spectrum_neb_extinct(vector<double> &lambda_out, 
					    vector<double> &L_lambda_out,
					    bool rest) { 
  set_spectrum(); 
  lambda_out = extinct->lambda_neb(rest); 
  L_lambda_out = L_lambda_neb_ext;
}


////////////////////////////////////////////////////////////////////////
// Routines to return the photometry
////////////////////////////////////////////////////////////////////////

const vector<double> &slug_cluster::get_photometry() { 
  set_photometry();
  return phot;
}

const vector<double> &slug_cluster::get_photometry_neb() {
  set_photometry();
  return phot_neb;
}

const vector<double> &slug_cluster::get_photometry_extinct() {
  set_photometry();
  return phot_ext;
}

const vector<double> &slug_cluster::get_photometry_neb_extinct() {
  set_photometry();
  return phot_neb_ext;
}


////////////////////////////////////////////////////////////////////////
// Routines to clear data
////////////////////////////////////////////////////////////////////////
void slug_cluster::clear_spectrum() {
  L_lambda.resize(0); 
  L_lambda_ext.resize(0); 
  L_lambda_neb.resize(0);
  L_lambda_neb_ext.resize(0);
  phot.resize(0);
  phot_neb.resize(0);
  phot_ext.resize(0);
  phot_neb_ext.resize(0);
  spec_set = false;
  phot_set = false;
}

////////////////////////////////////////////////////////////////////////
// Output physical properties
////////////////////////////////////////////////////////////////////////
void
slug_cluster::write_prop(ofstream& outfile, const outputMode out_mode,
			 const unsigned long trial,
			 bool cluster_only, const std::vector<double>& imfvp) const {
			 
  if (out_mode == ASCII) {
    outfile << setprecision(5) << scientific
	    << setw(11) << right << id << "   "
	    << setw(11) << right << curTime << "   "
	    << setw(11) << right << formationTime << "   "
	    << setw(11) << right << lifetime << "   "
	    << setw(11) << right << targetMass << "   "
	    << setw(11) << right << birthMass << "   "
	    << setw(11) << right << aliveMass << "   "
	    << setw(11) << right << stellarMass << "   "
	    << setw(11) << right << stars.size() << "   ";
    if (stars.size() > 0)
      outfile << setw(11) << right << stars[stars.size()-1];
    else
      outfile << setw(11) << right << 0.0;
    if (extinct != NULL)
      outfile << "   " << setw(11) << right << A_V;

    
    if (imfvp.size()>0)
    {
      //Loop over the variable parameters
      for (int p = 0; p<imfvp.size();p++)
      {
        outfile << "   " << setw(11) << right << imfvp[p];
      }
    
    }
    
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
    outfile.write((char *) &stellarMass, sizeof stellarMass);
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
      
    //Write out variable parameter values
    if (imfvp.size()>0)
    {
      //Loop over the variable parameters
      for (int p = 0; p<imfvp.size();p++)
      {
        outfile.write((char *) &imfvp[p], sizeof imfvp[p]);
      }
    
    }
   
    
  }
}


////////////////////////////////////////////////////////////////////////
// Output physical properties in FITS mode
////////////////////////////////////////////////////////////////////////
#ifdef ENABLE_FITS
void
slug_cluster::write_prop(fitsfile *out_fits, unsigned long trial,
                          const std::vector<double>& imfvp) {

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
  fits_write_col(out_fits, TDOUBLE, 9, nrows+1, 1, 1, &stellarMass,
		 &fits_status);
  vector<double>::size_type n = stars.size();
  fits_write_col(out_fits, TULONG, 10, nrows+1, 1, 1, &n,
		 &fits_status);
  double mstar;
  if (n>0) mstar = stars.back();
  else mstar = 0.0;
  fits_write_col(out_fits, TDOUBLE, 11, nrows+1, 1, 1, &mstar,
		 &fits_status);
  if (extinct != NULL)
    fits_write_col(out_fits, TDOUBLE, 12, nrows+1, 1, 1, &A_V,
		   &fits_status);
	if (imfvp.size()>0)
  {
    //Loop over the variable parameters
    int colnum=12;  //Current column number
    for (int p = 0; p<imfvp.size();p++)
    {
      colnum++;
      double vp_p=imfvp[p];
      fits_write_col(out_fits, TDOUBLE, colnum, nrows+1, 1, 1, &vp_p,
		  &fits_status);
    }
  
  }
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
    vector<double> lambda, L_lambda_star, L_lambda_star_ext;
    if (nebular == NULL) {
      lambda = specsyn->lambda();
      L_lambda_star = L_lambda;
      if (extinct != NULL) L_lambda_star_ext = L_lambda_ext;
    } else {
      // Need to output all data using nebular wavelength table in
      // this case, which means we need to interpolate the purely
      // stellar spectra onto it
      lambda = nebular->lambda();
      L_lambda_star = nebular->interp_stellar(L_lambda);
      if (extinct != NULL) 
	L_lambda_star_ext = nebular->interp_stellar(L_lambda_ext, 
						    extinct->off());
    }
    for (unsigned int i=0; i<lambda.size(); i++) {
      outfile << setprecision(5) << scientific 
	      << setw(11) << right << id << "   "
	      << setw(11) << right << curTime << "   "
	      << setw(11) << right << lambda[i] << "   "
	      << setw(11) << right << L_lambda_star[i];
      if (nebular != NULL)
	outfile << "   "
		<< setw(11) << right << L_lambda_neb[i];
      if (extinct != NULL) {
	int j;
	if (nebular == NULL) j = i - extinct->off();
	else j = i - extinct->off_neb();
	if ((j >= 0) && (j < L_lambda_star_ext.size())) {
	  outfile << "   "
		  << setw(11) << right << L_lambda_star_ext[j];
	  if (nebular != NULL)
	    outfile << "   "
		    << setw(11) << right << L_lambda_neb_ext[j];
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
    if (nebular != NULL)
      outfile.write((char *) L_lambda_neb.data(), 
		    L_lambda_neb.size()*sizeof(double));
    if (extinct != NULL) {
      outfile.write((char *) L_lambda_ext.data(), 
		    L_lambda_ext.size()*sizeof(double));
      if (nebular != NULL)
	outfile.write((char *) L_lambda_neb_ext.data(), 
		      L_lambda_neb_ext.size()*sizeof(double));
    }
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
  unsigned int colnum = 5;
  if (nebular != NULL) {
    fits_write_col(out_fits, TDOUBLE, colnum, nrows+1, 1, 
		   L_lambda_neb.size(), L_lambda_neb.data(), 
		   &fits_status);
    colnum++;
  }
  if (extinct != NULL) {
    fits_write_col(out_fits, TDOUBLE, colnum, nrows+1, 1, 
		   L_lambda_ext.size(), L_lambda_ext.data(), 
		   &fits_status);
    colnum++;
    if (nebular != NULL) {
      fits_write_col(out_fits, TDOUBLE, colnum, nrows+1, 1, 
		     L_lambda_neb_ext.size(), L_lambda_neb_ext.data(), 
		     &fits_status);
      colnum++;
    }
  }
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
    if (nebular != NULL)
      for (vector<double>::size_type i=0; i<phot_neb.size(); i++)
	outfile << "   " << setw(18) << right << phot_neb[i];
    if (extinct != NULL) {
      for (vector<double>::size_type i=0; i<phot_ext.size(); i++) {
	if (!std::isnan(phot_ext[i]))
	  outfile << "   " << setw(18) << right << phot_ext[i];
	else
	  outfile << "   " << setw(18) << right << " ";
      }
      if (nebular != NULL) {
	for (vector<double>::size_type i=0; i<phot_neb_ext.size(); i++) {
	  if (!std::isnan(phot_neb_ext[i]))
	    outfile << "   " << setw(18) << right << phot_neb_ext[i];
	  else
	    outfile << "   " << setw(18) << right << " ";
	}
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
    if (nebular != NULL)
      outfile.write((char *) phot_neb.data(), 
		    phot_neb.size()*sizeof(double));
    if (extinct != NULL) {
      outfile.write((char *) phot_ext.data(), 
		    phot_ext.size()*sizeof(double));
      if (nebular != NULL)
	outfile.write((char *) phot_neb_ext.data(), 
		      phot_neb_ext.size()*sizeof(double));
    }
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
  unsigned int colnum = 4;
  for (unsigned int i=0; i<phot.size(); i++) {
    fits_write_col(out_fits, TDOUBLE, colnum, nrows+1, 1, 1,
		   &(phot[i]), &fits_status);
    colnum++;
  }
  if (nebular != NULL) {
    for (unsigned int i=0; i<phot_neb.size(); i++) {
      fits_write_col(out_fits, TDOUBLE, colnum, nrows+1, 1, 1,
		     &(phot_neb[i]), &fits_status);
      colnum++;
    }
  }
  if (extinct != NULL) {
    for (unsigned int i=0; i<phot_ext.size(); i++) {
      fits_write_col(out_fits, TDOUBLE, colnum, nrows+1, 1, 1,
		     &(phot_ext[i]), &fits_status);
      colnum++;
    }
    if (nebular != NULL) {
      for (unsigned int i=0; i<phot_neb_ext.size(); i++) {
	fits_write_col(out_fits, TDOUBLE, colnum, nrows+1, 1, 1,
		       &(phot_neb_ext[i]), &fits_status);
	colnum++;
      }
    }
  }
}
#endif
