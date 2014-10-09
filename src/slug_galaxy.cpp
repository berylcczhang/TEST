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
#include "slug_cluster.H"
#include "slug_galaxy.H"
#include "slug_parmParser.H"
#include "slug_specsyn.H"
#include "slug_tracks.H"
#include <cassert>
#include <cmath>
#include <iomanip>

using namespace std;

////////////////////////////////////////////////////////////////////////
// A trivial little helper function, used below
////////////////////////////////////////////////////////////////////////
namespace galaxy {
  bool sort_death_time_decreasing(const slug_star star1, 
				  const slug_star star2) {
    return (star1.death_time > star2.death_time);
  }
}

////////////////////////////////////////////////////////////////////////
// The constructor
////////////////////////////////////////////////////////////////////////
slug_galaxy::slug_galaxy(const slug_parmParser& pp, 
			 const slug_PDF* my_imf,
			 const slug_PDF* my_cmf, 
			 const slug_PDF* my_clf, 
			 const slug_PDF* my_sfh, 
			 const slug_tracks* my_tracks, 
			 const slug_specsyn* my_specsyn,
			 const slug_filter_set* my_filters,
			 const slug_extinction* my_extinct) :
  imf(my_imf), 
  cmf(my_cmf), 
  clf(my_clf),
  sfh(my_sfh),
  tracks(my_tracks),
  specsyn(my_specsyn),
  filters(my_filters),
  extinct(my_extinct)
 {

  // Initialize mass and time
  curTime = 0.0;
  mass = 0.0;
  targetMass = 0.0;
  aliveMass = 0.0;
  clusterMass = 0.0;
  nonStochFieldMass = 0.0;

  // Get fc
  fc = pp.get_fClust();

  // Initialize the cluster ID pointer
  cluster_id = 0;

  // Initialize status flags
  Lbol_set = spec_set = field_data_set = phot_set = false;
}


////////////////////////////////////////////////////////////////////////
// Destructor
////////////////////////////////////////////////////////////////////////
slug_galaxy::~slug_galaxy() {

  // Destroy cluster lists
  while (disrupted_clusters.size() > 0) {
    if (disrupted_clusters.back() != nullptr)
      delete disrupted_clusters.back();
    disrupted_clusters.pop_back();
  }
  while (clusters.size() > 0) {
    if (clusters.back() != nullptr)
      delete clusters.back();
    clusters.pop_back();
  }
}


////////////////////////////////////////////////////////////////////////
// Reset function -- sets galaxy back to initial state
////////////////////////////////////////////////////////////////////////
void
slug_galaxy::reset(bool reset_cluster_id) {
  // N.B. By default we do NOT reset the cluster_id pointer here,
  // because if we're doing multiple trials, we want the cluster IDs
  // for different trials to be distinct.
  curTime = mass = targetMass = aliveMass = clusterMass 
    = nonStochFieldMass = 0.0;
  Lbol_set = spec_set = field_data_set = phot_set = false;
  field_stars.resize(0);
  while (disrupted_clusters.size() > 0) {
    if (disrupted_clusters.back() != nullptr)
      delete disrupted_clusters.back();
    disrupted_clusters.pop_back();
  }
  while (clusters.size() > 0) {
    if (clusters.back() != nullptr)
      delete clusters.back();
    clusters.pop_back();
  }
  L_lambda.resize(0);
  if (reset_cluster_id) cluster_id = 0;
}

////////////////////////////////////////////////////////////////////////
// Advance routine
////////////////////////////////////////////////////////////////////////
void
slug_galaxy::advance(double time) {

  // Make sure we're not trying to go back into the past
  assert(time >= curTime);

  // Compute mass of new stars to be created; this is equal to the
  // expected mass of stars in the next time interval, plus or minus
  // any deficit or over-production from the previous step. This
  // second part is added to ensure that, if we're doing stop nearest
  // or something like that, we get as close as possible at each
  // time.
  double new_mass = sfh->integral(curTime, time);
  double mass_to_draw = new_mass + targetMass - mass;
  targetMass += new_mass;


  // Create new clusters
  if (fc != 0) {

    // Get masses of new clusters
    vector<double> new_cluster_masses;
    cmf->drawPopulation(fc*mass_to_draw, new_cluster_masses);

    // Get birth times of new clusters
    vector<double> birth_times = sfh->draw(curTime, time,
					   new_cluster_masses.size());

    // Create clusters of chosen masses and birth times, and push them
    // onto the master cluster list
    for (unsigned int i=0; i<new_cluster_masses.size(); i++) {
      slug_cluster *new_cluster = 
	new slug_cluster(cluster_id++, new_cluster_masses[i],
			 birth_times[i], imf, tracks, specsyn, filters,
			 extinct, clf);
      clusters.push_back(new_cluster);
      mass += new_cluster->get_birth_mass();
      aliveMass += new_cluster->get_birth_mass();
      clusterMass += new_cluster->get_birth_mass();
    }
  }

  // Create new field stars
  if (fc != 1) {

    // Get masses of new field stars
    vector<double> new_star_masses;
    imf->drawPopulation((1.0-fc)*mass_to_draw, new_star_masses);

    // Get extinctions to field stars
    vector<double> AV;
    if (extinct != NULL) 
      AV = extinct->draw_AV(new_star_masses.size());

    // Push stars onto field star list; in the process, set the birth
    // time and death time for each of them
    for (unsigned int i=0; i<new_star_masses.size(); i++) {
      slug_star new_star;
      new_star.mass = new_star_masses[i];
      new_star.birth_time = sfh->draw(curTime, time);
      new_star.death_time = new_star.birth_time 
	+ tracks->star_lifetime(new_star.mass);
      field_stars.push_back(new_star);
      if (extinct != NULL) field_star_AV.push_back(AV[i]);
      mass += new_star.mass;
      aliveMass += new_star.mass;
    }

    // Sort field star list by death time, from largest to smallest
    sort(field_stars.begin(), field_stars.end(), 
	 galaxy::sort_death_time_decreasing);

    // Increase the non-stochastic field star mass by the mass of
    // field stars that should have formed below the stochstic limit
    if (imf->has_stoch_lim())
      nonStochFieldMass += 
	(1.0-fc)*new_mass*imf->integral_restricted();
  }

  // Advance all clusters to current time; track how the currently
  // alive star mass in the galaxy changes due to this evolution
  list<slug_cluster *>::iterator it;
  for (it = clusters.begin(); it != clusters.end(); it++) {
    aliveMass -= (*it)->get_alive_mass();
    clusterMass -= (*it)->get_alive_mass();
    (*it)->advance(time);
    aliveMass += (*it)->get_alive_mass();
    clusterMass += (*it)->get_alive_mass();
  }
  for (it = disrupted_clusters.begin(); 
       it != disrupted_clusters.end(); it++) {
    aliveMass -= (*it)->get_alive_mass();
    (*it)->advance(time);
    aliveMass += (*it)->get_alive_mass();
  }

  // See if any clusters were disrupted over the last time step, and,
  // if so, move them to the disrupted list
  it=clusters.begin();
  while (it != clusters.end()) {
    if ((*it)->disrupted()) {
      disrupted_clusters.push_back(*it);
      clusterMass -= (*it)->get_alive_mass();
      it = clusters.erase(it);
    } else {
      ++it;
    }
  }

  
  // Go through the field star list and remove any field stars that
  // have died
  while (field_stars.size() > 0) {
    if (field_stars.back().death_time < curTime) {
      aliveMass -= field_stars.back().mass;
      field_stars.pop_back();
      if (extinct != NULL) field_star_AV.pop_back();
    } else {
      break;
    }
  }

  // Flag that computed quantities are now out of date
  Lbol_set = spec_set = field_data_set = phot_set = false;

  // Store new time
  curTime = time;

}


////////////////////////////////////////////////////////////////////////
// Get stellar data on all field stars
////////////////////////////////////////////////////////////////////////
void
slug_galaxy::set_field_data() {

  // Do nothing if data is current
  if (field_data_set) return;

  // Initialize
  field_data.resize(0);

  // Get field star data
  for (unsigned int i=0; i<field_stars.size(); i++) {
    vector<double> m(1, field_stars[i].mass);
    const vector<slug_stardata> &stardata = 
      tracks->get_isochrone(curTime - field_stars[i].birth_time, m);
    if (stardata.size() > 0) field_data.push_back(stardata[0]);
  }

  // Set status flag
  field_data_set = true;
}


////////////////////////////////////////////////////////////////////////
// Compute Lbol
////////////////////////////////////////////////////////////////////////
void
slug_galaxy::set_Lbol() {

  // Do nothing if already set
  if (Lbol_set) return;

  // Initialize
  Lbol = 0.0;

  // First loop over non-disrupted clusters
  list<slug_cluster *>::iterator it;
  for (it = clusters.begin(); it != clusters.end(); it++)
    Lbol += (*it)->get_Lbol();

  // Now do disrupted clusters
  for (it = disrupted_clusters.begin(); 
       it != disrupted_clusters.end(); 
       it++)
    Lbol += (*it)->get_Lbol();

  // Now do stochastic field stars
  if (!field_data_set) set_field_data();
  for (unsigned int i=0; i<field_data.size(); i++) {
    Lbol += pow(10.0, field_data[0].logL);
  }

  // Now do non-stochastic field stars
  if (imf->has_stoch_lim())
    Lbol += specsyn->get_Lbol_cts_sfh(curTime);

  // Set flag
  Lbol_set = true;
}


////////////////////////////////////////////////////////////////////////
// Compute spectrum, getting Lbol in the process
////////////////////////////////////////////////////////////////////////
void
slug_galaxy::set_spectrum(const bool del_cluster) {

  // Do nothing if already set
  if (spec_set) return;

  // Initialize
  vector<double>::size_type nl = specsyn->n_lambda();
  L_lambda.assign(nl, 0.0);
  Lbol = 0.0;
  vector<double>::size_type nl_ext = 0;
  if (extinct != NULL) {
    nl_ext = extinct->n_lambda();
    L_lambda_ext.assign(nl_ext, 0.0);
    Lbol_ext = 0.0;
  }

  // Loop over non-disrupted clusters; for each one, get spectrum and
  // bolometric luminosity and add both to global sum
  list<slug_cluster *>::iterator it;
  for (it = clusters.begin(); it != clusters.end(); it++) {
    const vector<double>& spec = (*it)->get_spectrum();
    for (vector<double>::size_type i=0; i<nl; i++) 
      L_lambda[i] += spec[i];
    Lbol += (*it)->get_Lbol();
    if (extinct != NULL) {
      const vector<double>& spec_ext = (*it)->get_spectrum_extinct();
      for (vector<double>::size_type i=0; i<nl_ext; i++) 
	L_lambda_ext[i] += spec_ext[i];
      Lbol_ext += (*it)->get_Lbol_extinct();
    }
    if (del_cluster) {
      delete (*it);
      *it = nullptr;
    }
  }

  // Now do exactly the same thing for disrupted clusters
  for (it = disrupted_clusters.begin(); it != disrupted_clusters.end(); 
       it++) {
    const vector<double>& spec = (*it)->get_spectrum();
    for (vector<double>::size_type i=0; i<nl; i++) 
      L_lambda[i] += spec[i];
    Lbol += (*it)->get_Lbol();
    if (extinct != NULL) {
      const vector<double>& spec_ext = (*it)->get_spectrum_extinct();
      for (vector<double>::size_type i=0; i<nl_ext; i++) 
	L_lambda_ext[i] += spec_ext[i];
      Lbol_ext += (*it)->get_Lbol_extinct();
    }
    if (del_cluster) {
      delete (*it);
      *it = nullptr;
    }
  }

  // Now do stochastic field stars
  if (!field_data_set) set_field_data();
  for (vector<slug_stardata>::size_type i=0; i<field_data.size(); i++) {
    const vector<double> &spec = specsyn->get_spectrum(field_data[i]);
    for (vector<double>::size_type j=0; j<nl; j++) 
      L_lambda[j] += spec[j];
    Lbol += pow(10.0, field_data[i].logL);
    if (extinct != NULL) {
      vector<double> spec_ext = 
	extinct->spec_extinct(field_star_AV[i], spec);
      for (vector<double>::size_type j=0; j<nl_ext; j++) 
	L_lambda_ext[j] += spec_ext[j];
      Lbol_ext += int_tabulated::
	integrate(extinct->lambda(), spec_ext) / constants::Lsun;
    }
  }

  // Finally do non-stochastic field stars; note that non-stochastic
  // field stars get extincted by the expectation value of the AV
  // distribution.
  if (imf->has_stoch_lim()) {
    double Lbol_tmp;
    vector<double> spec;
    specsyn->get_spectrum_cts_sfh(curTime, spec, Lbol_tmp);
    for (vector<double>::size_type i=0; i<nl; i++) 
      L_lambda[i] += spec[i];
    Lbol += Lbol_tmp;
    if (extinct != NULL) {
      vector<double> spec_ext =
	extinct->spec_extinct(extinct->AV_expect(), spec);
      for (vector<double>::size_type i=0; i<nl_ext; i++) 
	L_lambda_ext[i] += spec_ext[i];
      Lbol_ext += int_tabulated::
	integrate(extinct->lambda(), spec_ext) / constants::Lsun;
    }
  }

  // Set flags
  Lbol_set = spec_set = true;
}


////////////////////////////////////////////////////////////////////////
// Compute photometry
////////////////////////////////////////////////////////////////////////
void
slug_galaxy::set_photometry(const bool del_cluster) {

  // Do nothing if already set
  if (phot_set) return;

  // Compute the spectrum
  set_spectrum(del_cluster);

  // Grab the wavelength table
  const vector<double>& lambda = specsyn->lambda();

  // Compute photometry
  phot = filters->compute_phot(lambda, L_lambda);

  // If any of the photometric values are -big, that indicates that we
  // want the bolometric luminosity, so insert that
  for (vector<double>::size_type i=0; i<phot.size(); i++)
    if (phot[i] == -constants::big) phot[i] = Lbol;

  // Repeat for extincted values
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
// Output integrated properties
////////////////////////////////////////////////////////////////////////
void
slug_galaxy::write_integrated_prop(ofstream& int_prop_file, 
				   const outputMode out_mode, 
				   const unsigned long trial) {

  if (out_mode == ASCII) {
    int_prop_file << setprecision(5) << scientific 
		  << setw(11) << right << curTime << "   "
		  << setw(11) << right << targetMass << "   "
		  << setw(11) << right << mass << "   "
		  << setw(11) << right << aliveMass << "   "
		  << setw(11) << right << clusterMass << "   "
		  << setw(11) << right << clusters.size() << "   "
		  << setw(11) << right << disrupted_clusters.size() << "   "
		  << setw(11) << right << field_stars.size()
		  << endl;
  } else {
    int_prop_file.write((char *) &trial, sizeof trial);
    int_prop_file.write((char *) &curTime, sizeof curTime);
    int_prop_file.write((char *) &targetMass, sizeof targetMass);
    int_prop_file.write((char *) &mass, sizeof mass);
    int_prop_file.write((char *) &aliveMass, sizeof aliveMass);
    int_prop_file.write((char *) &clusterMass, sizeof clusterMass);
    vector<slug_cluster *>::size_type n = clusters.size();
    int_prop_file.write((char *) &n, sizeof n);
    n = disrupted_clusters.size();
    int_prop_file.write((char *) &n, sizeof n);
    n = field_stars.size();
    int_prop_file.write((char *) &n, sizeof n);
  }
}


////////////////////////////////////////////////////////////////////////
// Output integrated properties to fits file
////////////////////////////////////////////////////////////////////////
#ifdef ENABLE_FITS
void
slug_galaxy::write_integrated_prop(fitsfile* int_prop_fits, 
				   unsigned long trial) {

  // Get current number of entries
  int fits_status = 0;
  long nrows = 0;
  fits_get_num_rows(int_prop_fits, &nrows, &fits_status);

  // Write a new entry
  fits_write_col(int_prop_fits, TULONG, 1, nrows+1, 1, 1, &trial, 
		 &fits_status);
  fits_write_col(int_prop_fits, TDOUBLE, 2, nrows+1, 1, 1, &curTime, 
		 &fits_status);
  fits_write_col(int_prop_fits, TDOUBLE, 3, nrows+1, 1, 1, &targetMass, 
		 &fits_status);
  fits_write_col(int_prop_fits, TDOUBLE, 4, nrows+1, 1, 1, &mass, 
		 &fits_status);
  fits_write_col(int_prop_fits, TDOUBLE, 5, nrows+1, 1, 1, &aliveMass, 
		 &fits_status);
  fits_write_col(int_prop_fits, TDOUBLE, 6, nrows+1, 1, 1, 
		 &clusterMass, &fits_status);
  vector<slug_cluster *>::size_type n = clusters.size();
  fits_write_col(int_prop_fits, TULONG, 7, nrows+1, 1, 1, &n,	 
		 &fits_status);
  n = disrupted_clusters.size();
  fits_write_col(int_prop_fits, TULONG, 8, nrows+1, 1, 1, &n,	 
		 &fits_status);
  n = field_stars.size();
  fits_write_col(int_prop_fits, TULONG, 9, nrows+1, 1, 1, &n,	 
		 &fits_status);
}
#endif

////////////////////////////////////////////////////////////////////////
// Output cluster properties
////////////////////////////////////////////////////////////////////////
void
slug_galaxy::write_cluster_prop(ofstream& cluster_prop_file, 
				const outputMode out_mode,
				const unsigned long trial) {

  // In binary mode, write out the time and the number of clusters
  // first, because individual clusters won't write this data
  if (out_mode == BINARY) {
    cluster_prop_file.write((char *) &trial, sizeof trial);
    cluster_prop_file.write((char *) &curTime, sizeof curTime);
    vector<double>::size_type n = clusters.size();
    cluster_prop_file.write((char *) &n, sizeof n);
  }

  // Now write out each cluster
  for (list<slug_cluster *>::iterator it = clusters.begin();
       it != clusters.end(); ++it)
    (*it)->write_prop(cluster_prop_file, out_mode, trial);
}

////////////////////////////////////////////////////////////////////////
// Output cluster properties in FITS mode
////////////////////////////////////////////////////////////////////////
#ifdef ENABLE_FITS
void
slug_galaxy::write_cluster_prop(fitsfile* cluster_prop_fits, 
				unsigned long trial) {
  for (list<slug_cluster *>::iterator it = clusters.begin();
       it != clusters.end(); ++it)
    (*it)->write_prop(cluster_prop_fits, trial);
}
#endif


////////////////////////////////////////////////////////////////////////
// Output integrated spectra
////////////////////////////////////////////////////////////////////////
void
slug_galaxy::write_integrated_spec(ofstream& int_spec_file, 
				   const outputMode out_mode,
				   const unsigned long trial,
				   const bool del_cluster) {

  // Make sure spectrum information is current. If not, compute it.
  if (!spec_set) set_spectrum(del_cluster);

  if (out_mode == ASCII) {
    vector<double> lambda = specsyn->lambda();
    for (vector<double>::size_type i=0; i<lambda.size(); i++) {
      int_spec_file << setprecision(5) << scientific 
		    << setw(11) << right << curTime << "   "
		    << setw(11) << right << lambda[i] << "   "
		    << setw(11) << right << L_lambda[i];
      if (extinct != NULL) {
	int j = i - extinct->off();
	if ((j >= 0) && (j < L_lambda_ext.size())) {
	  int_spec_file << "   "
			<< setw(11) << right << L_lambda_ext[j];
	}
      }
      int_spec_file << endl;
    }
  } else {
    int_spec_file.write((char *) &trial, sizeof trial);
    int_spec_file.write((char *) &curTime, sizeof curTime);
    int_spec_file.write((char *) L_lambda.data(), 
			sizeof(double)*L_lambda.size());
    if (extinct != NULL)
      int_spec_file.write((char *) L_lambda_ext.data(),
			  sizeof(double)*L_lambda_ext.size());
  }
}


////////////////////////////////////////////////////////////////////////
// Output integrated spectra in FITS mode
////////////////////////////////////////////////////////////////////////
#ifdef ENABLE_FITS
void
slug_galaxy::write_integrated_spec(fitsfile* int_spec_fits, 
				   unsigned long trial,
				   const bool del_cluster) {

  // Make sure spectrum information is current. If not, compute it.
  if (!spec_set) set_spectrum(del_cluster);

  // Get current number of entries
  int fits_status = 0;
  long nrows = 0;
  fits_get_num_rows(int_spec_fits, &nrows, &fits_status);

  // Write data
  fits_write_col(int_spec_fits, TULONG, 1, nrows+1, 1, 1, &trial, 
		 &fits_status);
  fits_write_col(int_spec_fits, TDOUBLE, 2, nrows+1, 1, 1, &curTime, 
		 &fits_status);
  fits_write_col(int_spec_fits, TDOUBLE, 3, nrows+1, 1, 
		 L_lambda.size(), L_lambda.data(), &fits_status);
  if (extinct != NULL)
    fits_write_col(int_spec_fits, TDOUBLE, 4, nrows+1, 1, 
		   L_lambda_ext.size(), L_lambda_ext.data(),
		   &fits_status);
}
#endif

////////////////////////////////////////////////////////////////////////
// Output cluster spectra
////////////////////////////////////////////////////////////////////////
void
slug_galaxy::write_cluster_spec(ofstream& cluster_spec_file, 
				const outputMode out_mode,
				const unsigned long trial) {

  // In binary mode, write out the time and the number of clusters
  // first, because individual clusters won't write this data
  if (out_mode == BINARY) {
    cluster_spec_file.write((char *) &trial, sizeof trial);
    cluster_spec_file.write((char *) &curTime, sizeof curTime);
    vector<double>::size_type n = clusters.size();
    cluster_spec_file.write((char *) &n, sizeof n);
  }

  // Now have each cluster write
  for (list<slug_cluster *>::iterator it = clusters.begin();
       it != clusters.end(); ++it)
    (*it)->write_spectrum(cluster_spec_file, out_mode, trial);
}


////////////////////////////////////////////////////////////////////////
// Output cluster spectra in FITS mode
////////////////////////////////////////////////////////////////////////
#ifdef ENABLE_FITS
void
slug_galaxy::write_cluster_spec(fitsfile* cluster_spec_fits, 
				unsigned long trial) {
  for (list<slug_cluster *>::iterator it = clusters.begin();
       it != clusters.end(); ++it)
    (*it)->write_spectrum(cluster_spec_fits, trial);
}
#endif

////////////////////////////////////////////////////////////////////////
// Output integrated photometry
////////////////////////////////////////////////////////////////////////
void
slug_galaxy::write_integrated_phot(ofstream& outfile, 
				   const outputMode out_mode,
				   const unsigned long trial,
				   const bool del_cluster) {

  // Make sure photometric information is current. If not, compute it.
  if (!phot_set) set_photometry(del_cluster);

  if (out_mode == ASCII) {
    outfile << setprecision(5) << scientific 
	    << setw(18) << right << curTime;
    for (vector<double>::size_type i=0; i<phot.size(); i++)
      outfile << "   " << setw(18) << right << phot[i];
    if (extinct != NULL) {
      for (vector<double>::size_type i=0; i<phot_ext.size(); i++) {
	if (!isnan(phot_ext[i]))
	  outfile << "   " << setw(18) << right << phot_ext[i];
	else
	  outfile << "   " << setw(18) << right << " ";
      }
    }
    outfile << endl;
  } else {
    outfile.write((char *) &trial, sizeof trial);
    outfile.write((char *) &curTime, sizeof curTime);
    outfile.write((char *) phot.data(), 
			sizeof(double)*phot.size());
    if (extinct != NULL)
      outfile.write((char *) phot_ext.data(), 
		    sizeof(double)*phot_ext.size());
  }
}

////////////////////////////////////////////////////////////////////////
// Output integrated photometry in FITS mode
////////////////////////////////////////////////////////////////////////
#ifdef ENABLE_FITS
void
slug_galaxy::write_integrated_phot(fitsfile* int_phot_fits, 
				   unsigned long trial,
				   const bool del_cluster) {

  // Make sure photometric information is current. If not, compute it.
  if (!phot_set) set_photometry(del_cluster);

  // Get current number of entries
  int fits_status = 0;
  long nrows = 0;
  fits_get_num_rows(int_phot_fits, &nrows, &fits_status);

  // Write data
  fits_write_col(int_phot_fits, TULONG, 1, nrows+1, 1, 1, &trial, 
		 &fits_status);
  fits_write_col(int_phot_fits, TDOUBLE, 2, nrows+1, 1, 1, &curTime, 
		 &fits_status);
  for (int i=0; i<phot.size(); i++) {
    int colnum = i+3;
    fits_write_col(int_phot_fits, TDOUBLE, colnum, nrows+1, 1, 1,
		   &(phot[i]), &fits_status);
  }
  if (extinct != NULL) {
    for (int i=0; i<phot_ext.size(); i++) {
      int colnum = i+3+phot.size();
      fits_write_col(int_phot_fits, TDOUBLE, colnum, nrows+1, 1, 1,
		     &(phot_ext[i]), &fits_status);
    }
  }

}
#endif

////////////////////////////////////////////////////////////////////////
// Output cluster photometry
////////////////////////////////////////////////////////////////////////
void
slug_galaxy::write_cluster_phot(ofstream& outfile, 
				const outputMode out_mode,
				const unsigned long trial) {

  // In binary mode, write out the time and the number of clusters
  // first, because individual clusters won't write this data
  if (out_mode == BINARY) {
    outfile.write((char *) &trial, sizeof trial);
    outfile.write((char *) &curTime, sizeof curTime);
    vector<double>::size_type n = clusters.size();
    outfile.write((char *) &n, sizeof n);
  }

  // Now have each cluster write
  for (list<slug_cluster *>::iterator it = clusters.begin();
       it != clusters.end(); ++it)
    (*it)->write_photometry(outfile, out_mode, trial);
}

////////////////////////////////////////////////////////////////////////
// Output cluster photometry in FITS mode
////////////////////////////////////////////////////////////////////////
#ifdef ENABLE_FITS
void
slug_galaxy::write_cluster_phot(fitsfile* cluster_phot_fits, 
				unsigned long trial) {
  for (list<slug_cluster *>::iterator it = clusters.begin();
       it != clusters.end(); ++it)
    (*it)->write_photometry(cluster_phot_fits, trial);
}
#endif

