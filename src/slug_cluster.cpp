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

#include "constants.H"
#include "slug_cluster.H"
#include <algorithm>
#include <cassert>
#include <iomanip>

////////////////////////////////////////////////////////////////////////
// The constructor
////////////////////////////////////////////////////////////////////////
slug_cluster::slug_cluster(const long my_id, const double my_mass, 
			   const double time, const slug_PDF *my_imf, 
			   const slug_tracks *my_tracks, 
			   const slug_specsyn *my_specsyn, 
			   const slug_PDF *my_clf) :
  targetMass(my_mass), imf(my_imf), clf(my_clf), tracks(my_tracks), 
  specsyn(my_specsyn), id(my_id), formationTime(time), curTime(time)
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

  // Initialize status flags for what data has been stored
  spec_set = Lbol_set = data_set = false;
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
  data_set = Lbol_set = spec_set = false;

  // Delete current stellar masses and data
  stars.resize(0);
  stardata.resize(0);

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
  int i = stars.size() - 1;
  if (i >= 0) {
    while (stars[i] > stellarDeathMass) {
      aliveMass -= stars.back();
      stars.pop_back();
      i--;
      if (i<0) break;
    }
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
  data_set = spec_set = Lbol_set = false;
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
  unsigned int nl = specsyn->n_lambda();
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

  // Flag that things are set
  spec_set = Lbol_set = true;
}


////////////////////////////////////////////////////////////////////////
// Output physical properties
////////////////////////////////////////////////////////////////////////
void
slug_cluster::write_prop(ofstream& outfile, const outputMode out_mode,
			 bool cluster_only) const {
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
    outfile << endl;
  } else if (out_mode == BINARY) {
    if (cluster_only) {
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
  }
}


////////////////////////////////////////////////////////////////////////
// Output spectrum properties
////////////////////////////////////////////////////////////////////////
void
slug_cluster::
write_spectrum(ofstream& outfile, const outputMode out_mode,
	       bool cluster_only) {

  // Make sure information is current
  if (!spec_set) set_spectrum();

  if (out_mode == ASCII) {
    vector<double> lambda = specsyn->lambda();
    for (unsigned int i=0; i<lambda.size(); i++) {
      outfile << setprecision(5) << scientific 
	      << setw(11) << right << id << "   "
	      << setw(11) << right << curTime << "   "
	      << setw(11) << right << lambda[i] << "   "
	      << setw(11) << right << L_lambda[i]
	      << endl;
    }
  } else {
    if (cluster_only) {
      outfile.write((char *) &curTime, sizeof curTime);
      vector<double>::size_type n = 1;
      outfile.write((char *) &n, sizeof n);
    }
    outfile.write((char *) &id, sizeof id);
    outfile.write((char *) L_lambda.data(), 
		  L_lambda.size()*sizeof(double));
  }
}
