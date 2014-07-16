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
#include "slug_cluster.H"
#include <algorithm>
#include <cassert>
#include <limits>
#include <iomanip>

#define BIG (numeric_limits<double>::max())

////////////////////////////////////////////////////////////////////////
// The constructor
////////////////////////////////////////////////////////////////////////
slug_cluster::slug_cluster(const long my_id, const double my_mass, 
			   const double time, slug_PDF *my_imf, 
			   slug_tracks *my_tracks, 
			   slug_specsyn *my_specsyn, slug_PDF *clf) :
  id(my_id), targetMass(my_mass), formationTime(time), curTime(time),
  imf(my_imf), tracks(my_tracks), specsyn(my_specsyn)
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
    lifetime = BIG;
  }

  // Initialize flags for the spectrum and Lbol
  spec_set = Lbol_set = false;
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

  // Mark that the spectrum and bolometric luminosity are not current
  spec_set = Lbol_set = false;
}

////////////////////////////////////////////////////////////////////////
// Get isochrone of log L, log Teff, log g values for all the stars in
// the cluster. Stars with masses above or below the limit of the
// tracks are omitted.
////////////////////////////////////////////////////////////////////////
void 
slug_cluster::get_isochrone(vector<double> &logL, 
			    vector<double> &logTeff,
			    vector<double> &logg,
			    vector<double> &logR) {
  tracks->get_isochrone(curTime-formationTime, stars, logL, logTeff, 
			logg, logR);
}


////////////////////////////////////////////////////////////////////////
// Routines to return the stored spectrum and bolometric luminosity
////////////////////////////////////////////////////////////////////////
vector<double>
slug_cluster::get_spectrum() {
  if (!spec_set) set_spectrum();
  return L_lambda;
}

void
slug_cluster::get_spectrum(vector<double>& lambda_out, 
			   vector<double>& L_lambda_out) {
  if (!spec_set) set_spectrum();
  lambda_out = specsyn->lambda();
  L_lambda_out = L_lambda;
}

double
slug_cluster::get_Lbol() {
  //if (!Lbol_set) set_Lbol();
  return Lbol;
}


////////////////////////////////////////////////////////////////////////
// Spectral synthesis routine
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

    // Do isochrone synthesis for stochastic stars
    vector<double> logL, logTeff, logg, logR;
    get_isochrone(logL, logTeff, logg, logR);

    // Get spectrum for stochastic stars
    specsyn->get_spectrum(logL, logTeff, logg, logR, L_lambda);

    // Add bolometric luminosity from stochastic stars
    for (unsigned int i=0; i<logL.size(); i++)
      Lbol += pow(10.0, logL[i]);
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
slug_cluster::write_prop(ofstream& outfile, const outputMode out_mode) {
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
write_spectrum(ofstream& outfile, const outputMode out_mode) {

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
    outfile.write((char *) &id, sizeof id);
    outfile.write((char *) L_lambda.data(), 
		  L_lambda.size()*sizeof(double));
  }
}
