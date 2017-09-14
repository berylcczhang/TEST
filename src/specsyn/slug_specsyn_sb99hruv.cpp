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
#include "slug_specsyn_sb99hruv.H"
#include "../constants.H"
#include <boost/filesystem.hpp>
#include <cstdlib>
#include <cmath>
#include <iostream>
#include <fstream>
#include <gsl/gsl_interp.h>
#include <gsl/gsl_integration.h>

using namespace std;
using namespace boost;
using namespace boost::filesystem;

////////////////////////////////////////////////////////////////////////
// The constructor
////////////////////////////////////////////////////////////////////////

slug_specsyn_sb99hruv::
slug_specsyn_sb99hruv(const char *dirname, const slug_tracks *my_tracks, 
                      const slug_PDF *my_imf, slug_PDF *my_sfh, 
                      slug_ostreams& ostreams_,	
                      const double z_in) :
                      slug_specsyn(my_tracks, my_imf, my_sfh, ostreams_, z_in),
                      powr(dirname, my_tracks, my_imf, my_sfh, ostreams_, z_in, false),    
                      hillier(dirname, my_tracks, my_imf, my_sfh, ostreams_, z_in, false),
                      kurucz(dirname, my_tracks, my_imf, my_sfh, ostreams_, z_in, false),
                      pauldrach(dirname, my_tracks, my_imf, my_sfh, ostreams_, z_in, false),
                      planck(kurucz.lambda(), my_tracks, my_imf, my_sfh, ostreams_, z_in)
{
  
  // LOAD IFA ATMOSPHERES
  
  // Choose the atmosphere file closest in metallicity to the
  // metallicity we are using for the tracks
  const string extensions[] = {"m13", "m07", "m04", "p00", "p03"};
  const double zrel[] = { 0.05, 0.2, 0.4, 1.0, 2.0 }; // Z relative to Solar
  const int nfile = 5;
  double zdiff = constants::big;
  int idx = -1;
  for (unsigned int i=0; i<nfile; i++) 
  {
    double zdiff1 = abs(log10(zrel[i]) - log10(my_tracks->get_metallicity()));
  
    if (zdiff > zdiff1) 
    {
      idx = (int) i;
      zdiff = zdiff1;
    }
  }

  // Print warning if this is a big extrapolation in metallicity
  if (zdiff > 0.1) 
  {
    ostreams.slug_warn_one << "stellar track metallicity is [Z/H] = "
                           << log10(my_tracks->get_metallicity())
                           << ", closest IFA atmosphere metallicity available is [Z/H] = "
                           << log10(zrel[idx])
                           << "; calculation will proceed" << endl;
  }

  // -----------------------------------------------------------------------

  // Construct the file name and try to open the files
  
  path dirname_path(dirname);  
  //Line
  string lfname = "ifa_line_"+extensions[idx]+".txt";
  std::ifstream atmos_file_l;
  
  path atmos_path_l = dirname_path / path(lfname.c_str());
  atmos_file_l.open(atmos_path_l.c_str());
  if (!atmos_file_l.is_open()) 
  {
    // Couldn't open file, so bail out
    ostreams.slug_err_one << "unable to open ifa line atmosphere file " 
                          << atmos_path_l.string() << endl;
    exit(1);
  }

  // Save file name
  file_name_l = atmos_path_l.string();
  
  //Continuum
  string cfname = "ifa_cont_"+extensions[idx]+".txt";
  std::ifstream atmos_file_c;

  path atmos_path_c = dirname_path / path(cfname.c_str());
  atmos_file_c.open(atmos_path_c.c_str());
  if (!atmos_file_c.is_open()) 
  {
    // Couldn't open file, so bail out
    ostreams.slug_err_one << "unable to open ifa continuum atmosphere file " 
                          << atmos_path_c.string() << endl;
    exit(1);
  }

  // Save file name
  file_name_c = atmos_path_c.string();


  // Wavelengths
  string wfname = "ifa_wave.txt";
  std::ifstream wave_file;
  
  path wave_path = dirname_path / path(wfname.c_str());
  wave_file.open(wave_path.c_str());
  if (!wave_file.is_open()) 
  {
    // Couldn't open file, so bail out
    ostreams.slug_err_one << "unable to open ifa wavelength file " 
                          << wave_path.string() << endl;
    exit(1);
  }

  // Save file name
  file_name_w = wave_path.string();

  // -----------------------------------------------------------------------

  // Load the high resolution wavelengths
  
  for (double cwl; wave_file >> cwl;)
  {
    lambda_rest_hr.push_back(cwl);
  }
  // Delete the duplicate wavelength
  lambda_rest_hr.erase(lambda_rest_hr.begin() + wlskip);
  
  // Compute observed frame wavelengths
  lambda_obs_hr.resize(lambda_rest_hr.size());
  for (unsigned int i=0; i<lambda_rest_hr.size(); i++)
  {
    lambda_obs_hr[i] = (1.0+z)*lambda_rest_hr[i];
  }
  
  // -----------------------------------------------------------------------

  // Load the atmosphere models for continuum and line

  // Note that the continuum produces a spectrum where the 900-3000A region
  // is continuum flux only from the IFA covered stars. The rest of the
  // spectrum and stellar types are as normal.

  // Resize the flux arrays
  array2d::extent_gen a2dext;                            // 2D array extent
  mflx_l_hr.resize(a2dext[86][lambda_rest_hr.size()]);   // Resize line flux array
  mflx_c_hr.resize(a2dext[86][lambda_rest_hr.size()]);   // Resize line flux array
  
  unsigned int lidx   = 0;                               // Line index
  unsigned int cidx   = 0;                               // Current index
  unsigned int modnum = 0;                               // Model number -1 (Leitherer 2010)

  
  // Line
  for (double lineitem; atmos_file_l >> lineitem;)
  {       
    if (lidx == 0)
    {
      teff_vector.push_back(lineitem);  // Store effective temperature
      lidx++;
    }
    else if (lidx ==1)
    {
      log_g_vector.push_back(lineitem); // Store log g
      lidx++;
    }
    else
    {
      if (lidx != 502) // Skip duplicate
      {
        // Store fluxes for this model
        mflx_l_hr[modnum][cidx] = lineitem; 
        
        //Reset indices if end of wave array
        if (cidx==4198) 
        {
          cidx=0;                         // Reset current index
          lidx=0;
          modnum++;                       // Increment model number
        }
        else
        {
          lidx++;
          cidx++;
        }
      }
      else
      {
        lidx++;
      }
    }
  }

  // Continuum
  lidx = 0;                               // Line index
  cidx = 0;                               // Current index
  modnum = 0;                             // Model number -1 (Leitherer 2010)  
  for (double lineitem; atmos_file_c >> lineitem;)
  {       
    if (lidx == 0 || lidx == 1)
    {
      lidx++;
    }
    else
    {
      if (lidx != 502)
      {
        // Store fluxes for this model
        mflx_c_hr[modnum][cidx] = lineitem; 
        
        //Reset indices if end of wave array
        if (cidx==4198) 
        {
          cidx=0;                         // Reset current index
          lidx=0;
          modnum++;                       // Increment model number
        }
        else
        {
          lidx++;
          cidx++;
        }
      }
      else
      {
        lidx++;
      }
    }
  }
  
  
  // -----------------------------------------------------------------------
   
  // Close files
  atmos_file_l.close();
  atmos_file_c.close();
  wave_file.close();

  // Handle low res construction to grab the wave for the low res case 
  lambda_obs_lr = kurucz.lambda();           //Low res observed spectrum
  lambda_rest_lr = kurucz.lambda(true);      //Low res rest spectrum
  
  
  // Now we make the final wavelength output, adding the high resolution
  // region to the 900-3000A region of the standard resolution spectrum.
  
  bool combined = false;    // Have we combined the High Res region?
  for (unsigned int i=0; i<lambda_rest_lr.size(); i++)
  {
    if (lambda_rest_lr[i] < lambda_rest_hr.front())
    {
      lambda_rest.push_back( lambda_rest_lr[i] );
      lambda_obs.push_back( lambda_obs_lr[i] );
    }
    else if ((lambda_rest_lr[i] < lambda_rest_hr.back()) && combined == false)
    {
      
      lambda_rest.insert(lambda_rest.end(),lambda_rest_hr.begin(),lambda_rest_hr.end());
      lambda_obs.insert(lambda_obs.end(),lambda_obs_hr.begin(),lambda_obs_hr.end());
      combined = true;
    }
    else if (lambda_rest_lr[i] > lambda_rest_hr.back())
    {
      lambda_rest.push_back( lambda_rest_lr[i] );
      lambda_obs.push_back( lambda_obs_lr[i] );
    }
  
  }
  
  v_integ.set_nvec(lambda_rest.size()+1);       //IMF integrator
  
  // Storage for rectified spectrum wavelengths
  recspec_lambda = lambda_rest_hr;              // Recspec WL set
  change_rectify(true);                         // Set flag for rest of code
}



////////////////////////////////////////////////////////////////////////
// Wrapper function that just decides which of type of atmosphere
// model to use for each star, calls the appropriate one for each, and
// returns the summed spectrum
// First implementation does nothing. The latter includes rectified spec
////////////////////////////////////////////////////////////////////////
vector<double>
slug_specsyn_sb99hruv::
get_spectrum(vector<slug_stardata>& stars) const 
{

  vector<double> L_lambda;

  cerr << "Not implemented" << endl;
  exit(1);
  
  return L_lambda;
}


////////////////////////////////////////////////////////////////////////
vector<double>
slug_specsyn_sb99hruv::
get_spectrum(vector<slug_stardata>& stars, vector<double>& recspec) const {

  // Initialize
  vector<double> L_lambda_lr(lambda_rest_lr.size());         // Low res flux
  vector<double> L_lambda_intpl(lambda_rest_lr.size());      // to be interpolated
  vector<double> L_lambda_hr(lambda_rest_hr.size());         // High res flux
  vector<double> L_lambda_hr_c(lambda_rest_hr.size());       // Continuum flux
  vector<double> L_lambda_hr_c_intpl(lambda_rest_hr.size()); // Interpolated cflux
  vector<double> L_lambda_hr_intpl(lambda_rest_hr.size());   // Interpolated flux
  vector<double> L_lambda;                                   // Full flux

  // -----------------------------------------------------------------------

  // Construct lists of stars to pass off the various model atmosphere
  // types
  vector<slug_stardata> stars_WR, stars_OB, stars_ku, stars_pl;
  for (unsigned int i=0; i<stars.size(); i++) 
  {
    if (stars[i].WR != NONE) 
    {
      // WR stars are tagged as such
      stars_WR.push_back(stars[i]);
    }
    else if ((stars[i].logTeff > kurucz.get_logTeff_max()) ||
             (stars[i].logTeff < kurucz.get_logTeff_min())) 
    {
      // Stars with temperature that are too high or low for the
      // models get treated as black bodies
      stars_pl.push_back(stars[i]);
    }
    else if ((stars[i].logTeff > pauldrach.get_logTeff_min()) &&
              (stars[i].logg >= pauldrach.get_logg_min()) &&
              (stars[i].logg <= pauldrach.get_logg_max(stars[i].logTeff))) 
    {
      // Stars in the right temperature and log g range get treated as
      // OB stars
      stars_OB.push_back(stars[i]);
    } 
    else 
    {
      // Otherwise use a Kurucz model
      stars_ku.push_back(stars[i]);
    }
  }

  // Now call each of the synthesizers we need to
  // Planck
  if (stars_pl.size() > 0) 
  {
    const vector<double> & L_lambda_tmp = planck.get_spectrum(stars_pl);
    for (unsigned int i=0; i<L_lambda_lr.size(); i++)
    {
      L_lambda_lr[i] += L_lambda_tmp[i];

    }
  }

  // Kurucz
  if (stars_ku.size() > 0) 
  {
    const vector<double> & L_lambda_tmp = kurucz.get_spectrum(stars_ku);
    for (unsigned int i=0; i<L_lambda_lr.size(); i++)
    {
      L_lambda_lr[i] += L_lambda_tmp[i];

    }
  }

  // Pauldrach / OB stars
  if (stars_OB.size() > 0) 
  {
    const vector<double> & L_lambda_tmp = pauldrach.get_spectrum(stars_OB);
    for (unsigned int i=0; i<L_lambda_lr.size(); i++)
    {

      L_lambda_lr[i] += L_lambda_tmp[i];

    }
  }
  
  // Hillier / WR stars
  if (stars_WR.size() > 0) 
  {
    const vector<double> & L_lambda_tmp = hillier.get_spectrum(stars_WR);
    for (unsigned int i=0; i<L_lambda_lr.size(); i++)
    {
      L_lambda_lr[i] += L_lambda_tmp[i];

    }
  }

  // -----------------------------------------------------------------------
  
  // Now we interpolate the spectra in the high resolution region,
  // and then calculate the IFA spectrum.

  
  vector<slug_stardata> stars_WR_hr, stars_pl_hr, stars_ku_hr_1, stars_ku_hr_2, 
                        stars_ifa_hr;
  
  for (unsigned int i=0; i<stars.size(); i++) 
  {
    if (stars[i].WR != NONE) 
    {
      // WR stars are tagged as such
      stars_WR_hr.push_back(stars[i]);
    }
    else if ((stars[i].logTeff > kurucz.get_logTeff_max()) ||
             (stars[i].logTeff < kurucz.get_logTeff_min()))
    {
      // Stars with temperature that are too high or low for the
      // models get treated as black bodies
      stars_pl_hr.push_back(stars[i]);
    }
    else if ((stars[i].logM > log10(10.0))
             && ( stars[i].logTeff < log10(22000.0)))
    {
      // First kurucz section in sb99
      stars_ku_hr_1.push_back(stars[i]);
    } 
    else if ( (stars[i].logM < log10(5.0)) || ((stars[i].logM >= log10(5.0))
             && ( stars[i].logTeff < log10(17000.0)))) 
    {
      // Second kurucz section in sb99
      stars_ku_hr_2.push_back(stars[i]);
    }
    else
    {
      // Otherwise use HRUV model
      stars_ifa_hr.push_back(stars[i]);
    }
  }  

  // Now we have lists of the correct stars we can call the appropriate
  // spectral synthesiser and interpolate the low resolution results

  // Regions to be interpolated go into a placeholder

  // Planck
  if (stars_pl_hr.size() > 0) 
  {
    const vector<double> & L_lambda_tmp = 
      planck.get_spectrum(stars_pl_hr);
    for (unsigned int i=0; i<L_lambda_intpl.size(); i++)
    {
      L_lambda_intpl[i] += L_lambda_tmp[i];
    }
  }

  // Kurucz
  if (stars_ku_hr_1.size() > 0) 
  {
    const vector<double> & L_lambda_tmp = 
      kurucz.get_spectrum(stars_ku_hr_1);
    for (unsigned int i=0; i<L_lambda_intpl.size(); i++)
    {
      L_lambda_intpl[i] += L_lambda_tmp[i];
    }
  }
  
  if (stars_ku_hr_2.size() > 0) 
  {
    const vector<double> & L_lambda_tmp = 
      kurucz.get_spectrum(stars_ku_hr_2);
    for (unsigned int i=0; i<L_lambda_intpl.size(); i++)
    {
      L_lambda_intpl[i] += L_lambda_tmp[i];
    }
  }
  
  // PoWR / WR stars
  if (stars_WR_hr.size() > 0) 
  {
    const vector<double> & L_lambda_tmp = 
      powr.get_spectrum(stars_WR_hr);
    for (unsigned int i=0; i<L_lambda_hr.size(); i++)
    {
      L_lambda_hr[i] += L_lambda_tmp[i];
      L_lambda_hr_c[i] += L_lambda_tmp[i];
    }
  }  
  
  
  // IFA HR Stars
  if (stars_ifa_hr.size() > 0) 
  {
    int modelpick = 0;                    // Nearest neighbour model
    double dimin;                         // Minimum difference
    double di;                            // Current difference
      
    for (unsigned int i=0; i < stars_ifa_hr.size(); i++)
    {

      dimin = 1.0e20;
      for (unsigned int m=0; m < 86; m++)
      {

        di = 1.0 * pow( ( stars_ifa_hr[i].logTeff - log10(teff_vector[m]) ),2.0)
        + 5.0 * pow( (stars_ifa_hr[i].logg - log_g_vector[m]),2.0);

        if (di < dimin)
        {
          dimin=di;
          modelpick=m;
        }

      }

      for (unsigned int wl=0; wl<lambda_rest_hr.size(); wl++)
      {
        double star_rad;      // Star radius
        star_rad = pow(10,stars_ifa_hr[i].logR) * constants::Rsun;
        L_lambda_hr[wl] += 4.0*M_PI*star_rad *star_rad* (mflx_l_hr[modelpick][wl]);
        L_lambda_hr_c[wl] += 4.0*M_PI*star_rad *star_rad* (mflx_c_hr[modelpick][wl]);

      }
    }
  }  
  
  // Now we interpolate the low res spectra to high res and combine 
  // with the low res.

  // Set up the interpolators
  gsl_interp_accel *acc = gsl_interp_accel_alloc();
  gsl_interp *hrinterp = gsl_interp_alloc(gsl_interp_linear,lambda_rest_lr.size() );
  gsl_interp_init(hrinterp, lambda_rest_lr.data(), L_lambda_intpl.data() , lambda_rest_lr.size() );  
  
  // Interpolate, adding on the high res uv
  
  for (unsigned int i=0; i<lambda_rest_hr.size(); i++) 
  {  
    double interpolatedvalue = gsl_interp_eval(hrinterp, lambda_rest_lr.data(), 
                               L_lambda_intpl.data(), lambda_rest_hr[i], acc); 

    L_lambda_hr_intpl[i] = L_lambda_hr[i] + interpolatedvalue;
    L_lambda_hr_c_intpl[i] = L_lambda_hr_c[i] + interpolatedvalue;
    //Rectify spectrum
    recspec[i] = L_lambda_hr_intpl[i]/L_lambda_hr_c_intpl[i];    
  }
  
  // Free the interpolators
  gsl_interp_free(hrinterp);
  gsl_interp_accel_free(acc); 



  // Rescale to offset the offset effect caused by the nearest neighbour method
  // This can be deactivated if using lines or rectified spectrum only.

  bool rescale = true;
  double rescale_ratio = 1;
  if (rescale)
  {
    
    double LRI = 0;
    double HRI = 0;
    
    double lim1 = 2400;
    double lim2 = 2900;
    
    // Integrate under low res
    for (unsigned int i=0; i<lambda_rest_lr.size(); i++)
    {
      if (lambda_rest_lr[i] >=lim1 && lambda_rest_lr[i] <=lim2)
      {
        LRI+=(lambda_rest_lr[i+1] - lambda_rest_lr[i])*(L_lambda_lr[i] + L_lambda_lr[i+1])/2.0 ;
      }
    }
    // Integrate under high res
    for (unsigned int i=0; i<lambda_rest_hr.size(); i++)
    {
      if (lambda_rest_hr[i] >=lim1 && lambda_rest_hr[i] <=lim2)
      {
        HRI+=(lambda_rest_hr[i+1] - lambda_rest_hr[i])*(L_lambda_hr_intpl[i] + L_lambda_hr_intpl[i+1])/2.0 ;
      }    


    }    
    
    // Get ratio
    rescale_ratio = LRI/HRI;

    //Rescale
    for (unsigned int i=0; i<lambda_rest_hr.size(); i++)
    {

      //Rescale spectrum
      L_lambda_hr_intpl[i] = L_lambda_hr_intpl[i]*rescale_ratio;    
    }
    
  }
  
  
  
  // Complete the high res region

  bool combined = false;    // Have we combined the High Res region?
  for (unsigned int i=0; i<lambda_rest_lr.size(); i++)
  {
  
    if (lambda_rest_lr[i] < lambda_rest_hr.front())
    {
      L_lambda.push_back( L_lambda_lr[i] );
    }
    else if ((lambda_rest_lr[i] < lambda_rest_hr.back()) && combined == false)
    {
      L_lambda.insert(L_lambda.end(),L_lambda_hr_intpl.begin(),L_lambda_hr_intpl.end());
      combined = true;       
    }
    else if (lambda_rest_lr[i] > lambda_rest_hr.back())
    {
      L_lambda.push_back( L_lambda_lr[i] );
    }
  
  }  


  return L_lambda;
}



////////////////////////////////////////////////////////////////////////
// Same as the previous function, but for one star at a time
////////////////////////////////////////////////////////////////////////
vector<double>
slug_specsyn_sb99hruv::
get_spectrum(const slug_stardata& stardata) const 
{

  vector<double> L_lambda;          // Final spectrum
  vector<double> L_lambda_hr;       // High resolution spectrum
  vector<double> L_lambda_intpl;    // Spectrum to be interpolated
  vector<double> L_lambda_lr;       // Low resolution spectrum
  // First we do the low resolution part of the spectrum
  
  // Choose the type of synthesizer, and return its output
  if (stardata.WR != NONE) 
  {
    // WR stars are tagged as such
    L_lambda_lr = hillier.get_spectrum(stardata);
  } 
  else if ((stardata.logTeff > kurucz.get_logTeff_max()) ||
	         (stardata.logTeff < kurucz.get_logTeff_min())) 
	{
    // Stars with temperature that are too high or low for the
    // models get treated as black bodies
    L_lambda_lr = planck.get_spectrum(stardata);
  } else if ((stardata.logTeff > pauldrach.get_logTeff_min()) &&
	           (stardata.logg >= pauldrach.get_logg_min()) &&
	           (stardata.logg <= pauldrach.get_logg_max(stardata.logTeff))) 
	{
    // Stars in the right temperature and log g range get treated as
    // OB stars
    L_lambda_lr = pauldrach.get_spectrum(stardata);
  } 
  else 
  {
    // Otherwise use a Kurucz model
    L_lambda_lr = kurucz.get_spectrum(stardata);
  }

  // -----------------------------------------------------------------------
  
  // Now we interpolate the spectrum in the high resolution region,
  // and then calculate the IFA spectrum.

  bool is_lr = false;     // Is this spectrum already at low resolution

  if (stardata.WR != NONE) 
  {
    // WR stars are tagged as such
    L_lambda_hr = powr.get_spectrum(stardata);
  }
  else if ((stardata.logTeff > kurucz.get_logTeff_max()) ||
     (stardata.logTeff < kurucz.get_logTeff_min()))
  {
    // Stars with temperature that are too high or low for the
    // models get treated as black bodies
    L_lambda_intpl = planck.get_spectrum(stardata);
    is_lr = true;
  }
  else if ((stardata.logM > log10(10.0))
           && ( stardata.logTeff < log10(22000.0)))
  {
    // First kurucz section in sb99
    L_lambda_intpl = kurucz.get_spectrum(stardata);
    is_lr = true;
  } 
  else if ( (stardata.logM < log10(5.0)) || ((stardata.logM >= log10(5.0))
           && ( stardata.logTeff < log10(17000.0)))) 
  {
    // Second kurucz section in sb99
    L_lambda_intpl = kurucz.get_spectrum(stardata);
    is_lr = true;
  }
  else
  {
    // Otherwise use HRUV model
    int modelpick = 0;                    // Nearest neighbour model
    double dimin;                         // Minimum difference
    double di;                            // Current difference
      

    dimin = 1.0e20;
    for (unsigned int m=0; m < 86; m++)
    {
        
      di = 1.0 * pow( ( stardata.logTeff - log10(teff_vector[m]) ),2.0)
           + 5.0 * pow( (stardata.logg - log_g_vector[m]),2.0);
      
      if (di < dimin)
      {
        dimin=di;
        modelpick=m;
      }
        
    }

    for (unsigned int wl=0; wl<lambda_rest_hr.size(); wl++)
    {
        double star_rad;      // Star radius
        star_rad = pow(10,stardata.logR) * constants::Rsun;
        L_lambda_hr[wl] += 4.0*M_PI*star_rad *star_rad* (mflx_l_hr[modelpick][wl]);;
        
    }
	  
  } 
  
  if (is_lr == true)
  {
    // Now we interpolate the low res spectra to high res and combine 
    // with the low res.

    // Set up the interpolators
    gsl_interp_accel *acc = gsl_interp_accel_alloc();
    gsl_interp *hrinterp = gsl_interp_alloc(gsl_interp_linear,lambda_rest_lr.size() );
    gsl_interp_init(hrinterp, lambda_rest_lr.data(), L_lambda_intpl.data() , lambda_rest_lr.size() );  
    
    // Interpolate, adding on the high res uv
    
    for (unsigned int i=0; i<lambda_rest_hr.size(); i++) 
    {  
      double interpolatedvalue = gsl_interp_eval(hrinterp, lambda_rest_lr.data(), 
                                 L_lambda_intpl.data(), lambda_rest_hr[i], acc); 

      L_lambda_hr[i] = interpolatedvalue;    
    }
    
    // Free the interpolators
    gsl_interp_free(hrinterp);
    gsl_interp_accel_free(acc); 

    /*
      This is where the rescaling code would go. It may not be necessary
      for the single-star case however.
    */
  }
  // Complete the high res region

  bool combined = false;    // Have we combined the High Res region?
  for (unsigned int i=0; i<lambda_rest_lr.size(); i++)
  {
  
    if (lambda_rest_lr[i] < lambda_rest_hr.front())
    {
      L_lambda.push_back( L_lambda_lr[i] );
    }
    else if ((lambda_rest_lr[i] < lambda_rest_hr.back()) && combined == false)
    {
      L_lambda.insert(L_lambda.end(),L_lambda_hr.begin(),L_lambda_hr.end());
      combined = true;       
    }
    else if (lambda_rest_lr[i] > lambda_rest_hr.back())
    {
      L_lambda.push_back( L_lambda_lr[i] );
    }
  
  }  


  return L_lambda;  
  
  
  
}
