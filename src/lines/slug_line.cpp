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

#include "slug_line.H"
#include "../constants.H"
#include "../utils/int_tabulated.H"
#include <cassert>
#include <cmath>
#include <gsl/gsl_interp.h>
#include <gsl/gsl_integration.h>

using namespace std;

////////////////////////////////////////////////////////////////////////
// Constructor
////////////////////////////////////////////////////////////////////////
slug_line::slug_line(const double line_wl_,
                     const vector<double>& integ_region_) :
                     line_wl(line_wl_),
                     integ_region(integ_region_)
{
  // Loaded
}

////////////////////////////////////////////////////////////////////////
// Routine to compute EW
////////////////////////////////////////////////////////////////////////
double
slug_line::compute_ew(const std::vector<double>& lambda,
                      const std::vector<double>& L_lambda) const 
{


  
  bool use_tab  = false; // Use SLUG's int_tabulated integration (takes priority
                         // over use_qags) Default is FALSE.

  bool use_qags = true;  // Use qags integration
                         // Added for testing purposes. Default to TRUE  
                         // If both are false, use eval_integ from gsl_interp
                         
  double ew     = 0.0;   // Equivalent width
  
  // Set up the interpolators
  gsl_interp_accel *acc = gsl_interp_accel_alloc();
  gsl_interp *spectrum_interp = gsl_interp_alloc(gsl_interp_linear,lambda.size() );
  gsl_interp_init(spectrum_interp, lambda.data(), L_lambda.data() , lambda.size() );    
  
  if (use_tab == false)
  {     
    if (use_qags == false)
    {

      // Use the gsl interpolator integrator.
      // This does not have the best results at the boundaries of the
      // integration region. It is slightly faster however. We assume the
      // bin edge is the midpoint between two bin centres.
      
      double gsl_I  = 0.0;  // Integral by gsl
      
      gsl_I =  gsl_interp_eval_integ(spectrum_interp, 
                                     lambda.data(), L_lambda.data(), 
                                     integ_region.front(), 
                                     integ_region.back(), 
                                     acc);
      
      ew = (integ_region.back() - integ_region.front()) - gsl_I;
      
    }
    else
    {
      
      // Use GSL integrator QAGS - need to tune tolerances
      // Also as we only have bin centres as points and these points
      // are not evenly spaced, this will also not be entirely correct.
      // Choosing between the qags, tab, and eval_integ is something of a 
      // choice between three imperfect methods.
      gsl_integration_workspace * interspace
                                = gsl_integration_workspace_alloc (1000);
      // We need a struct to hold the parameters to pass to the 
      // wrapper around the interpolator
      struct extra_params interfunc_params = { lambda, L_lambda, 
                                               acc, spectrum_interp};
                                               
      // Make a gsl function to pass to the integrator
      gsl_function interfunc;
      interfunc.function  = &slug_line::evaluate_interp;
      interfunc.params    = &interfunc_params;
      double gsl_I   = 0.0;  // Integral by gsl
      double gsl_err = 0.0;  // Error from gsl integration
      
      gsl_integration_qags(&interfunc, integ_region.front(),
                           integ_region.back(), 1e-4, 1e-3, 1000,
                           interspace, &gsl_I, &gsl_err);
                           
      ew = (integ_region.back() - integ_region.front()) - gsl_I;
      
      // Free the integration workspace
      gsl_integration_workspace_free(interspace);
    }
    

  }
  else
  {
    // Here we integrate using the int_tabulated method
    // used by slug for calculating the photometry
    
    // We need to find all the grid points within the 
    // integration region. Then we tack on the end
    // points of the integration region.
    // We then interpolate the value at the end points
    // using linear interpolation.
    // Then we pass this to the int_tabulated::integrate
    // function.

    std::vector<double> lambda_Integ;   // Wavelengths for integration
    std::vector<double> L_lambda_Integ; // L_lambda for integration
    
    // Loop through the grid to make a table of the wavelengths
    // and L_lambda that we need
    for (unsigned int wlidx = 0; wlidx < lambda.size(); wlidx++)
    {
      if (lambda[wlidx] >= integ_region.front() &&
          lambda[wlidx] <= integ_region.back())
      {
        lambda_Integ.push_back(lambda[wlidx]);
        L_lambda_Integ.push_back(L_lambda[wlidx]);
      }
    }
    // Handle the (likely) case where the integration region
    // edges are not present in the wavelength grid
    if (lambda_Integ.front() > integ_region.front())
    {
      // We dont have the lower limit
      
      double interpolatedvalue = gsl_interp_eval(spectrum_interp, lambda.data(), 
                                 L_lambda.data(), integ_region.front(), acc); 
      
      lambda_Integ.insert(lambda_Integ.begin(),integ_region.front());
      L_lambda_Integ.insert(L_lambda_Integ.begin(),interpolatedvalue);
    }
    if (lambda_Integ.back() < integ_region.back())
    {
      // We dont have the upper limit
      
      double interpolatedvalue = gsl_interp_eval(spectrum_interp, lambda.data(), 
                                 L_lambda.data(), integ_region.back(), acc); 
      
      lambda_Integ.push_back(integ_region.back());
      L_lambda_Integ.push_back(interpolatedvalue);
    }
    // Now we can call the integrator
    double tab_I = int_tabulated::integrate(lambda_Integ,L_lambda_Integ);
    ew = (integ_region.back() - integ_region.front()) - tab_I;
        
  }
  
  // Free the interpolators
  gsl_interp_free(spectrum_interp);
  gsl_interp_accel_free(acc);
  return ew;
}

// This is a wrapper around the interpolator so that it can be passed
// to QAGS as a gsl function
double
slug_line::evaluate_interp(double wl, void * input_params)
{
  double inter_val = 0.0;    // Interpolated value
  
  struct extra_params * params = (struct extra_params *)input_params;
  
  std::vector<double>lambda   = (params->lambda);     // Wavelengths
  std::vector<double>L_lambda = (params->L_lambda);   // Fluxes
  gsl_interp_accel  *acc      = (params->acc);        // Accelerator
  gsl_interp        *interpol = (params->interpol);   // Interpolator
  
  inter_val = gsl_interp_eval(interpol, lambda.data(), 
                              L_lambda.data(), wl, acc); 
  
  return inter_val;
}

