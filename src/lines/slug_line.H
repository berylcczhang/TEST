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
//
// class slug_line
//
// This class defines a single line, and can compute the
// equivalent width over this line
//
////////////////////////////////////////////////////////////////////////

#ifndef _slug_line_H_
#define _slug_line_H_

#include <string>
#include <vector>
#include "../constants.H"
#include "../utils/int_tabulated.H"
#include <gsl/gsl_interp.h>
class slug_line 
{

  public:
    
    // Constructor
    slug_line(const double line_wl_,
              const std::vector<double>& integ_region_);

    // Destructor
    ~slug_line() { }


    // Routines to compute equivalent width
    double compute_ew(const std::vector<double>& lambda_,
                      const std::vector<double>& L_lambda_) const;


    // Routines to return the wavelength and integration region for this
    // line; these are returned as const objects
    
    const double &get_wl() const 
    { 
      return line_wl; 
    }
    
    const std::vector<double> &get_integ_region() const 
    {
      return integ_region; 
    }   
    
    const double &get_wl_min() const 
    {
      return integ_region.front();
    }
    
    const double &get_wl_max() const
    {
      return integ_region.back();
    }
    
  private:
  
    // Data
    const double              line_wl;      // Line wavelength
    const std::vector<double> integ_region; // Integration region for EW
    
    // Needed for qags - evaluates the interpolated 
    // spectrum at wavelength wl
    static double evaluate_interp(double wl, void * input_params);    
    // Parameters for using GSL integration
    struct extra_params  
    {
      std::vector<double>   lambda;         // Wavelengths
      std::vector<double>   L_lambda;       // Fluxes
      gsl_interp_accel      *acc;           // Interpolator accelerator
      gsl_interp            *interpol;      // Interpolator object
    };
    
};

#endif
// _slug_line_H_

