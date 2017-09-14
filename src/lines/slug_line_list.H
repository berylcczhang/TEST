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
// class slug_line_list
//  Stores the line list for the equivalent width calculation
////////////////////////////////////////////////////////////////////////

#ifndef _slug_line_list_H_
#define _slug_line_list_H_

#include "../slug.H"
#include "../slug_IO.H"
#include "slug_line.H"
#include <vector>
#include <string>

class slug_line_list 
{


  public:

    // The constructor; takes the location of the line list file.
    
    slug_line_list(const std::vector<std::string>& line_names,
                   const char *line_dir, 
                   slug_ostreams& ostreams_);
    
    
    // The desctructor
    ~slug_line_list();
    
    // Return to return the names of the lines
    const std::vector<std::string>& get_line_names() const
    { 
      return line_names; 
    }

    // Routine to compute EW over a line list
    std::vector<double> 
    compute_ew(const std::vector<double>& lambda_,
	             const std::vector<double>& L_lambda) const;

    // Routine to provide read-only access to an individual line
    const slug_line *get_line(unsigned int i) const
    {
      return lines[i];
    }
    
  protected:
    
    // Classes to handle I/O
    slug_ostreams &ostreams;
    
  private:

    // Data
    std::vector<std::string>        line_names; // Line names requested
    std::vector<const slug_line *>  lines;      // Line properties 
};

#endif
// _slug_line_list_H_
