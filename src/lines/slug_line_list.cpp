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
#include "slug_line_list.H"
#include "../constants.H"
#include "../slug_MPI.H"
#include <cmath>
#include <iostream>
#include <fstream>
#include <boost/algorithm/string.hpp>
#include <boost/filesystem.hpp>
#include <boost/lexical_cast.hpp>


using namespace std;
using namespace boost;
using namespace boost::algorithm;
using namespace boost::filesystem;

////////////////////////////////////////////////////////////////////////
// The constructor
////////////////////////////////////////////////////////////////////////


slug_line_list::
slug_line_list(const std::vector<std::string>& line_names_,
               const char *line_dir,slug_ostreams &ostreams_) :
               ostreams(ostreams_),
               line_names(line_names_.size())
{
  
  // Try to open the master line list file
  //string fname = "LINE_LIST";
  string fname = "lines_all.dat";
  std::ifstream line_file;
  path dirname(line_dir);
  path line_path = dirname / path(fname.c_str());
  line_file.open(line_path.c_str());
  if (!line_file.is_open()) 
  {
    // Couldn't open file, so bail out
    ostreams.slug_err_one << "unable to open line file " 
                          << line_path.string() << endl;
    bailout(1);
  }
  
  // Read the list of available lines
  vector<string> avail_lines;
  vector<string> tokens;
  vector<double> line_wl, line_ir_lower, line_ir_upper;
  string curline;
  while (getline(line_file, curline)) 
  {
    // Split line into tokens; first token is index, second is line
    // name, third is line wavelength, fourth is lower bound
    // of integration region, fifth is upper bound of integration
    // region.
    trim(curline);
    split(tokens, curline, is_any_of("\t "), token_compress_on);
    avail_lines.push_back(tokens[1]);
    line_wl.push_back(lexical_cast<double>(tokens[2]));
    line_ir_lower.push_back(lexical_cast<double>(tokens[3]));
    line_ir_upper.push_back(lexical_cast<double>(tokens[4]));

  }
  line_file.close();
  ////////////////////////////////////////////////////////////////////
  
  // Find the lines we have been requested to read and store them
  vector<int> line_idx;
  for (unsigned int i = 0; i < line_names_.size(); i++)
  {
    string temp_name = line_names_[i];
    to_lower(temp_name);
    unsigned int j;
    for (j = 0; j < avail_lines.size(); j++)
    {
      string temp_name1 = avail_lines[j];
      to_lower(temp_name1);
      if (temp_name.compare(temp_name1) == 0) 
      {
        line_names[i] = avail_lines[j];
        vector<double> int_reg;                   // Integration region
        int_reg.push_back(line_ir_lower[j]);      // Lower limit
        int_reg.push_back(line_ir_upper[j]);      // Upper limit

        lines.push_back(new slug_line(line_wl[j],int_reg)); 
                                                  // Create the line
                                                  

        break;
      }
    }
    if (j == avail_lines.size())
    {
      ostreams.slug_err_one << "couldn't find line "
      << line_names_[i] << endl;
      bailout(1);
    }
  }

}
 

////////////////////////////////////////////////////////////////////////
// The destructor
////////////////////////////////////////////////////////////////////////
slug_line_list::
~slug_line_list() 
{
  
  for (vector<slug_line *>::size_type i = 0; i < lines.size(); i++)
    delete lines[i];
  
}
////////////////////////////////////////////////////////////////////////
// Compute equivalent width
////////////////////////////////////////////////////////////////////////
std::vector<double> 
slug_line_list::compute_ew(const std::vector<double>& lambda_,
                           const std::vector<double>& L_lambda_) const
{
  std::vector<double> ew;
  // Loop over the lines, calculating equivalent width for each
  for (unsigned int i=0; i<lines.size(); i++)
  {
    // If any of the lines are outside the available spectrum, bail
    if (lines[i]->get_wl_min() < lambda_.front() || 
        lines[i]->get_wl_max() > lambda_.back())
    {
      ostreams.slug_err_one << "line " << line_names[i]
                            << "is outside the available spectrum."
                            << endl;
      bailout(1);
    }
    else
    {
      // Calculate equivalent width for each line and store
      ew.push_back(lines[i]->compute_ew(lambda_,L_lambda_));
    }
  }

  
  return ew;
}
