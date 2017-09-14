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
#include "slug_PDF_segment.H"
#include "slug_PDF.H"
#include "../slug_MPI.H"
#include <boost/algorithm/string.hpp>
#include <boost/lexical_cast.hpp>

using namespace std;
using namespace boost;
using namespace boost::algorithm;

////////////////////////////////////////////////////////////////////////
// File parser
////////////////////////////////////////////////////////////////////////
parseStatus
slug_PDF_segment::parse(ifstream& file, int& lineCount, string &errMsg, 
			double *weight) {

  // First figure out which tokens we're supposed to get. Then allocate
  // a vector of doubles to hold them, and a vector of bools to track
  // if we've gotten them.
  const vector<string> tokList = tokenList();
  vector<double> tok_vals(tokList.size());
  vector<bool> have_tok(tokList.size(), false);
  
  //Set this segment's initial variability to false
  variable_seg = false;
		
  // Flags for min, max, and weight. Note that, if weight == NULL,
  // we've been called in basic mode, so we shouldn't expect to read a
  // min, max, or weight; in this case we set them to true to start.
  bool have_weight = (weight == NULL);
  bool have_min = (weight == NULL);
  bool have_max = (weight == NULL);

  // Set the error string, in case we need it
  string errStr = "Expected one of the following:";
  if (!have_weight) errStr += " 'min X', 'max X', 'weight X'";
  for (unsigned int i=0; i<tokList.size(); i++)
    errStr += ", '" + tokList[i] + " X'";

  // Read from file
  vector<string> tokens;
  string line, linecopy;
  while (!file.eof()) {

    // Get a line and trim whitespace
    getline(file, line);
    linecopy = line;
    lineCount++;
    trim(line);

    // Skip comment and blank lines
    if (line.length() == 0) continue;
    if (line.compare(0, 1, "#") == 0) continue;

    // Split line into tokens, and lowercase the first one
    split(tokens, line, is_any_of("\t ,"), token_compress_on);
    to_lower(tokens[0]);

    // Make sure there's no extraneous junk; if there is, bail out
    if (tokens.size() > 2) 
    {

      //Check if there is a variable parameter on this line first
      if (tokens[1].compare(0, 1, "V") == 0 || tokens[1].compare(0,1,"v") == 0)
      {

        //Check if there is extraneous junk on a variable line
        if (tokens.size() > 3)
        {
          if (tokens[3].compare(0, 1, "#") != 0)
          {
            errMsg = "Expected TYPE V PDFFilename.vpar";
            return PARSE_ERROR;
          }
        }
        //Set variable_seg to true
        //Note that this segment has a variable parameter
        variable_seg = true;

      }

      else if (tokens[2].compare(0, 1, "#") != 0) 
      {
        errMsg = errStr;
        return PARSE_ERROR;
      }

    }

    // Make sure we got two tokens
    if (tokens.size() == 1) {
      errMsg = errStr;
      return PARSE_ERROR;
    }

    // Did we get weight, min, or max?
    try {
      if ((tokens[0].compare("weight") == 0) && !have_weight) {
	*weight = lexical_cast<double>(tokens[1]);
	have_weight = true;
      } else if ((tokens[0].compare("min") == 0) && !have_min) {
	segMin = lexical_cast<double>(tokens[1]);
	have_min = true;
      } else if ((tokens[0].compare("max") == 0) && !have_max) {
	segMax = lexical_cast<double>(tokens[1]);
	have_max = true;
      } else {
      
      
	// Did we get one of the other tokens expected by this class?
	unsigned int i;
	for (i=0; i<tokList.size(); i++) {
	  if (tokens[0].compare(tokList[i]) == 0) {
	    
      if (tokens[1].compare(0,1,"V") == 0 || tokens[1].compare(0,1,"v") == 0)
      {

        //Initialise the parameter to a placeholder value
        tok_vals[i] = double(0.5);
        have_tok[i] = true;	    	

        //Store the index of the parameter and the path to the pdf     	    	
        variable_tok.push_back(i);
        variable_names.push_back(tokens[2]);

        alltoks.push_back(tok_vals[i]);	

        //Break
        break;

      }
      else
      {
        // Match found
        tok_vals[i] = lexical_cast<double>(tokens[1]);
        have_tok[i] = true;

        alltoks.push_back(tok_vals[i]);

        break;
      }
    }
  }

	// If we're here and i == tokList.size(), then we've compared
	// this token to min, max, weight, and all the tokens defined by
	// the class, and still not found a match. Bail out with an
	// error in this case
	if (i == tokList.size()) {
	  errMsg = errStr;
	  return PARSE_ERROR;
	}
      }
    } catch (const bad_lexical_cast& ia) {
      // If we're here, a type conversion failed
      (void) ia; // No-op to suppress compiler warning
      errMsg = errStr;
      return PARSE_ERROR;
    }

    // If we're read everything we need, call the initialize routine
    // and then return a success
    bool have_all_tok = have_weight && have_min && have_max;
    for (unsigned int i=0; i<have_tok.size(); i++)
      have_all_tok = have_all_tok && have_tok[i];
    if (have_all_tok) {
      initialize(tok_vals);
      initialised = true;
      return OK;
    }
  }

  // Special case: we're here because we got to EOF, but if we didn't
  // need to read any tokens, that's ok, and we should return
  // success.
  bool have_all_tok = have_weight && have_min && have_max;
  for (unsigned int i=0; i<have_tok.size(); i++)
    have_all_tok = have_all_tok && have_tok[i];
  if (have_all_tok) {
    initialize(tok_vals);
    initialised=true;
    return OK;
  }

  // If we got here, we've reached EOF without having the data we
  // need, so throw an EOF error
  errMsg = "Incomplete data on segment";
  return EOF_ERROR;
}

////////////////////////////////////////////////////////////////////////
// Update variable segments
////////////////////////////////////////////////////////////////////////
void slug_PDF_segment::update(const std::vector<double>& drawn_vals)
{

  //Make sure everything is the correct size
  if (variable_tok.size() != drawn_vals.size())
  {
    ostreams.slug_err_one
      << "Different number of variable tokens to drawn values" << endl;
    ostreams.slug_err_one
      << "slug: VT.size: " << variable_tok.size()
      << "  DV.size: " << drawn_vals.size() << endl;
    bailout(1);
  }

  //Loop over list of indices and modify the appropriate values
  //in the token list.	
  for (vector<double>::size_type i=0; i<variable_tok.size(); i++) {		

    //Assign the newly drawn value to the correct parameter.
    alltoks.at(variable_tok.at(i)) = drawn_vals.at(i);
    
  }

  //Initialize the segment again, now with new parameters
  initialize(alltoks);
   
}
////////////////////////////////////////////////////////////////////////
// Clean up variable segments
////////////////////////////////////////////////////////////////////////
void slug_PDF_segment::delete_v_pdfs()
{

  for (vector<double>::size_type i = 0; i < variable_param_pdfs.size(); ++i) 
  {
    delete variable_param_pdfs[i];
  }

  variable_param_pdfs.clear();

}

