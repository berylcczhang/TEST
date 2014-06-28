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
// This file defines the slug_PDF_segment class. This class represents
// a portion of a PDF defined by a particular functional form; it has
// a lower limit and an upper limit, and knows the value of the
// function evaluated at those limits (normalized so that the area
// under the entire function is 1). It also defines two pure virtual
// methods, one which reads data describing the segment from a file,
// and one which draws from the specified PDF.
// 
// slug_PDF_segment is an abstract base class, and particular
// functional forms (e.g. powerlaws, lognormals, Schecter functions,
// etc.) are implemented as derived classes. At present this file
// contains class implementing the following functional forms:
//
// slug_PDF_lognormal: lognormal PDF, characterized by mean and
// dispersion
// 
// slug_PDF_powerlaw: powerlaw PDF, characterized by slope
//
// slug_PDF_normal: normal PDF, characterized by mean and dispersion
//
// slug_PDF_schechter: schechter function PDF, with slope and cutoff
//
////////////////////////////////////////////////////////////////////////

#ifndef _slug_PDF_segment_H_
#define _slug_PDF_segment_H_

#include <boost/random/normal_distribution.hpp>
#include <boost/random/lognormal_distribution.hpp>
#include <boost/random/variate_generator.hpp>
#include <iostream>
#include <fstream>
#include "slug.H"

////////////////////////////////////////////////////////////////////////
// Useful enum
////////////////////////////////////////////////////////////////////////

enum parseStatus { OK, PARSE_ERROR, EOF_ERROR };

////////////////////////////////////////////////////////////////////////
// class slug_PDF_segment
////////////////////////////////////////////////////////////////////////
class slug_PDF_segment {

public:
  // Construct a segment with a specified min and max range
  slug_PDF_segment(double sMin, double sMax);

  // Construct an empty segment
  slug_PDF_segment();

  // Destructor
  virtual ~slug_PDF_segment();

  // Functions that return various internals
  double sMin() { return segMin; }     // Lower range limit
  double sMax() { return segMax; }     // Upper range limit
  double sMinVal() { return segMinVal; } // Value at lower range limit
  double sMaxVal() { return segMaxVal; } // Value at upper range limit

  // Function to draw a star from this IMF segment; pure virtual
  virtual double draw() = 0; 

  // Function to read data for a particular segment from a file; pure
  // virtual. Returns true on successful read, false on failure.
  virtual parseStatus parse(ifstream& file, int& lineCount,
			    string &errMsg, rng_type& rng,
			    double *weight = NULL) = 0;

protected:
  // Data on segment
  double segMin, segMax;  // Lower and upper limits of segment
  double segMinVal;       // Value of segment evaluated at segMin,
			  // before normalization is applied
  double segMaxVal;       // Value of segment evaluated at segMax,
			  // before normalization is applied
};


////////////////////////////////////////////////////////////////////////
// class slug_PDF_lognormal
////////////////////////////////////////////////////////////////////////
class slug_PDF_lognormal : public slug_PDF_segment {

public:
  // Construct a segment with no initialization
  slug_PDF_lognormal();

  // Construct a segment with a specified min and max range, and a
  // specfied center and dispersion; note that sMean is the mean mass,
  // not the mean log mass.
  slug_PDF_lognormal(double sMin, double sMax, double sMean, 
		     double sDisp, rng_type &rng);

  // Destructor
  ~slug_PDF_lognormal();

  // Draw function
  double draw();

  // Parser
  parseStatus parse(ifstream& file, int& lineCount, string &errMsg, 
		    rng_type& rng, double *weight = NULL);

private:
  // Initializer
  void initializer(double sMean, double sDisp, rng_type& rng);

  // The class that will generate lognormal deviates for us
  boost::variate_generator<rng_type&, 
			   boost::lognormal_distribution<> > *lndist;
};


////////////////////////////////////////////////////////////////////////
// class slug_PDF_normal
////////////////////////////////////////////////////////////////////////
class slug_PDF_normal : public slug_PDF_segment {

public:
  // Construct a segment with no initialization
  slug_PDF_normal();

  // Construct a segment with a specified min and max range, and a
  // specfied center and dispersion; note that sMean is the mean mass,
  // not the mean log mass.
  slug_PDF_normal(double sMin, double sMax, double sMean, 
		  double sDisp, rng_type &rng);

  // Destructor
  ~slug_PDF_normal();

  // Draw function
  double draw();

  // Parser
  parseStatus parse(ifstream& file, int& lineCount, string &errMsg, 
		    rng_type& rng, double *weight = NULL);

private:
  // Initializer
  void initializer(double sMean, double sDisp, rng_type& rng);

  // The class that will generate normal deviates for us
  boost::variate_generator<rng_type&, 
			   boost::normal_distribution<> > *ndist;
};


////////////////////////////////////////////////////////////////////////
// class slug_PDF_powerlaw
////////////////////////////////////////////////////////////////////////
class slug_PDF_powerlaw : public slug_PDF_segment {

public:
  // Construct a segment with no initialization
  slug_PDF_powerlaw();

  // Construct a segment with a specified min and max range, and a
  // specfied slope
  slug_PDF_powerlaw(double sMin, double sMax, double sSlope,
		    rng_type &rng);

  // Destructor
  ~slug_PDF_powerlaw();

  // Draw function
  double draw();

  // Parser
  parseStatus parse(ifstream& file, int& lineCount, string &errMsg, 
		    rng_type& rng, double *weight = NULL);

private:
  void initializer(double sSlope, rng_type& rng);

  double segSlope;              // Powerlaw slope

  // The class that will draw uniform deviates for us
  boost::variate_generator<rng_type&, 
			   boost::uniform_01<> > *unidist;
};


////////////////////////////////////////////////////////////////////////
// class slug_PDF_schechter
////////////////////////////////////////////////////////////////////////
class slug_PDF_schechter : public slug_PDF_segment {

public:
  // Construct a segment with no initialization
  slug_PDF_schechter();

  // Construct a segment with a specified min and max range, and a
  // specfied slope
  slug_PDF_schechter(double sMin, double sMax, double sSlope,
		     double sStar, rng_type &rng);

  // Destructor
  ~slug_PDF_schechter();

  // Draw function
  double draw();

  // Parser
  parseStatus parse(ifstream& file, int& lineCount, string &errMsg, 
		    rng_type& rng, double *weight = NULL);

private:
  void initializer(double sSlope, double sStar, rng_type& rng);

  double segSlope;              // Schechter function slope
  double segStar;               // Turnover value

  // The class that will draw uniform deviates for us
  boost::variate_generator<rng_type&, 
			   boost::uniform_01<> > *unidist;
};


#endif
// _slug_PDF_H_
