/*********************************************************************
Copyright (C) 2014-6 Robert da Silva, Michele Fumagalli, Mark
Krumholz, Evan Demers
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
// class slug_yields
//
// This class defines a set of yield tables to be used for the
// calculation of element yields from stellar nucleosynthesis. It can
// read yield tables from data files and interpolate on these tables
// using Steffen interpolation; we use Steffen rather than Akima to
// ensure that yields are never negative, because Steffen guarantees
// that no extrema appear between data points.
//
// In addition to the class itself, we also define two utility
// classes, element_data and isotope_data, to facilitate handling
// nuclear data.
//
////////////////////////////////////////////////////////////////////////

#ifndef _slug_yields_H_
#define _slug_yields_H_

#include "slug.H"
#include <vector>
#include <string>
#include <cmath>
#include <map>
#include <utility>
#include <boost/multi_array.hpp>
extern "C" {
#   include <gsl/gsl_errno.h>
#   include <gsl/gsl_interp.h>
#   include <gsl/gsl_spline.h>
}

typedef boost::multi_array<double, 2> array2d;

////////////////////////////////////////////////////////////////////////
// Utility classes for describing elements and isotopes, and to
// facilitate sorting them
////////////////////////////////////////////////////////////////////////
class element_data {
public:
  element_data(std::string name_, unsigned int number_) : 
    name(name_), number(number_) {};
  const std::string& symbol() const { return name; }
  unsigned int num() const { return number; }
  inline friend bool 
  operator<(const element_data& lhs, const element_data& rhs)
  { return lhs.number < rhs.number; }
  inline friend bool 
  operator> (const element_data& lhs, const element_data& rhs)
  { return rhs < lhs; }
  inline friend bool 
  operator<=(const element_data& lhs, const element_data& rhs)
  { return !(lhs > rhs); }
  inline friend bool 
  operator>=(const element_data& lhs, const element_data& rhs)
  { return !(lhs < rhs); }
  inline friend bool
  operator==(const element_data& lhs, const element_data& rhs)
  { return lhs.number == rhs.number; }
  inline friend bool
  operator!=(const element_data& lhs, const element_data& rhs)
  { return !(lhs == rhs); }
protected:
  std::string name;
  unsigned int number;
};

class isotope_data : public element_data {
public:
  isotope_data(std::string name_, unsigned int number_, 
	       unsigned int weight_, 
	       double lifetime_ = std::numeric_limits<double>::infinity()) :
    element_data(name_, number_), weight(weight_), lifetime(lifetime_) {};
  unsigned int wgt() const { return weight; }
  bool stable() const {
    return !(lifetime < std::numeric_limits<double>::infinity()); }
  double ltime() const { return lifetime; }
  double halflife() const { return lifetime/log(2.0); }
  void set_lifetime(double lifetime_) { this->lifetime = lifetime_; }
  inline friend bool
  operator<(const isotope_data& lhs, const isotope_data& rhs)
  { if (lhs.number == rhs.number) { return lhs.weight < rhs.weight; }
    else { return lhs.number < rhs.number; } }
  inline friend bool 
  operator> (const isotope_data& lhs, const isotope_data& rhs)
  { return rhs < lhs; }
  inline friend bool
  operator<=(const isotope_data& lhs, const isotope_data& rhs)
  { return !(lhs > rhs); }
  inline friend bool
  operator>=(const isotope_data& lhs, const isotope_data& rhs)
  { return !(lhs < rhs); }
  inline friend bool
  operator==(const isotope_data& lhs, const isotope_data& rhs)
  { return (lhs.number == rhs.number) && (lhs.weight == rhs.weight); }
  inline friend bool
  operator!=(const isotope_data& lhs, const isotope_data& rhs)
  { return !(lhs == rhs); }
private:
  unsigned int weight;
  double lifetime;
};


////////////////////////////////////////////////////////////////////////
// The slug_yields class
////////////////////////////////////////////////////////////////////////
class slug_yields {

public:
  // Construct a set of yield tables from files containing them
  slug_yields(const char *yield_dir, const bool no_decay_ = false);

  // Destructor
  ~slug_yields();

  // Returns the list of isotopes and how many of them there are
  std::vector<isotope_data>::size_type get_niso() const { return niso; }
  const std::vector<isotope_data>& get_isotopes() const { return isotopes; }
  const std::vector<isotope_data> get_stable_isotopes() const {
    std::vector<isotope_data> stable_isotopes;
    for (std::vector<isotope_data>::size_type i=0; i<isotopes.size(); i++)
      if (isotopes[i].stable()) stable_isotopes.push_back(isotopes[i]);
    return stable_isotopes;
  }
  const std::vector<isotope_data> get_unstable_isotopes() const {
    std::vector<isotope_data> unstable_isotopes;
    for (std::vector<isotope_data>::size_type i=0; i<isotopes.size(); i++)
      if (!isotopes[i].stable()) unstable_isotopes.push_back(isotopes[i]);
    return unstable_isotopes;
  }

  // Return the minimum and maximum masses in our table of yields
  double min_mass() const { return mass.front(); }
  double max_mass() const { return mass.back(); }

  // Return the yield from stars of specified initial mass; the vector
  // version returns the sum of the yields from a set of stars; the
  // alternative versions return the yields only from stable or
  // unstable isotopes; the optional second argument specifies that
  // the unstable isotopes should be decreased by an amount
  // corresponding to decay for a specified time
  std::vector<double> yield(const double m,
			    const double t_decay = 0.0) const;
  std::vector<double> yield(const std::vector<double>& m,
			    const std::vector<double>& t_decay
			    = std::vector<double>()) const;
  std::vector<double> yield_stable(const double m) const;
  std::vector<double> yield_stable(const std::vector<double>& m) const;
  std::vector<double> yield_unstable(const double m,
				     const double t_decay = 0.0) const;
  std::vector<double> yield_unstable(const std::vector<double>& m,
				     const std::vector<double>& t_decay
				     = std::vector<double>()) const;

  // Return the yield of a particular isotope; we have 4 versions of
  // this function, depending on whether they accept a single scalar
  // mass or a vector of masses, and on whether the accept an isotope
  // identified by index or by isotope_data
  double yield(const double m, const std::vector<double>::size_type i,
	       const double t_decay = 0.0) const;
  double yield(const double m, const isotope_data& iso,
	       const double t_decay = 0.0) const;
  double yield(const std::vector<double>& m,
	       const std::vector<double>::size_type i,
	       const std::vector<double>& t_decay
	       = std::vector<double>()) const;
  double yield(const std::vector<double>& m,
	       const isotope_data& iso,
	       const std::vector<double>& t_decay
	       = std::vector<double>()) const;

  // Return the range of masses that produce type II SNe; this is
  // produced in the form of a vector of masses giving the lower and
  // upper limits of mass ranges that produce SNe
  const std::vector<double>& sn_mass_range() const { return sn_mass; }

  // Function to determine if a given mass star produes a type II SN
  // when it dies
  bool produces_sn(const double m) const;

  // Function to return the index within the yield table for a
  // particular isotope, specified by its atomic number and atomic
  // weight
  std::vector<double>::size_type
  isotope_index(unsigned int num, unsigned int wgt) const {
    return isotope_map1[std::make_pair(num, wgt)];
  }
   
private:

  // Data
  std::vector<double>::size_type nmass;  // Number of masses in yield table
  std::vector<double>::size_type niso;   // Number of isotopes in table
  std::vector<double>::size_type nstable;   // Number of stable isotopes
  std::vector<double>::size_type nunstable; // Number of unstable isotopes
  std::vector<double> sn_mass;           // Mass range that produces SNe
  std::vector<double> mass;              // Mass of yield table entries
  std::vector<isotope_data> isotopes;    // Isotopes
  array2d sn_yield_tab;			 // Yield table values for SNe	
  array2d wind_yield_tab;		 // Yield table for pre-SN winds
  const bool no_decay;                   // Flag to turn off radioactive decay

  // Isotope indices
  mutable std::map<isotope_data, 
		   std::vector<double>::size_type> isotope_map;
  mutable std::map<std::pair<unsigned int, unsigned int>,
		   std::vector<double>::size_type> isotope_map1;

  // Interpolation machinery
  std::vector<gsl_spline *> sn_yield, wind_yield;
  std::vector<gsl_interp_accel *> sn_yield_accel, wind_yield_accel;
};
#endif
// _slug_yields_H_
