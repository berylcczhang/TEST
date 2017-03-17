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
// This class defines a generic interface to classes that do yield
// computations. It is pure virtual, since it just defines the
// interface.
//
////////////////////////////////////////////////////////////////////////

#ifndef _slug_yields_H_
#define _slug_yields_H_

#include "slug.H"
#include "slug_isotopes.H"
#include <map>
#include <utility>
#include <vector>

class slug_yields {

public:

  // Constructor
  slug_yields(const double metallicity_, const bool no_decay_ = false
#ifdef ENABLE_MPI
	      , const int rank_ = 0
#endif
	      )
    : no_decay(no_decay_), metallicity(metallicity_)
#ifdef ENABLE_MPI
    , rank(rank_)
#endif
  {};

  // Destructor
  virtual ~slug_yields() {};

  // Methods to access data on the isotopes we know about
  std::vector<isotope_data>::size_type get_niso() const { return niso; }
  const std::vector<isotope_data>& get_isotopes() const { return isotopes; }
  const std::vector<isotope_data> get_stable_isotopes() const;
  const std::vector<isotope_data> get_unstable_isotopes() const;

  // Return mass range from which yields are available
  double min_mass() const { return mmin; }
  double max_mass() const { return mmax; }

  // Return the yield from stars of specified initial mass, or from a
  // vector of stars
  std::vector<double> yield(const double m,
			    const double t_decay = 0.0) const;
  std::vector<double> yield(const std::vector<double>& m,
			    const std::vector<double>& t_decay
			    = std::vector<double>()) const;

  // Return the yield of a particular isotope, from a single star or a
  // vector of stars
  double yield(const double m,
	       const std::vector<double>::size_type i,
	       const double t_decay = 0.0) const;
  double yield(const std::vector<double>& m,
	       const std::vector<double>::size_type i,
	       const std::vector<double>& t_decay
	       = std::vector<double>()) const;

  // As a convenience for users, we provide an alternate interface to
  // the yield function where users can provide either isotope_data or
  // an atomic number and weight to specify the isotope whose yield
  // they want. These just wrap the yield functions above.
  double yield(const double m, const isotope_data& iso,
	       const double t_decay = 0.0) const {
    return yield(m, isotope_map.at(iso), t_decay);
  }
  double yield(const double m, const unsigned int num,
	       const unsigned int wgt,
	       const double t_decay = 0.0) const {
    return yield(m, isotope_map_za.at(std::make_pair(num, wgt)), t_decay);
  }
  double yield(const std::vector<double>& m,
	       const isotope_data& iso,
	       const std::vector<double>& t_decay
	       = std::vector<double>()) const {
    return yield(m, isotope_map.at(iso), t_decay);
  }
  double yield(const std::vector<double>& m,
	       const unsigned int num, const unsigned int wgt,
	       const std::vector<double>& t_decay
	       = std::vector<double>()) const {
    return yield(m, isotope_map_za.at(std::make_pair(num, wgt)), t_decay);
  }
  
  // Another convenience function for users; this one returns the
  // mapping between atomic number/weight and index, or between
  // isotope_data and index
  std::vector<double>::size_type
  isotope_index(const isotope_data& iso) const {
    return isotope_map.at(iso);
  }
  std::vector<double>::size_type
  isotope_index(unsigned int num, unsigned int wgt) const {
    return isotope_map_za.at(std::make_pair(num, wgt));
  }

  // Virtual methods to be implemented by derived classes. These
  // methods actually do the basic yield calculation, for all isotopes
  // or for a single isotope.
  virtual std::vector<double> get_yield(const double m) const = 0;
  virtual double get_yield(const double m,
			   const std::vector<double>::size_type i) const = 0;

  // Methods to return information on SNII and AGB stars. These are
  // defined here as returning no mass range and false, but are
  // declared virtual so they can be overridden by derived classes.
  virtual const std::vector<double>& sn_mass_range() const
  { return null_mass_range; }
  virtual const std::vector<double>& agb_mass_range() const
  { return null_mass_range; }
  virtual bool produces_sn(const double m) const { return false; }
  virtual bool produces_agb(const double m) const { return false; }

  // Flag to disable decay of unstable isotopes
  bool no_decay;

protected:

  // Data on isotopes; these are to be initialized by the derived class
  std::vector<double>::size_type niso;      // Number of isotopes
  std::vector<double>::size_type nstable;   // Number of stable isotopes
  std::vector<double>::size_type nunstable; // Number of unstable isotopes
  std::vector<isotope_data> isotopes;       // List of isotopes

  // Mappings between isotope data, or between atomic number and
  // weight, and index in the isotope table. Again, these are to be
  // instantiated by the derived class. They are declared as mutable
  // so that we can use the [] operator on them even if they are
  // attached to a const instance of the derived class.
  std::map<isotope_data, 
	   std::vector<double>::size_type> isotope_map;
  std::map<std::pair<unsigned int, unsigned int>,
	   std::vector<double>::size_type> isotope_map_za;

  // Mass range covered; to be set by derived class
  double mmin, mmax;

  // Metallicity
  double metallicity;

  // Null mass range, delared so we can return const pointers to an
  // empty vector
  std::vector<double> null_mass_range;

  // MPI rank
#ifdef ENABLE_MPI
  int rank;
#endif
};

#endif
// _slug_yields_H_
