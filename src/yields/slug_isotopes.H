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
// This module provides classes for storing and processing nuclear
// data -- tables of stable and unstable isotopes, lifetimes, etc.
//
////////////////////////////////////////////////////////////////////////

#ifndef _slug_isotopes_H_
#define _slug_isotopes_H_

#include <limits>
#include <string>
#include <vector>
#include <map>
#include <cmath>
#include "../slug_IO.H"
#include <boost/random/mersenne_twister.hpp>

////////////////////////////////////////////////////////////////////////
// The periodic table
////////////////////////////////////////////////////////////////////////
namespace periodic_table {
  // Element symbols
  const std::vector<std::string> symbols = {
    "h", "he", "li", "be", "b", "c", "n", "o", "f", "ne", "na", "mg", "al",
    "si", "p", "s", "cl", "ar", "k", "ca", "sc", "ti", "v", "cr", "mn", "fe", 
    "co", "ni", "cu", "zn", "ga", "ge", "as", "se", "br", "kr", "rb", "sr",
    "y", "zr", "nb", "mo", "tc", "ru", "rh", "pd", "ag", "cd", "in", "sn", 
    "sb", "te", "i", "xe", "cs", "ba", "la", "ce", "pr", "nd", "pm", "sm",
    "eu", "gd", "tb", "dy", "ho", "er", "tm", "yb", "lu", "hf", "ta", "w",
    "re", "os", "ir", "pt", "au", "hg", "tl", "pb", "bi", "po", "at", "rn",
    "fr", "ra", "ac", "th", "pa", "u", "np", "pu", "am", "cm", "bk", "cf",
    "es", "fm", "md", "no", "lr", "rf", "db", "sg", "bh", "hs", "mt", "ds",
    "rg", "cn", "uut", "fl", "uup", "lv", "uus", "uuo" };
}

////////////////////////////////////////////////////////////////////////
// A class that holds element symbol and number
////////////////////////////////////////////////////////////////////////
class element_data {
public:

  // Basic data
  element_data(std::string name_, unsigned int number_) : 
    name(name_), number(number_) {};
  const std::string& symbol() const { return name; }
  unsigned int num() const { return number; }

  // Define comparison operators to facilitate sorting
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

////////////////////////////////////////////////////////////////////////
// A class that holds information on isotopes; this class holds atomic
// symbol, atomic number, atomic mass number, lifetime, list of
// daughter nuclei, list of daughter stoichiometric factors
////////////////////////////////////////////////////////////////////////
class isotope_data : public element_data {

  friend class isotope_table;
  
public:

  // Constructors
  isotope_data(std::string name_, unsigned int number_, 
	       unsigned int weight_, 
	       double lifetime_ = -1.0,
	       const std::vector<unsigned int> &daughters_Z_
	       = std::vector<unsigned int>(),
	       const std::vector<unsigned int> &daughters_A_
	       = std::vector<unsigned int>(),
	       const std::vector<double> &daughters_branch_
	       = std::vector<double>()) :
    element_data(name_, number_), weight(weight_), lifetime(lifetime_),
    daughters_Z(daughters_Z_), daughters_A(daughters_A_),
    daughters_branch(daughters_branch_),
    daughters(daughters_Z_.size()) {};

  // Access to basic data
  unsigned int wgt() const { return weight; }
  bool stable() const { return lifetime < 0; }
  double ltime() const
  { return lifetime > 0 ? lifetime : std::numeric_limits<double>::infinity(); }
  double halflife() const { return ltime()/log(2.0); }
  const std::vector<const isotope_data *>& daught() const { return daughters; }
  const std::vector<double>& branch() const { return daughters_branch; }
  const std::vector<unsigned int>& daught_Z() const { return daughters_Z; }
  const std::vector<unsigned int>& daught_A() const { return daughters_A; }

  // Comparison operators to facilitate sorting
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

  // Private data
  const unsigned int weight;
  const double lifetime;
  const std::vector<unsigned int> daughters_Z, daughters_A;
  const std::vector<double> daughters_branch;

  // Pointers to daughters; this is marked as mutable due to a
  // sequencing issue: we won't have constructed all isotopes first,
  // before we can provide pointers to daughters, so we need to be
  // able to go back and modify the list of daughter pointers after
  // contsruction
  mutable std::vector<const isotope_data *> daughters;
};


////////////////////////////////////////////////////////////////////////
// A little utility method for sorting pointers to isotope_data
////////////////////////////////////////////////////////////////////////
namespace slug_isotopes {
  bool isotope_sort(const isotope_data *iso1,
		    const isotope_data *iso2);
}

////////////////////////////////////////////////////////////////////////
// A class that holds a table of isotopes. It reads the isotope data
// from a data file on construction, and has a variety of methods to
// expose that data to other classes.
////////////////////////////////////////////////////////////////////////

class isotope_table {

public:

  // Constructor
  isotope_table(const char *data_dir, slug_ostreams &ostreams);

  // Destructor
  ~isotope_table();

  // Methods to return lists of isotopes
  const std::vector<const isotope_data *>& all_isotopes() const
  { return all_isotopes_; };
  const std::vector<const isotope_data *>& stable_isotopes() const
  { return stable_isotopes_; };
  const std::vector<const isotope_data *>& unstable_isotopes() const
  { return unstable_isotopes_; };
  
  // Methods to return data on a single isotope. The methods and the
  // data they return are:
  // data -- returns full isotope_data
  // lfime -- returns lifetime
  // halflife -- returns half-life
  // stable -- returns whether isotope is stable
  // daught -- returns vector of pointers to daughter nuclei
  // daught_Z -- returns vector of daughter nuclei Z values
  // daught_A -- returns vector of daughter nuclei A values
  // branch -- returns vector of branching ratios for
  // creation of daughter nuclei; note that these may not add up to 1,
  // since a single decay channel may create multiple daughter nuclei,
  // and we have omitted certain decay channels where the nuclear data
  // are unavailable, or the number of possible branches is too large
  // to follow
  // Each of these methods comes in two versions, one that takes a
  // number and a weight, and one that takes a symbol of the form
  // symbol + number, e.g. c14 or fe60
  const isotope_data* data(unsigned int num, unsigned int wgt) const
  { return isotope_map_za.at(std::make_pair(num, wgt)); }
  const isotope_data* data(const std::string& symbol) const
  { return isotope_map_symbol.at(symbol); }
  double ltime(unsigned int num, unsigned int wgt) const
  { return isotope_map_za.at(std::make_pair(num, wgt))->ltime(); }
  double ltime(const std::string& symbol) const
  { return isotope_map_symbol.at(symbol)->ltime(); }
  double halflife(unsigned int num, unsigned int wgt) const
  { return isotope_map_za.at(std::make_pair(num, wgt))->halflife(); }
  double halflife(const std::string& symbol) const
  { return isotope_map_symbol.at(symbol)->halflife(); }
  bool stable(unsigned int num, unsigned int wgt) const
  { return isotope_map_za.at(std::make_pair(num, wgt))->stable(); }
  bool stable(const std::string& symbol) const
  { return isotope_map_symbol.at(symbol)->stable(); }
  const std::vector<const isotope_data *>& daught
  (unsigned int num, unsigned int wgt) const
  { return isotope_map_za.at(std::make_pair(num, wgt))->daught(); }
  const std::vector<const isotope_data *>& daught
  (const std::string& symbol) const
  { return isotope_map_symbol.at(symbol)->daught(); }
  const std::vector<double>& branch
  (unsigned int num, unsigned int wgt) const
  { return isotope_map_za.at(std::make_pair(num, wgt))->branch(); }
  const std::vector<double>& branch
  (const std::string& symbol) const
  { return isotope_map_symbol.at(symbol)->branch(); }
  const std::vector<unsigned int>& daught_Z
  (unsigned int num, unsigned int wgt) const
  { return isotope_map_za.at(std::make_pair(num, wgt))->daught_Z(); }
  const std::vector<unsigned int>& daught_Z
  (const std::string& symbol) const
  { return isotope_map_symbol.at(symbol)->daught_Z(); }
  const std::vector<unsigned int>& daught_A
  (unsigned int num, unsigned int wgt) const
  { return isotope_map_za.at(std::make_pair(num, wgt))->daught_A(); }
  const std::vector<unsigned int>& daught_A
  (const std::string& symbol) const
  { return isotope_map_symbol.at(symbol)->daught_A(); }

private:
  // Data
  std::vector<const isotope_data *> all_isotopes_;
  std::vector<const isotope_data *> stable_isotopes_;
  std::vector<const isotope_data *> unstable_isotopes_;
  std::map<std::string, const isotope_data *> isotope_map_symbol;
  std::map<std::pair<unsigned int, unsigned int>,
	   const isotope_data *> isotope_map_za;
};

#endif
// _slug_isotopes_H_
