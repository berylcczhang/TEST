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
// This module just provides little utility class that just provides a
// mechanism for storing information about isotopes and handling
// them. We also store some atomic data.
//
////////////////////////////////////////////////////////////////////////

#ifndef _slug_isotopes_H_
#define _slug_isotopes_H_

#include <limits>
#include <string>

// Information on elements
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

// Information on isotopes
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

// Namespace where we put some atomic data
namespace atomic_data {

  // Element symbols
  const std::vector<std::string> periodic_table = {
    "h", "he", "li", "be", "b", "c", "n", "o", "f", "ne", "na", "mg", "al",
    "si", "p", "s", "cl", "ar", "k", "ca", "sc", "ti", "v", "cr", "mn", "fe", 
    "co", "ni", "cu", "zn", "ga", "ge", "as", "se", "br", "kr", "rb", "sr",
    "y", "zr", "nb", "mo", "tc", "ru", "rh", "pd", "ag", "cd", "in", "sn", 
    "sb", "te", "i", "xe", "cs", "ba", "la", "ce", "pr", "nd", "pm", "sm",
    "eu", "gd", "tb", "dy", "ho", "er", "tm", "yb", "lu", "hf", "ta", "w",
    "re", "os", "ir", "pt", "au", "hg", "tl", "pb", "bi", "po", "at", "rn",
    "fr", "ra", "ac", "th", "pa", "u", "np", "pu", "am", "cm", "bk", "cf",
    "es", "fm", "md", "no", "lr", "rf" "db", "sg", "bh", "hs", "mt", "ds",
    "rg", "cn", "uut", "fl", "uup", "lv", "uus", "uuo" };

  // Lifetimes of some astrophysically-important radioactive isotopes
  const isotope_data c14("c", 6, 14, 1.81e11*log(2.0));
  const isotope_data na22("na", 11, 22, 8.21e7*log(2.0));
  const isotope_data al26("al", 13, 26, 2.26e13*log(2.0));
  const isotope_data si32("si", 14, 32, 5.4e9*log(2.0));
  const isotope_data cl36("cl", 17, 36, 9.5e12*log(2.0));
  const isotope_data ar39("ar", 18, 39, 8.5e9*log(2.0));
  const isotope_data k40("k", 19, 40, 4.03e16*log(2.0));
  const isotope_data ca41("ca", 20, 41, 3.3e12*log(2.0));
  const isotope_data ca45("ca", 20, 45, 1.41e7*log(2.0));
  const isotope_data ti44("ti", 22, 44, 2.0e9*log(2.0));
  const isotope_data v49("v", 23, 49, 2.9e7*log(2.0));
  const isotope_data mn53("mn", 25, 53, 1.2e13*log(2.0));
  const isotope_data mn54("mn", 25, 54, 2.70e7*log(2.0));
  const isotope_data fe55("fe", 26, 55, 8.64e7*log(2.0));
  const isotope_data fe60("fe", 26, 60, 4.7e13*log(2.0));
  const isotope_data co60("co", 27, 60, 6.08e10*log(2.0));
  const isotope_data ni56("ni", 28, 56, 1.87e7*log(2.0));
  const isotope_data ni57("ni", 28, 57, 3.08e6*log(2.0));
  const isotope_data ni59("ni", 28, 59, 2.41e12*log(2.0));
  const isotope_data ni63("ni", 28, 63, 3.16e9*log(2.0));
  const std::vector<const isotope_data> unstable_isotopes =
    { c14, na22, al26, si32, cl36, ar39, k40, ca41, ca45, ti44, v49, mn53, 
      mn54, fe55, fe60, co60, ni56, ni57, ni59, ni63 };
}

#endif
// _slug_isotopes_H_
