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

#include "slug_isotopes.H"
#include "../slug.H"
#include "../slug_IO.H"
#include <map>
#include <utility>
#include <vector>

// Enum of yield models
enum yieldMode { SNII_SUKHBOLD16, AGB_KARAKAS16_DOHERTY14,
		 SNII_SUKHBOLD16__AGB_KARAKAS16_DOHERTY14 };

////////////////////////////////////////////////////////////////////////
// A helper class that stores decay trees
////////////////////////////////////////////////////////////////////////

struct decay_tree_node {
  std::vector<double>::size_type idx;
  std::vector<struct decay_tree_node *> children;
  std::vector<double> creation_rates, decay_rates;
};

class slug_yield;

class decay_tree {

public:

  // Constructor and destructor
  decay_tree(const isotope_data *iso, const slug_yields *yld_) : yld(yld_)
  { root = create_node(iso, nullptr, 0.0); }
  ~decay_tree() { destroy_node(root); }

  // Function to evaluate tree and return amount of each nuclide
  // created or destroyed
  std::vector<double> decay_prod(const double t, const double m0) const;

private:

  // Functions to build and destroy nodes recursively
  struct decay_tree_node
  *create_node(const isotope_data *iso,
	       const struct decay_tree_node *parent,
	       const double creation_rate);
  void destroy_node(struct decay_tree_node *node);

  // Function to evaluate decay product on a node
  void decay_prod_node(const double t, const double m0,
		       std::vector<double> &m_prod,
		       const struct decay_tree_node *node) const;

  // Pointer to tree, and to parent yield object
  struct decay_tree_node *root;
  const slug_yields *yld;
};


////////////////////////////////////////////////////////////////////////
// The main slug_yields class
////////////////////////////////////////////////////////////////////////

class slug_yields {

public:

  // Constructor; we provide two implementations, one of which
  // initializes a new isotope table, and one of which gets a pointer
  // to an existing one
  slug_yields(const double metallicity_,
	      const char *iso_data_dir,
	      slug_ostreams &ostreams_,
	      const bool no_decay_ = false)
    : no_decay(no_decay_),
      ostreams(ostreams_),
      iso_table(new isotope_table(iso_data_dir, ostreams_)),
      own_isotope_table(true),
      metallicity(metallicity_)
  { };
  
  slug_yields(const double metallicity_,
	      const isotope_table *iso_table_,
	      slug_ostreams &ostreams_,
	      const bool no_decay_ = false)
    : no_decay(no_decay_),
      ostreams(ostreams_),
      iso_table(iso_table_),
      own_isotope_table(false),
      metallicity(metallicity_)
  { };

  // Destructor
  virtual ~slug_yields();

  // Methods to access data on the isotopes we know about
  std::vector<int>::size_type get_niso() const { return niso; }
  const std::vector<const isotope_data *>& get_isotopes() const
  { return isotopes; }
  const std::vector<const isotope_data *> get_stable_isotopes() const;
  const std::vector<const isotope_data *> get_unstable_isotopes() const;

  // Return mass range from which yields are available
  double min_mass() const { return mmin; }
  double max_mass() const { return mmax; }
  bool produces_yield(const double m) const {
    return (m >= mmin) && (m <= mmax);
  }

  // Return the yield from stars of specified initial mass, or from a
  // vector of stars
  std::vector<double> yield(const double m,
			    const double t_decay = 0.0) const;
  std::vector<double> yield(const std::vector<double>& m,
			    const std::vector<double>& t_decay
			    = std::vector<double>()) const;

  // Return the yield of a particular isotope, from a single star or a
  // vector of stars; note that these versions of yield do NOT compute
  // the full decay chain, and thus should only be used for isotopes
  // where the yield of that isotope contains non-negligible
  // contributions from the decay of other isotopes
  double yield(const double m,
	       const std::vector<double>::size_type i,
	       const double t_decay = 0.0) const;
  double yield(const std::vector<double>& m,
	       const std::vector<double>::size_type i,
	       const std::vector<double>& t_decay
	       = std::vector<double>()) const;

  // As a convenience for users, we provide an alternate interface to
  // the yield functions where the user provides either a pointer to
  // isotope_data, a string, or a (Z,A) pair
  double yield(const double m, const std::string &iso,
	       const double t_decay = 0.0) const {
    return yield(m, isotope_map_str.at(iso), t_decay);
  }
  double yield(const double m, const isotope_data *iso,
	       const double t_decay = 0.0) const {
    return yield(m, isotope_map_ptr.at(iso), t_decay);
  }
  double yield(const double m, const unsigned int num,
	       const unsigned int wgt,
	       const double t_decay = 0.0) const {
    return yield(m, isotope_map_za.at(std::make_pair(num, wgt)), t_decay);
  }
  double yield(const std::vector<double>& m,
	       const std::string &iso,
	       const std::vector<double>& t_decay
	       = std::vector<double>()) const {
    return yield(m, isotope_map_str.at(iso), t_decay);
  }
  double yield(const std::vector<double>& m,
	       const isotope_data *iso,
	       const std::vector<double>& t_decay
	       = std::vector<double>()) const {
    return yield(m, isotope_map_ptr.at(iso), t_decay);
  }
  double yield(const std::vector<double>& m,
	       const unsigned int num, const unsigned int wgt,
	       const std::vector<double>& t_decay
	       = std::vector<double>()) const {
    return yield(m, isotope_map_za.at(std::make_pair(num, wgt)), t_decay);
  }
  
  // Another convenience function for users; this one returns the
  // mapping between atomic number/weight, isotope_data ptr, or
  // string, and index
  std::vector<double>::size_type
  isotope_index(const std::string &iso) const {
    if (isotope_map_str.count(iso))
      return isotope_map_str.at(iso);
    else
      return std::numeric_limits<std::vector<double>::size_type>::max();
  }
  std::vector<double>::size_type
  isotope_index(const isotope_data *iso) const {
    if (isotope_map_ptr.count(iso))
      return isotope_map_ptr.at(iso);
    else
      return std::numeric_limits<std::vector<double>::size_type>::max();
  }
  std::vector<double>::size_type
  isotope_index(unsigned int num, unsigned int wgt) const {
    if (isotope_map_za.count(std::make_pair(num, wgt)))
      return isotope_map_za.at(std::make_pair(num, wgt));
    else
      return std::numeric_limits<std::vector<double>::size_type>::max();
  }

  // Function to compute changes in abundances due to radioactive
  // decay; m_init is a vector of initial masses for all isotopes
  // (must have the same size as isotopes), and returned vector is
  // masses after a time t passes
  std::vector<double> decay(const double t,
			    const std::vector<double> m_init) const;

  // Function to traverse a decay tree for a particular nuclide and
  // return the mass of decay products produced; t is time over which
  // decay occurs, m0 is initial mass of starting isotope, iso is
  // isotope data, return value is mass of all isotopes produced by
  // the decay, indexed per the stored list of isotopes
  std::vector<double> decay_prod(const double t, const double m0,
				 const isotope_data *iso) const;
  
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

  // This method does final initialization. It must be invokved by the
  // derived class's constructor. It records the number of stable and
  // unstable isotopes, sorts the isotope list, and builds all the
  // decay trees needed.
  void isotope_init();

  // IO handler
  slug_ostreams &ostreams;

  // A table of isotope properties
  const isotope_table *iso_table;
  bool own_isotope_table;

  // Data on isotopes; these are to be initialized by the derived class
  std::vector<double>::size_type niso;      // Number of isotopes
  std::vector<double>::size_type nstable;   // Number of stable isotopes
  std::vector<double>::size_type nunstable; // Number of unstable isotopes
  std::vector<const isotope_data *> isotopes; // List of isotopes

  // Decay trees for unstable isotopes
  std::vector<class decay_tree *> decay_trees;

  // Mappings between string, isotope_data pointer, or (Z,A) pair and
  // an index in the isotope table; these are created by the
  // isotope_init function.
  std::map<const std::string,
	   std::vector<double>::size_type> isotope_map_str;
  std::map<const isotope_data *, 
	   std::vector<double>::size_type> isotope_map_ptr;
  std::map<std::pair<unsigned int, unsigned int>,
	   std::vector<double>::size_type> isotope_map_za;

  // Mass range covered; to be set by derived class
  double mmin, mmax;

  // Metallicity
  double metallicity;

  // Null mass range, delared so we can return const pointers to an
  // empty vector
  std::vector<double> null_mass_range;
  
};

#endif
// _slug_yields_H_
