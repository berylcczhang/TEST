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
#include "../constants.H"
#include "../slug_MPI.H"
#include "slug_yields.H"
#include <cmath>
#include <gsl/gsl_odeiv2.h>
#include <gsl/gsl_matrix.h>
#include <algorithm>
using namespace std;

////////////////////////////////////////////////////////////////////////
// decay_tree class functions
////////////////////////////////////////////////////////////////////////

// Node builder
struct decay_tree_node *
decay_tree::create_node(const isotope_data *iso,
			const struct decay_tree_node *parent,
			const double creation_rate) {

  // Allocate memory
  struct decay_tree_node *node = new struct decay_tree_node;

  // Store my isotope index
  node->idx = yld->isotope_index(iso);

  // Copy decay rates and creation rates from parent; store my own
  // creation rate
  if (parent) {
    node->creation_rates = parent->creation_rates;
    node->decay_rates = parent->decay_rates;
    node->creation_rates.push_back(creation_rate);
  }

  // Store my own decay rate  
  if (iso->stable()) node->decay_rates.push_back(0.0);
  else node->decay_rates.push_back(constants::yr / iso->ltime());

  // Recursively build children
  node->children.resize(iso->daught().size());
  for (vector<int>::size_type i=0; i<iso->daught().size(); i++) {
    node->children[i] =
      create_node(iso->daught()[i], node,
		  node->decay_rates.back() * iso->branch()[i]);
  }

  // Return pointer to newly-created node
  return node;
}

// Node de-allocator
void decay_tree::destroy_node(struct decay_tree_node *node) {

  // Recursively de-allocate children
  for (vector<int>::size_type i=0; i<node->children.size(); i++)
    destroy_node(node->children[i]);

  // De-allocate this node
  delete node;
}

// Function to compute decay products
vector<double>
decay_tree::decay_prod(const double t, const double m0) const {

  // Allocate output holder
  vector<double> m_prod(yld->get_niso());

  // Begin recursive computation
  decay_prod_node(t, m0, m_prod, root);

  // Return result
  return m_prod;
}

// Function to compute decay products for a single node
void
decay_tree::decay_prod_node(const double t, const double m0,
			    vector<double> &m_prod,
			    const struct decay_tree_node *node) const {

  // Call decay recursively on our child nodes
  for (vector<int>::size_type i=0; i<node->children.size(); i++)
    decay_prod_node(t, m0, m_prod, node->children[i]);

  // If this isotope is not in our list, we're done
  if (node->idx == constants::sentinel) return;

  // Shorthands for convenience below
  const vector<double> &rd = node->decay_rates;
  const vector<double> &rc = node->creation_rates;

  // Compute decay for this node. This code looks sort of awful, but
  // it is the solution (computed with mathematica) for n_i(t) for the
  // ODE system
  // n_0' = -rd_0 n_0
  // n_1' = rc_0 n_0 - rd_1 n_1
  // n_2' = rc_1 n_1 - rd_2 n_2
  // n_3' = rc_2 n_2 - rd_3 n_3
  // ...
  // n_i' = rc_{i-1} n_{i-1} - rd_i n_i
  // with the initial conditions
  // n_0 = 1, n_j = 0 for all j != 0
  double rcprod = 1.0;
  for (vector<int>::size_type i=0; i<rc.size(); i++) rcprod *= rc[i];
  double denom = 1.0;
  for (vector<int>::size_type i=0; i<rd.size(); i++)
    for (vector<int>::size_type j=i+1; j<rd.size(); j++)
      denom *= (rd[i] - rd[j]);
  double num = 0.0;
  int sgn = 1;
  for (vector<int>::size_type i=0; i<rd.size(); i++) {
    double exfac = exp(-t*rd[i]);
    double fac = 1.0;
    for (vector<int>::size_type j=0; j<rd.size(); j++) {
      if (i == j) continue;
      for (vector<int>::size_type k=j+1; k<rd.size(); k++) {
	if (i == k) continue;
	fac *= (rd[j] - rd[k]);
      }
    }
    num += fac * exfac * sgn;
    sgn *= -1;
  }
  m_prod[node->idx] += m0 * rcprod * abs(num/denom);
}


////////////////////////////////////////////////////////////////////////
// slug_yields functions
////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////
// Destructor
////////////////////////////////////////////////////////////////////////
slug_yields::~slug_yields() {
  if (own_isotope_table) delete iso_table;
  for (std::vector<int>::size_type i=0; i<decay_trees.size(); i++)
    delete decay_trees[i];
}

////////////////////////////////////////////////////////////////////////
// Method to initialize isotope maps, etc.
////////////////////////////////////////////////////////////////////////
void slug_yields::isotope_init() {

  // Sort the isotopes, and return numbers of stable and unstable
  std::sort(isotopes.begin(), isotopes.end(), slug_isotopes::isotope_sort);
  niso = isotopes.size();
  nstable = 0;
  for (vector<double>::size_type i=0; i<niso; i++)
    if (isotopes[i]->stable()) nstable++;
  nunstable = niso - nstable;

  // Create maps that give different ways of accessing isotopes -- via
  // a string, an isotope_data pointer, or a (Z, A) pair
  for (vector<double>::size_type i=0; i<niso; i++) {
    stringstream ss;
    ss << isotopes[i]->symbol() << isotopes[i]->wgt();
    isotope_map_str[ss.str()] = i;
    isotope_map_ptr[isotopes[i]] = i;
    isotope_map_za[make_pair(isotopes[i]->num(), isotopes[i]->wgt())] = i;
  }
  
  // Initialize decay trees
  for (vector<double>::size_type i=0; i<isotopes.size(); i++) {
    if (isotopes[i]->stable()) {
      decay_trees.push_back(nullptr);
    } else {
      decay_tree *tree = new decay_tree(isotopes[i], this);
      decay_trees.push_back(tree);
    }
  }
}


////////////////////////////////////////////////////////////////////////
// Methods to return lists of only stable or only unstable isotopes
////////////////////////////////////////////////////////////////////////
inline const vector<const isotope_data *>
slug_yields::get_stable_isotopes() const {
  vector<const isotope_data *> stable_isotopes;
  for (vector<int>::size_type i=0; i<isotopes.size(); i++)
    if (isotopes[i]->stable()) stable_isotopes.push_back(isotopes[i]);
  return stable_isotopes;
}

inline const vector<const isotope_data *>
slug_yields::get_unstable_isotopes() const {
  vector<const isotope_data *> unstable_isotopes;
  for (vector<int>::size_type i=0; i<isotopes.size(); i++)
    if (!isotopes[i]->stable()) unstable_isotopes.push_back(isotopes[i]);
  return unstable_isotopes;
}

////////////////////////////////////////////////////////////////////////
// Methods to return yields for all elements with radioactive decay
// baked in
////////////////////////////////////////////////////////////////////////

// Version for a single star
vector<double>
slug_yields::yield(const double m,
		   const double t_decay) const {

  // Output holder
  vector<double> yld(niso);
  
  // Return 0 if outside our mass range
  if ((m < mmin) || (m > mmax)) return yld;

  // Get mass-dependent yield
  yld = get_yield(m);

  // Apply radioactive decay
  if (t_decay > 0 && !no_decay) yld = decay(t_decay, yld);

  // Return
  return yld;
}

// Version for a vector of stars
vector<double>
slug_yields::yield(const vector<double> &m,
		   const vector<double> &t_decay) const {
  // Output holder
  vector<double> yld(niso);

  // Loop over input stars
  for (vector<double>::size_type i=0; i<m.size(); i++) {

    // Get yield for this star
    vector<double> yld_tmp(niso);
    if (t_decay.size() > 0) yld_tmp = yield(m[i], t_decay[i]);
    else yld_tmp = yield(m[i]);

    // Add to running total
    for (vector<double>::size_type j=0; j<niso; j++) yld[j] += yld_tmp[j];
  }

  // Return
  return yld;
}

////////////////////////////////////////////////////////////////////////
// Methods to return yield of a single isotope
////////////////////////////////////////////////////////////////////////

// Version for a single star
double slug_yields::yield(const double m,
			  const std::vector<double>::size_type i,
			  const double t_decay) const {
  
  // Return 0 if outside our mass range
  if ((m < mmin) || (m > mmax)) return 0.0;

  // Get mass-dependent yield
  double yld = get_yield(m, i);

  // Apply radioactive decay
  if (t_decay > 0 && !no_decay && !isotopes[i]->stable())
    yld *= exp(-t_decay/isotopes[i]->ltime());

  // Return
  return yld;
}

// Version for a vector of stars
double slug_yields::yield(const std::vector<double>& m,
			  const std::vector<double>::size_type i,
			  const std::vector<double>& t_decay) const {
  // Output holder
  double yld = 0.0;

  // Loop over input stars
  for (vector<double>::size_type j=0; j<m.size(); j++) {

    // Add yield for this star
    if (t_decay.size() > 0) yld += yield(m[j], i, t_decay[j]);
    else yld += yield(m[j], i);
  }

  // Return
  return yld;
}

////////////////////////////////////////////////////////////////////////
// Functions to calculate radioactive decay chains
////////////////////////////////////////////////////////////////////////

vector<double>
slug_yields::decay(const double t,
		   const vector<double> m_init) const {

  // If decay is off or t == 0, just return a copy of the original masses
  if (no_decay || t == 0.0) {
    vector<double> m_fin = m_init;
    return m_fin;
  }

  // Initialize output holder
  assert(m_init.size() == isotopes.size());
  vector<double> m_fin(isotopes.size());

  // Compute final masses using decay trees
  for (vector<double>::size_type i=0; i<m_fin.size(); i++) {
    if (isotopes[i]->stable()) {
      m_fin[i] += m_init[i];
    } else {
      vector<double> m_prod = decay_trees[i]->decay_prod(t, m_init[i]);
      for (vector<double>::size_type j=0; j<m_fin.size(); j++) {
	m_fin[j] += m_prod[j];
      }
    }
  }

  // Return the final mass array
  return m_fin;
}

