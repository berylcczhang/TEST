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
// class slug_imf_integrator
//
// This is a helper class that can integrate various stellar
// quantities over an IMF. It is used for non-stochastic stellar
// populations, and can do things (for example) like computing the
// integrated bolometric luminosity of a stellar population. Formally,
// the integrate function here computes the integral
//
// <Q>(t) = (M_tot / <m>) int (dn/dm) Q(m, t) dm
//
// where M_tot is the total mass of the stellar population, <m> is the
// mean mass for the specified IMF, dn/dm is the IMF, and Q(m, t) is
// an arbitrary function of stellar mass m and age t. By default
// integration is only over the non-stochastic part of the IMF (i.e.,
// it contains only the contributions from stars that are not being
// account for stochastically), but it can optionally be over the full
// IMF.

// Similarly, integrate_sfh, which integrates over both the IMF
// and the SFH, computes the quantity
//
// <Q> = 1/<m> int_0^t SFR(t') int (dn/dm) Q(m, t-t') dm dt'
//
// where SFR(t') is the star formation rate at age t', and t is the
// time for which the star formation has been occurring.
//
// The algorithm here closely follows the GSL Gauss-Kronrod
// integration routines.
//
// Implementation note: this is a template class. It can be defined
// for any floating point type T, OR for a std::vector of any floating
// point type. In the vector case, integration is performed, and the
// tolerance is checked, on an element-wise basis, so that the results
// are identical (to within the specified tolerance) to what they
// would have been if the scalar version had been invoked on every
// element of the vector.
////////////////////////////////////////////////////////////////////////

#ifndef _slug_imf_integrator_H_
#define _slug_imf_integrator_H_

#include <vector>
#include <type_traits>
#include <boost/function.hpp>
#include "../constants.H"
#include "../slug_IO.H"
#include "../pdfs/slug_PDF.H"
#include "../tracks/slug_tracks.H"

////////////////////////////////////////////////////////////////////////
// Data for Gauss-Kronrod quadratures. These are copied directly from
// the GSL.
////////////////////////////////////////////////////////////////////////
namespace gkdata {
  const unsigned int gk_max_iter = 500; // Maximum number of iterations
  const unsigned int gknum = 61;  // Order of rule
  const unsigned int gknum1 = 31; // Number of independent points
  const double xgk[31] =          // Abcissae for 61-pt Kronrod rule
    {
      0.999484410050490637571325895705811,
      0.996893484074649540271630050918695,
      0.991630996870404594858628366109486,
      0.983668123279747209970032581605663,
      0.973116322501126268374693868423707,
      0.960021864968307512216871025581798,
      0.944374444748559979415831324037439,
      0.926200047429274325879324277080474,
      0.905573307699907798546522558925958,
      0.882560535792052681543116462530226,
      0.857205233546061098958658510658944,
      0.829565762382768397442898119732502,
      0.799727835821839083013668942322683,
      0.767777432104826194917977340974503,
      0.733790062453226804726171131369528,
      0.697850494793315796932292388026640,
      0.660061064126626961370053668149271,
      0.620526182989242861140477556431189,
      0.579345235826361691756024932172540,
      0.536624148142019899264169793311073,
      0.492480467861778574993693061207709,
      0.447033769538089176780609900322854,
      0.400401254830394392535476211542661,
      0.352704725530878113471037207089374,
      0.304073202273625077372677107199257,
      0.254636926167889846439805129817805,
      0.204525116682309891438957671002025,
      0.153869913608583546963794672743256,
      0.102806937966737030147096751318001,
      0.051471842555317695833025213166723,
      0.000000000000000000000000000000000
    };
  const double wg[15] =    /* weights of the 30-point gauss rule */
    {
      0.007968192496166605615465883474674,
      0.018466468311090959142302131912047,
      0.028784707883323369349719179611292,
      0.038799192569627049596801936446348,
      0.048402672830594052902938140422808,
      0.057493156217619066481721689402056,
      0.065974229882180495128128515115962,
      0.073755974737705206268243850022191,
      0.080755895229420215354694938460530,
      0.086899787201082979802387530715126,
      0.092122522237786128717632707087619,
      0.096368737174644259639468626351810,
      0.099593420586795267062780282103569,
      0.101762389748405504596428952168554,
      0.102852652893558840341285636705415
    };
  const double wgk[31] =   /* weights of the 61-point kronrod rule */
    {
      0.001389013698677007624551591226760,
      0.003890461127099884051267201844516,
      0.006630703915931292173319826369750,
      0.009273279659517763428441146892024,
      0.011823015253496341742232898853251,
      0.014369729507045804812451432443580,
      0.016920889189053272627572289420322,
      0.019414141193942381173408951050128,
      0.021828035821609192297167485738339,
      0.024191162078080601365686370725232,
      0.026509954882333101610601709335075,
      0.028754048765041292843978785354334,
      0.030907257562387762472884252943092,
      0.032981447057483726031814191016854,
      0.034979338028060024137499670731468,
      0.036882364651821229223911065617136,
      0.038678945624727592950348651532281,
      0.040374538951535959111995279752468,
      0.041969810215164246147147541285970,
      0.043452539701356069316831728117073,
      0.044814800133162663192355551616723,
      0.046059238271006988116271735559374,
      0.047185546569299153945261478181099,
      0.048185861757087129140779492298305,
      0.049055434555029778887528165367238,
      0.049795683427074206357811569379942,
      0.050405921402782346840893085653585,
      0.050881795898749606492297473049805,
      0.051221547849258772170656282604944,
      0.051426128537459025933862879215781,
      0.051494729429451567558340433647099
    };
}


////////////////////////////////////////////////////////////////////////
// Here we define a trivial little class whose sole job is to perform
// actions needed in the slug_imf_integrator class where the action we
// need to perform depends on whether the class argument is a
// vector or a scalar. We have to keep this out of the main
// slug_imf_integrator class, because we can't partially specialize in
// c++.
////////////////////////////////////////////////////////////////////////

// Version for scalar types
template<typename T>
class slug_imf_integrator_helper {
public:
  // Get size of vector
  typename std::vector<T>::size_type size(T &a) const { return 0; }
  // Initialize
  T init(typename std::vector<T>::size_type nvec = 0) 
    const { return 0.0; }
  // Do +=, *=, *, abs(a-b)
  void plusequal(T &a, const T &b) const { a += b; }
  void timesequal(T &a, const double &b) const { a *= b; }
  T times(const T &v1, const double &s1) const { return s1*v1; }  
  T absdiff(const T &a, const T &b) const { return std::abs(a-b); }
  // Do update needed in Gauss-Kronrod quadarture
  void gkupdate(T &lhs, const T &v1, const T &v2, const double &s1, 
		const double &s2, const double &fac) const
  { lhs += fac * (v1*s1 + v2*s2); }
  // Return abs(err/val)
  T rel_err(const T &err, const T &val) const 
  { return std::abs(err)/(std::abs(val)+constants::small); }
  // Update sum and error based on new values for an interval
  void update(T &val, const T &val_l, const T &val_r, 
	      const T &val_old) const {
    val += val_l + val_r - val_old;
  }
};

// Version for vector types
template<typename T, typename A>
class slug_imf_integrator_helper<std::vector<T, A>> {
public:
  // Get size of vector
  typename std::vector<T, A>::size_type size(std::vector<T, A> &a) 
  const { return a.size(); }
  // Initialize vector of zeros
  std::vector<T, A> init(typename std::vector<T, A>::size_type nvec = 0) 
    const { std::vector<T, A> vec(nvec, 0.0); return vec; }
  // Do +=, *=, *, abs(a-b) elementwise
  void plusequal(std::vector<T, A> &a, const std::vector<T, A> &b) const { 
    for (typename std::vector<T, A>::size_type i=0; i<a.size(); i++) 
      a[i] += b[i]; 
  }
  void timesequal(std::vector<T, A> &a, const double &b) const { 
    for (typename std::vector<T, A>::size_type i=0; i<a.size(); i++) 
      a[i] *= b; 
  }
  std::vector<T, A> times(const std::vector<T, A> &v1, 
			  const double &s1) const {
    std::vector<T, A> vec(v1.size(), 0.0);
    for (typename std::vector<T, A>::size_type i=0; i<v1.size(); i++) 
      vec[i] = s1*v1[i];
    return vec;
  }
  std::vector<T, A> absdiff(const std::vector<T, A> &a, 
			    const std::vector<T, A> &b) const {
    std::vector<T, A> vec(a.size(), 0.0);
    for (typename std::vector<T, A>::size_type i=0; i<a.size(); i++) 
      vec[i] = std::abs(a[i]-b[i]);
    return vec;
  }
  // Do update needed in Gauss-Kronrod quadarture
  void gkupdate(std::vector<T, A> &lhs, const std::vector<T, A> &v1,
		const std::vector<T, A> &v2, const double &s1, 
		const double &s2, const double &fac) const {
    for (typename std::vector<T, A>::size_type i=0; i<lhs.size(); i++) 
      lhs[i] += fac * (s1*v1[i] + s2*v2[i]);
  }
  // Return element-wise max of abs(err/val)
  double rel_err(const std::vector<T, A> &err, 
		 const std::vector<T, A> &val) const {
    double re = 0.0;
    for (typename std::vector<T, A>::size_type i=0; i<err.size(); i++) 
      re = std::max(re, std::abs(err[i])/(std::abs(val[i])+constants::small));
    return re;
  }
  void update(std::vector<T, A> &val, const std::vector<T, A> &val_l,
	      const std::vector<T, A> &val_r, 
	      const std::vector<T, A> &val_old) const {
    for (typename std::vector<T, A>::size_type i=0; i<val.size(); i++) 
      val[i] += val_l[i] + val_r[i] - val_old[i];
  }
};

////////////////////////////////////////////////////////////////////////
// Main IMF integrator class
////////////////////////////////////////////////////////////////////////

// Scalar case
template <typename T>
class slug_imf_integrator {

public:

  // Constructor; this sets some default options
  slug_imf_integrator(const slug_tracks *tracks_, 
		      const slug_PDF *imf_,
		      const slug_PDF *sfh_,
		      slug_ostreams &ostreams_,
		      typename std::vector<T>::size_type nvec_ = 0,
		      const double tol_ = 1.0e-3,
		      const bool include_stoch_ = false) :
    ostreams(ostreams_), tracks(tracks_), imf(imf_), sfh(sfh_),
    help(), nvec(nvec_), tol(tol_), include_stoch(include_stoch_) { }

  // Destructor
  ~slug_imf_integrator() { }

  // Routines to get and set various options:
  // -- nvec specifies the length of the vector returned for vector
  //    functions
  // -- tol specifies the relative accuracy tolerance of the
  //    integration
  // -- include_stoch specifies whether the integration should just be
  //    over the non-stochastic portion of the IMF (if false) or the
  //    entire IMF (if true)
  typename std::vector<T>::size_type get_nvec() const { return nvec; }
  void set_nvec(typename std::vector<T>::size_type nvec_) { nvec = nvec_; }
  double get_tol() const { return tol; }
  void set_tol(double tol_) { tol = tol_; }
  bool get_include_stoch() const { return include_stoch; }
  void set_include_stoch(bool include_stoch_) 
  { include_stoch = include_stoch_; }

  // Routine to integrate a specified quantity over IMF; the
  // integrate routine takes a total mass for the population, an age,
  // and a function that takes stardata as input and returns the desired
  // quantity. Note that, because this routine makes
  // use of the tracks, it treats the integrand as zero for any parts
  // of the IMF that are outside the mass range that the tracks
  // cover.
  T integrate(const double m_tot, const double age,
	      boost::function<T(const slug_stardata &)> func_in = 0) const;

  // Version of the integrate function that works over a limited mass range
  T integrate_lim(const double m_tot, const double age,
		  const double m_min, const double m_max,
		  boost::function<T(const slug_stardata &)>
		  func_in = 0) const;
  
  
  // These two functions are simlar to the integrate function, but
  // they take as arguments functions that just take the mass, or take
  // just the mass and age (or just mass), as inputs, rather than
  // taking stardata as input. Because they do not need stardata,
  // these functions do not call the track routines to generate 
  // isochrones (hence the _nt extension), and so can integrate over
  // mass ranges that the tracks do not cover.
  T integrate_nt(const double m_tot, const double age,
		 boost::function<T(const double &, const double &)> 
		 func_ = 0) const;
  T integrate_nt(const double m_tot,
		 boost::function<T(const double &)> func_ = 0) const;

  // Versions of the previous two functions that accept a limited mass range
  T integrate_nt_lim(const double m_tot, const double age,
		     const double m_min, const double m_max,
		     boost::function<T(const double &, const double &)> 
		     func_ = 0) const;
  T integrate_nt_lim(const double m_tot,
		     const double m_min, const double m_max,
		     boost::function<T(const double &)> func_ = 0) const;
  

  // Routines to integrate a specified quantity over both the IMF and
  // the star formation history. This has two of the three forms that
  // appear for the integrate routines. It is missing the form for a
  // function that only takes mass (as opposed to mass and age) as an
  // argument, because in that case the integration is trivial, since
  // the mass and time parts are separable.
  T integrate_sfh(const double t,
		  boost::function<T(const slug_stardata &)> func_) const;
  T integrate_sfh_nt(const double t,
		     boost::function<T(const double &, const double &)> 
		     func_) const;

protected:
  // IO handler
  slug_ostreams& ostreams;

private:

  // Helper function to do integral over a particular mass range; the
  // IMF must be defined over this range
  T integrate_range(const double m_tot, const double age,
		    const double m_min, const double m_max) const;

  // Helper function to do Gauss-Kronrod integration on a single mass
  // interval
  void integrate_gk(const double m_min, const double m_max,
		    const double age, 
		    T &result, T &err) const;
  void integrate_sfh_gk(const double t_min, const double t_max,
			const double t, T &result, T &err) const;

  // Data on the tracks, IMF, SFH
  const slug_tracks *tracks;    // Tracks
  const slug_PDF *imf;          // IMF
  const slug_PDF *sfh;          // Star formation history

  // Helper class to handle stuff that depends on whether we're
  // operating on scalar or vector data
  const slug_imf_integrator_helper<T> help;

  // Variables that control internal behavior
  mutable boost::function<T(const slug_stardata &)> *func = 0;
  mutable boost::function<T(const double &)> *func_mass = 0;
  mutable boost::function<T(const double &, const double &)> *func_mass_age = 0;
  mutable typename std::vector<T>::size_type nvec;
  mutable double tol;
  bool include_stoch;
};

#endif // _slug_imf_integrator_H_
