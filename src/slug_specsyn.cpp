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

#include <algorithm>
#include <cassert>
#include <cmath>
#include <iterator>
#include <iostream>
#include "constants.H"
#include "slug_specsyn.H"

using namespace std;

#define GK_MAX_ITER 500

////////////////////////////////////////////////////////////////////////
// Gauss-Kronrod abcissae and weights, copied directly from GSL. Used
// for the Gauss-Kronrod integrations done in this class.
////////////////////////////////////////////////////////////////////////

#if 1
// Order of kronrod rule, and number of independent points
static unsigned int gknum = 15;
static unsigned int gknum1 = 8;

static const double xgk[8] =    /* abscissae of the 15-point kronrod rule */
{
  0.991455371120812639206854697526329,
  0.949107912342758524526189684047851,
  0.864864423359769072789712788640926,
  0.741531185599394439863864773280788,
  0.586087235467691130294144838258730,
  0.405845151377397166906606412076961,
  0.207784955007898467600689403773245,
  0.000000000000000000000000000000000
};

/* xgk[1], xgk[3], ... abscissae of the 7-point gauss rule. 
   xgk[0], xgk[2], ... abscissae to optimally extend the 7-point gauss rule */

static const double wg[4] =     /* weights of the 7-point gauss rule */
{
  0.129484966168869693270611432679082,
  0.279705391489276667901467771423780,
  0.381830050505118944950369775488975,
  0.417959183673469387755102040816327
};

static const double wgk[8] =    /* weights of the 15-point kronrod rule */
{
  0.022935322010529224963732008058970,
  0.063092092629978553290700663189204,
  0.104790010322250183839876322541518,
  0.140653259715525918745189590510238,
  0.169004726639267902826583426598550,
  0.190350578064785409913256402421014,
  0.204432940075298892414161999234649,
  0.209482141084727828012999174891714
};
#endif

#if 0
// Order of kronrod rule, and number of independent points
static unsigned int gknum = 61;
static unsigned int gknum1 = 31;

static const double xgk[31] =   /* abscissae of the 61-point kronrod rule */
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

/* xgk[1], xgk[3], ... abscissae of the 30-point gauss rule. 
   xgk[0], xgk[2], ... abscissae to optimally extend the 30-point gauss rule */

static const double wg[15] =    /* weights of the 30-point gauss rule */
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

static const double wgk[31] =   /* weights of the 61-point kronrod rule */
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
#endif


#if 0
// Order of kronrod rule, and number of independent points
static unsigned int gknum = 21;
static unsigned int gknum1 = 11;

static const double xgk[11] =   /* abscissae of the 21-point kronrod rule */
{
  0.995657163025808080735527280689003,
  0.973906528517171720077964012084452,
  0.930157491355708226001207180059508,
  0.865063366688984510732096688423493,
  0.780817726586416897063717578345042,
  0.679409568299024406234327365114874,
  0.562757134668604683339000099272694,
  0.433395394129247190799265943165784,
  0.294392862701460198131126603103866,
  0.148874338981631210884826001129720,
  0.000000000000000000000000000000000
};

/* xgk[1], xgk[3], ... abscissae of the 10-point gauss rule. 
   xgk[0], xgk[2], ... abscissae to optimally extend the 10-point gauss rule */

static const double wg[5] =     /* weights of the 10-point gauss rule */
{
  0.066671344308688137593568809893332,
  0.149451349150580593145776339657697,
  0.219086362515982043995534934228163,
  0.269266719309996355091226921569469,
  0.295524224714752870173892994651338
};

static const double wgk[11] =   /* weights of the 21-point kronrod rule */
{
  0.011694638867371874278064396062192,
  0.032558162307964727478818972459390,
  0.054755896574351996031381300244580,
  0.075039674810919952767043140916190,
  0.093125454583697605535065465083366,
  0.109387158802297641899210590325805,
  0.123491976262065851077958109831074,
  0.134709217311473325928054001771707,
  0.142775938577060080797094273138717,
  0.147739104901338491374841515972068,
  0.149445554002916905664936468389821
};
#endif


////////////////////////////////////////////////////////////////////////
// Constructor
////////////////////////////////////////////////////////////////////////
slug_specsyn::slug_specsyn(const slug_tracks *my_tracks, 
			   const slug_PDF *my_imf,
			   const slug_PDF *my_sfh, const double z_in) :
  z(z_in), tracks(my_tracks), imf(my_imf), sfh(my_sfh)
{ }


////////////////////////////////////////////////////////////////////////
// Destructor
////////////////////////////////////////////////////////////////////////
slug_specsyn::~slug_specsyn() { }


////////////////////////////////////////////////////////////////////////
// Spectral synthesis function for a continuous IMF
////////////////////////////////////////////////////////////////////////

// The get_spectrum_cts function returns the specific luminosity
// L_lambda and bolometric luminosity L_bol for a mono-age stellar
// population of a specified mass, including only the contribution
// from stars that are being treated non-stochastically. Formally,
// this function computes the integrals
//
// L_lambda = (M_tot / <M>) .
//    int_{ln M_min}^{ln M_min,stoch} (dN/dln M) L_lambda(M, t) dlnM
// L_bol = (M_tot / <M>) 
//    int_{ln M_min}^{ln M_min,stoch} (dN/dln M) L_bol(M, t) dlnM
//
// where dN/dM is the IMF (normalized to have an integral of unity
// over its full range), L_lambda(M,t) is the luminosity at
// wavelength lambda for a star of mass M and age t, L_bol(M,t) is the
// bolometric luminosity of a star of mass M and age t, M_tot is the
// total mass in stars over the full mass range, M_min is the minimum
// mass in the IMF, M_min,stoch is the minimm mass being treated
// stochastically, and <M> is the expectation mass of a star computed
// over the full IMF.
//
// The algorithm implemented here is an adaptive Gauss-Kronrod method;
// the code is based on the GSL implementation (significantly
// simplified and ported to c++), adapted to the particular case we
// have where the output is a vector rather than a single number, and
// where it is much more efficient to do interpolation on the tracks for a
// vector of masses all at once, rather than doing them one at a time.

void
slug_specsyn::get_spectrum_cts(const double m_tot, const double age,
			       vector<double>& L_lambda, double& L_bol,
			       const double tol) const {

  // Initialize L_lambda and L_bol
  L_lambda.assign(lambda_rest.size(), 0.0);
  L_bol = 0.0;

  // Allocate workspace
  qag_wksp q(lambda_rest.size(), gknum);

  // Get the range of integration from the IMF and the stellar tracks:
  // minimum mass is the larger of the smallest mass in the IMF and
  // the lowest mass track maximum mass is the smallest of the edge of
  // the non-stochastic range, the largest mass track, and the death
  // mass at this age
  double m_min = max(imf->get_xMin(), tracks->min_mass());
  double m_max = min(min(imf->get_xStochMin(), tracks->max_mass()),
			    tracks->death_mass(age));

  // If m_min > m_max (can happen for strange IMFs if all stars have
  // died), just leave L_lambda and L_bol as 0
  if (m_min >= m_max) return;

  // Do initial integration with Gauss-Kronrod
  double bol_err;
  get_spectrum_cts_gk(m_min, m_max, age, L_lambda, L_bol, 
		      q.errsum, bol_err, q);

  // Get error estimates
  double spec_err = 0.0;
  for (vector<double>::size_type i=0; i<L_lambda.size(); i++)
    spec_err = max(spec_err, q.errsum[i]/L_lambda[i]);
  double err = max(spec_err, bol_err/L_bol);

  // If error is not below tolerance, begin recursive bisection
  if (err > tol) {

    // Initialize the interval, result, and error pointers
    q.a.assign(1, m_min);
    q.b.assign(1, m_max);
    q.r.assign(1, L_lambda);
    q.e.assign(1, q.errsum);
    q.me.assign(1, err);
    q.rbol.assign(1, L_bol);
    q.ebol.assign(1, bol_err);
    vector<double>::size_type intervalptr = 0;
    unsigned int itCounter = 1;

    // Begin iterating
    while (1) {

      // Figure out which interval to work on
      double m_left = q.a[intervalptr];
      double m_right = q.b[intervalptr];
      double m_cen = 0.5 * (m_left + m_right);

      // Compute integrals on two bisected sub-sections
      double L_bol1, L_bol2, bol_err1, bol_err2;
      get_spectrum_cts_gk(m_left, m_cen, age, q.L_out1, L_bol1, 
			  q.err1, bol_err1, q);
      get_spectrum_cts_gk(m_cen, m_right, age, q.L_out2, L_bol2, 
			  q.err2, bol_err2, q);

      // Update result and the error estimate
      for (unsigned int j=0; j<lambda_rest.size(); j++) {
	L_lambda[j] += q.L_out1[j] + q.L_out2[j] - q.r[intervalptr][j];
	q.errsum[j] += q.err1[j] + q.err2[j] - q.e[intervalptr][j];
      }
      L_bol += L_bol1 + L_bol2 - q.rbol[intervalptr];
      bol_err += bol_err1 + bol_err2 - q.ebol[intervalptr];

      // Have we converged? If so, stop iterating
      spec_err = 0.0;
      for (vector<double>::size_type i=0; i<L_lambda.size(); i++)
	spec_err = max(spec_err, q.errsum[i]/L_lambda[i]);
      err = max(spec_err, bol_err/L_bol);
      if (err < tol) break;

      // If we're here, we haven't converged. Replace the current
      // interval with the left half, then push the right half onto
      // the list.
      q.b[intervalptr] = m_cen;
      q.r[intervalptr] = q.L_out1;
      q.e[intervalptr] = q.err1;
      q.rbol[intervalptr] = L_bol1;
      q.ebol[intervalptr] = bol_err1;
      spec_err = 0.0;
      for (vector<double>::size_type i=0; i<L_lambda.size(); i++)
	spec_err = max(spec_err, q.err1[i]/L_lambda[i]);
      q.me[intervalptr] = max(spec_err, bol_err1/L_bol);
      q.a.push_back(m_cen);
      q.b.push_back(m_right);
      q.r.push_back(q.L_out2);
      q.e.push_back(q.err2);
      q.rbol.push_back(L_bol2);
      q.ebol.push_back(bol_err2);
      spec_err = 0.0;
      for (vector<double>::size_type i=0; i<L_lambda.size(); i++)
	spec_err = max(spec_err, q.err2[i]/L_lambda[i]);
      q.me.push_back(max(spec_err, bol_err2/L_bol));

      // Traverse the list of intervals to decide which to work on next
      intervalptr = (vector<double>::size_type) 
	distance(q.me.begin(), 
		 max_element(q.me.begin(), q.me.end()));

      // Update the iteration counter, and check against maximum
      itCounter++;
      if (itCounter > GK_MAX_ITER) {
	cerr << "Error: non-convergence in non-stochastic "
	     << "spectral integration over mass!" << endl;
	exit(1);
      }
    }
  }

  // Apply final normalization
  for (unsigned int j=0; j<lambda_rest.size(); j++) {
    L_lambda[j] *= m_tot / imf->expectationVal();
  }
  L_bol *= m_tot / imf->expectationVal();
}


////////////////////////////////////////////////////////////////////////
// Exactly the same as the previous function, except that this one
// only computes L_bol
////////////////////////////////////////////////////////////////////////
double
slug_specsyn::get_Lbol_cts(const double m_tot, const double age,
			   const double tol) const {

  // Initialize L_bol
  double L_bol = 0.0;

  // Allocate workspace
  qag_wksp q(1, gknum);

  // Get the range of integration from the IMF and the stellar tracks:
  // minimum mass is the larger of the smallest mass in the IMF and
  // the lowest mass track
  // maximum mass is the smallest of the edge of the non-stochastic
  // range, the largest mass track, and the death mass at this age
  double m_min = max(imf->get_xMin(), tracks->min_mass());
  double m_max = min(min(imf->get_xStochMin(), tracks->max_mass()),
			    tracks->death_mass(age));

  // Do initial integration with Gauss-Kronrod
  double bol_err;
  get_Lbol_cts_gk(m_min, m_max, age, L_bol, bol_err);

  // Get error estimate
  double err = bol_err/L_bol;

  // If error is not below tolerance, begin recursive bisection
  if (err > tol) {

    // Initialize the interval, result, and error pointers
    q.a.assign(1, m_min);
    q.b.assign(1, m_max);
    q.rbol.assign(1, L_bol);
    q.ebol.assign(1, bol_err);
    vector<double>::size_type intervalptr = 0;
    unsigned int itCounter = 1;

    // Begin iterating
    while (1) {

      // Figure out which interval to work on
      double m_left = q.a[intervalptr];
      double m_right = q.b[intervalptr];
      double m_cen = 0.5 * (m_left + m_right);

      // Compute integrals on two bisected sub-sections
      double L_bol1, L_bol2, bol_err1, bol_err2;
      get_Lbol_cts_gk(m_left, m_cen, age, L_bol1, bol_err1);
      get_Lbol_cts_gk(m_cen, m_right, age, L_bol2, bol_err2);

      // Update result and the error estimate
      L_bol += L_bol1 + L_bol2 - q.rbol[intervalptr];
      bol_err += bol_err1 + bol_err2 - q.ebol[intervalptr];

      // Have we converged? If so, stop iterating
      err = bol_err/L_bol;
      if (err < tol) break;

      // If we're here, we haven't converged. Replace the current
      // interval with the left half, then push the right half onto
      // the list.
      q.b[intervalptr] = m_cen;
      q.rbol[intervalptr] = L_bol1;
      q.ebol[intervalptr] = bol_err1;
      q.a.push_back(m_cen);
      q.b.push_back(m_right);
      q.rbol.push_back(L_bol2);
      q.ebol.push_back(bol_err2);

      // Traverse the list of intervals to decide which to work on next
      intervalptr = (vector<double>::size_type)
	distance(q.ebol.begin(), 
		 max_element(q.ebol.begin(), q.ebol.end()));

      // Update the iteration counter, and check against maximum
      itCounter++;
      if (itCounter > GK_MAX_ITER) {
	cerr << "Error: non-convergence in non-stochastic "
	     << "spectral integration over time!" << endl;
	exit(1);
      }
    }
  }

  // Apply final normalization
  L_bol *= m_tot / imf->expectationVal();
  return L_bol;
}


////////////////////////////////////////////////////////////////////////
// Spectral synthesis function for a continuous IMF and continuous SFH
////////////////////////////////////////////////////////////////////////

// The get_spectrum_cts_sfh function returns the specific luminosity
// L_lambda and bolometric luminosity L_bol for a population of stars
// over a particular time range at a particular time, including only
// the stars being treated non-stochastically. Formally, this function
// computes the two double integrals:
//
// L_lambda = 1 / <M> .
//    int_{0}^{t} SFR(t') .
//    int_{M_min}^{M_min,stoch} (dN/dM) L_lambda(M, t-t') dM dt'
// L_bol = 1 / <M> .
//    int_{0}^{t} SFR(t') .
//    int_{M_min}^{M_min,stoch} (dN/dM) L_bol(M, t-t') dM dt'
//
// where dN/dM is the IMF (normalized to have an integral of unity
// over its full range), L_lambda(M, t) and L_bol(M, t) are the
// specific luminosity at wavelength lambda and bolometric luminosity
// for a star of mass M and age t, M_tot is the total mass in stars
// over the full mass range, M_min is the minimum mass in the IMF,
// M_min,stoch is the minimm mass being treated stochastically, <M> is
// the expectation mass of a star computed over the full IMF, SFR(t)
// is the star formation rate at age t, and t is the time for which
// star formation has been going on.
//
// The algorithm implemented here is an adaptive Gauss-Kronrod method,
// heavily cribbed from the GSL.

void
slug_specsyn::
get_spectrum_cts_sfh(const double t, vector<double>& L_lambda, 
		     double& L_bol, const double tol) const {

  // Initialize L_lambda
  L_lambda.assign(lambda_rest.size(), 0.0);

  // Allocate workspace
  qag_wksp q(lambda_rest.size(), gknum);

  // Do initial integration with Gauss-Kronrod
  double err_bol;
  get_spectrum_cts_sfh_gk(0, t, t, L_lambda, L_bol,
			  q.errsum, err_bol, q, tol);

  // Get error estimates
  double max_L = *max_element(L_lambda.begin(), L_lambda.end());
  double max_err = *max_element(q.errsum.begin(), q.errsum.end());
  double err = max(max_err/max_L, err_bol/L_bol);

  // If error is not below tolerance, begin recursive bisection
  if (err > tol) {

    // Initialize the interval, result, and error pointers
    q.a.assign(1, 0);
    q.b.assign(1, t);
    q.r.assign(1, L_lambda);
    q.e.assign(1, q.errsum);
    q.me.assign(1, err);
    q.rbol.assign(1, L_bol);
    q.ebol.assign(1, err_bol);
    vector<double>::size_type intervalptr = 0;
    unsigned int itCounter = 1;

    // Begin iterating
    while (1) {

      // Figure out which interval to work on
      double t_left = q.a[intervalptr];
      double t_right = q.b[intervalptr];
      double t_cen = 0.5 * (t_left + t_right);

      // Compute integrals on two bisected sub-sections
      double L_bol1, L_bol2, err_bol1, err_bol2;
      get_spectrum_cts_sfh_gk(t_left, t_cen, t, q.L_out1, L_bol1,
			      q.err1, err_bol1, q, tol);
      get_spectrum_cts_sfh_gk(t_cen, t_right, t, q.L_out2, L_bol2,
			      q.err2, err_bol2, q, tol);

      // Update result and the error estimate
      for (unsigned int j=0; j<lambda_rest.size(); j++) {
	L_lambda[j] += q.L_out1[j] + q.L_out2[j] - q.r[intervalptr][j];
	q.errsum[j] += q.err1[j] + q.err2[j] - q.e[intervalptr][j];
      }
      L_bol += L_bol1 + L_bol2 - q.rbol[intervalptr];
      err_bol += err_bol1 + err_bol2 - q.ebol[intervalptr];

      // Have we converged? If so, stop iterating
      max_L = *max_element(L_lambda.begin(), L_lambda.end());
      max_err = *max_element(q.errsum.begin(), q.errsum.end());
      err = max(max_err/max_L, err_bol/L_bol);
      if (err < tol) break;

      // If we're here, we haven't converged. Replace the current
      // interval with the left half, then push the right half onto
      // the list.
      q.b[intervalptr] = t_cen;
      q.r[intervalptr] = q.L_out1;
      q.e[intervalptr] = q.err1;
      q.rbol[intervalptr] = L_bol1;
      q.ebol[intervalptr] = err_bol1;
      q.me[intervalptr] = 
	max((*max_element(q.err1.begin(), q.err1.end())) / max_L,
	    err_bol1 / L_bol);
      q.a.push_back(t_cen);
      q.b.push_back(t_right);
      q.r.push_back(q.L_out2);
      q.e.push_back(q.err2);
      q.rbol.push_back(L_bol2);
      q.ebol.push_back(err_bol2);
      q.me.push_back(
	max((*max_element(q.err2.begin(), q.err2.end())) / max_L,
	    err_bol2 / L_bol));

      // Traverse the list of intervals to decide which to work on next
      intervalptr = (vector<double>::size_type)
	distance(q.me.begin(), 
		 max_element(q.me.begin(), q.me.end()));

      // Update the iteration counter, and check against maximum
      itCounter++;
      if (itCounter > GK_MAX_ITER) {
	cerr << "Error: non-convergence in non-stochastic "
	     << "spectral integration for SFH!" << endl;
	exit(1);
      }
    }
  }

  // Apply final normalization
  for (unsigned int j=0; j<lambda_rest.size(); j++) {
    L_lambda[j] /= imf->expectationVal();
  }
  L_bol /= imf->expectationVal();
}


////////////////////////////////////////////////////////////////////////
// Exactly the same as the previous function, except that this one
// only computes L_bol
////////////////////////////////////////////////////////////////////////
double
slug_specsyn::
get_Lbol_cts_sfh(const double t, const double tol) const {

  // Allocate workspace
  double L_bol, err_bol;
  qag_wksp q(1, gknum);

  // Do initial integration with Gauss-Kronrod
  get_Lbol_cts_sfh_gk(0, t, t, L_bol, err_bol, tol);

  // Get error estimate
  double err = err_bol/L_bol;

  // If error is not below tolerance, begin recursive bisection
  if (err > tol) {

    // Initialize the interval, result, and error pointers
    q.a.assign(1, 0);
    q.b.assign(1, t);
    q.rbol.assign(1, L_bol);
    q.ebol.assign(1, err_bol);
    vector<double>::size_type intervalptr = 0;
    unsigned int itCounter = 1;

    // Begin iterating
    while (1) {

      // Figure out which interval to work on
      double t_left = q.a[intervalptr];
      double t_right = q.b[intervalptr];
      double t_cen = 0.5 * (t_left + t_right);

      // Compute integrals on two bisected sub-sections
      double L_bol1, L_bol2, err_bol1, err_bol2;
      get_Lbol_cts_sfh_gk(t_left, t_cen, t, L_bol1, err_bol1, tol);
      get_Lbol_cts_sfh_gk(t_cen, t_right, t, L_bol2, err_bol2, tol);

      // Update result and the error estimate
      L_bol += L_bol1 + L_bol2 - q.rbol[intervalptr];
      err_bol += err_bol1 + err_bol2 - q.ebol[intervalptr];

      // Have we converged? If so, stop iterating
      err = err_bol/L_bol;
      if (err < tol) break;

      // If we're here, we haven't converged. Replace the current
      // interval with the left half, then push the right half onto
      // the list.
      q.b[intervalptr] = t_cen;
      q.rbol[intervalptr] = L_bol1;
      q.ebol[intervalptr] = err_bol1;
      q.a.push_back(t_cen);
      q.b.push_back(t_right);
      q.rbol.push_back(L_bol2);
      q.ebol.push_back(err_bol2);

      // Traverse the list of intervals to decide which to work on next
      intervalptr = (vector<double>::size_type)
	distance(q.ebol.begin(), 
		 max_element(q.ebol.begin(), q.ebol.end()));

      // Update the iteration counter, and check against maximum
      itCounter++;
      if (itCounter > GK_MAX_ITER) {
	cerr << "Error: non-convergence in non-stochastic "
	     << "spectral integration!" << endl;
	exit(1);
      }
    }
  }

  // Apply final normalization
  L_bol /= imf->expectationVal();
  return L_bol;
}


////////////////////////////////////////////////////////////////////////
// Helper function to evaluate GK rule on a mass interval. This
// structure of this code closely follows the GSL qk routine.
////////////////////////////////////////////////////////////////////////
void
slug_specsyn::
get_spectrum_cts_gk(const double m_min, const double m_max,
		    const double age, vector<double>& L_lambda, 
		    double& L_bol, vector<double>& err,
		    double& err_bol, qag_wksp& q) const {

  // Initialize the Gauss sum accumulators zero
  double L_bol_gauss = 0.0;
  q.gaussQuad.assign(lambda_rest.size(), 0.0);

  // Construct grid of mass points
  double m_cen = 0.5 * (m_min + m_max);
  double half_length = 0.5 * (m_max - m_min);
  for (unsigned int i=0; i<gknum/2; i++) {
    q.x_k[i] = m_cen - half_length * xgk[i];
    q.x_k[gknum-i-1] = m_cen + half_length * xgk[i];
  }
  q.x_k[gknum/2] = m_cen;

  // Get stellar data for the mass grid
  const vector<slug_stardata> &stardata = 
    tracks->get_isochrone(age, q.x_k);

  // Now form the Gauss and Kronrod sums

  // Central mass point
  unsigned int ptr1 = gknum/2;
  unsigned int ptr2;

  // Get spectrum at this mass
  q.L_tmp1 = get_spectrum(stardata[ptr1]);

  // Get IMF at this mass
  double imf_val1 = (*imf)(q.x_k[ptr1]);
  double imf_val2;

  // Get bolometric luminosity at this mass
  double L1 = pow(10.0, stardata[ptr1].logL);
  double L2;

  // Add to Kronrod sum, and to Gauss sum if central point appears in
  // it
  for (unsigned int j=0; j<lambda_rest.size(); j++)
    L_lambda[j] = q.L_tmp1[j] * imf_val1 * wgk[gknum1-1];
  L_bol = L1 * imf_val1 * wgk[gknum1-1];
  if (gknum1 % 2 == 0) {
    for (unsigned int j=0; j<lambda_rest.size(); j++)
      q.gaussQuad[j] = q.L_tmp1[j] * imf_val1 * wg[gknum1 / 2 - 1];
    L_bol_gauss = L1 * imf_val1 * wg[gknum1 / 2 - 1];
  }
  
  // Compute terms that are common to both Gauss and Kronrod sum
  for (unsigned int i=0; i<(gknum1-1)/2; i++) {

    // Point on the left side of the mass interval
    ptr1 = 2*i+1;
    q.L_tmp1 = get_spectrum(stardata[ptr1]);
    imf_val1 = (*imf)(q.x_k[ptr1]);
    L1 = pow(10.0, stardata[ptr1].logL);

    // Point on the right side of the mass interval
    ptr2 = gknum - 2*i - 2;
    q.L_tmp2 = get_spectrum(stardata[ptr2]);
    imf_val2 = (*imf)(q.x_k[ptr2]);
    L2 = pow(10.0, stardata[ptr2].logL);

    // Compute the contribution to the Gaussian and Kronrod quadratures
    for (unsigned int j=0; j<lambda_rest.size(); j++) {
      q.gaussQuad[j] += wg[i] * 
	(q.L_tmp1[j]*imf_val1 + q.L_tmp2[j]*imf_val2);
      L_lambda[j] += wgk[ptr1] * 
	(q.L_tmp1[j]*imf_val1 + q.L_tmp2[j]*imf_val2);
    }
    L_bol += wgk[ptr1] * (L1*imf_val1 + L2*imf_val2);
    L_bol_gauss += wg[i] * (L1*imf_val1 + L2*imf_val2);
  }

  // Compute terms that appear only in the Kronrod sum
  for (unsigned int i=0; i<gknum1/2; i++) {

    // Point on left half of interval
    ptr1 = 2*i;
    q.L_tmp1 = get_spectrum(stardata[ptr1]);
    imf_val1 = (*imf)(q.x_k[ptr1]);
    L1 = pow(10.0, stardata[ptr1].logL);

    // Point on right half of interval
    ptr2 = gknum - 2*i - 1;
    q.L_tmp2 = get_spectrum(stardata[ptr2]);
    imf_val2 = (*imf)(q.x_k[ptr2]);
    L2 = pow(10.0, stardata[ptr2].logL);

    // Add to Kronrod sum
    for (unsigned int j=0; j<lambda_rest.size(); j++)
      L_lambda[j] += wgk[ptr1] * 
	(q.L_tmp1[j]*imf_val1 + q.L_tmp2[j]*imf_val2);
    L_bol += wgk[ptr1] * (L1*imf_val1 + L2*imf_val2);
  }

  // Scale results by length of mass interval to properly normalize
  for (unsigned int j=0; j<lambda_rest.size(); j++) {
    L_lambda[j] *= half_length;
    q.gaussQuad[j] *= half_length;
  }
  L_bol *= half_length;
  L_bol_gauss *= half_length;

  // Compute error
  for (unsigned int j=0; j<lambda_rest.size(); j++) {
    err[j] = abs(L_lambda[j] - q.gaussQuad[j]);
  }
  err_bol = abs(L_bol - L_bol_gauss);
}

////////////////////////////////////////////////////////////////////////
// Same as previous function, but only for Lbol; all steps related to
// the spectrum are omitted
////////////////////////////////////////////////////////////////////////
void
slug_specsyn::
get_Lbol_cts_gk(const double m_min, const double m_max,
		const double age, double& L_bol, 
		double& err_bol) const {

  // Initialize the Gauss sum accumulator zero
  double L_bol_gauss = 0.0;

  // Construct grid of mass points
  vector<double> x_k(gknum);
  double m_cen = 0.5 * (m_min + m_max);
  double half_length = 0.5 * (m_max - m_min);
  for (unsigned int i=0; i<gknum/2; i++) {
    x_k[i] = m_cen - half_length * xgk[i];
    x_k[gknum-i-1] = m_cen + half_length * xgk[i];
  }
  x_k[gknum/2] = m_cen;

  // Get stellar data for the mass grid
  const vector<slug_stardata> &stardata = 
    tracks->get_isochrone(age, x_k);

  // Now form the Gauss and Kronrod sums

  // Central mass point
  unsigned int ptr1 = gknum/2;
  unsigned int ptr2;

  // Get IMF at this mass
  double imf_val1 = (*imf)(x_k[ptr1]);
  double imf_val2;

  // Get bolometric luminosity at this mass
  double L1 = pow(10.0, stardata[ptr1].logL);
  double L2;

  // Add to Kronrod sum, and to Gauss sum if central point appears in
  // it
  L_bol = L1 * imf_val1 * wgk[gknum1-1];
  if (gknum1 % 2 == 0) {
    L_bol_gauss = L1 * imf_val1 * wg[gknum1 / 2 - 1];
  }
  
  // Compute terms that are common to both Gauss and Kronrod sum
  for (unsigned int i=0; i<(gknum1-1)/2; i++) {

    // Point on the left side of the mass interval
    ptr1 = 2*i+1;
    imf_val1 = (*imf)(x_k[ptr1]);
    L1 = pow(10.0, stardata[ptr1].logL);

    // Point on the right side of the mass interval
    ptr2 = gknum - 2*i - 2;
    imf_val2 = (*imf)(x_k[ptr2]);
    L2 = pow(10.0, stardata[ptr2].logL);

    // Compute the contribution to the Gauss and Kronrod quadratures
    L_bol += wgk[ptr1] * (L1*imf_val1 + L2*imf_val2);
    L_bol_gauss += wg[i] * (L1*imf_val1 + L2*imf_val2);
  }

  // Compute terms that appear only in the Kronrod sum
  for (unsigned int i=0; i<gknum1/2; i++) {

    // Point on left half of interval
    ptr1 = 2*i;
    imf_val1 = (*imf)(x_k[ptr1]);
    L1 = pow(10.0, stardata[ptr1].logL);

    // Point on right half of interval
    ptr2 = gknum - 2*i - 1;
    imf_val2 = (*imf)(x_k[ptr2]);
    L2 = pow(10.0, stardata[ptr2].logL);

    // Add to Kronrod sum
    L_bol += wgk[ptr1] * (L1*imf_val1 + L2*imf_val2);
  }

  // Scale results by length of mass interval to properly normalize
  L_bol *= half_length;
  L_bol_gauss *= half_length;

  // Compute error
  err_bol = abs(L_bol - L_bol_gauss);
}


////////////////////////////////////////////////////////////////////////
// Helper function to evaluate GK rule on a time interval. This
// structure of this code closely follows the GSL qk routine.
////////////////////////////////////////////////////////////////////////
void
slug_specsyn::
get_spectrum_cts_sfh_gk(const double t_min, const double t_max, 
			const double t, vector<double>& L_lambda, 
			double& L_bol, vector<double>& err, 
			double& err_bol, qag_wksp& q,
			const double tol) const {

  // Initialize the accumulator for the Gauss sum to zero
  double L_bol_gauss = 0.0;
  q.gaussQuad.assign(lambda_rest.size(), 0.0);

  // Construct grid of time points
  double t_cen = 0.5 * (t_min + t_max);
  double half_length = 0.5 * (t_max - t_min);
  for (unsigned int i=0; i<gknum/2; i++) {
    q.x_k[i] = t_cen - half_length * xgk[i];
    q.x_k[gknum-i-1] = t_cen + half_length * xgk[i];
  }
  q.x_k[gknum/2] = t_cen;

  // Now form the Gauss and Kronrod sums

  // Central mass point
  unsigned int ptr1 = gknum/2;
  unsigned int ptr2;

  // Get spectrum at the central time, leaving a fair margin of error
  // in the tolerance; note that we normalize to 1 Msun here, and fix
  // the normalization later
  double L_bol1, L_bol2;
  get_spectrum_cts(1.0, t-q.x_k[gknum/2],  q.L_tmp1, L_bol1, tol/10.0);

  // Get SFR at central time
  double sfh_val1 = (*sfh)(q.x_k[ptr1]);
  double sfh_val2;

  // Add to Kronrod sum, and to Gauss sum if central point appears in
  // it
  for (unsigned int j=0; j<lambda_rest.size(); j++)
    L_lambda[j] = q.L_tmp1[j] * sfh_val1 * wgk[gknum1-1];
  L_bol = L_bol1 * sfh_val1 * wgk[gknum1-1];
  if (gknum1 % 2 == 0) {
    for (unsigned int j=0; j<lambda_rest.size(); j++)
      q.gaussQuad[j] = q.L_tmp1[j] * sfh_val1 * wg[gknum1 / 2 - 1];
    L_bol_gauss = L_bol1 * sfh_val1 * wg[gknum1 / 2 - 1];
  }
  
  // Compute terms that are common to both Gauss and Kronrod sum
  for (unsigned int i=0; i<(gknum1-1)/2; i++) {

    // Point on the left side of the mass interval
    ptr1 = 2*i+1;
    get_spectrum_cts(1.0, t-q.x_k[ptr1], q.L_tmp1, L_bol1, tol/10.0);
    sfh_val1 = (*sfh)(q.x_k[ptr1]);

    // Point on the right side of the mass interval
    ptr2 = gknum - 2*i - 2;
    get_spectrum_cts(1.0, t-q.x_k[ptr2], q.L_tmp2, L_bol2, tol/10.0);
    sfh_val2 = (*sfh)(q.x_k[ptr2]);

    // Compute the contribution to the Gaussian and Kronrod quadratures
    for (unsigned int j=0; j<lambda_rest.size(); j++) {
      q.gaussQuad[j] += wg[i] * 
	(q.L_tmp1[j]*sfh_val1 + q.L_tmp2[j]*sfh_val2);
      L_lambda[j] += wgk[ptr1] * 
	(q.L_tmp1[j]*sfh_val1 + q.L_tmp2[j]*sfh_val2);
    }
    L_bol += wgk[ptr1] * (L_bol1*sfh_val1 + L_bol2*sfh_val2);
    L_bol_gauss += wg[i] * (L_bol1*sfh_val1 + L_bol2*sfh_val2);
  }

  // Compute terms that appear only in the Kronrod sum
  for (unsigned int i=0; i<gknum1/2; i++) {

    // Point on left half of interval
    ptr1 = 2*i;
    get_spectrum_cts(1.0, t-q.x_k[ptr1], q.L_tmp1, L_bol1, tol/10.0);
    sfh_val1 = (*sfh)(q.x_k[ptr1]);

    // Point on right half of interval
    ptr2 = gknum - 2*i - 1;
    get_spectrum_cts(1.0, t-q.x_k[ptr2], q.L_tmp2, L_bol2, tol/10.0);
    sfh_val2 = (*sfh)(q.x_k[ptr2]);

    // Add to Kronrod sum
    for (unsigned int j=0; j<lambda_rest.size(); j++)
      L_lambda[j] += wgk[ptr1] * 
	(q.L_tmp1[j]*sfh_val1 + q.L_tmp2[j]*sfh_val2);
    L_bol += wgk[ptr1] * (L_bol1*sfh_val1 + L_bol2*sfh_val2);
  }

  // Scale results by length of mass interval to properly normalize
  for (unsigned int j=0; j<lambda_rest.size(); j++) {
    L_lambda[j] *= half_length;
    q.gaussQuad[j] *= half_length;
  }
  L_bol *= half_length;
  L_bol_gauss *= half_length;

  // Compute error
  for (unsigned int j=0; j<lambda_rest.size(); j++) {
    err[j] = abs(L_lambda[j] - q.gaussQuad[j]);
  }
  err_bol = abs(L_bol - L_bol_gauss);
}


////////////////////////////////////////////////////////////////////////
// Same as previous function, but only for Lbol; all steps related to
// the spectrum are omitted
////////////////////////////////////////////////////////////////////////
void
slug_specsyn::
get_Lbol_cts_sfh_gk(const double t_min, const double t_max, 
		    const double t, double& L_bol, double& err_bol,
		    const double tol) const {

  // Initialize the accumulator for the Gauss sum to zero
  double L_bol_gauss = 0.0;

  // Construct grid of time points
  vector<double> x_k(gknum);
  double t_cen = 0.5 * (t_min + t_max);
  double half_length = 0.5 * (t_max - t_min);
  for (unsigned int i=0; i<gknum/2; i++) {
    x_k[i] = t_cen - half_length * xgk[i];
    x_k[gknum-i-1] = t_cen + half_length * xgk[i];
  }
  x_k[gknum/2] = t_cen;

  // Now form the Gauss and Kronrod sums

  // Central mass point
  unsigned int ptr1 = gknum/2;
  unsigned int ptr2;

  // Get Lbol at the central time, leaving a fair margin of error
  // in the tolerance; note that we normalize to 1 Msun here, and fix
  // the normalization later
  double L_bol1 = get_Lbol_cts(1.0, t-x_k[gknum/2], tol/10.0);
  double L_bol2;

  // Get SFR at central time
  double sfh_val1 = (*sfh)(x_k[ptr1]);
  double sfh_val2;

  // Add to Kronrod sum, and to Gauss sum if central point appears in
  // it
  L_bol = L_bol1 * sfh_val1 * wgk[gknum1-1];
  if (gknum1 % 2 == 0) {
    L_bol_gauss = L_bol1 * sfh_val1 * wg[gknum1 / 2 - 1];
  }
  
  // Compute terms that are common to both Gauss and Kronrod sum
  for (unsigned int i=0; i<(gknum1-1)/2; i++) {

    // Point on the left side of the mass interval
    ptr1 = 2*i+1;
    L_bol1 = get_Lbol_cts(1.0, t-x_k[ptr1], tol/10.0);
    sfh_val1 = (*sfh)(x_k[ptr1]);

    // Point on the right side of the mass interval
    ptr2 = gknum - 2*i - 2;
    L_bol2 = get_Lbol_cts(1.0, t-x_k[ptr2], tol/10.0);
    sfh_val2 = (*sfh)(x_k[ptr2]);

    // Compute the contribution to the Gaussian and Kronrod quadratures
    L_bol += wgk[ptr1] * (L_bol1*sfh_val1 + L_bol2*sfh_val2);
    L_bol_gauss += wg[i] * (L_bol1*sfh_val1 + L_bol2*sfh_val2);
  }

  // Compute terms that appear only in the Kronrod sum
  for (unsigned int i=0; i<gknum1/2; i++) {

    // Point on left half of interval
    ptr1 = 2*i;
    L_bol1 = get_Lbol_cts(1.0, t-x_k[ptr1], tol/10.0);
    sfh_val1 = (*sfh)(x_k[ptr1]);

    // Point on right half of interval
    ptr2 = gknum - 2*i - 1;
    L_bol2 = get_Lbol_cts(1.0, t-x_k[ptr2], tol/10.0);
    sfh_val2 = (*sfh)(x_k[ptr2]);

    // Add to Kronrod sum
    L_bol += wgk[ptr1] * (L_bol1*sfh_val1 + L_bol2*sfh_val2);
  }

  // Scale results by length of mass interval to properly normalize
  L_bol *= half_length;
  L_bol_gauss *= half_length;

  // Compute error
  err_bol = abs(L_bol - L_bol_gauss);
}
