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

// Constants, unit conversions, and related stuff used in various places
// throughout the code.

#ifndef _constants_H_
#define _constants_H_

#include <limits>
#include <vector>

namespace constants {

  // Physical constants -- 2010 CODATA values, or derived therefrom,
  // all in CGS units
  const double G = 6.67384e-8;
  const double c = 2.99792458e10;
  const double h = 6.62606957e-27;
  const double hbar = 1.054571725e-27;
  const double kB = 1.3806488e-16;
  const double sigmaSB = 5.670373e-5;
  const double hcOverkB = 1.4387769599838156;
  const double hc = 1.986445683269303e-16;
  const double hc2 = 5.955214340707938e-06;
  const double Ryd = 2.1798721703680373e-11;
  const double echarge = 4.80320425e-10;
  const double melectron = 9.10938291e-28;

  // Ionization threshhold wavelengths for HI, HeI, HeII, in Angstrom
  const double lambdaHI = 911.2670505509151;
  const double lambdaHeI = 503.9888;
  const double lambdaHeII = 227.2979;

  // Number ratio of He to H
  const double xHe = 0.1;

  // Unit conversions to CGS
  const double Msun = 1.989e33;
  const double Rsun = 6.955e10;
  const double Lsun = 3.846e33;
  const double Angstrom = 1.0e-8;
  const double pc = 3.0856775814671918e18;
  const double kpc = 3.0856775814671918e21;
  const double yr = 31557600.0;

  // Log base 10 of constants and conversions
  const double logsigmaSB = -4.2463884;
  const double logG = -7.1756242;
  const double logMsun = 33.2986566;
  const double logLsun = 33.5850093;
  const double logRsun = 10.8422971;

  // Log base 10 of e
  const double loge = 0.43429448190325182;

  // Numeric limits and sentinels
  const double big = std::numeric_limits<double>::max();
  const double small = std::numeric_limits<double>::min();
  const std::vector<int>::size_type sentinel =
    std::numeric_limits<std::vector<int>::size_type>::max();
}

#endif
// _constants_H_
