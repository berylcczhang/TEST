/*********************************************************************
Copyright (C) 2014-7 Robert da Silva, Michele Fumagalli, Mark
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

#include "slug_yields_snii.H"

using namespace std;

bool slug_yields_snii::produces_sn(const double m) const {
  bool sn_flag = false;
  for (vector<double>::size_type i=0; i<sn_mass.size(); i++) {
    if (m < sn_mass[i]) break;
    sn_flag = !sn_flag;
  }
  return sn_flag;
}
