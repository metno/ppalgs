#ifndef UTILS_HPP_
#define UTILS_HPP_

/*
 ppalgs - Utils

 $Id: utils.hpp 1205 2014-12-11 16:03:33Z martinls $

 Copyright (C) 2008 met.no

 Contact information:
 Norwegian Meteorological Institute
 Box 43 Blindern
 0313 OSLO
 NORWAY

 This is free software; you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation; either version 2 of the License, or
 (at your option) any later version.
 
 You should have received a copy of the GNU General Public License
 along with this software; if not, write to the Free Software
 Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
 */

#include <string>
#include <sstream>
#include <cmath>

using namespace std;

/**
 * Find height in meters.
 * @param t Temperature in Kelvin
 * @param p Pressure
 * @param ps Surface pressure
 * @return Height in meters
 */
inline float z(float t, float p, float ps)
{
  // constants
  float r = 287.05;;
  float g = 9.80665;

  return -(r*t/g) * (log(p/ps));
}

#endif /* UTILS_HPP_ */
