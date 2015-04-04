/**
 Contrails
 @author Martin Lilleeng Saetra <martinls@met.no>

 $Id: Contrails.h 1203 2014-11-27 17:05:25Z martinls $

 Copyright (C) 2009 met.no

 Contact information:
 Norwegian Meteorological Institute
 Box 43 Blindern
 0313 OSLO
 NORWAY
 email: diana@met.no

 This is free software; you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation; either version 2 of the License, or
 (at your option) any later version.

 You should have received a copy of the GNU General Public License
 along with this software; if not, write to the Free Software
 Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
 */

#ifndef CONTRAILS_
#define CONTRAILS_

//#define DEBUG

#include <cmath>
#include <cstddef>

using namespace std;

class Contrails {
private:
  /// domain size
  int nx, ny;

  /// initializing constants
  float p0;
  float cp;
  float r;
  float rcp;
  float epsi;
  float T0;

  float a;
  float b;
  float c;
  float d;
  float f;
  float h;
  float est; // mb
  float Ts;  // K

  /// The sum of all vertical fields.
  float* sum;

  /// special fields
  float* contrailsBottom;
  float* contrailsTop;

  /**
   * Set domain size.
   * @param nx Size in the x-dimension
   * @param ny Size in the y-dimension
   */
  void initialize(int nx, int ny);

  /**
   * Saturation specific humidity.
   * @param esi Saturation vapor pressure over ICE
   * @param p Pressure
   */
  float qs(float esi, float p);

  /**
   * Saturation vapor pressure over WATER.
   * (esi can be found in several ways, here we have used the
   * Goff-Gratch equation:
   * http://hurri.kean.edu/~yoh/calculations/satvap/satvap.html
   * (see Saturation Vapor Pressure over WATER))
   * @param t Temperature
   */
  float esi(float t);

  /**
   * Relative humidity.
   * @param q Specific humidity
   * @param qs Saturation specific humidity
   */
  float RH(float q, float qs);

  /**
   * The Schmidt-Appleman criterion.
   * @param Tc Critical temperature
   * @param t Environment temperature
   */
  int SchmidtAppleman(float Tc, float t);

  /**
   * Check if domain size have been set.
   * @return True if domain size is set, false otherwise
   */
  bool initialized();

public:
  /**
   * Constructor.
   */
  Contrails();

  /**
   * Destructor.
   */
  virtual ~Contrails();

  /**
   * Returns the sum of all vertical fields.
   */
  float* getSumField();

  float* getContrailsBottomField();
  float* getContrailsTopField();

  /**
   * Reset integration sum between timesteps and sum-field.
   */
  void reset();

  /**
   * Calculates the areas in a two-dimensional (spatial) field
   * where condensation trails can occur.
   * @param nx width of 2D fields
   * @param ny height of 2D fields
   * @param ap hybrid coordinate
   * @param b hybrid coordinate
   * @param data Results are written to this array
   * @param ps Surface pressure
   * @param t Air temperature
   * @param q Specific humidity
   * @param p Pressure. (If not given it will be calculated from ps.)
   */
  void calcContrails2D(int nx, int ny, float ap, float b, float *data,
		    float *ps, float *t, float *q, float *p=NULL);

};

/**
 * Calculate critical temperature, used in the
 * Schmidt-Appleman criterion.
 * @param pamb Pressure
 * @param rhamb Relative humidity
 */
float Tcrit(float pamb, float rhamb);

inline float Contrails::qs(float esi, float p)
{
  return epsi*esi/p;
}

inline float Contrails::RH(float q, float qs)
{
  return q/qs;
}

inline int Contrails::SchmidtAppleman(float Tc, float t)
{
#ifdef DEBUG
  cout << "SchmidtAppleman: Tc: " << Tc << " - t: " << t << " = " << Tc-t << endl;
#endif

  if (Tc > t) {
    return 1;
  }
  return 0;
}

#endif /*CONTRAILS_*/
