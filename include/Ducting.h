/**
  Ducting
  @author Martin Lilleeng Saetra <martinls@met.no>
  (Based on version in fortran77 by DNMI/FoU  07.06.1994  Laila Sidselrud)
  
  $Id: Ducting.h 1196 2014-10-01 14:14:18Z martinls $

  Copyright (C) 2008 met.no

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

#ifndef DUCTING_H_
#define DUCTING_H_

#include <cmath>
#include <cstddef>

using namespace std;

class Ducting
{
private:
	/// domain size
	int nx, ny;
	
	/// initializing constants
	float p0;
	float cp;
	float r;
	float rcp;
	float eps;
	
	/// last z (z minus) and M values
	float* zm;
	float* Mm;

	// 1.0 (true) indicates that we are inside a duct for this vertical column
	float* inDuct;
	float* currentElevatedDuctBottom;
	float* currentElevatedDuctTop;

	/// single-level parameters
	float* dMdzMin;
	float* surfaceDuctBottom;
	float* surfaceDuctTop;
	float* elevatedDuctBottom;
	float* elevatedDuctTop;
	float* noOfElevatedDucts;
	
	/**
	 * Set domain size.
	 * @param nx Size in the x-dimension
	 * @param ny Size in the y-dimension
	 */
	void initialize(int nx, int ny);
	
	/**
	 * Check if domain size have been set.
	 * @return True if domain size is set, false otherwise
	 */
	bool initialized();
	
	/**
	 * Find the z-differential.
	 * @param z Height in meters
	 * @param i Index in x-dimension
	 * @param j Index in y-dimension
	 * @return z-differential
	 */
	float dz(float z, int i, int j);
	
	/**
	 * The Exner function.
	 * @param p Pressure
	 * @return The pi of the Exner function
	 */
	float exner(float p);
	
	/**
	 * Calculate absolute temperature from potential temperature and the
	 * Exner function result (pi).
	 * @param th Potential temperature
	 * @param pi The pi of the Exner function
	 * @return Absolute temperature in Kelvin
	 */
	float t(float th, float pi);
	
	/**
	 * Calculate potential temperature from absolute temperature and the
	 * Exner function result (pi).
	 * @param t Absolute temperature
	 * @param pi The pi of the Exner function
	 * @return Potential temperature in Kelvin
	 */
	float th(float t, float pi);

	/**
	 * Calculate M.
	 * @param N A value to describe ducting; N-units
	 * @param z Height in meters
	 * @return M-units; N-units with the curve of the earth surface taken into account
	 */ 
	float M(float N, float z);
	
	/**
	 * Calculate N.
	 * @param t Temperature in Kelvin
	 * @param q Specific humidity
	 * @param p Pressure
	 * @return N-units describing ducting conditions
	 */
	float N(float t, float q, float p);
	
	/**
	 * Calculate the ducting gradient for (i,j)
	 * @param i Index in x-dimension
	 * @param j Index in y-dimension
	 * @param data_m M-values are written to this array
     * @param data_z Heights in meters are written to this array
	 * @param ps Pointer to surface pressure data (hPa)
	 * @param p Pointer to pressure data (hPa)
	 * @param th Pointer to potential temperature data
	 * @param q Pointer to specific humidity data
	 * @return dM/dz for (i,j)
	 */
	float dMdz(int i, int j, float* data_m, float* data_z, float ps, float p, float th, float q);
	
public:
	/**
	 * Constructor.
	 */
	Ducting();
	
	/**
	 * Destructor.
	 */
	virtual ~Ducting();
	
	/**
	 * Get the minimum dMdz field.
	 */
	float* getdMdzMinField();
	
	/**
	 * Lowest height of surface duct (below 1000 ft)
	 *
	 * NOTE:
	 * - dM/dz < 0 indicates a duct
	 * - if multiple ducts are present, the thickest is chosen
	 * - if multiple ducts have the same thickness, the lowest is chosen
	 */
	float* getSurfaceDuctBottom();

	/**
	 * Highest height of surface duct (below 1000 ft)
	 *
	 * NOTE:
	 * - dM/dz < 0 indicates a duct
	 * - if multiple ducts are present, the thickest is chosen
	 * - if multiple ducts have the same thickness, the lowest is chosen
	 */
	float* getSurfaceDuctTop();

	/**
	 * Lowest height of elevated duct (1000 ft or higher)
	 *
	 * NOTE:
	 * - dM/dz < 0 indicates a duct
	 * - if multiple ducts are present, the thickest is chosen
	 * - if multiple ducts have the same thickness, the lowest is chosen
	 */
	float* getElevatedDuctBottom();

	/**
	 * Highest height of elevated duct (1000 ft or higher)
	 *
	 * NOTE:
	 * - dM/dz < 0 indicates a duct
	 * - if multiple ducts are present, the thickest is chosen
	 * - if multiple ducts have the same thickness, the lowest is chosen
	 */
	float* getElevatedDuctTop();

	/**
	 * Total number of elevated ducts (in the vertical dimension)
	 */
	float* getNoOfElevatedDucts();

	/**
	 * Reset all max/min fields.
	 */
	void reset();
	
	/**
	 * Calculates the ducting gradients (dM/dz) of a 
	 * two-dimensional (spatial) field.
	 * @param nx width of 2D fields
	 * @param ny height of 2D fields
	 * @param ap hybrid coordinate
	 * @param b hybrid coordinate
	 * @param data Ducting gradients dM/dz are written to this array
	 * @param data_m M-values are written to this array
	 * @param data_z Heights in meters are written to this array
	 * @param ps Surface pressure
	 * @param t Air temperature
	 * @param q Specific humidity
	 * @param p Pressure. (If not given it will be calculated from ps.)
	 */
	void calcDuctingGrads2D(int nx, int ny, float ap, float b,
			float* data, float* data_m, float* data_z, float* ps, float* t, float* q, float* p=NULL);
};

inline float Ducting::exner(float p)
{
	return pow((p / p0), rcp);
}

inline float Ducting::t(float th, float pi)
{
	return pi * th;
}

inline float Ducting::th(float t, float pi)
{
	return pi / t;
}

inline float Ducting::N(float t, float q, float p)
{
	return 77.6 * (p/t) + 373000. * (q*p)/(eps*t*t);
}

inline float Ducting::M(float N, float z)
{
	return N + 157. * z * 0.001;
}

#endif /*DUCTING_H_*/
