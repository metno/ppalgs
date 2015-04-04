/**
 Icing
 @author Martin Lilleeng Saetra <martinls@met.no>

 $Id: Icing.h 1200 2014-11-17 22:08:14Z martinls $

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

#ifndef ICING_
#define ICING_

#include <cmath>
#include <cstddef>

using namespace std;

class Icing {
private:
	/// domain size
	int nx, ny;
		
	/// initializing constants
	float p0;
	float r;
	float cp;
	float rcp;
	float rho;
	///  g, acceleration of gravity
	float g;
	
	/// integration sum
	float* sum;
	
	/// special fields
	float* icingindexMax;
	float* icingindexMaxHeight;
	float* icingindexBottomGt4;
	float* icingindexTopGt4;
	
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
	 * Icing index A.
	 * @param cw Total amount of cloud water (g/kg)
	 * @param t Temperature (C)
	 */
	int A(float cw, float t);
	
	/**
	 * Icing index B.
	 * @param A Icing index A
	 * @param w Vertical velocity (cm/s)
	 */
	int B(int A, float w);
	
	/**
	 * Compute vertical velocity w from variables at model levels.
	 * NOTE: zsSum and sum must be initialized to 0.
	 * (Ported from a fortran77 routine written by Jan Erik Haugen.)
	 * @param k Vertical level index
	 * @param boundarySouth South boundary of model (degrees)
	 * @param dx Grid distance in x-direction (degrees)
	 * @param dy Grid distance in y-direction (degrees)
	 * @param af A-hybrid
	 * @param bf B-hybrid
	 * @param zsSum Used for vertical integration of z
	 * @param sum Used for vertical integration of w
	 * @param zs Surface geopotential (meters)
	 * @param ps Surface pressure (Pa)
	 * @param u Velocity component in x-direction (m/s)
	 * @param v Velocity component in y-direction (m/s)
	 * @param t Absolute temperature (K)
	 */
  float* initW(int k, float boundarySouth, float dx, float dy,
      float af, float bf, float *ah, float *bh, float *sum,
      float *z, float *ps, float *u, float *v, float *t);
	
public:
	/**
	 * Constructor.
	 */
	Icing();

	/**
	 * Destructor.
	 */
	virtual ~Icing();
	
	float* getIcingindexMaxField();
	float* getIcingindexMaxHeightField();
	float* getIcingindexBottomGt4Field();
	float* getIcingindexTopGt4Field();

	/**
	 * Reset integration sum between timesteps.
	 */
	void reset();
	
	/**
	 * Computes halflevel values for some parameter.
	 * @param fullLevelValues The existing fullevel values
	 * @param size Size of fullLevelValues array
	 * @param first Specify first element of halflevel array
	 * @param last Specify last element of halflevel array
	 * @return The computed halflevel values
	 */
	template<class T>
	T* calcHalfLevelValues(T* fullLevelValues, int size, T first=0, T last=0)
	{
	  T* halfLevelValues = new T[size];

	  halfLevelValues[0] = first;
	  halfLevelValues[size - 1] = last;

	  for (int i = size - 2; i > 0; --i) {
	    halfLevelValues[i] = 2. * fullLevelValues[i] - halfLevelValues[i + 1];
	  }

	  return halfLevelValues;
	}

	/**
	 * Calculates the icing indices of a 
	 * two-dimensional (spatial) field.
	 * @param nx width of 2D fields
	 * @param ny height of 2D fields
	 * @param ap hybrid coordinate
	 * @param b hybrid coordinate
	 * @param data Results are written to this array
	 * @param ps Surface pressure
	 * @param t Air temperature
	 * @param cw Cloud water
	 * @param w Upward air velocity
	 * @param p Pressure. (If not given it will be calculated from ps.)
	 */
	void calcIcingIndices2D(int nx, int ny, float ap, float b, float *data,
		    float *ps, float *t, float *cw, float *w,
		    float *p=NULL);
};

#endif /*ICING_*/
