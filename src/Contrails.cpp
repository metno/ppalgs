/**
 Contrails

 $Id: Contrails.cpp 1203 2014-11-27 17:05:25Z martinls $

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

#include "Contrails.h"
#include "utils.hpp"
#include "NetCDFHandler.h"
#include "GribHandler.h"

#define ONE_FEET 3.2808399
//#define NDEBUG

#include <assert.h>
#include <memory>
#include <set>
#include <vector>
#include <iostream>
#include <algorithm>
#include <fstream>
#include <pthread.h>
#include <time.h>
#include <netcdf>

#include <diField/diField.h>
#include <diField/diFieldManager.h>
#include <diField/diProjection.h>
#include <diField/diFieldSource.h>
#include <diField/diMetnoFieldFile.h>
#include <diField/diMetConstants.h>

#include <puTools/miCommandLine.h>
#include <puTools/miString.h>
#include <puTools/miTime.h>

#include <milib/milib.h>

using namespace std;

struct threadArg_t {
  int n;
  float *data;
  float *p;
  float *rh;
};

Contrails::Contrails()
{
  /// initializing constants
  p0 = 1000.;
  cp = 1004.;
  r = 287.05;
  rcp = r / cp;
  epsi = 0.6220;
  T0 = 273.16;

  a = -7.90298;
  b = 5.02808;
  c = -1.3816 * pow(10., -7);
  d = 11.344;
  f = 8.1328 * pow(10., -3);
  h = -3.49149;
  est = 1013.246;
  Ts = 373.16;

  nx = 0;
  ny = 0;
}

Contrails::~Contrails()
{ }

bool Contrails::initialized()
{
  return nx && ny;
}

void Contrails::initialize(int nx, int ny)
{
  this->nx = nx;
  this->ny = ny;

  int elements = nx * ny;

  sum = new float[elements];

  contrailsBottom = new float[elements];
  contrailsTop = new float[elements];

  for (int i = 0; i < elements; ++i) {
    sum[i] = 0.;
    contrailsBottom[i] = 0.;
    contrailsTop[i] = 0.;
  }
}

float* Contrails::getSumField()
{
  return sum;
}

float* Contrails::getContrailsTopField()
{
  return contrailsTop;
}

float* Contrails::getContrailsBottomField()
{
  return contrailsBottom;
}

void Contrails::reset()
{
  int elements = nx * ny;

  for (int i = 0; i < elements; ++i) {
    sum[i] = 0.;
    contrailsBottom[i] = 0.;
    contrailsTop[i] = 0.;
  }
}

float Contrails::esi(float t)
{
  float Z = a * ((Ts / t) - 1) + b * log10(Ts / t) +
      c * (pow(10, d*(1-(t/Ts))) - 1) + f * (pow(10, h*((Ts/t)-1)) - 1);

  return est * pow(10, Z);
}

float Tcrit(float pamb, float rhamb)
{
  // Physical constants used

  // Absolute T
  float t0 = 273.16;
  // Gas constant for dry air
  float rd = 287.05;
  // Gas constant for H2O
  float rh2o = 461.5;
  // Specific  heat capacity of air
  float cp = 1.005;
  // Boltzman's constant

  // liquid saturation
  float bal = 6.1121;
  float bbl = 18.729;
  float bcl = 257.87;
  float bdl = 227.3;

  // Heat liberated per gram fuel (Karcher [1994]??)
  float delh = 42.;
  // Water vapor mixing ratio in exhaust
  float wexh = 1.25;
  // dbg: adjust wexh
  // wexh = 1.3
  // fraction of combustion heat converted to propulsion
  float effic = 0.3;

  // Initialisation
  float tamax = -32.0;
  float tamin = -62.0;
  float tamb = 0.5*(tamin+tamax);
  float dTmin = 0.0;
  float dTmax = 100.0;
  float dT = 10.0;
  // No need to compute this in the inner loops
  float C1 = 1.e-3*pamb*wexh*cp*rh2o / (bal*delh*(1. - effic)*rd);
  float C2 = 1.0e3;  // scaling of equations, might be removed
  float rha1 = rhamb;
  if (rhamb >= 0.9999) { rha1 = 0.9999; }
  float bdli = 1.0 / bdl;
  // iteration control
  int mxiter = 10;
  int iter;
  float eps1 = 1.0e-8;
  float eps2 = 1.0e-4;

  // Find tamb and deltaT simultaneously by Newton iteration
  // by solving slplume(dT,tamb) - 1 = 0 and d(slplume)/d(dT) = 0 (for max.)
  // The iterations don't always converge
  // It works OK if we restrict to tamin < tamb < tamax and dTmin < dT < dTmax
  //  cout << "s1     s2      ta       dT     a11    a22    d1     d2" << endl;
  for (iter = 0; iter <= mxiter; ++iter) {
    float tp = tamb + dT;   // tplume
    float ga = (bbl - tamb * bdli) * tamb / (tamb + bcl);
    float gp = (bbl - tp * bdli) * tp / (tp + bcl);
    float fa = exp(ga) / (tamb + t0);
    float fp = exp(gp) / (tp + t0);
    float gda = (bbl*bcl - bdli*tamb*(tamb+2.*bcl))/((tamb+bcl)*(tamb+bcl));
    float gdp = (bbl*bcl - bdli*tp*(tp+2.*bcl))/((tp+bcl)*(tp+bcl));
    float ta0i = 1.0/(tamb+t0);
    float tp0i = 1.0/(tp+t0);
    float fda = fa*(gda - ta0i);
    float fdp = fp*(gdp - tp0i);
    float gddp = -2.0*(bdli + gdp)/(tp + bcl);
    float s2 = rha1*fa + C1*dT*ta0i - fp;  // s2 is slplume(dT,tamb) - 1
    float s1 = C2*(C1*ta0i - (s2 + fp)*(gdp - tp0i));  // s1 is d(slplume)/d(dT)
    float atmp = C2*(s2 + fp)*(gddp + tp0i*tp0i);
    float a11 = -C2*C1*ta0i*(gdp - tp0i) - atmp; // a11 is d�(slplume)/d(dT)�
    float a22 = rha1*fda - C1*dT*ta0i*ta0i - fdp;
    float a21 = C1*ta0i - fdp;
    float a12 = -C2*C1*ta0i*ta0i - C2*(a22 + fdp)*(gdp - tp0i) - atmp;
    // if a11 becomes >= 0, there is no max for deltaT in the plume
    if ( a11 > -eps1 ) {
      cerr << "a11 too big: " << a11 << " Abort!" << endl;
      cerr << "iter,a21,a12,a22,s1,s2=" << iter <<" "<< a21 <<" "<< a12 <<" "<< a22 <<" "<< s1 <<" "<< s2 << endl;
      cerr << "pamb: " << pamb << endl;
      cerr << "rhamb: " << rhamb << endl;
      cerr << "determinant: " << a11*a22 - a12*a21;
      exit(1);
    }
    // Solve equation system by gaussian elimination on the 2x2 Jacobian aij
    float qq = a21/a11;
    float dta = (qq*s1 - s2)/(a22 - qq*a12); // check determinant ??
    float ddT = (-s1 - a12*dta)/a11;
    //    float d1 = a11*ddT + a12*dta + s1;  //for debug, d1 should be 0
    //    float d2 = a21*ddT + a22*dta + s2;  //for debug, d2 should be 0
    // Restrict step length to be inside min,max box
    float alfa = 1.0;
    if ( tamb+alfa*dta < tamin ) {
      alfa = (tamin-tamb)/dta;
    } else if ( tamb+alfa*dta > tamax ) {
      alfa = (tamax-tamb)/dta;
    }
    if ( dT+alfa*ddT < dTmin ) {
      alfa = (dTmin-dT)/ddT;
    } else if ( dT+alfa*ddT > dTmax ) {
      alfa = (dTmax-dT)/ddT;
    }
    //    if ( alfa < 0.0 ) {  //debug only, should not happen
    //      cerr << "alfa became negative: " << alfa << " Abort!" << endl;
    //      exit(1);
    //    }
    // Update variables and check for convergence
    tamb = tamb + alfa*dta;
    dT = dT + alfa*ddT;
    //    cout << s1 <<" "<< s2 <<" "<< tamb <<" "<< dT <<" "<< a11 <<" "<< a22 <<" "<< d1 <<" "<< d2 << endl;
    if ( (dta*dta + ddT*ddT) < eps2 ) {
      break;
    }
  }
  float Tlc;
  if (iter < mxiter) {
    Tlc = tamb;
  } else {
    //cout << "pamb: " << pamb << endl;
    //cout << "rhamb: " << rhamb << endl;
    //    cout << "max iterations reached, tamb: " << tamb << endl;
    Tlc = tamin;
  }

  return Tlc;
}

void Contrails::calcContrails2D(int nx, int ny, float ap, float b, float *data,
	    float *ps, float *t, float *q, float *p)
{
  if (!initialized()) {
    initialize(nx, ny);
  }

  float *tCritValues = new float[nx * ny];

  for (int j = 0; j < ny; ++j) {
    for (int i = 0; i < nx; ++i) {
      float pVal = 0;
	  if (p == NULL) {
		pVal = ap + b * (ps[i + (nx * j)]);
	  } else {
		pVal = p[i + (nx * j)];
	  }

      float esiValue = esi(t[i + (nx * j)]);
      // (scale pressure Pa -> hPa)
      float qsValue = qs(esiValue, pVal * 0.01);
      float rhValue = RH(q[i + (nx * j)], qsValue);
      // (scale pressure Pa -> hPa)
      tCritValues[i + (nx * j)] = Tcrit(pVal * 0.01, rhValue);
    }
  }

  for (int j = 0; j < ny; ++j) {
    for (int i = 0; i < nx; ++i) {
      float pVal = 0;
      if (p == NULL) {
       	pVal = ap + b * (ps[i + (nx * j)]);
      } else {
       	pVal = p[i + (nx * j)];
      }

      // (scale pressure Pa -> hPa)
      float zValue = ::z(t[i + (nx * j)], pVal * 0.01, ps[i + (nx * j)] * 0.01);

      // K -> deg. C
      float tVal = t[i + (nx * j)] - 273.15;

      data[i + (nx * j)] = SchmidtAppleman(tCritValues[i + (nx * j)], tVal);

      // save this point/cell in the vertical sum-field as well
      if(data[i + (nx * j)] != 0) {
        sum[i + (nx * j)] = data[i + (nx * j)];

        // save maximum and minimum heights where contrails can occur
        // NOTE: Assumes layer 1 is TOA
        if (contrailsTop[i + (nx * j)] == 0) {
          contrailsTop[i + (nx * j)] = zValue*ONE_FEET; ///< convert to ft;
          // make sure bottom value is initialized (in case there are contrails in only one level/height)
          contrailsBottom[i + (nx * j)] = zValue*ONE_FEET; ///< convert to ft
        } else {
          contrailsBottom[i + (nx * j)] = zValue*ONE_FEET; ///< convert to ft
        }
      }

    }
  }

  delete[] tCritValues;
}
