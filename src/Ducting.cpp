/**
 Ducting

 $Id: Ducting.cpp 1207 2015-02-24 20:54:56Z martinls $

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

#include "Ducting.h"
#include "utils.hpp"
#include "NetCDFHandler.h"
#include "GribHandler.h"

#define ONE_FEET 3.2808399
#define THOUSAND_FEET_IN_METERS 304.8
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

using namespace std;


Ducting::Ducting()
{
  /// initializing constants
  p0 = 1000.;
  cp = 1004.;
  r = 287.;
  rcp = r / cp;
  eps = 0.622;

  nx = 0;
  ny = 0;
}

Ducting::~Ducting()
{
  delete [] zm;
  delete [] Mm;
  delete [] inDuct;
  delete [] currentElevatedDuctBottom;
  delete [] currentElevatedDuctTop;

  delete [] dMdzMin;
  delete [] surfaceDuctBottom;
  delete [] surfaceDuctTop;
  delete [] elevatedDuctBottom;
  delete [] elevatedDuctTop;
  delete [] noOfElevatedDucts;
}

bool Ducting::initialized()
{
  return nx && ny;
}

void Ducting::initialize(int nx, int ny)
{
  this->nx = nx;
  this->ny = ny;

  int elements = nx * ny;

  zm = new float[elements];
  Mm = new float[elements];
  inDuct = new float[elements];
  currentElevatedDuctBottom = new float[elements];
  currentElevatedDuctTop = new float[elements];
  dMdzMin = new float[elements];
  surfaceDuctBottom = new float[elements];
  surfaceDuctTop = new float[elements];
  elevatedDuctBottom = new float[elements];
  elevatedDuctTop = new float[elements];
  noOfElevatedDucts = new float[elements];

  for (int i = 0; i < elements; ++i) {
    zm[i] = 0.;
    Mm[i] = 0.;
    inDuct[i] = 0.;
    currentElevatedDuctBottom[i] = numeric_limits<float>::max();
    currentElevatedDuctTop[i] = -numeric_limits<float>::max();
    dMdzMin[i] = numeric_limits<float>::max();
    surfaceDuctBottom[i] = numeric_limits<float>::max();
    surfaceDuctTop[i] = -numeric_limits<float>::max();
    elevatedDuctBottom[i] = numeric_limits<float>::max();
    elevatedDuctTop[i] = -numeric_limits<float>::max();
    noOfElevatedDucts[i] = 0.0f;
  }
}

float* Ducting::getdMdzMinField()
{
  return dMdzMin;
}

float* Ducting::getSurfaceDuctBottom()
{
  return surfaceDuctBottom;
}

float* Ducting::getSurfaceDuctTop()
{
  return surfaceDuctTop;
}

float* Ducting::getElevatedDuctBottom()
{
  return elevatedDuctBottom;
}

float* Ducting::getElevatedDuctTop()
{
  return elevatedDuctTop;
}

float* Ducting::getNoOfElevatedDucts()
{
  return noOfElevatedDucts;
}

void Ducting::reset()
{
  int elements = nx * ny;

  for (int i = 0; i < elements; ++i) {
    zm[i] = 0.;
    Mm[i] = 0.;
    inDuct[i] = 0.;
    currentElevatedDuctBottom[i] = numeric_limits<float>::max();
    currentElevatedDuctTop[i] = -numeric_limits<float>::max();
    dMdzMin[i] = numeric_limits<float>::max();
    surfaceDuctBottom[i] = numeric_limits<float>::max();
    surfaceDuctTop[i] = -numeric_limits<float>::max();
    elevatedDuctBottom[i] = numeric_limits<float>::max();
    elevatedDuctTop[i] = -numeric_limits<float>::max();
    noOfElevatedDucts[i] = 0.0f;
  }
}

float Ducting::dz(float z, int i, int j)
{
  return -z * 0.001 + zm[i + (nx * j)] * 0.001;
}

float Ducting::dMdz(int i, int j, float* data_m, float* data_z, float ps, float p, float t, float q)
{
  int ind = i + (nx * j);

  //DEBUG
  //if(i==j && i==100) cout << "Sfc pressure: " << ps << " Pressure: " << p << endl;

  float piValue = exner(p);
  float zValue = ::z(t, p, ps);
  data_z[ind] = zValue;
  float dzValue = dz(zValue, i, j);

  if (dzValue != 0) {
    float MValue = M(N(t, q, p), zValue);
    data_m[ind] = MValue;
    float dMdz = (-MValue + Mm[ind]) / dzValue;

    /// save lowest dMdz
    if (dMdz < dMdzMin[ind]) {
      dMdzMin[ind] = dMdz;
    }

    /// surface duct
    /// WARNING: Assumes TOA as first level
    if (zValue < THOUSAND_FEET_IN_METERS) {
    	if(dMdz < 0 && zValue < surfaceDuctBottom[ind]) {
    		surfaceDuctBottom[ind] = zValue;
    	}

    	if(dMdz < 0 && zValue > surfaceDuctTop[ind]) {
    		surfaceDuctTop[ind] = zValue;
    	}

	/// elevated ducts
	/// WARNING: Assumes TOA as first level
	/// entering duct or already inside duct
    } else { // zValue >= THOUSAND_FEET_IN_METERS
		if (dMdz < 0) {
			inDuct[ind] = 1.0;

			if(zValue > currentElevatedDuctTop[ind])
				currentElevatedDuctTop[ind] = zValue;
		/// exiting duct
		} else if (inDuct[ind]) {
			inDuct[ind] = 0.0;
			noOfElevatedDucts[ind]++;

			//if(zValue < currentElevatedDuctBottom[ind])
				currentElevatedDuctBottom[ind] = zValue;

			// set new elevatedDuct if currentElevatedDuct is thicker
			if(currentElevatedDuctTop[ind] - currentElevatedDuctBottom[ind] >
				elevatedDuctTop[ind] - elevatedDuctBottom[ind]) {
				elevatedDuctBottom[ind] = currentElevatedDuctBottom[ind];
				elevatedDuctTop[ind] = currentElevatedDuctTop[ind];
			}

			// reset currentElevatedDuct
			currentElevatedDuctBottom[ind] = numeric_limits<float>::max();
			currentElevatedDuctTop[ind] = -numeric_limits<float>::max();
		}
    }

    /// save current level as new previous level
    Mm[ind] = MValue;
    zm[ind] = zValue;

    return dMdz;
  } else {
    /// save current level as new previous level
    Mm[ind] = 0.;
    zm[ind] = 0.;

    data_m[ind] = 0.;

    return 0.;
  }
}

void Ducting::calcDuctingGrads2D(int nx, int ny, float ap, float b,
		float* data, float* data_m, float* data_z, float* ps, float* t, float* q, float* p)
{
  if (!initialized()) {
    initialize(nx, ny);
  }

#pragma omp parallel for
  for (int j = 0; j < ny; ++j) {
    for (int i = 0; i < nx; ++i) {
      int ind = i + (nx * j);
      float pVal = 0;
      if (p == NULL) {
        pVal = ap + b * (ps[ind]);
      } else {
        pVal = p[ind];
      }

      // (scale pressure Pa -> hPa)
      data[ind] = dMdz(i, j, data_m, data_z, ps[ind] * 0.01, pVal * 0.01, t[ind], q[ind]);
    }
  }
}
