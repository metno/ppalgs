/**
 Icing

 $Id: Icing.cpp 1200 2014-11-17 22:08:14Z martinls $

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

#include "Icing.h"
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

using namespace std;

Icing::Icing()
{
  /// initializing constants
  p0 = 1000.;
  cp = 1004.;
  r = 287.;
  rcp = r / cp;
  rho = 1.293;

  ///  g, acceleration of gravity
  g = 9.8;

  nx = 0;
  ny = 0;
}

Icing::~Icing()
{
  delete[] sum;
  delete[] icingindexMax;
  delete[] icingindexMaxHeight;
  delete[] icingindexBottomGt4;
  delete[] icingindexTopGt4;
}

bool Icing::initialized()
{
  return nx && ny;
}

void Icing::initialize(int nx, int ny)
{
  this->nx = nx;
  this->ny = ny;

  sum = new float[nx * ny]; ///< used in initW
  icingindexMax = new float[nx * ny];
  icingindexMaxHeight = new float[nx * ny];
  icingindexBottomGt4 = new float[nx * ny];
  icingindexTopGt4 = new float[nx * ny];

  for (int i = 0; i < nx * ny; ++i) {
    sum[i] = 0.;
    icingindexMax[i] = 0.;
    icingindexMaxHeight[i] = 0.;
    icingindexBottomGt4[i] = 0.;
    icingindexTopGt4[i] = 0.;
  }
}

float* Icing::getIcingindexMaxField()
{
  return icingindexMax;
}

float* Icing::getIcingindexMaxHeightField()
{
  return icingindexMaxHeight;
}

float* Icing::getIcingindexBottomGt4Field()
{
  return icingindexBottomGt4;
}

float* Icing::getIcingindexTopGt4Field()
{
  return icingindexTopGt4;
}

void Icing::reset()
{
  int elements = nx * ny;

  for (int i = 0; i < elements; ++i) {
    sum[i] = 0.;
    icingindexMax[i] = 0.;
    icingindexMaxHeight[i] = 0.;
    icingindexBottomGt4[i] = 0.;
    icingindexTopGt4[i] = 0.;
  }
}

template <typename T>
bool isCloseToZero(T x)
{
    return std::abs(x) < std::numeric_limits<T>::epsilon();
}

int Icing::A(float cw, float t)
{
  // avoid log(0) for very small values of cw
  if(isCloseToZero(cw))
	  cw = 0.f;

  float A;
  float C = t - 273.15;

  if (cw != 0 && C < 0) {
    if (C >= -15) {
      A = 5 + log(cw);
    } else { /// t < -15
      A = 5 + log(cw * (40 + C) / 25);
    }
  } else {
    A = -99;
  }

  return static_cast<int> (round(A));
}

int Icing::B(int A, float w)
{
  int B;

  if (A != -99) {
    if (w > 30) {
      B = A + 4;
    } else if (w > 20) {
      B = A + 3;
    } else if (w > 10) {
      B = A + 2;
    } else if (w >= 0) {
      B = A + 1;
    } else { ///< w < 0
      B = A - 1;
    }
  } else {
    B = 0;
  }

  if (B < 0) {
    return 0;
  } else {
    return B;
  }
}

float* Icing::initW(int k, float boundarySouth, float dx, float dy, float af,
    float bf, float *ah, float *bh, float *sum, float *z, float *ps, float *u,
    float *v, float *t)
{
  /*cout << endl << "k: " << k << " boundarySouth: " << boundarySouth << " dx: " << dx << " dy: " << dy
   << " af: " << af << " bf: " << bf << " ah: " << ah[10] << " bh: " << bh[10] << " zs: " << zs[10]
   << " ps: " << ps[10] << " u: " << u[10] << " v: " << v[10] << " t: " << t[10] << endl;*/

  /// allocate arrays
  int nelements = nx * ny;

  float *w = new float[nelements];

  float *hx = new float[nelements];
  float *rhx = new float[nelements];
  float *hy = new float[nelements];
  float *rhy = new float[nelements];
  float *rhxy = new float[nelements];
  float *dp = new float[nelements];
  float *dlnp = new float[nelements];
  float *alfa = new float[nelements];
  float *uu = new float[nelements];
  float *vv = new float[nelements];

  /// earth radius
  float a = 6.37e+6;

  float pir180 = 4. * atanf(1.) / 180.;
  float y1r = boundarySouth * pir180;
  float dxr = dx * pir180;
  float dyr = dy * pir180;
  float zrdx2 = 1. / (2. * dxr);
  float zrdy2 = 1. / (2. * dyr);

  /// compute mapfactors for computation of derivatives
  for (int j = 0; j < ny; ++j) {
    //float cosy=cos(static_cast<float>(j)*dyr); ///< no displacement (y1r) of grid
    float cosy = cos(y1r + static_cast<float> (j) * dyr);

    for (int i = 0; i < nx; ++i) {
      int ind = i + (nx * j);

      hx[ind] = a * cosy;
      rhx[ind] = 1. / hx[ind];
      hy[ind] = a;
      rhy[ind] = 1. / hy[ind];
      rhxy[ind] = rhx[ind] * rhy[ind];
    }
  }

  /// compute pressure variables needed only once
  float ln2 = logf(2.);

  float da = ah[1] - ah[0];
  float db = bh[1] - bh[0];

  if (k == 0) {
    for (int j = 0; j < ny; ++j) {
      for (int i = 0; i < nx; ++i) {
        int ind = i + (nx * j);

        dp[ind] = da + db * ps[ind];

        dlnp[ind] = 0.;
        alfa[ind] = ln2;
      }
    }
  } else {
    da = ah[k + 1] - ah[k];
    db = bh[k + 1] - bh[k];

    for (int j = 0; j < ny; ++j) {
      for (int i = 0; i < nx; ++i) {
        int ind = i + (nx * j);

        float pm = ah[k] + bh[k] * ps[ind];
        float pp = ah[k + 1] + bh[k + 1] * ps[ind];

        dp[ind] = da + db * ps[ind];
        dlnp[ind] = logf(pp / pm);
        alfa[ind] = 1. - pm * dlnp[ind] / dp[ind];
      }
    }
  }

  /// vertical integral of divergence, compute w
  for (int j = 0; j < ny; ++j) {
    for (int i = 0; i < nx; ++i) {
      int ind = i + (nx * j);

      uu[ind] = hy[ind] * u[ind] * dp[ind];
      vv[ind] = hx[ind] * v[ind] * dp[ind];
    }
  }

  for (int j = 1; j < ny - 1; ++j) {
    for (int i = 1; i < nx - 1; ++i) {
      int ind = i + (nx * j);

      float div = rhxy[ind] * (zrdx2 * (uu[ind + 1] - uu[ind - 1])
	  + zrdy2 * (vv[ind + nx] - vv[ind - nx]));

      float w1 = r * t[ind] * (dlnp[ind] * sum[ind] + alfa[ind] * div)
          / dp[ind];

      float w2 = rhx[ind] * zrdx2 * (z[ind + 1] - z[ind - 1])
          + rhy[ind] * zrdy2 * (z[ind + nx] - z[ind - nx]);
      w[ind] = (w1 + w2) / g;
      sum[ind] = sum[ind] + div;
    }
  }

  for (int i = 1; i < nx - 1; ++i) {
    w[i] = w[i + nx];
    w[i + (nx * (ny - 1))] = w[i + (nx * (ny - 2))];
  }
  for (int j = 0; j < ny; ++j) {
    w[nx * j] = w[1 + (nx * j)];
    w[nx - 1 + (nx * j)] = w[nx - 2 + (nx * j)];
  }

  /// clean up
  delete[] hx;
  delete[] rhx;
  delete[] hy;
  delete[] rhy;
  delete[] rhxy;
  delete[] dp;
  delete[] dlnp;
  delete[] alfa;
  delete[] uu;
  delete[] vv;

  return w;
}

void Icing::calcIcingIndices2D(int nx, int ny, float ap, float b, float *data,
    float *ps, float *t, float *cw, float *w, float *p)
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
      float zValue = ::z(t[ind], pVal * 0.01, ps[ind] * 0.01);
      ///FIXME: double-check conversions!
      int AValue = A(1000 * cw[ind], t[ind]); ///< convert kg/kg -> g/kg

      data[ind] = B(AValue, w[ind] * 100); ///< convert m/s -> cm/s

      /// save highest icingindex and current height
      if (data[ind] > icingindexMax[ind]) {
        icingindexMax[ind] = data[ind];
        icingindexMaxHeight[ind] = zValue*ONE_FEET; ///< convert to ft
      }
      /// save bottom and top layer with icingindex > 4
      /// NOTE: Assumes layer 1 is TOA
      if (data[ind] > 4) {
        if (icingindexTopGt4[ind] == 0) {
          icingindexTopGt4[ind] = zValue*ONE_FEET; ///< convert to ft
          // make sure bottom value is initialized (in case there are icing only in one level/height)
          icingindexBottomGt4[ind] = zValue*ONE_FEET; ///< convert to ft
        } else {
          icingindexBottomGt4[ind] = zValue*ONE_FEET; ///< convert to ft
        }
      }
    }
  }
}
