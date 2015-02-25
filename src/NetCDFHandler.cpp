/**
 NetCDFHandler

 $Id: NetCDFHandler.cpp 1205 2014-12-11 16:03:33Z martinls $

 Copyright (C) 2014 MET Norway

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

#include "NetCDFHandler.h"
#include <sstream>

using namespace std;
using namespace netCDF;

NetCDFHandler::NetCDFHandler(string _filename, NcFile::FileMode filemode) :
		filemode(filemode)
{
	filename = _filename;
	type = FileType::NETCDF;
}

NetCDFHandler::~NetCDFHandler() { }

void NetCDFHandler::init(pt::ptime referenceTime, const vector<pt::ptime> times)
{
	NcFile dataFile(filename, filemode);
	if (dataFile.isNull()) {
		cout << "Couldn't open file!\n";
		exit(-1);
	}

	// adjusting locale (to fit template file)
	pt::time_facet *min_time_facet = new pt::time_facet("%Y-%m-%d %H:%M:%SZ");
	pt::time_facet *max_time_and_expires_facet = new pt::time_facet("%Y-%m-%d");
	stringstream minTimeAtt, maxTimeAtt, expireTimeAtt;
	minTimeAtt.imbue(locale(minTimeAtt.getloc(), min_time_facet));
	maxTimeAtt.imbue(locale(maxTimeAtt.getloc(), max_time_and_expires_facet));
	expireTimeAtt.imbue(locale(expireTimeAtt.getloc(), max_time_and_expires_facet));

	// prepare the attributes that need to be changed in the output nc-file
	minTimeAtt << times[0];
	maxTimeAtt << times[times.size()-1];

	// set expire time one month after last time step in file
	expireTimeAtt << times[times.size()-1] + boost::gregorian::months(1);

	try {
		// change the attributes, the reftime var, and add times (using netcdf c++ api)
		dataFile.putAtt("min_time", minTimeAtt.str());
		dataFile.putAtt("max_time", maxTimeAtt.str());
		dataFile.putAtt("Expires", expireTimeAtt.str());
	} catch (exceptions::NcException & e) {
		cout << "unknown error" << endl;
		cout << e.what() << endl;
	}

	// write reference time variable
	netCDF::NcVar refTimeVar = dataFile.getVar("forecast_reference_time");
	pt::ptime epoch(date(1970,Jan,1));
	pt::time_duration timeFromEpoch = referenceTime - epoch;
	double timeAsDouble = static_cast<double> (timeFromEpoch.total_seconds());
	refTimeVar.putVar(&timeAsDouble);

	// write time variable
	netCDF::NcVar timeVar = dataFile.getVar("time");
	int index = 0;
	vector<size_t> index_vec;
	index_vec.resize(1);
	for(pt::ptime time : times) {
		pt::time_duration timeFromEpoch = time - epoch;
		double timeAsDouble = static_cast<double> (timeFromEpoch.total_seconds());

		index_vec[0] = 0;
		timeVar.putVar(index_vec, timeAsDouble);

		++index;
	}
}

shared_ptr<vector<double> > NetCDFHandler::getTimes(string variableName) {
	NcFile dataFile(filename, filemode);
	if (dataFile.isNull()) {
		cout << "Couldn't open file!\n";
		exit(-1);
	}

	shared_ptr<vector<double> > times(new vector<double>());

	netCDF::NcVar var = dataFile.getVar(variableName);

	int dims = var.getDimCount();

	int dim = 0;
	for(int i=0; i < dims; ++i) {
		// NOTE: Assumes time dimension is unlimited!
		if (var.getDim(i).isUnlimited()) {
			dim = i;
			break;
		}
	}

	times->resize(var.getDim(dim).getSize());
	var.getVar(&times->at(0));

	return times;
}

shared_ptr<vector<double> > NetCDFHandler::getLevels(string variableName) {
	NcFile dataFile(filename, filemode);
	if (dataFile.isNull()) {
		cout << "Couldn't open file!\n";
		exit(-1);
	}

	shared_ptr<vector<double> > levels(new vector<double>());

	netCDF::NcVar var = dataFile.getVar(variableName);

	int dims = var.getDimCount();

	int dim = 0;
	for(int i=0; i < dims; ++i) {
		// NOTE: Assumes height/z dimension is named "hybrid0"!
		if (var.getDim(i).getName().compare("hybrid0") == 0) {
			dim = i;
			break;
		}
	}

	levels->resize(var.getDim(dim).getSize());
	var.getVar(&levels->at(0));

	return levels;
}
