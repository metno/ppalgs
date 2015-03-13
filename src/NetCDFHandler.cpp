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

#include "fimex/CoordinateSystemSliceBuilder.h"
#include "fimex/CDMFileReaderFactory.h"
#include "fimex/CDMVariable.h"
#include "fimex/CDMReaderUtils.h"
#include "fimex/XMLInput.h"
#include "fimex/DataDecl.h"
#include "fimex/Data.h"

#include "fimex/CDMReaderWriter.h"
#include "fimex/NetCDF_CDMReader.h"

using namespace std;

using namespace netCDF;
using namespace MetNoFimex;

NetCDFHandler::NetCDFHandler(string _filename, NcFile::FileMode filemode) :
		filemode(filemode)
{
	filename = _filename;

	reader = CDMFileReaderFactory::create("netcdf", filename);
	cdm = reader->getCDM();

	// get all coordinate systems from file, usually one, but may be a few (theoretical limit: # of variables)
	coordSys = listCoordinateSystems(reader);

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

bool NetCDFHandler::writeSpatialGriddedLevel(string variableName, double _time, double _level,
		size_t size, vector<float>& data)
{
	// refresh
	reader = CDMFileReaderFactory::create("netcdf", filename);
	cdm = reader->getCDM();

	// get all coordinate systems from file, usually one, but may be a few (theoretical limit: # of variables)
	coordSys = listCoordinateSystems(reader);

	int time_offset = 0;
	int level_offset = 0;

	// XXX: do some sanity checks (the variable must exist, the dimensions x,y,z, and time should exist, and there should be only one refTime)

	// find an appropriate coordinate system for a variable
	vector<boost::shared_ptr<const CoordinateSystem> >::iterator varSysIt =
			find_if(coordSys.begin(), coordSys.end(), CompleteCoordinateSystemForComparator(variableName));

	if (varSysIt != coordSys.end()) {
		if ((*varSysIt)->isSimpleSpatialGridded()) {
			CoordinateSystem::ConstAxisPtr xAxis = (*varSysIt)->getGeoXAxis(); // X or Lon
			CoordinateSystem::ConstAxisPtr yAxis = (*varSysIt)->getGeoYAxis(); // Y or Lat
			CoordinateSystem::ConstAxisPtr zAxis = (*varSysIt)->getGeoZAxis(); // Z or height (ml, pl, etc.)
			CoordinateSystem::ConstAxisPtr tAxis = (*varSysIt)->getTimeAxis(); // time
			CoordinateSystemSliceBuilder sb(cdm, *varSysIt);

			// handling of time
			if (tAxis.get() != 0) {
				DataPtr times = reader->getDataSlice(tAxis->getName(),
						sb.getTimeVariableSliceBuilder());

				// find time offset
				boost::shared_array<double> times_data = times->asDouble();
				for(int i=0; i < times->size(); ++i, time_offset++)
					if(times_data[i] == _time) {
						break;
					}

			}

			sb.setTimeStartAndSize(time_offset, 1);

			// further selection of data
			sb.setAll(xAxis);
			sb.setAll(yAxis);

			// find level offset
			DataPtr levels = reader->getData(zAxis->getName());
			boost::shared_array<double> levels_data = levels->asDouble();
			for(int i=0; i < levels->size(); ++i, level_offset++)
				if (isClose(levels_data[i], _level)) {
					break;
				}

			sb.setStartAndSize(zAxis, level_offset, 1);

			// write the data
			//boost::shared_ptr<CDMReaderWriter> writer = boost::dynamic_pointer_cast<CDMReaderWriter>(reader);
			boost::shared_ptr<CDMReaderWriter> writer = boost::shared_ptr<CDMReaderWriter>(new NetCDF_CDMReader(filename, true));

			DataPtr write_data = createData(CDMDataType::CDM_FLOAT, data.begin(), data.end());
			writer->putDataSlice(variableName, sb, write_data);
			writer->sync();

			return true;
		}
	}

	return false;
}
