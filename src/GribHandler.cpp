/**
 GribHandler

 $Id: GribHandler.cpp 1205 2014-12-11 16:03:33Z martinls $

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

#include "GribHandler.h"
#include <sstream>

#include "fimex/CoordinateSystemSliceBuilder.h"
#include "fimex/CDMFileReaderFactory.h"
#include "fimex/CDMVariable.h"
#include "fimex/CDMReaderUtils.h"
#include "fimex/XMLInput.h"
#include "fimex/DataDecl.h"
#include "fimex/Data.h"

using namespace MetNoFimex;
using namespace std;

GribHandler::GribHandler(string _filename, string config_filename) {
	filename = _filename;

	reader = CDMFileReaderFactory::create("grib", filename, XMLInputFile(config_filename));
	cdm = reader->getCDM();

	// get all coordinate systems from file, usually one, but may be a few (theoretical limit: # of variables)
	coordSys = listCoordinateSystems(reader);

	type = FileType::GRIB;
}

GribHandler::~GribHandler() { }

shared_ptr<vector<string> > GribHandler::getVariables() {
	shared_ptr<vector<string> > variables(new vector<string>());

	const vector<CDMVariable>& cdm_variables = cdm.getVariables();

	for(CDMVariable variable : cdm_variables)
		variables->push_back(variable. getName());

	return variables;
}

shared_ptr<vector<string> > GribHandler::getDimensions() {
	shared_ptr<vector<string> > dimensions(new vector<string>());

	const vector<CDMDimension>& cdm_dimensions = cdm.getDimensions();
	for(CDMDimension dimension : cdm_dimensions)
		dimensions->push_back(dimension.getName());

	return dimensions;
}

shared_ptr<vector<double> > GribHandler::getTimes(string variableName) {
	shared_ptr<vector<double> > times(new vector<double>());
	DataPtr cdm_times;

	// XXX: do some sanity checks (the variable must exist, time variable must exist, etc.)

	// find an appropriate coordinate system for a variable
	vector<boost::shared_ptr<const CoordinateSystem> >::iterator varSysIt =
			find_if(coordSys.begin(), coordSys.end(), CompleteCoordinateSystemForComparator(variableName));

	CoordinateSystem::ConstAxisPtr tAxis = (*varSysIt)->getTimeAxis(); // time

	if (tAxis.get() != 0) {
		cdm_times = tAxis->getData();
	}

	boost::shared_array<double> times_data = cdm_times->asDouble();

	for(int i=0; i < cdm_times->size(); ++i)
		times->push_back(times_data[i]);

	return times;
}

int GribHandler::getNx(string variableName)
{
	DataPtr nx;

	// XXX: do some sanity checks (the variable must exist, time variable must exist, etc.)

	// find an appropriate coordinate system for a variable
	vector<boost::shared_ptr<const CoordinateSystem> >::iterator varSysIt =
			find_if(coordSys.begin(), coordSys.end(), CompleteCoordinateSystemForComparator(variableName));

	CoordinateSystem::ConstAxisPtr xAxis = (*varSysIt)->getGeoXAxis(); // X or Lon

	if (xAxis.get() != 0) {
		nx = xAxis->getData();
	}

	return nx->size();
}

int GribHandler::getNy(string variableName)
{
	DataPtr ny;

	// XXX: do some sanity checks (the variable must exist, time variable must exist, etc.)

	// find an appropriate coordinate system for a variable
	vector<boost::shared_ptr<const CoordinateSystem> >::iterator varSysIt =
			find_if(coordSys.begin(), coordSys.end(), CompleteCoordinateSystemForComparator(variableName));

	CoordinateSystem::ConstAxisPtr yAxis = (*varSysIt)->getGeoYAxis(); // Y or Lat

	if (yAxis.get() != 0) {
		ny = yAxis->getData();
	}

	return ny->size();
}

pt::ptime GribHandler::getRefTime() {
	return getUniqueForecastReferenceTime(reader);
}

shared_ptr<vector<double> > GribHandler::getLevels(string variableName) {
	shared_ptr<vector<double> > levels(new vector<double>());
	DataPtr cdm_levels;

	// XXX: do some sanity checks (the variable must exist, time variable must exist, etc.)

	// find an appropriate coordinate system for a variable
	vector<boost::shared_ptr<const CoordinateSystem> >::iterator varSysIt =
			find_if(coordSys.begin(), coordSys.end(), CompleteCoordinateSystemForComparator(variableName));

	CoordinateSystem::ConstAxisPtr zAxis = (*varSysIt)->getGeoZAxis(); // Z or height (ml, pl, etc.)

	if (zAxis.get() != 0) {
		cdm_levels = zAxis->getData();
	}

	boost::shared_array<double> levels_data = cdm_levels->asDouble();

	for(int i=0; i < cdm_levels->size(); ++i)
		levels->push_back(levels_data[i]);

	return levels;
}

float GribHandler::readSingleValueLevel(std::string variableName, double _time, double _level)
{
	DataPtr data;
	int time_offset = 0;
	int level_offset = 0;

	// XXX: do some sanity checks (the variable must exist, the dimensions x,y,z, and time should exist, and there should be only one refTime)

	// find an appropriate coordinate system for a variable
	vector<boost::shared_ptr<const CoordinateSystem> >::iterator varSysIt =
			find_if(coordSys.begin(), coordSys.end(), CompleteCoordinateSystemForComparator(variableName));

	if (varSysIt != coordSys.end()) {
			CoordinateSystem::ConstAxisPtr zAxis = (*varSysIt)->getGeoZAxis(); // Z or height (ml, pl, etc.)
			CoordinateSystem::ConstAxisPtr tAxis = (*varSysIt)->getTimeAxis(); // time
			CoordinateSystemSliceBuilder sb(cdm, *varSysIt);

			// handling of time
			if (tAxis.get() != 0) {
				DataPtr times = reader->getDataSlice(tAxis->getName(),
						sb.getTimeVariableSliceBuilder());

				// find time offset
				boost::shared_array<double> times_data = times->asDouble();
				for(int i=0; i < times->size(); ++i)
					if(times_data[i] == _time)
						break;
					else
						++time_offset;

				sb.setTimeStartAndSize(time_offset, 1);
			}

			// find level offset
			DataPtr levels = zAxis->getData();
			boost::shared_array<double> levels_data = levels->asDouble();
			for(int i=0; i < levels->size(); ++i)
				if (levels_data[i] == _level)
					break;
				else
					++level_offset;

			sb.setStartAndSize(zAxis, level_offset, 1);

			// fetch the data
			data = reader->getDataSlice(variableName, sb);
	}

	return data->asFloat()[0];
}

boost::shared_array<float> GribHandler::readSpatialGriddedLevel(string variableName, double _time, double _level)
{
	DataPtr data;
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
					if(times_data[i] == _time)
						break;

				sb.setTimeStartAndSize(time_offset, 1);
			}
			// further selection of data
			sb.setAll(xAxis);
			sb.setAll(yAxis);

			// find level offset
			if(_level != -1) {
				DataPtr levels = zAxis->getData();
				boost::shared_array<double> levels_data = levels->asDouble();
				for(int i=0; i < levels->size(); ++i, level_offset++)
					if (levels_data[i] == _level)
						break;

				sb.setStartAndSize(zAxis, level_offset, 1);
			}

			// fetch the data
			data = reader->getDataSlice(variableName, sb);
		}
	}

	return data->asFloat();
}
