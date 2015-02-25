/**
 GribHandler

 $Id$

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

#ifndef GRIBHANDLER_
#define GRIBHANDLER_

#include "FileHandler.h"

#include <string>
#include <memory>
#include <vector>

#include <boost/shared_ptr.hpp>
#include <boost/shared_array.hpp>
#include <boost/date_time/posix_time/posix_time.hpp>

#include "fimex/CDMReader.h"
#include "fimex/coordSys/CoordinateSystem.h"
#include "fimex/CDM.h"

namespace pt = boost::posix_time;

/**
 * Class for reading Grib files using Fimex (no support for writing yet)
 */
class GribHandler : public FileHandler {
	public:
	GribHandler(std::string filename, std::string config_filename);
	virtual ~GribHandler();

	//PLACEHOLDERS
	virtual void init(pt::ptime referenceTime, const std::vector<pt::ptime> times) { };

	/**
	 * Fetch all variable names
	 */
	virtual std::shared_ptr<std::vector<std::string>  > getVariables();

	/**
	 * Fetch all dimension names
	 */
	virtual std::shared_ptr<std::vector<std::string> > getDimensions();

	/**
	 * Get width of variableName 2D field
	 */
	virtual int getNx(std::string variableName);

	/**
	 * Get height of variableName 2D field
	 * */
	virtual int getNy(std::string variableName);

	/**
	 * Get reference time
	 */
	virtual pt::ptime getRefTime();

	/**
	 * Fetch all available times (in epoch format)
	 */
	virtual std::shared_ptr<std::vector<double> > getTimes(std::string variableName);

	/**
	 * Fetch all available levels
	 */
	virtual std::shared_ptr<std::vector<double> > getLevels(std::string variableName);

	/**
	 * Read single-value level from grib file
	 * NOTE: Assumes one refTime!
	 * @param variableName The variable to read from
	 * @param time The time to read
	 * @param level The level to read
	 * @return level from variableName as a single value
	 */
	virtual float readSingleValueLevel(std::string variableName, double _time, double _level);

	/**
	 * Read horizontal 2D level from grib file
	 * NOTE: Assumes one refTime!
	 * @param variableName The variable to read from
	 * @param time The time to read
	 * @param level The level to read
	 * @return level from variableName as a 2D field
	 */
	virtual boost::shared_array<float> readSpatialGriddedLevel(std::string variableName, double _time, double _level);

	private:
	boost::shared_ptr<MetNoFimex::CDMReader> reader;
	MetNoFimex::CDM cdm;
	std::vector<boost::shared_ptr<const MetNoFimex::CoordinateSystem> > coordSys;
};

#endif /*GRIBHANDLER_*/
