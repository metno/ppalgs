/**
 FileHandler

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

#ifndef FILEHANDLER_
#define FILEHANDLER_

#include <string>
#include <memory>
#include <vector>

#include <boost/shared_ptr.hpp>
#include <boost/shared_array.hpp>
#include <boost/date_time/posix_time/posix_time.hpp>

enum class FileType {
	GRIB,
	NETCDF
};

namespace pt = boost::posix_time;

/**
 * Superclass for file reading using Fimex
 */
class FileHandler {
	public:
	virtual ~FileHandler() { };

	std::string getFilename() { return filename; }
	FileType getFileType() { return type; }

	/**
	 * Initialize the file with correct attributes
	 */
	virtual void init(pt::ptime referenceTime, const std::vector<pt::ptime> times) = 0;

	/**
	 * Fetch all variable names
	 */
	virtual std::shared_ptr<std::vector<std::string>  > getVariables() = 0;

	/**
	 * Fetch all dimension names
	 */
	virtual std::shared_ptr<std::vector<std::string> > getDimensions() = 0;

	/**
	 * Get width of variableName 2D field
	 */
	virtual int getNx(std::string variableName) = 0;

	/**
	 * Get height of variableName 2D field
	 * */
	virtual int getNy(std::string variableName) = 0;

	/**
	 * Get reference time
	 */
	virtual pt::ptime getRefTime() = 0;

	/**
	 * Fetch all available times (in epoch format)
	 */
	virtual std::shared_ptr<std::vector<double> > getTimes(std::string variableName) = 0;

	/**
	 * Fetch all available levels
	 */
	virtual std::shared_ptr<std::vector<double> > getLevels(std::string variableName) = 0;

	/**
	 * Read single-value level from grib file
	 * NOTE: Assumes one refTime!
	 * @param variableName The variable to read from
	 * @param time The time to read
	 * @param level The level to read
	 * @return level from variableName as a single value
	 */
	virtual float readSingleValueLevel(std::string variableName, double _time, double _level) = 0;

	/**
	 * Read horizontal 2D level from grib file
	 * NOTE: Assumes one refTime!
	 * @param variableName The variable to read from
	 * @param time The time to read
	 * @param level The level to read
	 * @return level from variableName as a 2D field
	 */
	virtual boost::shared_array<float> readSpatialGriddedLevel(std::string variableName, double _time, double _level) = 0;

	protected:
	FileType type;
	std::string filename;
};

#endif /*FILEHANDLER_*/
