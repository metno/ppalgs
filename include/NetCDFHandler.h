/**
 NetCDFHandler

 $Id: NetCDFHandler.h 1205 2014-12-11 16:03:33Z martinls $

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

#ifndef NETCDFHANDLER_
#define NETCDFHANDLER_

#include "FileHandler.h"

#include <string>
#include <memory>
#include <vector>
#include <netcdf>
#include <set>
#include <iterator>

#include <boost/shared_ptr.hpp>
#include <boost/shared_array.hpp>
#include <boost/date_time/posix_time/posix_time.hpp>

#include "fimex/CDMReader.h"
#include "fimex/coordSys/CoordinateSystem.h"
#include "fimex/CDM.h"

using namespace boost::gregorian;
namespace pt = boost::posix_time;

/**
 * Class for handling netCDF-files, with both reading and writing capabilities,
 * using netCDF API.
 *
 * NOTE: For some reason dataFile cannot be a pointer. This results in silent
 * errors, in which data is simply not written to the open netCDF-file.
 */
class NetCDFHandler : public FileHandler {
	public:
	NetCDFHandler(std::string _filename, netCDF::NcFile::FileMode filemode=netCDF::NcFile::read);
	virtual ~NetCDFHandler();

	virtual std::string getFilename() { return filename; }

	/**
	 * Initialize the file with correct attributes
	 */
	virtual void init(pt::ptime referenceTime, const std::vector<pt::ptime> times);

	/**
	 * Fetch all available times (in epoch format)
	 */
	virtual std::shared_ptr<std::vector<double> > getTimes(std::string variableName);

	/**
	 * Fetch all available levels
	 */
	virtual std::shared_ptr<std::vector<double>  > getLevels(std::string variableName);

	/**
	 * Write horizontal 2D level to netCDF file
	 * NOTE: Assumes one refTime!
	 * @param variableName The variable to write to
	 * @param size Size of data
	 * @param data The data to write
	 * @param _time The time to write
	 * @param _level The level to write
	 * @return True on success
	 */
	virtual bool writeSpatialGriddedLevel(std::string variableName, size_t size, std::vector<float>& data,
			double _time, double _level=-1);

	// PLACEHOLDERS
	virtual std::shared_ptr<std::vector<std::string>  > getVariables() { return std::shared_ptr<std::vector<std::string> >(); };
	virtual std::shared_ptr<std::vector<std::string> > getDimensions() { return std::shared_ptr<std::vector<std::string> >(); };
	virtual int getNx(std::string variableName) { return -1; };
	virtual int getNy(std::string variableName) { return -1; };
	virtual pt::ptime getRefTime() { return pt::ptime(); };
	virtual float readSingleValueLevel(std::string variableName, double _time, double _level) { return 0.0f; };
	virtual boost::shared_array<float> readSpatialGriddedLevel(std::string variableName, double _time, double _level=-1) { return boost::shared_array<float>(); };

	private:
	inline bool isClose(double a, double b, double epsilon = 1e-5) { return std::fabs(a - b) < epsilon;	};

	netCDF::NcFile::FileMode filemode;

	boost::shared_ptr<MetNoFimex::CDMReader> reader;
	MetNoFimex::CDM cdm;
	std::vector<boost::shared_ptr<const MetNoFimex::CoordinateSystem> > coordSys;
};

#endif /*NETCDFHANDLER_*/
