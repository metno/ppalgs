/**
 ppalgs

 $Id: ppalgs.cpp 1207 2015-02-24 20:54:56Z martinls $

 Copyright (C) 2014 met.no

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

#include "FileHandler.h"
#include "GribHandler.h"
#include "NetCDFHandler.h"
#include "Ducting.h"
#include "Icing.h"
#include "Contrails.h"

#include <boost/date_time/posix_time/posix_time.hpp>
#include <boost/date_time/posix_time/posix_time_io.hpp>
#include <boost/date_time/gregorian_calendar.hpp>
#include <boost/program_options.hpp>

using namespace std;
using namespace boost::gregorian;

namespace po = boost::program_options;
namespace pt = boost::posix_time;

/**
 * Copy a spatial gridded field from one file to another.
 * @param input FileHandler for input file
 * @param output FileHandler for output file
 * @param variableName Name of field to copy
 * @param size Size (nx*ny) of field to copy
 * @param time Time step
 * @param level Level (defaults to -1; used for single-level variables)
 */
void copySpatialGriddedLevel(unique_ptr<FileHandler>& input,
		unique_ptr<FileHandler>& output, string variableName, int size,
		double time, double level=-1)
{
	boost::shared_array<float> tmp_data = input->readSpatialGriddedLevel(variableName, time, level);
	vector<float> data (tmp_data.get(), tmp_data.get() + size);

	output->writeSpatialGriddedLevel(variableName, size, data, time, level);
}

/**
 * Compute ducting gradients.
 *
 * Variables used:
 * 18/theta_potensiell_temp - air_potential_temperature - use air_temperature_ml, and compute potential temperature
 * 9/q_spesifikk_fuktighet - specific_humidity - specific_humidity_ml
 * 8/p_lufttrykk - ps (surface_air_pressure) - surface_air_pressure
 */
void executeDucting(unique_ptr<FileHandler>& input,
		unique_ptr<FileHandler>& output,
		const pt::ptime & startTime = pt::ptime(),
		const pt::ptime & stopTime = pt::ptime(),
		const int _startLevel = -1, const int _stopLevel = -1)
{
	// adjusting locale (to fit FieldManager (diField) and miutil (miTime) requirements)
	pt::time_facet *facet = new pt::time_facet("%Y-%m-%dT%H:%M:%S");
	cout.imbue(locale(cout.getloc(), facet));

	bool initializedOutputFile = false;
	bool verbose = false;
	Ducting ducting;

	shared_ptr<vector<string> > variables = input->getVariables();
	if (verbose) {
	  cout << "Available variables: " << endl;
	  for(string variable : *variables) {
	    cout << variable << endl;
	  }
	  cout << endl;
	}
	shared_ptr<vector<string> > dimensions = input->getDimensions();
	if (verbose) {
	  cout << "Available dimensions: " << endl;
	  for(string dimension : *dimensions) {
	    cout << dimension << endl;
	  }
	  cout << endl;
	}

	shared_ptr<vector<double> > times = input->getTimes("specific_humidity_ml");
	vector<pt::ptime> ptimes;
	if (verbose) cout << "Available times:" << endl;
	for(int time : *times) {
		time_t tmp_time(time);
		pt::ptime boost_time = pt::from_time_t(tmp_time);
		ptimes.push_back(boost_time);
		if (verbose) cout << boost_time << ": ";
	}
	if (verbose) cout << endl << endl;

	shared_ptr<vector<double> > levels = input->getLevels("specific_humidity_ml");
	if (verbose) {
	  cout << "Available levels:" << endl;
	  for (double level : *levels) {
	    cout << level << ": ";
	  }
	  cout << endl << endl;
	}

	cout << "nx: " << input->getNx("specific_humidity_ml") << endl
	     << "ny: " << input->getNy("specific_humidity_ml") << endl << endl;

	pt::ptime refTime = input->getRefTime();
	cout << "Reference time: " << refTime << endl;

	//FIXME: check that the number of vertical levels matches in input and output, and
	// that the horizontal dimensions match for all input and output 2D fields
	// fetch netcdf output levels (k)
	//shared_ptr<vector<double> > k(output.getLevels());

	/// iterate times
	int time_offset=0;
	for(double time : *times) {
		cout << "Processing time " << ptimes[time_offset] << endl;

		/// read ps
		boost::shared_array<float> ps = input->readSpatialGriddedLevel("surface_air_pressure", time, 0.0);

		int nx = input->getNx("surface_air_pressure");
		int ny = input->getNy("surface_air_pressure");
		int size = nx*ny;

		/// iterate levels
		int level_offset=0;
		for(double level : *levels) {
		        if (verbose) cout << "Processing level " << level << endl;

			/// read fields/data
			boost::shared_array<float> t = input->readSpatialGriddedLevel("air_temperature_ml", time, level);
			boost::shared_array<float> q = input->readSpatialGriddedLevel("specific_humidity_ml", time, level);

			string hnum = input->getHybridNumber();
			float ap = input->readSingleValueLevel("ap"+hnum, time, level);
			float b = input->readSingleValueLevel("b"+hnum, time, level);

			/// calculate ducting gradients and write gradients to file
			vector<float> data = vector<float>(size);
			vector<float> data_m = vector<float>(size);
			vector<float> data_z = vector<float>(size);

			ducting.calcDuctingGrads2D(nx, ny, ap, b, data.data(), data_m.data(), data_z.data(), ps.get(), t.get(), q.get());

			// write current time and level to file
			if (!initializedOutputFile) {
				output->init(refTime, ptimes);
				initializedOutputFile = true;
			}

			output->writeSpatialGriddedLevel("ducting_ml", size, data, time, level);
			output->writeSpatialGriddedLevel("ducting_m_ml", size, data_m, time, level);
			output->writeSpatialGriddedLevel("z_ml", size, data_z, time, level);

			++level_offset;
		}

		{
		float* dMdzMinField = ducting.getdMdzMinField();
		vector<float> data(&dMdzMinField[0], &dMdzMinField[size]);

		output->writeSpatialGriddedLevel("ducting_sum", size, data, time);
		}

		{
		float* surfaceDuctBottom = ducting.getSurfaceDuctBottom();
		vector<float> data(&surfaceDuctBottom[0], &surfaceDuctBottom[size]);

		output->writeSpatialGriddedLevel("ducting_surface_bottom", size, data, time);
		}

		{
		float* surfaceDuctTop = ducting.getSurfaceDuctTop();
		vector<float> data(&surfaceDuctTop[0], &surfaceDuctTop[size]);

		output->writeSpatialGriddedLevel("ducting_surface_top", size, data, time);
		}

		{
		float* elevatedDuctBottom = ducting.getElevatedDuctBottom();
		vector<float> data(&elevatedDuctBottom[0], &elevatedDuctBottom[size]);

		output->writeSpatialGriddedLevel("ducting_elevated_bottom", size, data, time);
		}

		{
		float* elevatedDuctTop = ducting.getElevatedDuctTop();
		vector<float> data(&elevatedDuctTop[0], &elevatedDuctTop[size]);

		output->writeSpatialGriddedLevel("ducting_elevated_top", size, data, time);
		}

		{
		float* noOfElevatedDucts = ducting.getNoOfElevatedDucts();
		vector<float> data(&noOfElevatedDucts[0], &noOfElevatedDucts[size]);

		output->writeSpatialGriddedLevel("ducting_no_of_elevated_ducts", size, data, time);
		}

		copySpatialGriddedLevel(input, output, "surface_air_pressure", size, time);


		ducting.reset();
		++time_offset;
	}
}

/**
 * Compute icing indices.
 *
 * Variables used:
 * 1 - geopotential_height - not needed, as vertical velocity w is available in AROME MetCoOp (as upward_air_velocity_ml)
 * 18 - air_potential_temperature - use air_temperature_ml, and compute potential temperature
 * 13 - omega (lagrangian_tendency_of_air_pressure) - not needed, as vertical velocity w is available in AROME MetCoOp (as upward_air_velocity_ml)
 * 22 - cloud_liquid_water_content_of_atmosphere_layer - atmosphere_cloud_condensed_water_content_ml
 * 8 - ps (surface_air_pressure) - surface_air_pressure
 */
void executeIcing(unique_ptr<FileHandler>&  input,
		unique_ptr<FileHandler>& output,
		const pt::ptime & startTime = pt::ptime(),
		const pt::ptime & stopTime = pt::ptime(),
		const int _startLevel = -1, const int _stopLevel = -1)
{
	// adjusting locale (to fit FieldManager (diField) and miutil (miTime) requirements)
	pt::time_facet *facet = new pt::time_facet("%Y-%m-%dT%H:%M:%S");
	cout.imbue(locale(cout.getloc(), facet));

	bool initializedOutputFile = false;
	bool verbose = false;
	Icing icing;

	shared_ptr<vector<double> > times = input->getTimes("air_temperature_ml");
	vector<pt::ptime> ptimes;
	for(int time : *times) {
		time_t tmp_time(time);
		pt::ptime boost_time = pt::from_time_t(tmp_time);
		ptimes.push_back(boost_time);
	}

	shared_ptr<vector<double> > levels = input->getLevels("air_temperature_ml");

	pt::ptime refTime = input->getRefTime();

	//FIXME: check that the number of vertical levels matches in input and output, and
	// that the horizontal dimensions match for all input and output 2D fields
	// fetch netcdf output levels (k)
	//shared_ptr<vector<double> > k(output.getLevels());

	int time_offset=0;
	for(double time : *times) {
		cout << "Processing time " << ptimes[time_offset] << endl;

		/// read ps
		boost::shared_array<float> ps = input->readSpatialGriddedLevel("surface_air_pressure", time, 0.0);

		int nx = input->getNx("surface_air_pressure");
		int ny = input->getNy("surface_air_pressure");
		int size = nx*ny;

		/// iterate levels
		int level_offset=0;
		for(double level : *levels) {
		        if (verbose) cout << "Processing level " << level << endl;

			/// read fields/data
			boost::shared_array<float> t = input->readSpatialGriddedLevel("air_temperature_ml", time, level);
			boost::shared_array<float> cw;
			try {
			  if (verbose) cout << "Reading atmosphere_cloud_condensed_water_content_ml" << endl;
			  cw = input->readSpatialGriddedLevel("atmosphere_cloud_condensed_water_content_ml", time, level);
			} catch(...) {
			  if (verbose) cout << "Reading mass_fraction_of_cloud_condensed_water_in_air_ml" << endl;
			  cw = input->readSpatialGriddedLevel("mass_fraction_of_cloud_condensed_water_in_air_ml", time, level);
			}
			boost::shared_array<float> w = input->readSpatialGriddedLevel("upward_air_velocity_ml", time, level);

			string hnum = input->getHybridNumber();
			float ap = input->readSingleValueLevel("ap"+hnum, time, level);
			float b = input->readSingleValueLevel("b"+hnum, time, level);

			/// calculate ducting gradients and write gradients to file
			vector<float> data = vector<float>(size);

			icing.calcIcingIndices2D(nx, ny, ap, b, data.data(), ps.get(), t.get(), cw.get(), w.get());

			// write current time and level to file
			if (!initializedOutputFile) {
				output->init(refTime, ptimes);
				initializedOutputFile = true;
			}

			output->writeSpatialGriddedLevel("icing_index_ml", size, data, time, level);

			++level_offset;
		}

		{
		float* icingindexMaxField = icing.getIcingindexMaxField();
		vector<float> data(&icingindexMaxField[0], &icingindexMaxField[size]);

		output->writeSpatialGriddedLevel("icing_sum", size, data, time);
		}

		{
		float* icingindexMaxHeightField = icing.getIcingindexMaxHeightField();
		vector<float> data(&icingindexMaxHeightField[0], &icingindexMaxHeightField[size]);

		output->writeSpatialGriddedLevel("icing_height", size, data, time);
		}

		{
		float* icingindexBottomGt4Field = icing.getIcingindexBottomGt4Field();
		vector<float> data(&icingindexBottomGt4Field[0], &icingindexBottomGt4Field[size]);

		output->writeSpatialGriddedLevel("icing_height_bottom", size, data, time);
		}

		{
		float* icingindexTopGt4Field = icing.getIcingindexTopGt4Field();
		vector<float> data(&icingindexTopGt4Field[0], &icingindexTopGt4Field[size]);

		output->writeSpatialGriddedLevel("icing_height_top", size, data, time);
		}

		icing.reset();
		++time_offset;
	}
}

/**
 * Compute contrails forecast.
 *
 * Variables used:
 * 18/theta_potensiell_temp - air_potential_temperature - use air_temperature_ml, and compute potential temperature
 * 9/q_spesifikk_fuktighet - specific_humidity - specific_humidity_ml
 * 8/p_lufttrykk - ps (surface_air_pressure) - surface_air_pressure
 */
void executeContrails(unique_ptr<FileHandler>& input,
		unique_ptr<FileHandler>& output,
		const pt::ptime & startTime = pt::ptime(),
		const pt::ptime & stopTime = pt::ptime(),
		const int _startLevel = -1, const int _stopLevel = -1)
{
	// adjusting locale (to fit FieldManager (diField) and miutil (miTime) requirements)
	pt::time_facet *facet = new pt::time_facet("%Y-%m-%dT%H:%M:%S");
	cout.imbue(locale(cout.getloc(), facet));

	bool initializedOutputFile = false;
	bool verbose = false;
	Contrails contrails;

	shared_ptr<vector<double> > times = input->getTimes("specific_humidity_ml");
	vector<pt::ptime> ptimes;
	for(int time : *times) {
		time_t tmp_time(time);
		pt::ptime boost_time = pt::from_time_t(tmp_time);
		ptimes.push_back(boost_time);
	}

	shared_ptr<vector<double> > levels = input->getLevels("specific_humidity_ml");

	pt::ptime refTime = input->getRefTime();

	//FIXME: check that the number of vertical levels matches in input and output, and
	// that the horizontal dimensions match for all input and output 2D fields
	// fetch netcdf output levels (k)
	//shared_ptr<vector<double> > k(output.getLevels());

	int time_offset=0;
	for(double time : *times) {
		cout << "Processing time " << ptimes[time_offset] << endl;

		/// read ps
		boost::shared_array<float> ps = input->readSpatialGriddedLevel("surface_air_pressure", time, 0.0);

		int nx = input->getNx("surface_air_pressure");
		int ny = input->getNy("surface_air_pressure");
		int size = nx*ny;

		/// iterate levels
		int level_offset=0;
		for(double level : *levels) {
		        if (verbose) cout << "Processing level " << level << endl;

			/// read fields/data
			boost::shared_array<float> t = input->readSpatialGriddedLevel("air_temperature_ml", time, level);
			boost::shared_array<float> q = input->readSpatialGriddedLevel("specific_humidity_ml", time, level);

			string hnum = input->getHybridNumber();
			float ap = input->readSingleValueLevel("ap"+hnum, time, level);
			float b = input->readSingleValueLevel("b"+hnum, time, level);

			/// calculate ducting gradients and write gradients to file
			vector<float> data = vector<float>(size);

			contrails.calcContrails2D(nx, ny, ap, b, data.data(), ps.get(), t.get(), q.get());

			// write current time and level to file
			if (!initializedOutputFile) {
				output->init(refTime, ptimes);
				initializedOutputFile = true;
			}

			output->writeSpatialGriddedLevel("contrails_ml", size, data, time, level);

			++level_offset;
		}

		{
		float* sumField = contrails.getSumField();
		vector<float> data(&sumField[0], &sumField[size]);

		output->writeSpatialGriddedLevel("contrails_sum", size, data, time);
		}

		{
		float* contrailsBottomField = contrails.getContrailsBottomField();
		vector<float> data(&contrailsBottomField[0], &contrailsBottomField[size]);

		output->writeSpatialGriddedLevel("contrails_bottom", size, data, time);
		}

		{
		float* contrailsTopField = contrails.getContrailsTopField();
		vector<float> data(&contrailsTopField[0], &contrailsTopField[size]);

		output->writeSpatialGriddedLevel("contrails_top", size, data, time);
		}

		contrails.reset();
		++time_offset;
	}
}

int main(int argc, char* argv[]) {
	string version_string = "1.1";
	int startLevel, stopLevel;
	startLevel = -1;
	stopLevel = -1;
	pt::ptime startTime, stopTime;
	string gribConfigFileName;

	// Declare the supported options.
	po::options_description cmdline_options;

	po::options_description generic("Options");
	generic.add_options()("help,h", "produce help message")
			("version,v", "produce version message")("ducting,d", "compute ducting gradients")
			("icing,i", "compute icing indices")
			("contrails,c",	"compute contrails indication")
			("levels,l", po::value<vector<int> >(), "compute only given interval of vertical levels (as integers): <from> <to>")
			("times,t", po::value<vector<string> >(), "compute only given time interval (as ISO datetime, e.g., 2014-08-18T17:50:11): <from> <to>")
			("input-type", po::value<string>(), "file format of input-file")
			("output-type", po::value<string>(), "file format of output-file")
			("extra,e", po::value<vector<string> >(), "additional input-file")
			("config", po::value<string>(), "grib config file");

	po::options_description hidden("Hidden options");
	hidden.add_options()
		("input-file", po::value<string>(), "input data file")
		("output-file", po::value<string>(), "output data file");

	po::positional_options_description p;
	p.add("input-file", 1);
	p.add("output-file", 1);

	cmdline_options.add(generic).add(hidden);

	po::variables_map vm;
	po::store(po::command_line_parser(argc, argv).options(cmdline_options).positional(p).run(), vm);
	po::notify(vm);

	if (vm.count("version")) {
		cout << version_string << endl;
		return 0;
	}
	if (vm.count("help")) {
		cout << "Usage: ppalgs <input-file> <output-file> --input-type=<grib|nc> --output-type=<nc> <option> [-e <additional input-file>]"
				<< endl << generic << endl;
		return 0;
	}
	if (argc < 5) {
		cout << "ERROR: Too few arguments!" << endl
				<< "Usage: ppalgs <input-file> <output-file> --input-type=<grib|nc> --output-type=<nc> <option> [-e <additional input-file>]"
				<< endl	<< generic << endl;
		return -1;
	}

	if (vm.count("time")) {
		vector<string> time_interval = vm["time"].as<vector<string> >();
		startTime = pt::from_iso_string(time_interval[0]);
		if (time_interval.size() > 1)
			stopTime = pt::from_iso_string(time_interval[1]);
		else
			stopTime = pt::from_iso_string(time_interval[0]);
	}
	if (vm.count("level")) {
		vector<int> level_interval = vm["level"].as<vector<int> >();
		startLevel = level_interval[0];
		if (level_interval.size() > 1)
			stopLevel = level_interval[1];
		else
			stopLevel = level_interval[0];
	}

	if (!vm.count("input-file") || !vm.count("output-file")) {
		cerr << "Input- and/or output filename(s) not specified, exiting."
				<< endl;
		return -2;
	}
	if (!vm.count("input-type") || !vm.count("output-type")) {
		cerr << "Input- and/or output file type(s) not specified, exiting."
				<< endl;
		return -3;
	}
	if (vm["input-type"].as<string>().compare("grib") == 0) {
		if(vm.count("config")) {
			gribConfigFileName = vm["config"].as<string>();
		} else {
			cerr << "Grib config file not given, exiting."
				<< endl;
			return -4;
		}
	}

	// handle input/output
	std::unique_ptr<FileHandler> inputHandler;
	std::unique_ptr<FileHandler> outputHandler;

	vector<string> inputfilenamev;
	inputfilenamev.push_back(vm["input-file"].as<string>());

	if (vm.count("extra")) { ///< if given, include extra inputfile
		cout << "Using extra inputfile: " << vm["extra"].as<string>() << endl;
		inputfilenamev.push_back(vm["extra"].as<string>());
	}

	if(vm["input-type"].as<string>().compare("grib") == 0) {
		inputHandler.reset(new GribHandler(inputfilenamev[0], gribConfigFileName));
	} else if(vm["input-type"].as<string>().compare("nc") == 0) {
		inputHandler.reset(new NetCDFHandler(inputfilenamev[0], netCDF::NcFile::read));
	} else {
		cout << "ERROR: Unsupported input type! Supported types: Grib and netCDF" << endl
				<< "Usage: ppalgs <input-file> <output-file> --input-type=<grib|nc> --output-type=<nc> <option> [-e <additional input-file>]"
				<< endl	<< generic << endl;
		return -4;
	}

	if(vm["output-type"].as<string>().compare("nc") == 0) {
		outputHandler.reset(new NetCDFHandler(vm["output-file"].as<string>(), netCDF::NcFile::write));
	} else {
		cout << "ERROR: Unsupported output type! Supported types: NetCDF" << endl
				<< "Usage: ppalgs <input-file> <output-file> --input-type=<grib|nc> --output-type=<nc> <option> [-e <additional input-file>]"
				<< endl	<< generic << endl;
		return -4;
	}

	// run chosen algorithm
	if (vm.count("ducting")) {
		executeDucting(inputHandler, outputHandler, startTime, stopTime, startLevel, stopLevel);
	}
	if (vm.count("icing")) {
		executeIcing(inputHandler, outputHandler, startTime, stopTime, startLevel, stopLevel);
	}
	if (vm.count("contrails")) {
		executeContrails(inputHandler, outputHandler, startTime, stopTime, startLevel, stopLevel);
	}
}
