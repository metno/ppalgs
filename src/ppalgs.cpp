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
 * Compute ducting gradients.
 *
 * Variables used:
 * 18/theta_potensiell_temp - air_potential_temperature - use air_temperature_ml, and compute potential temperature
 * 9/q_spesifikk_fuktighet - specific_humidity - specific_humidity_ml
 * 8/p_lufttrykk - ps (air_pressure_at_sea_level) - air_pressure_at_sea_level
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
	Ducting ducting;

	shared_ptr<vector<string> > variables = input->getVariables();
	cout << "Available variables: " << endl;
	for(string variable : *variables) {
		cout << variable << endl;
	}
	cout << endl;
	shared_ptr<vector<string> > dimensions = input->getDimensions();
	cout << "Available dimensions: " << endl;
	for(string dimension : *dimensions) {
		cout << dimension << endl;
	}
	cout << endl;

	shared_ptr<vector<double> > times = input->getTimes("specific_humidity_ml");
	vector<pt::ptime> ptimes;
	cout << "Available times:" << endl;
	for(int time : *times) {
		time_t tmp_time(time);
		pt::ptime boost_time = pt::from_time_t(tmp_time);
		ptimes.push_back(boost_time);
		cout << boost_time << ": ";
	}
	cout << endl << endl;

	shared_ptr<vector<double> > levels = input->getLevels("specific_humidity_ml");
	cout << "Available levels:" << endl;
	for (double level : *levels) {
		cout << level << ": ";
	}
	cout << endl << endl;

	cout << "nx: " << input->getNx("specific_humidity_ml") << endl
			<< "ny: " << input->getNy("specific_humidity_ml") << endl << endl;

	cout << "Reference time:" << endl;
	pt::ptime refTime = input->getRefTime();
	cout << refTime << endl;

	//FIXME: check that the number of vertical levels matches in input and output, and
	// that the horizontal dimensions match for all input and output 2D fields
	// fetch netcdf output levels (k)
	//shared_ptr<vector<double> > k(output.getLevels());

	/// iterate times
	int time_offset=0;
	for(double time : *times) {
		cout << "Processing time " << ptimes[time_offset] << endl;

		/// read ps
		boost::shared_array<float> ps = input->readSpatialGriddedLevel("air_pressure_at_sea_level", time, 0.0);

		int nx = input->getNx("air_pressure_at_sea_level");
		int ny = input->getNy("air_pressure_at_sea_level");
		int size = nx*ny;

		/// iterate levels
		int level_offset=0;
		for(double level : *levels) {
			cout << "Processing level " << level << endl;

			/// read fields/data
			boost::shared_array<float> t = input->readSpatialGriddedLevel("air_temperature_ml", time, level);
			boost::shared_array<float> q = input->readSpatialGriddedLevel("specific_humidity_ml", time, level);

			float ap = input->readSingleValueLevel("ap0", time, level);
			float b = input->readSingleValueLevel("b0", time, level);

			/// calculate ducting gradients and write gradients to file
			vector<float> data = vector<float>(size);

			ducting.calcDuctingGrads2D(nx, ny, ap, b, data.data(), ps.get(), t.get(), q.get());

			// write current time and level to file
			if (!initializedOutputFile) {
				output->init(refTime, ptimes);
				initializedOutputFile = true;
			}

			output->writeSpatialGriddedLevel("ducting_ml", time, level, size, data);

			++level_offset;

			break;
		}

		/// write minimum dMdz field
		// build fieldrequest for current time and level
		/*
		{
		FieldRequest fieldrequest;
		fieldrequest.modelName = "AROME-MetCoOp";
		fieldrequest.paramName = "ducting_sum";
		stringstream tmp_time;
		tmp_time.imbue(locale(cout.getloc(), facet));
		tmp_time << ptimes[time_offset];
		///FIXME: Remove miTime badness (e.g., by using fimex instead of diField(FieldManager))
		fieldrequest.ptime = miutil::miTime(tmp_time.str());
		stringstream ref_time;
		ref_time.imbue(locale(cout.getloc(), facet));
		ref_time << refTime;
		fieldrequest.refTime = ref_time.str();
		fieldrequest.taxis = "time";

		std::vector<std::string> modelConfigInfo;
		modelConfigInfo.push_back(
				"model=" + fieldrequest.modelName
						+ " t=fimex sourcetype=netcdf file=" + output->getFilename()
						+ " writeable=true");
		fieldManager->addModels(modelConfigInfo);

		Field* fieldW = 0;
		fieldManager->makeField(fieldW, fieldrequest);
		fieldW->nx = nx;
		fieldW->ny = ny;
		float *data = new float[size]; ///< no memleak! deleted by Field::cleanup()
		memcpy(data, ducting.getdMdzMinField(), size * sizeof(float));
		fieldW->data = data;
		fieldManager->writeField(fieldrequest, fieldW);
		delete fieldW;
		}
		*/

		/// write air_pressure_at_sea_level field


		// Will enable if needed...
		//    /// write heights of minimum dMdz field
		//    dMdzMinHeightField.shallowMemberCopy(*psField); ///< use metadata from psField
		//    dMdzMinHeightField.data = ducting.getdMdzMinHeightField();
		//
		//    output.writeFieldMetnoFormat(&dMdzMinHeightField, outputfilename, funit, 125);

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
 * 22 - cloud_liquid_water_content_of_atmosphere_layer - ga_76_253_251_109 (will be atmosphere_cloud_condensed_water_content_ml???)
 * 8 - ps (air_pressure_at_sea_level) - air_pressure_at_sea_level
 */
/*
void executeIcing(unique_ptr<FileHandler>&  input,
		unique_ptr<FileHandler>& output,
		const pt::ptime & startTime = pt::ptime(),
		const pt::ptime & stopTime = pt::ptime(),
		const int _startLevel = -1, const int _stopLevel = -1)
{
	// adjusting locale (to fit FieldManager (diField) and miutil (miTime) requirements)
	pt::time_facet *facet = new pt::time_facet("%Y-%m-%dT%H:%M:%S");
	cout.imbue(locale(cout.getloc(), facet));

	std::unique_ptr<FieldManager> fieldManager(new FieldManager());
	bool initializedOutputFile = false;
	Icing icing;

	shared_ptr<vector<double> > times = input->getTimes("specific_humidity_ml");
	vector<pt::ptime> ptimes;
	cout << "Available times:" << endl;
	for(int time : *times) {
		time_t tmp_time(time);
		pt::ptime boost_time = pt::from_time_t(tmp_time);
		ptimes.push_back(boost_time);
		cout << boost_time << ": ";
	}
	cout << endl << endl;

	shared_ptr<vector<double> > levels = input->getLevels("specific_humidity_ml");
	cout << "Available levels:" << endl;
	for (double level : *levels) {
		cout << level << ": ";
	}
	cout << endl << endl;

	cout << "nx: " << input->getNx("specific_humidity_ml") << endl
			<< "ny: " << input->getNy("specific_humidity_ml") << endl << endl;

	cout << "Reference time:" << endl;
	pt::ptime refTime = input->getRefTime();
	cout << refTime << endl;

	//FIXME: check that the number of vertical levels matches in input and output, and
	// that the horizontal dimensions match for all input and output 2D fields
	// fetch netcdf output levels (k)
	//shared_ptr<vector<double> > k(output.getLevels());

	int time_offset=0;
	for(double time : *times) {
		cout << "Processing time " << ptimes[time_offset] << endl;

		/// read ps
		boost::shared_array<float> ps = input->readSpatialGriddedLevel("air_pressure_at_sea_level", time, 0.0);

		int nx = input->getNx("air_pressure_at_sea_level");
		int ny = input->getNy("air_pressure_at_sea_level");
		int size = nx*ny;

		/// iterate levels
		int level_offset=0;
		for(double level : *levels) {
			cout << "Processing level " << level << endl;

			/// read fields/data
			boost::shared_array<float> t = input->readSpatialGriddedLevel("air_temperature_ml", time, level);
			boost::shared_array<float> cw = input->readSpatialGriddedLevel("ga_76_253_251_109", time, level);
			boost::shared_array<float> w = input->readSpatialGriddedLevel("upward_air_velocity_ml", time, level);

			float ap = input->readSingleValueLevel("ap0", time, level);
			float b = input->readSingleValueLevel("b0", time, level);

			/// calculate ducting gradients and write gradients to file
			vector<float> data = vector<float>(size);

			icing.calcIcingIndices2D(nx, ny, ap, b, data.data(), ps.get(), t.get(), cw.get(), w.get());

			// write current time and level to file
			FieldRequest fieldrequest;
			fieldrequest.modelName = "AROME-MetCoOp";
			fieldrequest.paramName = "icing_index_ml";
			stringstream tmp_time;
			tmp_time.imbue(locale(cout.getloc(), facet));
			tmp_time << ptimes[time_offset];
			///FIXME: Remove miTime badness (e.g., by using fimex instead of diField(FieldManager))
			fieldrequest.ptime = miutil::miTime(tmp_time.str());
			stringstream ref_time;
			ref_time.imbue(locale(cout.getloc(), facet));
			ref_time << refTime;
			fieldrequest.refTime = ref_time.str();
			fieldrequest.taxis = "time";
			fieldrequest.zaxis = "hybrid0";
			stringstream plevelSS;
			plevelSS << level;
			fieldrequest.plevel = plevelSS.str();

			if (!initializedOutputFile) {
				output->init(refTime, ptimes);
				initializedOutputFile = true;
			}
			std::vector<std::string> modelConfigInfo;
			modelConfigInfo.push_back(
					"model=" + fieldrequest.modelName
							+ " t=fimex sourcetype=netcdf file="
							+ output->getFilename() + " writeable=true");
			fieldManager->addModels(modelConfigInfo);

			Field* fieldW = 0;
			fieldManager->makeField(fieldW, fieldrequest);
			fieldW->nx = nx;
			fieldW->ny = ny;
			fieldW->data = data.data();
			fieldManager->writeField(fieldrequest, fieldW);
			//delete fieldW;

			++level_offset;
		}

		/// write icing_sum field
		{
		// build fieldrequest for current time and level
		FieldRequest fieldrequest;
		fieldrequest.modelName = "AROME-MetCoOp";
		fieldrequest.paramName = "icing_sum";
		stringstream tmp_time;
		tmp_time.imbue(locale(cout.getloc(), facet));
		tmp_time << ptimes[time_offset];
		///FIXME: Remove miTime badness (e.g., by using fimex instead of diField(FieldManager))
		fieldrequest.ptime = miutil::miTime(tmp_time.str());
		stringstream ref_time;
		ref_time.imbue(locale(cout.getloc(), facet));
		ref_time << refTime;
		fieldrequest.refTime = ref_time.str();
		fieldrequest.taxis = "time";

		std::vector<std::string> modelConfigInfo;
		modelConfigInfo.push_back(
				"model=" + fieldrequest.modelName
						+ " t=fimex sourcetype=netcdf file=" + output->getFilename()
						+ " writeable=true");
		fieldManager->addModels(modelConfigInfo);

		Field* fieldW = 0;
		fieldManager->makeField(fieldW, fieldrequest);
		fieldW->nx = nx;
		fieldW->ny = ny;
		float *data = new float[size]; ///< no memleak! deleted by Field::cleanup()
		memcpy(data, icing.getIcingindexMaxField(), size * sizeof(float));
		fieldW->data = data;
		fieldManager->writeField(fieldrequest, fieldW);
		delete fieldW;
		}

		/// write icing_height field
		{
		// build fieldrequest for current time and level
		FieldRequest fieldrequest;
		fieldrequest.modelName = "AROME-MetCoOp";
		fieldrequest.paramName = "icing_height";
		stringstream tmp_time;
		tmp_time.imbue(locale(cout.getloc(), facet));
		tmp_time << ptimes[time_offset];
		///FIXME: Remove miTime badness (e.g., by using fimex instead of diField(FieldManager))
		fieldrequest.ptime = miutil::miTime(tmp_time.str());
		stringstream ref_time;
		ref_time.imbue(locale(cout.getloc(), facet));
		ref_time << refTime;
		fieldrequest.refTime = ref_time.str();
		fieldrequest.taxis = "time";

		std::vector<std::string> modelConfigInfo;
		modelConfigInfo.push_back(
				"model=" + fieldrequest.modelName
						+ " t=fimex sourcetype=netcdf file=" + output->getFilename()
						+ " writeable=true");
		fieldManager->addModels(modelConfigInfo);

		Field* fieldW = 0;
		fieldManager->makeField(fieldW, fieldrequest);
		fieldW->nx = nx;
		fieldW->ny = ny;
		float *data = new float[size]; ///< no memleak! deleted by Field::cleanup()
		memcpy(data, icing.getIcingindexMaxHeightField(), size * sizeof(float));
		fieldW->data = data;
		fieldManager->writeField(fieldrequest, fieldW);
		delete fieldW;
		}

		/// write icing_height_bottom field
		// build fieldrequest for current time and level
		{
		FieldRequest fieldrequest;
		fieldrequest.modelName = "AROME-MetCoOp";
		fieldrequest.paramName = "icing_height_bottom";
		stringstream tmp_time;
		tmp_time.imbue(locale(cout.getloc(), facet));
		tmp_time << ptimes[time_offset];
		///FIXME: Remove miTime badness (e.g., by using fimex instead of diField(FieldManager))
		fieldrequest.ptime = miutil::miTime(tmp_time.str());
		stringstream ref_time;
		ref_time.imbue(locale(cout.getloc(), facet));
		ref_time << refTime;
		fieldrequest.refTime = ref_time.str();
		fieldrequest.taxis = "time";

		std::vector<std::string> modelConfigInfo;
		modelConfigInfo.push_back(
				"model=" + fieldrequest.modelName
						+ " t=fimex sourcetype=netcdf file=" + output->getFilename()
						+ " writeable=true");
		fieldManager->addModels(modelConfigInfo);

		Field* fieldW = 0;
		fieldManager->makeField(fieldW, fieldrequest);
		fieldW->nx = nx;
		fieldW->ny = ny;
		float *data = new float[size]; ///< no memleak! deleted by Field::cleanup()
		memcpy(data, icing.getIcingindexBottomGt4Field(), size * sizeof(float));
		fieldW->data = data;
		fieldManager->writeField(fieldrequest, fieldW);
		delete fieldW;
		}

		/// write icing_height_top field
		{
		// build fieldrequest for current time and level
		FieldRequest fieldrequest;
		fieldrequest.modelName = "AROME-MetCoOp";
		fieldrequest.paramName = "icing_height_top";
		stringstream tmp_time;
		tmp_time.imbue(locale(cout.getloc(), facet));
		tmp_time << ptimes[time_offset];
		///FIXME: Remove miTime badness (e.g., by using fimex instead of diField(FieldManager))
		fieldrequest.ptime = miutil::miTime(tmp_time.str());
		stringstream ref_time;
		ref_time.imbue(locale(cout.getloc(), facet));
		ref_time << refTime;
		fieldrequest.refTime = ref_time.str();
		fieldrequest.taxis = "time";

		std::vector<std::string> modelConfigInfo;
		modelConfigInfo.push_back(
				"model=" + fieldrequest.modelName
						+ " t=fimex sourcetype=netcdf file=" + output->getFilename()
						+ " writeable=true");
		fieldManager->addModels(modelConfigInfo);

		Field* fieldW = 0;
		fieldManager->makeField(fieldW, fieldrequest);
		fieldW->nx = nx;
		fieldW->ny = ny;
		float *data = new float[size]; ///< no memleak! deleted by Field::cleanup()
		memcpy(data, icing.getIcingindexTopGt4Field(), size * sizeof(float));
		fieldW->data = data;
		fieldManager->writeField(fieldrequest, fieldW);
		delete fieldW;
		}

		icing.reset();
		++time_offset;
	}
}
*/

/**
 * Compute contrails forecast.
 *
 * Variables used:
 * 18/theta_potensiell_temp - air_potential_temperature - use air_temperature_ml, and compute potential temperature
 * 9/q_spesifikk_fuktighet - specific_humidity - specific_humidity_ml
 * 8/p_lufttrykk - ps (air_pressure_at_sea_level) - air_pressure_at_sea_level
 */
/*
void executeContrails(unique_ptr<FileHandler>& input,
		unique_ptr<FileHandler>& output,
		const pt::ptime & startTime = pt::ptime(),
		const pt::ptime & stopTime = pt::ptime(),
		const int _startLevel = -1, const int _stopLevel = -1)
{
	// adjusting locale (to fit FieldManager (diField) and miutil (miTime) requirements)
	pt::time_facet *facet = new pt::time_facet("%Y-%m-%dT%H:%M:%S");
	cout.imbue(locale(cout.getloc(), facet));

	std::unique_ptr<FieldManager> fieldManager(new FieldManager());
	bool initializedOutputFile = false;
	Contrails contrails;

	shared_ptr<vector<double> > times = input->getTimes("specific_humidity_ml");
	vector<pt::ptime> ptimes;
	cout << "Available times:" << endl;
	for(int time : *times) {
		time_t tmp_time(time);
		pt::ptime boost_time = pt::from_time_t(tmp_time);
		ptimes.push_back(boost_time);
		cout << boost_time << ": ";
	}
	cout << endl << endl;

	shared_ptr<vector<double> > levels = input->getLevels("specific_humidity_ml");
	cout << "Available levels:" << endl;
	for (double level : *levels) {
		cout << level << ": ";
	}
	cout << endl << endl;

	cout << "nx: " << input->getNx("specific_humidity_ml") << endl
			<< "ny: " << input->getNy("specific_humidity_ml") << endl << endl;

	cout << "Reference time:" << endl;
	pt::ptime refTime = input->getRefTime();
	cout << refTime << endl;

	//FIXME: check that the number of vertical levels matches in input and output, and
	// that the horizontal dimensions match for all input and output 2D fields
	// fetch netcdf output levels (k)
	//shared_ptr<vector<double> > k(output.getLevels());

	int time_offset=0;
	for(double time : *times) {
		cout << "Processing time " << ptimes[time_offset] << endl;

		/// read ps
		boost::shared_array<float> ps = input->readSpatialGriddedLevel("air_pressure_at_sea_level", time, 0.0);

		int nx = input->getNx("air_pressure_at_sea_level");
		int ny = input->getNy("air_pressure_at_sea_level");
		int size = nx*ny;

		/// iterate levels
		int level_offset=0;
		for(double level : *levels) {
			cout << "Processing level " << level << endl;

			/// read fields/data
			boost::shared_array<float> t = input->readSpatialGriddedLevel("air_temperature_ml", time, level);
			boost::shared_array<float> q = input->readSpatialGriddedLevel("specific_humidity_ml", time, level);

			float ap = input->readSingleValueLevel("ap0", time, level);
			float b = input->readSingleValueLevel("b0", time, level);

			/// calculate ducting gradients and write gradients to file
			vector<float> data = vector<float>(size);

			contrails.calcContrails2D(nx, ny, ap, b, data.data(), ps.get(), t.get(), q.get());

			// write current time and level to file
			FieldRequest fieldrequest;
			fieldrequest.modelName = "AROME-MetCoOp";
			fieldrequest.paramName = "contrails_ml";
			stringstream tmp_time;
			tmp_time.imbue(locale(cout.getloc(), facet));
			tmp_time << ptimes[time_offset];
			///FIXME: Remove miTime badness (e.g., by using fimex instead of diField(FieldManager))
			fieldrequest.ptime = miutil::miTime(tmp_time.str());
			stringstream ref_time;
			ref_time.imbue(locale(cout.getloc(), facet));
			ref_time << refTime;
			fieldrequest.refTime = ref_time.str();
			fieldrequest.taxis = "time";
			fieldrequest.zaxis = "hybrid0";
			stringstream plevelSS;
			plevelSS << level;
			fieldrequest.plevel = plevelSS.str();

			if (!initializedOutputFile) {
				output->init(refTime, ptimes);
				initializedOutputFile = true;
			}
			std::vector<std::string> modelConfigInfo;
			modelConfigInfo.push_back(
					"model=" + fieldrequest.modelName
							+ " t=fimex sourcetype=netcdf file="
							+ output->getFilename() + " writeable=true");
			fieldManager->addModels(modelConfigInfo);

			Field* fieldW = 0;
			fieldManager->makeField(fieldW, fieldrequest);
			fieldW->nx = nx;
			fieldW->ny = ny;
			fieldW->data = data.data();
			fieldManager->writeField(fieldrequest, fieldW);
			//delete fieldW;

			++level_offset;
		}

		/// write contrails_sum field
		{
		// build fieldrequest for current time and level
		FieldRequest fieldrequest;
		fieldrequest.modelName = "AROME-MetCoOp";
		fieldrequest.paramName = "contrails_sum";
		stringstream tmp_time;
		tmp_time.imbue(locale(cout.getloc(), facet));
		tmp_time << ptimes[time_offset];
		///FIXME: Remove miTime badness (e.g., by using fimex instead of diField(FieldManager))
		fieldrequest.ptime = miutil::miTime(tmp_time.str());
		stringstream ref_time;
		ref_time.imbue(locale(cout.getloc(), facet));
		ref_time << refTime;
		fieldrequest.refTime = ref_time.str();
		fieldrequest.taxis = "time";

		std::vector<std::string> modelConfigInfo;
		modelConfigInfo.push_back(
				"model=" + fieldrequest.modelName
						+ " t=fimex sourcetype=netcdf file=" + output->getFilename()
						+ " writeable=true");
		fieldManager->addModels(modelConfigInfo);

		Field* fieldW = 0;
		fieldManager->makeField(fieldW, fieldrequest);
		fieldW->nx = nx;
		fieldW->ny = ny;
		float *data = new float[size]; ///< no memleak! deleted by Field::cleanup()
		memcpy(data, contrails.getSumField(), size * sizeof(float));
		fieldW->data = data;
		fieldManager->writeField(fieldrequest, fieldW);
		delete fieldW;
		}

		/// write contrails_bottom field
		// build fieldrequest for current time and level
		{
		FieldRequest fieldrequest;
		fieldrequest.modelName = "AROME-MetCoOp";
		fieldrequest.paramName = "contrails_bottom";
		stringstream tmp_time;
		tmp_time.imbue(locale(cout.getloc(), facet));
		tmp_time << ptimes[time_offset];
		///FIXME: Remove miTime badness (e.g., by using fimex instead of diField(FieldManager))
		fieldrequest.ptime = miutil::miTime(tmp_time.str());
		stringstream ref_time;
		ref_time.imbue(locale(cout.getloc(), facet));
		ref_time << refTime;
		fieldrequest.refTime = ref_time.str();
		fieldrequest.taxis = "time";

		std::vector<std::string> modelConfigInfo;
		modelConfigInfo.push_back(
				"model=" + fieldrequest.modelName
						+ " t=fimex sourcetype=netcdf file=" + output->getFilename()
						+ " writeable=true");
		fieldManager->addModels(modelConfigInfo);

		Field* fieldW = 0;
		fieldManager->makeField(fieldW, fieldrequest);
		fieldW->nx = nx;
		fieldW->ny = ny;
		float *data = new float[size]; ///< no memleak! deleted by Field::cleanup()
		memcpy(data, contrails.getContrailsBottomField(), size * sizeof(float));
		fieldW->data = data;
		fieldManager->writeField(fieldrequest, fieldW);
		delete fieldW;
		}

		/// write contrails_top field
		{
		// build fieldrequest for current time and level
		FieldRequest fieldrequest;
		fieldrequest.modelName = "AROME-MetCoOp";
		fieldrequest.paramName = "contrails_top";
		stringstream tmp_time;
		tmp_time.imbue(locale(cout.getloc(), facet));
		tmp_time << ptimes[time_offset];
		///FIXME: Remove miTime badness (e.g., by using fimex instead of diField(FieldManager))
		fieldrequest.ptime = miutil::miTime(tmp_time.str());
		stringstream ref_time;
		ref_time.imbue(locale(cout.getloc(), facet));
		ref_time << refTime;
		fieldrequest.refTime = ref_time.str();
		fieldrequest.taxis = "time";

		std::vector<std::string> modelConfigInfo;
		modelConfigInfo.push_back(
				"model=" + fieldrequest.modelName
						+ " t=fimex sourcetype=netcdf file=" + output->getFilename()
						+ " writeable=true");
		fieldManager->addModels(modelConfigInfo);

		Field* fieldW = 0;
		fieldManager->makeField(fieldW, fieldrequest);
		fieldW->nx = nx;
		fieldW->ny = ny;
		float *data = new float[size]; ///< no memleak! deleted by Field::cleanup()
		memcpy(data, contrails.getContrailsTopField(), size * sizeof(float));
		fieldW->data = data;
		fieldManager->writeField(fieldrequest, fieldW);
		delete fieldW;
		}

		contrails.reset();
		++time_offset;
	}
}
*/

int main(int argc, char* argv[]) {
	//FIXME: this should probably be a commandline option
	const std::string gribConfigFileName("/home/martinsa/tmp_metno/AromeMetCoOpGribReaderConfig.xml");

	string version_string = "1.0";
	int startLevel, stopLevel;
	startLevel = -1;
	stopLevel = -1;
	pt::ptime startTime, stopTime;

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
			("extra,e", po::value<vector<string> >(), "additional input-file");

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
	//	executeIcing(inputHandler, outputHandler, startTime, stopTime, startLevel, stopLevel);
	}
	if (vm.count("contrails")) {
	//	executeContrails(inputHandler, outputHandler, startTime, stopTime, startLevel, stopLevel);
	}
}
