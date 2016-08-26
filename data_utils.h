#ifndef DATA_UTILS_H
#define DATA_UTILS_H
//Include database dependencies
//QT
#include "qtextstream.h"
#include "qstring.h"
#include "qvector.h"
#include "qevent.h"
#include "qinputdialog.h"
#include "qmessagebox.h"
#include "qfiledialog.h"
//STD
#include <vector>
#include <fstream>
#include <sstream>
#include <iostream>
#include <iomanip>
#include <time.h>       /* time */
//Boost
#include <boost/make_shared.hpp>
#include <boost/scoped_ptr.hpp>
#include <boost/thread.hpp>
#include <boost/log/sources/logger.hpp>
#include <boost/log/sources/record_ostream.hpp>
#include <boost/math/special_functions/sinc.hpp>
#include <boost/date_time/posix_time/posix_time.hpp>
#include "boost/random.hpp"
#include "boost/generator_iterator.hpp"
//Self-Written
#include "data_def.h"

namespace saxs
{
	/***************************************************************************************************************/
	// GUI Messaging
	/***************************************************************************************************************/

	//	File dialog box
	QString filedialog(QString dialogTitle, QString FileFilter)
	{
		QFileDialog dialog;
		dialog.setNameFilter(FileFilter);
		dialog.setFileMode(QFileDialog::ExistingFile);
		dialog.setWindowTitle(dialogTitle);
		dialog.setViewMode(QFileDialog::Detail);
		dialog.setOption(QFileDialog::DontUseNativeDialog, true);
		QStringList fileNames;
		if (dialog.exec())
		{
			fileNames = dialog.selectedFiles();
		}
		else
		{
			fileNames.append("none");
		}
		return fileNames[0];
	}

	//	Multiple Files dialog box
	QStringList multifiledialog(QString dialogTitle, QString FileFilter)
	{
		QFileDialog dialog;
		dialog.setNameFilter(FileFilter);
		dialog.setFileMode(QFileDialog::ExistingFiles);
		dialog.setWindowTitle(dialogTitle);
		dialog.setViewMode(QFileDialog::Detail);
		dialog.setOption(QFileDialog::DontUseNativeDialog, true);
		QStringList fileNames;
		if (dialog.exec())
		{
			fileNames = dialog.selectedFiles();
		}
		else
		{
			fileNames.append("none");
		}
		return fileNames;
	}

	//	Error Message box
	void error_message(QString errormessage)
	{
		//Predefine Messagebox for error warnings
		QMessageBox msgBox;
		msgBox.setText("Error");
		msgBox.setStandardButtons(QMessageBox::Ok);
		msgBox.setIcon(QMessageBox::Critical);

		msgBox.setInformativeText(errormessage);
		int ret = msgBox.exec();

	}

	//	Integer inputbox
	int input_message(QString text, int curr, int min, int max)
	{
		//Predefine Messagebox for error warnings
		bool ok;
		int input = -1;
		input = QInputDialog::getInt(NULL, "Specify..", text, curr, min, max, 1, &ok);

		if (ok && input > 0) return input;
		else return curr;
	}

	//	Double inputbox
	double input_message_doub(QString text, double curr, double min, double max)
	{
		//Predefine Messagebox for error warnings
		bool ok;
		double input = -1;
		input = QInputDialog::getDouble(NULL, "Specify..", text, curr, min, max, 1, &ok);

		if (ok && input > 0) return input;
		else return curr;
	}




	/***************************************************************************************************************/
	// Data Handling functions
	/***************************************************************************************************************/

	//------------------------------------------------------------------------------
	//	Data Import
	//------------------------------------------------------------------------------

	//	Import scattering data from file
	bool import_scatteringdata(const QString file_name, std::vector<scatteringdata_sp>& data)
	{
		//transform qstring into string
		std::string file_name_string = file_name.toStdString();
		std::string file_ext = file_name_string.substr(file_name_string.find_last_of(".") + 1);

		//Predefine Messagebox for error warnings
		QMessageBox msgBox;
		msgBox.setText("Error");
		msgBox.setStandardButtons(QMessageBox::Ok);
		msgBox.setIcon(QMessageBox::Critical);

		//Temporary line used during import
		std::string line;

		//Open file
		std::ifstream infile(file_name_string.c_str(), std::ifstream::in);
		//Check if correctly open
		if (!infile.is_open())
		{
			saxs::error_message("Could not open file " + file_name + "!");
			return false;
		}

		//Rule for qI
		if (file_ext=="qI")
		{
			//Check number of columns to read error
			getline(infile, line);
			std::istringstream iss(line);
			int columns = 0;
			do
			{
				std::string sub;
				iss >> sub;
				if (sub.length())
					++columns;
			} while (iss);
			double error_scale = 0;
			bool error_specified;
			switch (columns)
			{
			case 0: saxs::error_message("No data in file:" + file_name + "!");
				infile.close();
				return false;
			case 1: saxs::error_message("No scattering data in file:" + file_name + "!");
				infile.close();
				return false;
			case 2: //error_scale = saxs::input_message_doub("No error given in file "+ file_name +"\nSpecify artificial error weight in %:",3,0.0001,100);
				error_specified = false;
				break;
			case 3: error_specified = true;
				break;
			}

			//Go back to line 0
			infile.clear();
			infile.seekg(0, std::ios::beg);

			//Temporary variable used during import
			double q = 0;
			double I = 0;
			double e = 0;


			//Read and parse file line by line
			int counter = 0;
			while (getline(infile, line))
			{
				counter += 1;

				std::istringstream iss(line);
				if (error_specified) iss >> q >> I >> e;
				else
				{
					iss >> q >> I;
					e = 0;
				}

				//Check for errors
				if (iss.fail())
				{
					saxs::error_message("Error in " + file_name + " at line " + QString::number(counter));
					infile.close();
					return false;
				}

				//Store data
				data.push_back(scatteringdata_sp(new scatteringdata(q, I, e)));
			}
		}

		//Rule for chifile
		if (file_ext == "chi")
		{
			//Skip first 4 lines
			for (int i = 0; i < 4;i++) getline(infile, line);

			//Check number of columns to read error
			getline(infile, line);
			std::istringstream iss(line);
			int columns = 0;
			do
			{
				std::string sub;
				iss >> sub;
				if (sub.length())
					++columns;
			} while (iss);
			double error_scale = 0;
			bool error_specified;
			switch (columns)
			{
			case 0: saxs::error_message("No data in file:" + file_name + "!");
				infile.close();
				return false;
			case 1: saxs::error_message("No scattering data in file:" + file_name + "!");
				infile.close();
				return false;
			case 2: //error_scale = saxs::input_message_doub("No error given in file "+ file_name +"\nSpecify artificial error weight in %:",3,0.0001,100);
				error_specified = false;
				break;
			case 3: error_specified = true;
				break;
			}

			//Go back to line 0
			infile.clear();
			infile.seekg(0, std::ios::beg);
			for (int i = 0; i < 4; i++) getline(infile, line);

			//Temporary variable used during import
			double q = 0;
			double I = 0;
			double e = 0;


			//Read and parse file line by line
			int counter = 3;
			while (getline(infile, line))
			{
				counter += 1;

				//Check for nan
				if (line.find("nan") != std::string::npos) continue;
				//Check for empty lines
				if (line.length()<3) continue;

				std::istringstream iss(line);
				if (error_specified) iss >> q >> I >> e;
				else
				{
					iss >> q >> I;
					e = 0;
				}

				//Check for errors
				if (iss.fail())
				{
					saxs::error_message("Error in " + file_name + " at line " + QString::number(counter));
					infile.close();
					return false;
				}

				//Store data
				data.push_back(scatteringdata_sp(new scatteringdata(q, I, e)));
			}
		}


		infile.close();
		return true;
	}

	//	Import 3d model data from xyz file
	bool import_xyz_model(const QString file_name, std::vector<coordinate_sp>& data)
	{
		//transform qstring into string
		std::string file_name_string = file_name.toStdString();

		//Predefine Messagebox for error warnings
		QMessageBox msgBox;
		msgBox.setText("Error");
		msgBox.setStandardButtons(QMessageBox::Ok);
		msgBox.setIcon(QMessageBox::Critical);

		//Open file 
		std::ifstream infile(file_name_string.c_str(), std::ifstream::in);

		//Check if correctly open
		if (!infile.is_open())
		{
			std::stringstream error_stream;
			msgBox.setInformativeText("Could not open file " + file_name + "!");
			int ret = msgBox.exec();
			return false;
		}
		//Temporary variable used during import
		double x = 0;
		double y = 0;
		double z = 0;

		//Temporary line used during import
		std::string line;

		//Read and parse file line by line
		int counter = 0;
		while (getline(infile, line))
		{
			counter += 1;

			std::istringstream iss(line);
			iss >> x >> y >> z;

			//Check for errors
			if (iss.fail())
			{
				std::stringstream error_stream;
				msgBox.setInformativeText("Error in " + file_name + " at line " + QString::number(counter));
				int ret = msgBox.exec();
				infile.close();
				return false;
			}

			//Store coordinates
			data.push_back(coordinate_sp(new coordinate(x, y, z, 0)));
		}

		infile.close();
		return true;
	}

	//	Import 3d model data from pdb file
	bool import_pdb_model(const QString file_name, std::vector<coordinate_sp>& data)
	{
		//transform qstring into string
		std::string file_name_string = file_name.toStdString();

		//Predefine Messagebox for error warnings
		QMessageBox msgBox;
		msgBox.setText("Error");
		msgBox.setStandardButtons(QMessageBox::Ok);
		msgBox.setIcon(QMessageBox::Critical);

		//Open file 
		std::ifstream infile(file_name_string.c_str(), std::ifstream::in);

		//Check if correctly open
		if (!infile.is_open())
		{
			std::stringstream error_stream;
			msgBox.setInformativeText("Could not open file " + file_name + "!");
			int ret = msgBox.exec();
			return false;
		}
		//Temporary variable used during import
		double x = 0;
		double y = 0;
		double z = 0;
		std::string catchstring;

		//Temporary line used during import
		std::string line;

		//Read and parse file line by line
		int counter = 0;
		while (getline(infile, line))
		{
			counter += 1;

			//Check if ATOM is first in lin
			if (line.substr(0, 4) != "ATOM") continue;
			else line = line.substr(30, 80);

			std::istringstream iss(line);
			iss >> x >> y >> z;

			//Check for errors
			if (iss.fail())
			{
				std::stringstream error_stream;
				msgBox.setInformativeText("Error in " + file_name + " at line " + QString::number(counter));
				int ret = msgBox.exec();
				infile.close();
				return false;
			}

			//Store coordinates
			data.push_back(coordinate_sp(new coordinate(x, y, z, 0)));
		}

		infile.close();
		return true;
	}

	//	Generate artifical voxel grid
	void generate_art_grid(std::vector<coordinate_sp>& model, const double rmax, 
					const std::pair<double,double>z_boundaries, const double resolution)
	{
		double zmin = z_boundaries.first;
		double zmax = z_boundaries.second;

		double xstart = -rmax + resolution / 2.;
		double ystart = xstart;
		double zstart = zmin + resolution / 2.;

		double x = xstart, y=ystart, z=zstart;

		while (z < zmax)
		{
			y = ystart;
			while (y < rmax)
			{
				x = xstart;
				while (x < rmax)
				{
					model.push_back(coordinate_sp(new coordinate(x, y, z, 0)));
					x += resolution;
				}
				y += resolution;
			}
			z += resolution;
		}
	}

	//	Superimposes merged model on artifical grid
	void impose_model_on_grid(std::vector<coordinate_sp>& model, std::vector<coordinate_sp>& grid)
	{
		int max_count = 0;
		double r_min;
		int min_pntr;
		double r;

		for (int i = 0; i < model.size(); i++)
		{
			r_min = 10000000000;
			for (int k = 0; k < grid.size(); k++)
			{
				r = std::pow(grid[k]->m_x - model[i]->m_x, 2) +
					std::pow(grid[k]->m_y - model[i]->m_y, 2) +
					std::pow(grid[k]->m_z - model[i]->m_z, 2);
				if (r < r_min) 
				{
					r_min = r;
					min_pntr = k;
				}
			}
			grid[min_pntr]->m_nr_contacts += 1;
			if (grid[min_pntr]->m_nr_contacts>max_count) max_count = grid[min_pntr]->m_nr_contacts;
		}

		//Now delete zeros
		int count = 0;
		do
		{
			if(grid[count]->m_nr_contacts == 0)grid.erase(grid.begin() + count);
			else 
			{
				count++;
			}
		} while (count < grid.size());
	}


	//------------------------------------------------------------------------------
	//	Data Type conversions
	//------------------------------------------------------------------------------

	//	Extract values from scatteringdataobject into globals 
	void convertScatteringToVector(std::vector<saxs::scatteringdata_sp> scatterobject, QVector<double>& x, QVector<double>& y, QVector<double>& err)
	{
		int n = scatterobject.size();
		int i = 0;
		x.resize(n);
		y.resize(n);
		err.resize(n);

		for (int i = 0; i<n; ++i)
		{
			x[i] = scatterobject[i]->m_q;
			y[i] = scatterobject[i]->m_I;
			err[i] = scatterobject[i]->m_e;
		}
	}

	//	Rescales scattering data so smallest value is 1
	void rescaleScatteringData(std::vector<saxs::scatteringdata_sp> scatterobject, QVector<double>& y, QVector<double>& err)
	{
		int n = scatterobject.size();
		double ymin = *std::min_element(y.constBegin(), y.constEnd());
		int i = 0;
		bool calcerr = false;
		if (err[0] == 0) calcerr = true;

		for (int i = 0; i<n; ++i)
		{
			y[i] = scatterobject[i]->m_I / ymin;
			if (!calcerr) err[i] = scatterobject[i]->m_e / ymin;
			else err[i] = 0.1*std::pow(y[i], 0.5);
			scatterobject[i]->m_I = y[i];
			scatterobject[i]->m_e = err[i];
		}
	}

	//	Checks the number of points in the scattering data and reduces it
	void reduce_scatteringdata_size(std::vector<saxs::scatteringdata_sp>& scatterobject)
	{
		//Firt check size
		int num_points,num_newpoints,rest_newpoints;
		num_points = scatterobject.size();
		if (num_points > 128)
		{
			std::string str;
			str = "The loaded data contains " + boost::lexical_cast<std::string>(num_points) + " points." +
				"\n\nSpecify divisor or abort: ";
			int divisor = input_message(QString::fromStdString(str), 1,1,20);
			//Check if "abort" was pressed.
			if (divisor == 1) return;
			num_newpoints = num_points / divisor;

			//Make temporary vectors/variables
			std::vector<double> q_vec, int_vec, err_vec;
			int current_idnex;

			//averaging into temporary vectors
			for (int i = 0; i < num_newpoints; i++)
			{
				q_vec.push_back(scatterobject[i*divisor]->m_q);
				int_vec.push_back(scatterobject[i*divisor]->m_I);
				err_vec.push_back(scatterobject[i*divisor]->m_e);
			}

			//now write to scatterobject
			scatterobject.clear();

			for (int i = 0; i < num_newpoints; i++)
			{
				scatterobject.push_back(scatteringdata_sp(new scatteringdata(
					q_vec[i],
					int_vec[i],
					err_vec[i])));
			}
			return;
		}
		else
		{
			return;
		}
	}

	//	Extract values from Coordinateobject into globals 
	void convertCoordinateToVector(std::vector<saxs::coordinate_sp> coordinate, QVector<double>& x, QVector<double>& y, QVector<double>& z)
	{
		int n = coordinate.size();
		int i = 0;
		x.resize(n);
		y.resize(n);
		z.resize(n);

		for (int i = 0; i<n; ++i)
		{
			x[i] = coordinate[i]->m_x;
			y[i] = coordinate[i]->m_y;
			z[i] = coordinate[i]->m_z;
		}
	}

	//	Convert QVector values into Coordinateobject
	void convertQvectorsToCoordinate(std::vector<saxs::coordinate_sp>& loc_coordinate, QVector<double>& x, QVector<double>& y, QVector<double>& z)
	{
		int n = x.size();
		int i = 0;
		loc_coordinate.clear();
		for (int i = 0; i<n; ++i)
		{
			loc_coordinate.push_back(coordinate_sp(new coordinate(x[i], y[i], z[i], 0)));
		}
	}

	//	Convert values from Coordinateobject into cylindrical cartesian system 
	std::vector<double> convertCoordinateToCylinder(const std::vector<saxs::coordinate_sp>& cartesian, std::vector<saxs::coordinate_sp>& cylinder)
	{
		cylinder.clear();
		int n = cartesian.size();
		int i = 0;
		std::vector<double> boundaries(3);
		double r, phi, z;
		double rmax = 0;
		double zmax = 0;
		double zmin = 100;
		for (int i = 0; i<n; ++i)
		{
			r = std::pow(std::pow(cartesian[i]->m_x, 2) + std::pow(cartesian[i]->m_y, 2), 0.5);
			phi = std::atan2(cartesian[i]->m_y, cartesian[i]->m_x);
			z = cartesian[i]->m_z;
			if (rmax < r) { rmax = r; }
			if (zmax < z) { zmax = z; }
			if (zmin > z) { zmin = z; }
			cylinder.push_back(coordinate_sp(new coordinate(r, phi, z, 0)));
		}
		boundaries[0] = rmax;
		boundaries[1] = zmax;
		boundaries[2] = zmin;
		return boundaries;
	}

	//	Convert values from cylindrical object into plotting QVectors 
	void convertCylinderToQVectors(const std::vector<saxs::coordinate_sp>& cylindrical,
		QVector<double>& x, QVector<double>& y, const double phi)
	{
		int n = cylindrical.size();
		int i = 0;
		for (int i = 0; i<n; ++i)
		{
			x[i] = cylindrical[i]->m_x * std::cos(cylindrical[i]->m_y + phi);
			y[i] = cylindrical[i]->m_x * std::sin(cylindrical[i]->m_y + phi);
		}
	}

	//	Convert values from cylindrical object into std vector 
	void convertCylinderToStdVectors(const std::vector<saxs::coordinate_sp>& cylindrical,
		std::vector<double>& x, std::vector<double>& y, const double phi)
	{
		int n = cylindrical.size();
		int i = 0;
		for (int i = 0; i<n; ++i)
		{
			x[i] = cylindrical[i]->m_x * std::cos(cylindrical[i]->m_y + phi);
			y[i] = cylindrical[i]->m_x * std::sin(cylindrical[i]->m_y + phi);
		}
	}

	//	Generate random cylindrical coordinates into coordinate structure
	void gen_cyl_randomseed(std::vector<coordinate_sp>& cylindrical_coord, double height, double diameter, const int nr_atoms)
	{
		//pseudo coords: m_x = r, m_y = phi, m_z = z

		typedef boost::mt19937 RNGType;
		RNGType rng(std::time(0));
		boost::uniform_real<> zero_to_one(0, 1);
		boost::variate_generator<RNGType, boost::uniform_real<>> dice(rng, zero_to_one);

		for (int i = 0; i < nr_atoms; i++)
			cylindrical_coord.push_back(coordinate_sp(new coordinate(
				double(std::pow(dice(), 0.5) * diameter / 2),
				double(dice() * 2 * 3.14159265359),
				double(dice() * height),
				0
				)));
	}

	// Expand current fittingobject by num_expand
	void expand_fittingobject(fittingobject_sp& reftofittingobject, int num_expand)
	{
		typedef boost::mt19937 RNGType;
		RNGType rng(std::time(0));
		boost::uniform_real<> zero_to_one(0, 1);
		boost::variate_generator<RNGType, boost::uniform_real<>> dice(rng, zero_to_one);
		int index;
		int i=0;

		do
		{
			index = int(dice()*reftofittingobject->m_model.size());
			if (index > (reftofittingobject->m_model.size() - 1)) continue;
			reftofittingobject->m_model.push_back(coordinate_sp(new coordinate(
				double(reftofittingobject->m_model[index]->m_x + (double(dice()) - 0.5) * reftofittingobject->m_diameter*0.2),
				double(reftofittingobject->m_model[index]->m_y + (double(dice()) - 0.5) * reftofittingobject->m_diameter*0.2),
				double(reftofittingobject->m_model[index]->m_z + (double(dice()) - 0.5) * reftofittingobject->m_diameter*0.2),
				0
				)));
			i++;
		} while (i < num_expand);
	}



	/***************************************************************************************************************/
	// Output functions/File Writers
	/***************************************************************************************************************/
	
	//Returns current time as string
	std::string now_str()
	{
		// Get current time from the clock, using microseconds resolution
		const boost::posix_time::ptime now =
			boost::posix_time::microsec_clock::local_time();

		// Get the time offset in current day
		const boost::posix_time::time_duration td = now.time_of_day();

		//
		// Extract hours, minutes, seconds and milliseconds.
		//
		// Since there is no direct accessor ".milliseconds()",
		// milliseconds are computed _by difference_ between total milliseconds
		// (for which there is an accessor), and the hours/minutes/seconds
		// values previously fetched.
		//
		const long hours = td.hours();
		const long minutes = td.minutes();
		const long seconds = td.seconds();
		const long milliseconds = td.total_milliseconds() -
			((hours * 3600 + minutes * 60 + seconds) * 1000);

		//
		// Format like this:
		//
		//      hh:mm:ss.SSS
		//
		// e.g. 02:15:40:321
		//
		//      ^          ^
		//      |          |
		//      123456789*12
		//      ---------10-     --> 12 chars + \0 --> 13 chars should suffice
		//  
		// 
		char buf[40];
		sprintf(buf, "%02ld:%02ld:%02ld.%03ld",
			hours, minutes, seconds, milliseconds);

		return buf;
	}

	//------------------------------------------------------------------------------
	//	Base functions to save model as *.pdb file and scattering data as *.chi
	//------------------------------------------------------------------------------
	//Helper function: prints double as string
	std::string prd(double x, int decDigits, const int width)
	{
		std::stringstream ss;
		ss << std::right;
		ss.fill(' ');        // fill space around displayed #
		ss.width(width);     // set  width around displayed #
		ss.precision(decDigits); // set # places after decimal
		ss << std::fixed;
		ss << x;
		return ss.str();
	}

	//Helper function: prints integer as string
	std::string pri(int x, const int width)
	{
		std::stringstream ss;
		ss << std::right;
		ss.fill(' ');        // fill space around displayed #
		ss.width(width);     // set  width around displayed #
		ss.precision(0);	// set # places after decimal
		ss << std::fixed;
		ss << x;
		return ss.str();
	}

	//Save Single Model as PDB file
	void writemodeltopdb(std::string file_name, std::vector<coordinate_sp>& m_model)
	{
		//add fileextension
		file_name += ".pdb";

		//Open file 
		std::ofstream outfile(file_name.c_str(), std::ifstream::out);

		//Check if correctly open
		if (!outfile.is_open())
		{
			std::stringstream error_stream;
			error_message(QString::fromStdString("Cannot open " + file_name));
			return;
		}

		//Store results to file
		for (int i = 0; i<m_model.size(); ++i)
		{
			outfile << "ATOM  " << pri((i + 1), 5) << "  CA  ASP" << pri((1 + int(i / 10)), 5) << "    "
				<< prd(m_model[i]->m_x, 3, 8)
				<< prd(m_model[i]->m_y, 3, 8)
				<< prd(m_model[i]->m_z, 3, 8)
				<< prd(1.0, 2, 6)
				<< prd(20.0, 2, 6)
				<< "           C"
				<< std::endl;
		}
		//Close file
		outfile.close();
	}

	//Save Single Model as PDB file with occupancy
	void writeoccmodeltopdb(std::string file_name, std::vector<coordinate_sp>& m_model)
	{
		//add fileextension
		file_name += ".pdb";

		//Open file 
		std::ofstream outfile(file_name.c_str(), std::ifstream::out);

		//Check if correctly open
		if (!outfile.is_open())
		{
			std::stringstream error_stream;
			error_message(QString::fromStdString("Cannot open " + file_name));
			return;
		}

		//Find Maximum Occupancy
		int occ_max=0;
		for (int i = 0; i<m_model.size(); ++i)
			if (m_model[i]->m_nr_contacts>occ_max)occ_max = m_model[i]->m_nr_contacts;

		//Store results to file
		for (int i = 0; i<m_model.size(); ++i)
		{
			outfile << "ATOM  " << pri((i + 1), 5) << "  CA  ASP" << pri((1 + int(i / 10)), 5) << "    "
				<< prd(m_model[i]->m_x, 3, 8)
				<< prd(m_model[i]->m_y, 3, 8)
				<< prd(m_model[i]->m_z, 3, 8)
				<< prd(double(m_model[i]->m_nr_contacts)/double(occ_max), 2, 6)
				<< prd(20.0, 2, 6)
				<< "           C"
				<< std::endl;
		}
		//Close file
		outfile.close();
	}

	//Save Stacked Model as PDB file
	void writestackedmodeltopdb(std::string file_name, std::vector<coordinate_sp>& m_model, int num_stacks, double stack_spacing)
	{
		//add fileextension
		file_name += "_stck.pdb";

		//Open file 
		std::ofstream outfile(file_name.c_str(), std::ifstream::out);

		//Check if correctly open
		if (!outfile.is_open())
		{
			std::stringstream error_stream;
			error_message(QString::fromStdString("Cannot open " + file_name));
			return;
		}

		//Store results to file
		for (int j = 0; j < num_stacks; j++)
		{
			for (int i = 0; i<m_model.size(); ++i)
			{
				outfile << "ATOM  " << pri((j + num_stacks + i + 1), 5)
					<< "  CA  ASP" 
					<< pri((1 + int( (j + num_stacks + i) / 10)), 5)
					<< "    "
					<< prd(m_model[i]->m_x, 3, 8)
					<< prd(m_model[i]->m_y, 3, 8)
					<< prd(m_model[i]->m_z + double(j) * double(stack_spacing), 3, 8)
					<< prd(1.0, 2, 6)
					<< prd(20.0, 2, 6)
					<< "           C"
					<< std::endl;
			}
		}
		//Close file
		outfile.close();
	}

	//Save Stacked Model as PDB file with occupancy
	void writestackedoccmodeltopdb(std::string file_name, std::vector<coordinate_sp>& m_model, int num_stacks, double stack_spacing)
	{
		//add fileextension
		file_name += "_stck.pdb";

		//Open file 
		std::ofstream outfile(file_name.c_str(), std::ifstream::out);

		//Check if correctly open
		if (!outfile.is_open())
		{
			std::stringstream error_stream;
			error_message(QString::fromStdString("Cannot open " + file_name));
			return;
		}

		//Find Maximum Occupancy
		int occ_max = 0;
		for (int i = 0; i<m_model.size(); ++i)
			if (m_model[i]->m_nr_contacts>occ_max)occ_max = m_model[i]->m_nr_contacts;

		//Store results to file
		for (int j = 0; j < num_stacks; j++)
		{
			for (int i = 0; i<m_model.size(); ++i)
			{
				outfile << "ATOM  " << pri((j + num_stacks + i + 1), 5)
					<< "  CA  ASP"
					<< pri((1 + int((j + num_stacks + i) / 10)), 5)
					<< "    "
					<< prd(m_model[i]->m_x, 3, 8)
					<< prd(m_model[i]->m_y, 3, 8)
					<< prd(m_model[i]->m_z + double(j) * double(stack_spacing), 3, 8)
					<< prd(double(m_model[i]->m_nr_contacts) / double(occ_max), 2, 6)
					<< prd(20.0, 2, 6)
					<< "           C"
					<< std::endl;
			}
		}
		//Close file
		outfile.close();
	}

	//Save Scattering Data as CHI file
	void writedatatochi(std::string file_name, const std::vector<double>& m_data_q, const std::vector<double>& m_data_I,
		const std::vector<double>& m_fitted_I)
	{
		//add fileextension
		file_name += ".chi";

		//Open file 
		std::ofstream outfile(file_name.c_str(), std::ifstream::out);

		//Check if correctly open
		if (!outfile.is_open())
		{
			std::stringstream error_stream;
			error_message(QString::fromStdString("Cannot open " + file_name));
			return;
		}

		//Set precision and notation
		outfile.precision(16);
		outfile << std::scientific;

		//Store results to file
		for (int i = 0; i<m_data_q.size(); ++i)
		{
			outfile << " " << m_data_q[i]
				<< " " << m_data_I[i]
				<< " " << m_fitted_I[i]
				<< std::endl;
		}
		//Close file
		outfile.close();
	}

	//Saves Tracer values in log file
	void writechisquretofile(std::string file_name, int coreid, std::vector<saxs::model_sp>& result_tracer)
	{
		//add fileextension
		file_name += ".log";

		//Open file 
		std::ofstream outfile(file_name.c_str(), std::ifstream::out);

		//Check if correctly open
		if (!outfile.is_open())
		{
			std::stringstream error_stream;
			error_message(QString::fromStdString("Cannot open " + file_name));
			return;
		}

		//Set precision and notation
		
		outfile << std::scientific;

		//Make File Header
		outfile << "run";
		outfile << "\t" << "Target function";
		outfile << "\t" << "ChiSquare";
		outfile << "\t" << "Connectivity";
		outfile << "\t" << "Compactness";
		outfile << std::endl;

		//Store results to file
		for (int j = 0; j < result_tracer[coreid]->m_chi_tracer.size(); j++)
		{
			outfile.precision(0);
			outfile << j;
			outfile.precision(8);
			outfile << std::scientific;
			outfile << "\t" << result_tracer[coreid]->m_target_f_tracer[j];
			outfile << "\t" << result_tracer[coreid]->m_chi_tracer[j];
			outfile << "\t" << result_tracer[coreid]->m_connectivity_tracer[j];
			outfile << "\t" << result_tracer[coreid]->m_compactness_tracer[j];
			outfile  << std::endl;
		}
		//Close file
		outfile.close();
	}

}//End of SAXS namespace

#endif /*!DATA_UTILS_H*/