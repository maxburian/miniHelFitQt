#ifndef THREAD_UTILS_H
#define THREAD_UTILS_H
//Include database dependencies
//QT
#include <QtCore>
#include <QApplication>
//STD
#include <vector>
#include <random>
#include <stdint.h>
#include <iostream>
//BOOST
#include <boost\shared_ptr.hpp>
#include <boost/make_shared.hpp>
#include <boost/log/sources/logger.hpp>
#include <boost/log/sources/record_ostream.hpp>
#include <boost/thread.hpp>
#include <boost/math/special_functions/sinc.hpp>
#include "boost/random.hpp"
#include "boost/generator_iterator.hpp"
#include <boost/date_time/posix_time/posix_time.hpp>
#include <boost/lexical_cast.hpp>
//Self-Written
#include "data_def.h"
#include "data_utils.h"

namespace saxs
{
	//helper function
	uint32_t randseed = 7;  // 100% random seed value

	//------------------------------------------------------------------------------
	//	Pseudo-Random number generation
	//------------------------------------------------------------------------------
	double xor128(void) {
		randseed ^= randseed << 13;
		randseed ^= randseed >> 17;
		randseed ^= randseed << 5;
		return double(randseed) / 2 / double(2147483647);
	}



	/***************************************************************************************************************/
	// Thread dedicated functions
	/***************************************************************************************************************/
	//------------------------------------------------------------------------------
	//	Returning dimensional parameters of Model
	//------------------------------------------------------------------------------

	//Returns max. sphere radius of close-packed DA in given cylnder Volume
	double get_critradius_from_coordinates(double cyl_diameter, double stack_spacing, int N)
	{
		return std::pow((3. / 16.*0.74*stack_spacing*cyl_diameter*cyl_diameter / double(N)), 0.3333);
	}

	//Returns the diameter of a cylinder with the same Rg than the model
	double get_cyldiameter_from_coordinates(std::vector<coordinate_sp>& model)
	{
		double rsq = 0;
		for (int i = 0; i < model.size(); i++)
		{
			rsq += model[i]->m_x*model[i]->m_x + model[i]->m_y*model[i]->m_y;
		}
		return 2. * std::sqrt(2 * rsq / double(model.size()));
	}

	//Returns R-max of model
	double get_rmax_from_coordinates(std::vector<coordinate_sp>& model)
	{
		double rmax = 0;
		double r;
		for (int i = 0; i < model.size(); i++)
		{
			r = model[i]->m_x*model[i]->m_x + model[i]->m_y*model[i]->m_y;
			if (r > rmax)rmax = r;
		}
		return std::sqrt(r);
	}

	//Returns pari<z_min,z_max> of model
	std::pair<double,double> get_zboundaries_from_coordinates(std::vector<coordinate_sp>& model)
	{
		double zmin = 100000;
		double zmax = -100000;
		std::pair<double, double> returnvals;
		for (int i = 0; i < model.size(); i++)
		{
			if (model[i]->m_z > zmax)zmax = model[i]->m_z;
			if (model[i]->m_z < zmin)zmin = model[i]->m_z;
		}
		returnvals.first = zmin;
		returnvals.second = zmax;

		return returnvals;
	}


	//------------------------------------------------------------------------------
	//	Returning Compactness/Connectivity parameters of Model
	//------------------------------------------------------------------------------

	//Calculates formfactor of fittingobject
	void calc_ff(fittingobject_sp& reftofittingobject)
	{
		reftofittingobject->m_ff_I.clear();
		for (int i = 0; i < reftofittingobject->m_data_q.size(); i++)
		{
			reftofittingobject->m_ff_I.push_back(
				std::pow(
					std::exp(
						-(reftofittingobject->m_data_q[i]* reftofittingobject->m_data_q[i])/
						(reftofittingobject->m_sigma_ff*reftofittingobject->m_sigma_ff))
				, 2)
			);
		}
	}

	//Calculates the number of contacts for each DA
	void calc_contacts_of_model(fittingobject_sp& reftofittingobject)
	{
		for (int i = 0; i < reftofittingobject->m_model.size(); i++) reftofittingobject->m_model[i]->m_nr_contacts = 0;

		double dx, dy, dz, rsqu;
		for (int i = 0; i < reftofittingobject->m_model.size() - 1; i++)
		{
			for (int j = i; j < reftofittingobject->m_model.size(); j++)
			{
				dx = reftofittingobject->m_model[j]->m_x - reftofittingobject->m_model[i]->m_x;
				dy = reftofittingobject->m_model[j]->m_y - reftofittingobject->m_model[i]->m_y;
				dz = reftofittingobject->m_model[j]->m_z - reftofittingobject->m_model[i]->m_z;
				if ((dx*dx + dy*dy + dz*dz) < reftofittingobject->m_contact_d_sq)
				{
					reftofittingobject->m_model[i]->m_nr_contacts++;
					reftofittingobject->m_model[j]->m_nr_contacts++;
				}
			}
			//Account for atoms at upper z-border of BB:
			if (reftofittingobject->m_model[i]->m_z > (reftofittingobject->m_stack_spacing - std::sqrt(reftofittingobject->m_contact_d_sq)))
			{
				for (int j = i; j < reftofittingobject->m_model.size(); j++)
				{
					dx = reftofittingobject->m_model[j]->m_x - reftofittingobject->m_model[i]->m_x;
					dx = reftofittingobject->m_model[j]->m_y - reftofittingobject->m_model[i]->m_y;
					dx = reftofittingobject->m_model[j]->m_z + reftofittingobject->m_stack_spacing - reftofittingobject->m_model[i]->m_z;
					if ((dx*dx + dy*dy + dz*dz) < reftofittingobject->m_contact_d_sq)
					{
						reftofittingobject->m_model[i]->m_nr_contacts++;
						reftofittingobject->m_model[j]->m_nr_contacts++;
					}
				}
			}
			//Account for atoms at lower z-border of BB:
			if (reftofittingobject->m_model[i]->m_z < (std::sqrt(reftofittingobject->m_contact_d_sq)))
			{
				for (int j = i; j < reftofittingobject->m_model.size(); j++)
				{
					dx = reftofittingobject->m_model[j]->m_x - reftofittingobject->m_model[i]->m_x;
					dx = reftofittingobject->m_model[j]->m_y - reftofittingobject->m_model[i]->m_y;
					dx = reftofittingobject->m_model[j]->m_z - reftofittingobject->m_stack_spacing - reftofittingobject->m_model[i]->m_z;
					if ((dx*dx + dy*dy + dz*dz) < reftofittingobject->m_contact_d_sq)
					{
						reftofittingobject->m_model[i]->m_nr_contacts++;
						reftofittingobject->m_model[j]->m_nr_contacts++;
					}
				}
			}
		}
	}

	//returns the mean number of contacts of the model
	double return_mean_contacts_of_model(std::vector<coordinate_sp>& model)
	{
		int contactsum = 0;
		for (int i = 0; i < model.size(); i++)
		{
			if (model[i]->m_nr_contacts<12)	contactsum += model[i]->m_nr_contacts;
			else contactsum += 12;
			 
		}
		return double(contactsum) / double(model.size());
	}

	//returns the mean number of contacts of contactvector
	double return_mean_contacts_of_contacvector(std::vector<int>& contacvector)
	{
		int contactsum = 0;
		for (int i = 0; i < contacvector.size(); i++)
		{
			contactsum += contacvector[i];
		}
		return double(contactsum) / double(contacvector.size());
	}

	//returns the mean conenctivity of the fittingobject
	double return_mean_connectivity_of_fittingobject(fittingobject_sp& reftofittingobject)
	{
		int contacts = 0;
		double connectivity = 0;
		for (int i = 0; i < reftofittingobject->m_model.size(); i++)
		{
			if (reftofittingobject->m_model[i]->m_nr_contacts<13) connectivity +=
				(std::exp(-0.5 / reftofittingobject->m_tauConn*double(reftofittingobject->m_model[i]->m_nr_contacts)) - 
					std::exp(-0.5 / reftofittingobject->m_tauConn*12.));
			else connectivity += (1- (1 - (std::exp(-0.5 / reftofittingobject->m_tauConn*double(24 - reftofittingobject->m_model[i]->m_nr_contacts)))))-
				std::exp(-0.5 / reftofittingobject->m_tauConn*12.);
		}
		return double(connectivity) / double(reftofittingobject->m_model.size());
	}

	//returns the mean connectivity of contactvector
	double return_mean_connectivity_of_contactvector(std::vector<int>& contacvector, fittingobject_sp& reftofittingobject)
	{
		int contacts = 0;
		double connectivity = 0;
		for (int i = 0; i < contacvector.size(); i++)
		{
			if (contacvector[i]<13) connectivity += (std::exp(-0.5/reftofittingobject->m_tauConn*double(contacvector[i])) - 
				std::exp(-0.5 / reftofittingobject->m_tauConn*12.));
			else connectivity += (1-(1 - (std::exp(-0.5 /reftofittingobject->m_tauConn*double(24 - contacvector[i]))))) - 
				std::exp(-0.5 / reftofittingobject->m_tauConn*12.);
			/*
			if (contacvector[i]<12) contacts = contacvector[i];
			else contacts = 12;
			connectivity += (1 - (std::exp(-0.25*double(contacts)) - std::exp(-6.0)));
			*/
		}
		return double(connectivity) / double(contacvector.size());
	}

	//returns the fraction of saturated DAs
	double return_fraction_contactsaturated_das(std::vector<int>& contacvector)
	{
		int contacts = 0;
		double connectivity = 0;
		for (int i = 0; i < contacvector.size(); i++)
		{
			if (contacvector[i] > 11) contacts += 1;
		}
		return double(contacts) / double(contacvector.size());
	}

	//returns the compactness of fittingobject
	double return_compactness_of_fittingobject(fittingobject_sp& reftofittingobject)
	{
		double compactness = 0;
		double x;
		double y;
		for (int i = 0; i < reftofittingobject->m_model.size(); i++)
		{
			x = reftofittingobject->m_model[i]->m_x;
			y = reftofittingobject->m_model[i]->m_y;
			compactness += (1 - std::erfc(0.5 - (1. - std::sqrt(x*x + y*y)/ (reftofittingobject->m_diameter/2.))*reftofittingobject->m_sigmaComp*3.1415) / 2.0);
		}
		return double(compactness) / double(reftofittingobject->m_model.size());
	}

	//recenter fittingobject to xy=(0,0)
	void recenter_fittingobject(fittingobject_sp& reftofittingobject)
	{
		int nr_da_inside = 0;
		double x_com = 0;
		double y_com = 0;
		for (int i = 0; i < reftofittingobject->m_model.size(); i++)
		{
			x_com += reftofittingobject->m_model[i]->m_x;
			y_com += reftofittingobject->m_model[i]->m_y;
		}

		x_com = x_com / double(reftofittingobject->m_model.size());
		y_com = y_com / double(reftofittingobject->m_model.size());

		for (int i = 0; i < reftofittingobject->m_model.size(); i++)
		{
			reftofittingobject->m_model[i]->m_x-=x_com;
			reftofittingobject->m_model[i]->m_y-=y_com;
		}
	}


	//------------------------------------------------------------------------------
	//	Theoretical/Experimental data fitting
	//------------------------------------------------------------------------------
	
	//Returns ChiSquare Value of two vectors
	double chisquare(std::vector<double>& exp, std::vector<double>& err, std::vector<double>& model)
	{
		//Chisquare between data1(exp) and data2(model)
		double chi = 0;
		double weight = 0;
		for (int i = 0; i < exp.size(); i++)
		{
			chi += std::pow((exp[i] - model[i]) / err[i], 2);
			weight += 1;
		}
		return chi / weight;
	}

	//Linear Regression fit of two vectors (data is weighted)
	std::pair<double, double> lindatafit(std::vector<double>& data, std::vector<double>& fit, std::vector<double>& weight)
	{
		double f2sum = 0;
		double fsum = 0;
		double dfsum = 0;
		double dsum = 0;
		double weightsum = 0;

		for (int i = 0; i < data.size(); i++)
		{
			f2sum = f2sum + fit[i] * fit[i] * weight[i];
			fsum = fsum + fit[i] * weight[i];
			dfsum = dfsum + data[i] * fit[i] * weight[i];
			dsum = dsum + data[i] * weight[i];
			weightsum = weightsum + weight[i];
		}
		double k = (weightsum*dfsum - dsum*fsum) / (weightsum*f2sum - fsum*fsum);
		double d = (dfsum - f2sum*k) / fsum;

		std::pair<double, double> fitresults(d, k);
		return fitresults;
	}

	//Fits model_i to data_i of fittingobject - result in fitted_i
	void fit_modeltoexp(fittingobject_sp& reftofittingobject)
	{
		//Curve weighting: 0...q,1...q^2, 2...none
		std::vector<double> exp;
		std::vector<double> model;
		std::vector<double> weight;
		std::pair<double, double> fitparams;

		reftofittingobject->m_norm_offset = *std::min_element(reftofittingobject->m_model_I.begin(), reftofittingobject->m_model_I.end()) + 1;
		reftofittingobject->m_norm = (*std::max_element(reftofittingobject->m_model_I.begin(), reftofittingobject->m_model_I.end()) - reftofittingobject->m_norm_offset);

		for (int i = 0; i < reftofittingobject->m_model_I.size(); i++)
		{
			reftofittingobject->m_fitted_I[i] = reftofittingobject->m_ff_I[i]*(reftofittingobject->m_model_I[i] -
				reftofittingobject->m_norm_offset) / (reftofittingobject->m_norm)*reftofittingobject->m_data_I[0];
		}

		//Loading data in standard vectors and finding minimum value
		for (int i = 0; i < reftofittingobject->m_data_q.size(); i++)
		{
			exp.push_back(reftofittingobject->m_data_I[i]);
			model.push_back(reftofittingobject->m_fitted_I[i]);
			weight.push_back(1. / reftofittingobject->m_data_e[i]);
		}
		weight.resize(model.size());
		//Apply weighting conditions
		for (int i = 0; i < reftofittingobject->m_data_q.size(); i++)
		{
			if (reftofittingobject->m_weighing == 0)
			{
				weight[i] *= reftofittingobject->m_data_q[i];
			}
			if (reftofittingobject->m_weighing == 1)
			{
				weight[i] *= std::pow(reftofittingobject->m_data_q[i], 2);
			}
		}

		//GetFitparams
		fitparams = lindatafit(exp, model, weight);
		reftofittingobject->m_linreg_results = fitparams;
		for (int i = 0; i < reftofittingobject->m_data_q.size(); i++)
		{
			reftofittingobject->m_fitted_I[i] = fitparams.first + fitparams.second * reftofittingobject->m_fitted_I[i];
		}

		reftofittingobject->m_chi = chisquare(reftofittingobject->m_data_I, reftofittingobject->m_data_e,
			reftofittingobject->m_fitted_I);
	}

	//Recalculates fitted_i from model_i using fitparams
	void recalc_fitted_i(fittingobject_sp& localfittingobject_sp)
	{
		std::pair<double, double> fitparams = localfittingobject_sp->m_linreg_results;
		for (int i = 0; i < localfittingobject_sp->m_data_q.size(); i++)
		{
			localfittingobject_sp->m_fitted_I[i] = localfittingobject_sp->m_ff_I[i] * (localfittingobject_sp->m_model_I[i] -
				localfittingobject_sp->m_norm_offset) / (localfittingobject_sp->m_norm)*localfittingobject_sp->m_data_I[0];
			localfittingobject_sp->m_fitted_I[i] = localfittingobject_sp->m_linreg_results.first + localfittingobject_sp->m_linreg_results.second * localfittingobject_sp->m_fitted_I[i];
		}
		localfittingobject_sp->m_chi = chisquare(localfittingobject_sp->m_data_I, localfittingobject_sp->m_data_e,
			localfittingobject_sp->m_fitted_I);
	}

	//Returns target function that is minimized during the fitting procedure
	double tartget_function(double chisq, double conn, double conn_w,
		double comp, double comp_w, fittingobject_sp& reftofittingobject)
	{
		return chisq + conn_w*(conn) +  comp_w*0.5*comp;
	}



	/***************************************************************************************************************/
	/*								MAIN FITTING AND CALCULATION FUNCTTIONS										   */
	/***************************************************************************************************************/


	/***************************************************************/
	// DebyeCalc Functions
	/***************************************************************/
	//------------------------------------------------------------------------------
	//	Thread: task to calculate static debye function of fittingobject in specified q-regime
	//------------------------------------------------------------------------------
	void calculate_debye(boost::mutex& io_mutex, int start_index, int end_index,fittingobject_sp& fittingobject)
	{
		//Thread local temporary vector
		std::vector<double> local_results(fittingobject->m_data_q.size(), 0.0f);

		//Keep size of lookup table
		std::size_t sinc_lookup_size = saxs::sinc_lookup.size();

		//Temporary variables used in Debye function
		double r = 0;
		double qr = 0;

		//local coordinate vector
		std::vector<coordinate_sp> coordinates = fittingobject->m_model;

		//Compute Debye function [Max Code]
		//Loop around every atom in BuildingBlock
		for (std::size_t i = 0; i < (coordinates.size()); ++i)
		{
			//BuildingBlockCodnition
			if (i < (coordinates.size() - 1))
			{
				//loop over j points
				for (std::size_t j = i + 1; j < (coordinates.size()); ++j)
				{
					r = std::pow(
						std::pow(coordinates[i]->m_x - coordinates[j]->m_x, 2) +
						std::pow(coordinates[i]->m_y - coordinates[j]->m_y, 2) +
						std::pow(coordinates[i]->m_z - coordinates[j]->m_z, 2),
						0.5);//loop over q-points

					for (int k = start_index; k < end_index; k++)
					{
						qr = fittingobject->m_data_q[k] * r * 1000 + 0.5;
						//Lookup Condition
						if (qr < sinc_lookup_size)
						{
							local_results[k] += saxs::sinc_lookup[(int)qr] * 2.0*double(fittingobject->m_num_stacks);
						}
						else
						{
							local_results[k] += boost::math::sinc_pi(double(qr) / 1000)*2.0*double(fittingobject->m_num_stacks);
						}
					}//k index loop
				}//j index loop
			}//if Building block Condition

			 //Loop over stacks
			for (int l = 1; l < fittingobject->m_num_stacks; ++l)
			{
				//loop over j - radial pairs
				for (std::size_t j = 0; j < (coordinates.size()); ++j)
				{
					r = std::pow(
						std::pow(coordinates[i]->m_x - coordinates[j]->m_x, 2) +
						std::pow(coordinates[i]->m_y - coordinates[j]->m_y, 2) +
						std::pow(coordinates[i]->m_z - (coordinates[j]->m_z + double(l*fittingobject->m_stack_spacing)), 2),
						0.5);
					//loop over q-points
					for (int k = start_index; k < end_index; k++)
					{
						qr = fittingobject->m_data_q[k] * r * 1000 + 0.5;
						//Lookup Condition
						if (qr < sinc_lookup_size)
						{
							local_results[k] += 2 * saxs::sinc_lookup[(int)qr] * double(fittingobject->m_num_stacks - l);
						}
						else
						{
							local_results[k] += 2 * boost::math::sinc_pi(double(qr) / 1000)* double(fittingobject->m_num_stacks - l);
						}
					}// k index
				}//j index loop
			}//l index loop
			 //Notify user to current status


		}// i index loop

		 //Lock mutex
		boost::mutex::scoped_lock lock(io_mutex);

		//Copy data from local temporary vector to results
		for (unsigned int k = start_index; k < end_index; k++)
		{
			fittingobject->m_model_I[k] = local_results[k] + double(fittingobject->m_num_stacks*fittingobject->m_model.size());
		}

		//Unlock mutex
		lock.unlock();
	}

	//------------------------------------------------------------------------------
	//	Core: Function to calculate static debye function of fittingobect
	//------------------------------------------------------------------------------
	void calcDebyeStackCurrentModel(fittingobject_sp& reftofittingobject)
	{
		//Thread container
		boost::scoped_ptr<boost::thread_group> threadGroup_sp(new boost::thread_group);

		//Mutex to protect IO operations
		boost::mutex io_mutex;

		//Initialize thread parameters
		int thread_step = static_cast<int>(reftofittingobject->m_data_q.size() / reftofittingobject->m_num_cores);
		int index_delta = reftofittingobject->m_data_q.size() - thread_step*reftofittingobject->m_num_cores;
		//Taking into account index offset due to core number
		//Putting indexboundaries into boundaries vector of size cores+1
		std::vector<int> boundaries;
		boundaries.push_back(0);
		int offset = 1;
		for (int i = 0; i < reftofittingobject->m_num_cores; i++)
		{
			if (i < index_delta) offset = 1;
			else offset = 0;
			boundaries.push_back(boundaries[i] + thread_step + offset);
		}

		//Instantiate and start threads
		for (unsigned int i = 0; i<reftofittingobject->m_num_cores; i++)
		{
			threadGroup_sp->add_thread(new boost::thread(calculate_debye, boost::ref(io_mutex), boundaries[i], boundaries[i + 1], boost::ref(reftofittingobject)));
		}

		//Wait until the end of their jobs
		if (threadGroup_sp)
			threadGroup_sp->join_all();

		//normfittingdata
		if (reftofittingobject->m_data_I.size() == 0)
		{
			double norm_offset = *std::min_element(reftofittingobject->m_model_I.begin(), reftofittingobject->m_model_I.end()) - 1;
			double norm = (*std::max_element(reftofittingobject->m_model_I.begin(), reftofittingobject->m_model_I.end()) - norm_offset);

			for (int i = 0; i < reftofittingobject->m_model_I.size(); i++)
			{
				reftofittingobject->m_fitted_I[i] = reftofittingobject->m_ff_I[i]*(reftofittingobject->m_model_I[i] - norm_offset) / (norm);
			}
			reftofittingobject->m_chi = NAN;
			reftofittingobject->m_linreg_results.first = NAN;
			reftofittingobject->m_linreg_results.second = NAN;
		}
		else
		{
			fit_modeltoexp(reftofittingobject);
		}
		//Make linear regression
	}


	/***************************************************************/
	// DebyeFit Functions
	/***************************************************************/

	//------------------------------------------------------------------------------
	// Recalculates the scattering pattern if changecood is moved by movement - Result will be put into references new_model_i
	//------------------------------------------------------------------------------
	void thread_recalc_change(fittingobject_sp& reftofittingobject, int& changedcoord, coordinate_sp& movement, 
													std::vector<int>& contact_vector)
	{
		//tempvars
		double r_old = 0;
		double qr_old = 0;
		double sinc_qr_old = 0;
		double r_new = 0;
		double qr_new = 0;
		double sinc_qr_new = 0;
		//tempvars for stack
		double r_old_st = 0;
		double qr_old_st = 0;
		double sinc_qr_old_st = 0;
		double r_new_st = 0;
		double qr_new_st = 0;
		double sinc_qr_new_st = 0;
		int sinc_lookup_size = saxs::sinc_lookup.size();
		double x, y;
		
		//reset contactcount on changed atom
		contact_vector[changedcoord] = 0;

		for (int i = 0; i < reftofittingobject->m_model.size(); i++)
		{
			if (i != changedcoord)
			{
				//Building block contribution
				r_old = std::pow(
					std::pow(reftofittingobject->m_model[i]->m_x - reftofittingobject->m_model[changedcoord]->m_x, 2) +
					std::pow(reftofittingobject->m_model[i]->m_y - reftofittingobject->m_model[changedcoord]->m_y, 2) +
					std::pow(reftofittingobject->m_model[i]->m_z - reftofittingobject->m_model[changedcoord]->m_z, 2),
					0.5);
				r_new = std::pow(
					std::pow(reftofittingobject->m_model[i]->m_x - (reftofittingobject->m_model[changedcoord]->m_x + movement->m_x), 2) +
					std::pow(reftofittingobject->m_model[i]->m_y - (reftofittingobject->m_model[changedcoord]->m_y + movement->m_y), 2) +
					std::pow(reftofittingobject->m_model[i]->m_z - (reftofittingobject->m_model[changedcoord]->m_z + movement->m_z), 2),
					0.5);
				//loop over q-points
				for (int q = 0; q < reftofittingobject->m_data_q.size(); q++)
				{
					//Lookup Condition
					qr_old = reftofittingobject->m_data_q[q] * r_old * 1000 + 0.5;
					if (qr_old < sinc_lookup_size) sinc_qr_old = saxs::sinc_lookup[(int)qr_old];
					else sinc_qr_old = boost::math::sinc_pi(double(qr_old) / 1000);

					qr_new = reftofittingobject->m_data_q[q] * r_new * 1000 + 0.5;
					if (qr_new < sinc_lookup_size) sinc_qr_new = saxs::sinc_lookup[(int)qr_new];
					else sinc_qr_new = boost::math::sinc_pi(double(qr_new) / 1000);
					//Lookup Condition
					reftofittingobject->m_model_I[q] += 2 * (sinc_qr_new - sinc_qr_old)*double(reftofittingobject->m_num_stacks);
				}//enfor q
				//now adjust contacts
				if (r_old*r_old < reftofittingobject->m_contact_d_sq)
				{
					contact_vector[i] -= 1;
				}
				if (r_new*r_new < reftofittingobject->m_contact_d_sq)
				{
					contact_vector[i] += 1;
					contact_vector[changedcoord] += 1;
				}	
				 //Stack contribution
				if (reftofittingobject->m_num_stacks>1)
				{
					for (int j = 1; j < reftofittingobject->m_num_stacks; j++)
					{
						r_old = std::pow(
							std::pow(reftofittingobject->m_model[i]->m_x - reftofittingobject->m_model[changedcoord]->m_x, 2) +
							std::pow(reftofittingobject->m_model[i]->m_y - reftofittingobject->m_model[changedcoord]->m_y, 2) +
							std::pow(reftofittingobject->m_model[i]->m_z + double(j) * reftofittingobject->m_stack_spacing 
																				- reftofittingobject->m_model[changedcoord]->m_z, 2),
							0.5);
						r_new = std::pow(
							std::pow(reftofittingobject->m_model[i]->m_x - (reftofittingobject->m_model[changedcoord]->m_x + movement->m_x), 2) +
							std::pow(reftofittingobject->m_model[i]->m_y - (reftofittingobject->m_model[changedcoord]->m_y + movement->m_y), 2) +
							std::pow(reftofittingobject->m_model[i]->m_z + double(j) * reftofittingobject->m_stack_spacing 
																			- (reftofittingobject->m_model[changedcoord]->m_z + movement->m_z), 2),
							0.5);
						r_old_st = std::pow(
							std::pow(reftofittingobject->m_model[i]->m_x - reftofittingobject->m_model[changedcoord]->m_x, 2) +
							std::pow(reftofittingobject->m_model[i]->m_y - reftofittingobject->m_model[changedcoord]->m_y, 2) +
							std::pow(reftofittingobject->m_model[i]->m_z - double(j) * reftofittingobject->m_stack_spacing 
																				- reftofittingobject->m_model[changedcoord]->m_z, 2),
							0.5);
						r_new_st = std::pow(
							std::pow(reftofittingobject->m_model[i]->m_x - (reftofittingobject->m_model[changedcoord]->m_x + movement->m_x), 2) +
							std::pow(reftofittingobject->m_model[i]->m_y - (reftofittingobject->m_model[changedcoord]->m_y + movement->m_y), 2) +
							std::pow(reftofittingobject->m_model[i]->m_z - double(j) * reftofittingobject->m_stack_spacing 
																				- (reftofittingobject->m_model[changedcoord]->m_z + movement->m_z), 2),
							0.5);
						//loop over q-points
						for (int q = 0; q < reftofittingobject->m_data_q.size(); q++)
						{
							//Lookup Condition
							qr_old = reftofittingobject->m_data_q[q] * r_old * 1000 + 0.5;
							if (qr_old < sinc_lookup_size) sinc_qr_old = saxs::sinc_lookup[(int)qr_old];
							else sinc_qr_old = boost::math::sinc_pi(double(qr_old) / 1000);

							qr_new = reftofittingobject->m_data_q[q] * r_new * 1000 + 0.5;
							if (qr_new < sinc_lookup_size) sinc_qr_new = saxs::sinc_lookup[(int)qr_new];
							else sinc_qr_new = boost::math::sinc_pi(double(qr_new) / 1000);

							qr_old_st = reftofittingobject->m_data_q[q] * r_old_st * 1000 + 0.5;
							if (qr_old_st < sinc_lookup_size) sinc_qr_old_st = saxs::sinc_lookup[(int)qr_old_st];
							else sinc_qr_old_st = boost::math::sinc_pi(double(qr_old_st) / 1000);

							qr_new_st = reftofittingobject->m_data_q[q] * r_new_st * 1000 + 0.5;
							if (qr_new_st < sinc_lookup_size) sinc_qr_new_st = saxs::sinc_lookup[(int)qr_new_st];
							else sinc_qr_new_st = boost::math::sinc_pi(double(qr_new_st) / 1000);
							//Lookup Condition

							reftofittingobject->m_model_I[q] += 2 * (sinc_qr_new - sinc_qr_old)*double(reftofittingobject->m_num_stacks - j);

							reftofittingobject->m_model_I[q] += 2 * (sinc_qr_new_st - sinc_qr_old_st)*double(reftofittingobject->m_num_stacks - j);
						}//enfor q
					}//endfor stacks
				}//endif stack condition
			}//endif i!=changedcoord
		}//endfor i

		//Also recalc contacts for upper boundary atoms
		if (reftofittingobject->m_model[changedcoord]->m_z > (reftofittingobject->m_stack_spacing-std::sqrt(reftofittingobject->m_contact_d_sq))
			|| reftofittingobject->m_model[changedcoord]->m_z + movement->m_z > (reftofittingobject->m_stack_spacing - std::sqrt(reftofittingobject->m_contact_d_sq)))
		{
			for (int i = 0; i < reftofittingobject->m_model.size(); i++)
			{
				if (i != changedcoord)
				{
					r_old = std::pow(reftofittingobject->m_model[i]->m_x - reftofittingobject->m_model[changedcoord]->m_x, 2) +
						std::pow(reftofittingobject->m_model[i]->m_y - reftofittingobject->m_model[changedcoord]->m_y, 2) +
						std::pow(reftofittingobject->m_model[i]->m_z+ reftofittingobject->m_stack_spacing - reftofittingobject->m_model[changedcoord]->m_z, 2);
					r_new = std::pow(reftofittingobject->m_model[i]->m_x - (reftofittingobject->m_model[changedcoord]->m_x + movement->m_x), 2) +
						std::pow(reftofittingobject->m_model[i]->m_y - (reftofittingobject->m_model[changedcoord]->m_y + movement->m_y), 2) +
						std::pow(reftofittingobject->m_model[i]->m_z+ reftofittingobject->m_stack_spacing - (reftofittingobject->m_model[changedcoord]->m_z + movement->m_z), 2);
				if (r_old < reftofittingobject->m_contact_d_sq)
					{
						contact_vector[i] -= 1;
					}
					if (r_new < reftofittingobject->m_contact_d_sq)
					{
						contact_vector[i] += 1;
						contact_vector[changedcoord] += 1;
					}
				}
			}
		}
		//Also recalc contacts for lower boundary atoms
		if (reftofittingobject->m_model[changedcoord]->m_z <std::sqrt(reftofittingobject->m_contact_d_sq)
			|| reftofittingobject->m_model[changedcoord]->m_z + movement->m_z < (std::sqrt(reftofittingobject->m_contact_d_sq)))
		{
			for (int i = 0; i < reftofittingobject->m_model.size(); i++)
			{
				if (i != changedcoord)
				{
					r_old = std::pow(reftofittingobject->m_model[i]->m_x - reftofittingobject->m_model[changedcoord]->m_x, 2) +
						std::pow(reftofittingobject->m_model[i]->m_y - reftofittingobject->m_model[changedcoord]->m_y, 2) +
						std::pow(reftofittingobject->m_model[i]->m_z - reftofittingobject->m_stack_spacing - reftofittingobject->m_model[changedcoord]->m_z, 2);
					r_new = std::pow(reftofittingobject->m_model[i]->m_x - (reftofittingobject->m_model[changedcoord]->m_x + movement->m_x), 2) +
						std::pow(reftofittingobject->m_model[i]->m_y - (reftofittingobject->m_model[changedcoord]->m_y + movement->m_y), 2) +
						std::pow(reftofittingobject->m_model[i]->m_z - reftofittingobject->m_stack_spacing - (reftofittingobject->m_model[changedcoord]->m_z + movement->m_z), 2);
					if (r_old < reftofittingobject->m_contact_d_sq)
					{
						contact_vector[i] -= 1;
					}
					if (r_new < reftofittingobject->m_contact_d_sq)
					{
						contact_vector[i] += 1;
						contact_vector[changedcoord] += 1;
					}
				}
			}
		}

	
		//Now recalc compactness
		x = reftofittingobject->m_model[changedcoord]->m_x;
		y = reftofittingobject->m_model[changedcoord]->m_y;
		r_old = std::sqrt(x*x + y*y);
		reftofittingobject->m_compactness += -(1.0 - std::erfc(
			-0.5 - (1. - r_old / reftofittingobject->m_diameter / 2.)*reftofittingobject->m_sigmaComp*3.1415) / 2.0)/ double(reftofittingobject->m_model.size())
											+ (1.0 - std::erfc(
			-0.5 - (1. - r_new / reftofittingobject->m_diameter / 2.)*reftofittingobject->m_sigmaComp*3.1415) / 2.0)/double(reftofittingobject->m_model.size());
	
	}

}	//End of namespace

	//------------------------------------------------------------------------------
	//	Thread: Stop Signal handler
	//------------------------------------------------------------------------------
	void DebyeFitThread::stop()
	{
		stop_thread = true;
	}

	//------------------------------------------------------------------------------
	//	Core: Node function using qthread worker
	//------------------------------------------------------------------------------
	void HelFitQt::startDebyeFitCurrentModel()
	{
		//Preparing tracer object: Clear and initilaizes
		int num_threads;
		if (globalfittingobject->m_multicore) num_threads = 1;
		else num_threads = globalfittingobject->m_num_cores;
		result_tracer.clear();
		for (unsigned int i = 0; i < num_threads; i++)
		{
			result_tracer.push_back(saxs::model_sp(new saxs::model));

			result_tracer[i]->m_chi = globalfittingobject->m_chi;
			result_tracer[i]->m_chi_tracer.push_back(globalfittingobject->m_chi);

			result_tracer[i]->m_mean_nr_contacts = globalfittingobject->m_mean_nr_contacts;
			result_tracer[i]->m_contact_tracer.push_back(globalfittingobject->m_mean_nr_contacts);

			result_tracer[i]->m_mean_connectivity = globalfittingobject->m_mean_connectivity;
			result_tracer[i]->m_connectivity_tracer.push_back(globalfittingobject->m_mean_connectivity);

			result_tracer[i]->m_compactness = globalfittingobject->m_compactness;
			result_tracer[i]->m_compactness_tracer.push_back(globalfittingobject->m_compactness);

			result_tracer[i]->m_targetf = globalfittingobject->m_chi*globalfittingobject->m_chi;
			result_tracer[i]->m_target_f_tracer.push_back(globalfittingobject->m_compactness);

			result_tracer[i]->m_thread_status = false;
			for (int j = 0; j < globalfittingobject->m_model.size(); j++)
			{
				result_tracer[i]->m_coordinate.push_back(saxs::coordinate_sp(new saxs::coordinate(
					globalfittingobject->m_model[j]->m_x,
					globalfittingobject->m_model[j]->m_y,
					globalfittingobject->m_model[j]->m_z,
					globalfittingobject->m_model[j]->m_nr_contacts)));
			}
			for (int j = 0; j < globalfittingobject->m_model_I.size(); j++)
			{
				result_tracer[i]->m_model_I.push_back(0.0f);
				result_tracer[i]->m_fitted_I.push_back(0.0f);
			}
		}

		//Now everything is ready for calculation
		//Instantiate threadvector
		std::vector<DebyeFitThread*> thread_vector;

		//Registering SignalHandling datatypes
		qRegisterMetaType<std::string>();
		qRegisterMetaType<QVector<double>>();



		//Writing variables into thread-object
		for (int i = 0; i < num_threads; i++)
		{
			thread_vector.push_back(new DebyeFitThread);
			thread_vector[i]->reftofittingobject = globalfittingobject;
			thread_vector[i]->result_tracer = result_tracer;
			thread_vector[i]->n_runs = globalfittingobject->m_num_runs;
			thread_vector[i]->starttemp = globalfittingobject->m_start_temp;
			thread_vector[i]->endtemp = 0;
			thread_vector[i]->deltatemp = globalfittingobject->m_delta_temp;
			thread_vector[i]->coreid = i;
			//Connecting QT Signals
			QObject::connect(thread_vector[i], SIGNAL(sig_writelogext(std::string)), this, SLOT(writetologext(std::string)));
			QObject::connect(thread_vector[i], SIGNAL(sig_plot_current_thread_Data(QVector<double>)), this, SLOT(plot_current_thread_Data(QVector<double>)));
			QObject::connect(thread_vector[i], SIGNAL(sig_fittingthread_done(int)), this, SLOT(fittingthread_done(int)));
			QObject::connect(thread_vector[i], SIGNAL(sig_update_statusbar(double)), this, SLOT(update_statusbar(double)));
			QObject::connect(this, SIGNAL(sig_stop_thread()), thread_vector[i], SLOT(stop()));
			QObject::connect(thread_vector[i], SIGNAL(sig_live_update_model(QVector<double>, QVector<double>, QVector<double>)),
								this, SLOT(live_update_model(QVector<double>, QVector<double>, QVector<double>)));
			//Connect core zero to plotting thread
			if (i == 0) QObject::connect(this, SIGNAL(sig_ui_ready()), thread_vector[i], SLOT(receive_ui_ready()));
		}

		//Start Threads
		for (int i = 0; i < num_threads; i++)
		{
			thread_vector[i]->start();
		}
	}

	//------------------------------------------------------------------------------
	//	Thread: ui ready for graphing
	//------------------------------------------------------------------------------
	void DebyeFitThread::receive_ui_ready()
	{
		ui_ready = true;
	}




	/***************************************************************/
	// Multicore DebyeFit Functions
	/***************************************************************/
	namespace saxs {

	//------------------------------------------------------------------------------
	// Recalculates the scattering pattern if changecood is moved by movement - Result will be put into references new_model_i
	//------------------------------------------------------------------------------
	void mc_recalc_change_I(fittingobject_sp& reftofittingobject, int& changedcoord, coordinate_sp& movement,
		const int start_index, const int end_index)
	{
		//tempvars
		double r_old = 0;
		double qr_old = 0;
		double sinc_qr_old = 0;
		double r_new = 0;
		double qr_new = 0;
		double sinc_qr_new = 0;
		//tempvars for stack
		double r_old_st = 0;
		double qr_old_st = 0;
		double sinc_qr_old_st = 0;
		double r_new_st = 0;
		double qr_new_st = 0;
		double sinc_qr_new_st = 0;
		int sinc_lookup_size = saxs::sinc_lookup.size();
		double x, y;
		double coord_dsq = std::pow(reftofittingobject->m_diameter / 4., 2);

		for (int i = 0; i < reftofittingobject->m_model.size(); i++)
		{
			if (i != changedcoord)
			{
				//Building block contribution
				r_old = std::pow(
					std::pow(reftofittingobject->m_model[i]->m_x - reftofittingobject->m_model[changedcoord]->m_x, 2) +
					std::pow(reftofittingobject->m_model[i]->m_y - reftofittingobject->m_model[changedcoord]->m_y, 2) +
					std::pow(reftofittingobject->m_model[i]->m_z - reftofittingobject->m_model[changedcoord]->m_z, 2),
					0.5);
				r_new = std::pow(
					std::pow(reftofittingobject->m_model[i]->m_x - (reftofittingobject->m_model[changedcoord]->m_x + movement->m_x), 2) +
					std::pow(reftofittingobject->m_model[i]->m_y - (reftofittingobject->m_model[changedcoord]->m_y + movement->m_y), 2) +
					std::pow(reftofittingobject->m_model[i]->m_z - (reftofittingobject->m_model[changedcoord]->m_z + movement->m_z), 2),
					0.5);
				//loop over q-points
				for (int q = start_index; q < end_index; q++)
				{
					//Lookup Condition
					qr_old = reftofittingobject->m_data_q[q] * r_old * 1000 + 0.5;
					if (qr_old < sinc_lookup_size) sinc_qr_old = saxs::sinc_lookup[(int)qr_old];
					else sinc_qr_old = boost::math::sinc_pi(double(qr_old) / 1000);

					qr_new = reftofittingobject->m_data_q[q] * r_new * 1000 + 0.5;
					if (qr_new < sinc_lookup_size) sinc_qr_new = saxs::sinc_lookup[(int)qr_new];
					else sinc_qr_new = boost::math::sinc_pi(double(qr_new) / 1000);
					//Lookup Condition
					reftofittingobject->m_model_I[q] += 2 * (sinc_qr_new - sinc_qr_old)*double(reftofittingobject->m_num_stacks);
				}//enfor q
				//Stack contribution
				if (reftofittingobject->m_num_stacks>1)
				{
					for (int j = 1; j < reftofittingobject->m_num_stacks; j++)
					{
						r_old = std::pow(
							std::pow(reftofittingobject->m_model[i]->m_x - reftofittingobject->m_model[changedcoord]->m_x, 2) +
							std::pow(reftofittingobject->m_model[i]->m_y - reftofittingobject->m_model[changedcoord]->m_y, 2) +
							std::pow(reftofittingobject->m_model[i]->m_z + double(j) * reftofittingobject->m_stack_spacing
								- reftofittingobject->m_model[changedcoord]->m_z, 2),
							0.5);
						r_new = std::pow(
							std::pow(reftofittingobject->m_model[i]->m_x - (reftofittingobject->m_model[changedcoord]->m_x + movement->m_x), 2) +
							std::pow(reftofittingobject->m_model[i]->m_y - (reftofittingobject->m_model[changedcoord]->m_y + movement->m_y), 2) +
							std::pow(reftofittingobject->m_model[i]->m_z + double(j) * reftofittingobject->m_stack_spacing
								- (reftofittingobject->m_model[changedcoord]->m_z + movement->m_z), 2),
							0.5);
						r_old_st = std::pow(
							std::pow(reftofittingobject->m_model[i]->m_x - reftofittingobject->m_model[changedcoord]->m_x, 2) +
							std::pow(reftofittingobject->m_model[i]->m_y - reftofittingobject->m_model[changedcoord]->m_y, 2) +
							std::pow(reftofittingobject->m_model[i]->m_z - double(j) * reftofittingobject->m_stack_spacing
								- reftofittingobject->m_model[changedcoord]->m_z, 2),
							0.5);
						r_new_st = std::pow(
							std::pow(reftofittingobject->m_model[i]->m_x - (reftofittingobject->m_model[changedcoord]->m_x + movement->m_x), 2) +
							std::pow(reftofittingobject->m_model[i]->m_y - (reftofittingobject->m_model[changedcoord]->m_y + movement->m_y), 2) +
							std::pow(reftofittingobject->m_model[i]->m_z - double(j) * reftofittingobject->m_stack_spacing
								- (reftofittingobject->m_model[changedcoord]->m_z + movement->m_z), 2),
							0.5);
						//loop over q-points
						for (int q = start_index; q < end_index; q++)
						{
							//Lookup Condition
							qr_old = reftofittingobject->m_data_q[q] * r_old * 1000 + 0.5;
							if (qr_old < sinc_lookup_size) sinc_qr_old = saxs::sinc_lookup[(int)qr_old];
							else sinc_qr_old = boost::math::sinc_pi(double(qr_old) / 1000);

							qr_new = reftofittingobject->m_data_q[q] * r_new * 1000 + 0.5;
							if (qr_new < sinc_lookup_size) sinc_qr_new = saxs::sinc_lookup[(int)qr_new];
							else sinc_qr_new = boost::math::sinc_pi(double(qr_new) / 1000);

							qr_old_st = reftofittingobject->m_data_q[q] * r_old_st * 1000 + 0.5;
							if (qr_old_st < sinc_lookup_size) sinc_qr_old_st = saxs::sinc_lookup[(int)qr_old_st];
							else sinc_qr_old_st = boost::math::sinc_pi(double(qr_old_st) / 1000);

							qr_new_st = reftofittingobject->m_data_q[q] * r_new_st * 1000 + 0.5;
							if (qr_new_st < sinc_lookup_size) sinc_qr_new_st = saxs::sinc_lookup[(int)qr_new_st];
							else sinc_qr_new_st = boost::math::sinc_pi(double(qr_new_st) / 1000);
							//Lookup Condition

							reftofittingobject->m_model_I[q] += 2 * (sinc_qr_new - sinc_qr_old)*double(reftofittingobject->m_num_stacks - j);

							reftofittingobject->m_model_I[q] += 2 * (sinc_qr_new_st - sinc_qr_old_st)*double(reftofittingobject->m_num_stacks - j);
						}//enfor q
					}//endfor stacks
				}//endif stack condition
			}//endif i!=changedcoord
		}//endfor i
	}

	void mc_recalc_change_regparams(fittingobject_sp& reftofittingobject, int& changedcoord, coordinate_sp& movement,
		std::vector<int>& contact_vector)
	{
		//tempvars
		double r_old = 0;
		double qr_old = 0;
		double sinc_qr_old = 0;
		double r_new = 0;
		double qr_new = 0;
		double sinc_qr_new = 0;
		//tempvars for stack
		double r_old_st = 0;
		double qr_old_st = 0;
		double sinc_qr_old_st = 0;
		double r_new_st = 0;
		double qr_new_st = 0;
		double sinc_qr_new_st = 0;
		int sinc_lookup_size = saxs::sinc_lookup.size();
		double x, y;

		//reset contactcount on changed atom
		contact_vector[changedcoord] = 0;

		for (int i = 0; i < reftofittingobject->m_model.size(); i++)
		{
			if (i != changedcoord)
			{
				//Building block contribution
				r_old = std::pow(
					std::pow(reftofittingobject->m_model[i]->m_x - reftofittingobject->m_model[changedcoord]->m_x, 2) +
					std::pow(reftofittingobject->m_model[i]->m_y - reftofittingobject->m_model[changedcoord]->m_y, 2) +
					std::pow(reftofittingobject->m_model[i]->m_z - reftofittingobject->m_model[changedcoord]->m_z, 2),
					0.5);
				r_new = std::pow(
					std::pow(reftofittingobject->m_model[i]->m_x - (reftofittingobject->m_model[changedcoord]->m_x + movement->m_x), 2) +
					std::pow(reftofittingobject->m_model[i]->m_y - (reftofittingobject->m_model[changedcoord]->m_y + movement->m_y), 2) +
					std::pow(reftofittingobject->m_model[i]->m_z - (reftofittingobject->m_model[changedcoord]->m_z + movement->m_z), 2),
					0.5);
				 //now adjust contacts
				if (r_old*r_old < reftofittingobject->m_contact_d_sq)
				{
					contact_vector[i] -= 1;
				}
				if (r_new*r_new < reftofittingobject->m_contact_d_sq)
				{
					contact_vector[i] += 1;
					contact_vector[changedcoord] += 1;
				}
				//Stack contribution
				if (reftofittingobject->m_num_stacks>1)
				{
					for (int j = 1; j < reftofittingobject->m_num_stacks; j++)
					{
						r_old = std::pow(
							std::pow(reftofittingobject->m_model[i]->m_x - reftofittingobject->m_model[changedcoord]->m_x, 2) +
							std::pow(reftofittingobject->m_model[i]->m_y - reftofittingobject->m_model[changedcoord]->m_y, 2) +
							std::pow(reftofittingobject->m_model[i]->m_z + double(j) * reftofittingobject->m_stack_spacing
								- reftofittingobject->m_model[changedcoord]->m_z, 2),
							0.5);
						r_new = std::pow(
							std::pow(reftofittingobject->m_model[i]->m_x - (reftofittingobject->m_model[changedcoord]->m_x + movement->m_x), 2) +
							std::pow(reftofittingobject->m_model[i]->m_y - (reftofittingobject->m_model[changedcoord]->m_y + movement->m_y), 2) +
							std::pow(reftofittingobject->m_model[i]->m_z + double(j) * reftofittingobject->m_stack_spacing
								- (reftofittingobject->m_model[changedcoord]->m_z + movement->m_z), 2),
							0.5);
						r_old_st = std::pow(
							std::pow(reftofittingobject->m_model[i]->m_x - reftofittingobject->m_model[changedcoord]->m_x, 2) +
							std::pow(reftofittingobject->m_model[i]->m_y - reftofittingobject->m_model[changedcoord]->m_y, 2) +
							std::pow(reftofittingobject->m_model[i]->m_z - double(j) * reftofittingobject->m_stack_spacing
								- reftofittingobject->m_model[changedcoord]->m_z, 2),
							0.5);
						r_new_st = std::pow(
							std::pow(reftofittingobject->m_model[i]->m_x - (reftofittingobject->m_model[changedcoord]->m_x + movement->m_x), 2) +
							std::pow(reftofittingobject->m_model[i]->m_y - (reftofittingobject->m_model[changedcoord]->m_y + movement->m_y), 2) +
							std::pow(reftofittingobject->m_model[i]->m_z - double(j) * reftofittingobject->m_stack_spacing
								- (reftofittingobject->m_model[changedcoord]->m_z + movement->m_z), 2),
							0.5);
					}//endfor stacks
				}//endif stack condition
			}//endif i!=changedcoord
		}//endfor i

		 //Also recalc contacts for upper boundary atoms
		if (reftofittingobject->m_model[changedcoord]->m_z > (reftofittingobject->m_stack_spacing - std::sqrt(reftofittingobject->m_contact_d_sq))
			|| reftofittingobject->m_model[changedcoord]->m_z + movement->m_z > (reftofittingobject->m_stack_spacing - std::sqrt(reftofittingobject->m_contact_d_sq)))
		{
			for (int i = 0; i < reftofittingobject->m_model.size(); i++)
			{
				if (i != changedcoord)
				{
					r_old = std::pow(reftofittingobject->m_model[i]->m_x - reftofittingobject->m_model[changedcoord]->m_x, 2) +
						std::pow(reftofittingobject->m_model[i]->m_y - reftofittingobject->m_model[changedcoord]->m_y, 2) +
						std::pow(reftofittingobject->m_model[i]->m_z + reftofittingobject->m_stack_spacing - reftofittingobject->m_model[changedcoord]->m_z, 2);
					r_new = std::pow(reftofittingobject->m_model[i]->m_x - (reftofittingobject->m_model[changedcoord]->m_x + movement->m_x), 2) +
						std::pow(reftofittingobject->m_model[i]->m_y - (reftofittingobject->m_model[changedcoord]->m_y + movement->m_y), 2) +
						std::pow(reftofittingobject->m_model[i]->m_z + reftofittingobject->m_stack_spacing - (reftofittingobject->m_model[changedcoord]->m_z + movement->m_z), 2);
					if (r_old < reftofittingobject->m_contact_d_sq)
					{
						contact_vector[i] -= 1;
					}
					if (r_new < reftofittingobject->m_contact_d_sq)
					{
						contact_vector[i] += 1;
						contact_vector[changedcoord] += 1;
					}
				}
			}
		}
		//Also recalc contacts for lower boundary atoms
		if (reftofittingobject->m_model[changedcoord]->m_z <std::sqrt(reftofittingobject->m_contact_d_sq)
			|| reftofittingobject->m_model[changedcoord]->m_z + movement->m_z < (std::sqrt(reftofittingobject->m_contact_d_sq)))
		{
			for (int i = 0; i < reftofittingobject->m_model.size(); i++)
			{
				if (i != changedcoord)
				{
					r_old = std::pow(reftofittingobject->m_model[i]->m_x - reftofittingobject->m_model[changedcoord]->m_x, 2) +
						std::pow(reftofittingobject->m_model[i]->m_y - reftofittingobject->m_model[changedcoord]->m_y, 2) +
						std::pow(reftofittingobject->m_model[i]->m_z - reftofittingobject->m_stack_spacing - reftofittingobject->m_model[changedcoord]->m_z, 2);
					r_new = std::pow(reftofittingobject->m_model[i]->m_x - (reftofittingobject->m_model[changedcoord]->m_x + movement->m_x), 2) +
						std::pow(reftofittingobject->m_model[i]->m_y - (reftofittingobject->m_model[changedcoord]->m_y + movement->m_y), 2) +
						std::pow(reftofittingobject->m_model[i]->m_z - reftofittingobject->m_stack_spacing - (reftofittingobject->m_model[changedcoord]->m_z + movement->m_z), 2);
					if (r_old < reftofittingobject->m_contact_d_sq)
					{
						contact_vector[i] -= 1;
					}
					if (r_new < reftofittingobject->m_contact_d_sq)
					{
						contact_vector[i] += 1;
						contact_vector[changedcoord] += 1;
					}
				}
			}
		}
		//Now recalc compactness
		x = reftofittingobject->m_model[changedcoord]->m_x;
		y = reftofittingobject->m_model[changedcoord]->m_y;
		r_old = std::sqrt(x*x + y*y);
		reftofittingobject->m_compactness += -(1.0 - std::erfc(
			-0.5 - (1. - r_old / reftofittingobject->m_diameter / 2.)*reftofittingobject->m_sigmaComp*3.1415) / 2.0) / double(reftofittingobject->m_model.size())
			+ (1.0 - std::erfc(
				-0.5 - (1. - r_new / reftofittingobject->m_diameter / 2.)*reftofittingobject->m_sigmaComp*3.1415) / 2.0) / double(reftofittingobject->m_model.size());

	}

	}

	//------------------------------------------------------------------------------
	//	Thread: Fitting function  function
	//------------------------------------------------------------------------------
	void DebyeFitThread::run()
	{
		//Make local copy of fitting object so changes don't affect global object
		saxs::fittingobject localfittingobject = *reftofittingobject;
		boost::shared_ptr<saxs::fittingobject> localfittingobject_sp = boost::make_shared<saxs::fittingobject>(localfittingobject);

		//Deepcopy model
		localfittingobject_sp->m_model.clear();
		for (int i = 0; i < reftofittingobject->m_model.size(); i++)
			localfittingobject_sp->m_model.push_back(saxs::coordinate_sp(new saxs::coordinate(reftofittingobject->m_model[i]->m_x,
				reftofittingobject->m_model[i]->m_y,
				reftofittingobject->m_model[i]->m_z,
				reftofittingobject->m_model[i]->m_nr_contacts)));

		//Deppcopy data
		int datasize = localfittingobject_sp->m_model_I.size();
		localfittingobject_sp->m_model_I.clear();
		localfittingobject_sp->m_fitted_I.clear();
		localfittingobject_sp->m_data_I.clear();
		for (int i = 0; i < datasize; i++)
		{
			localfittingobject_sp->m_model_I.push_back(reftofittingobject->m_model_I[i]);
			localfittingobject_sp->m_fitted_I.push_back(reftofittingobject->m_fitted_I[i]);
			localfittingobject_sp->m_data_I.push_back(reftofittingobject->m_data_I[i]);
		}

		//Thread container
		boost::scoped_ptr<boost::thread_group> threadGroup_sp(new boost::thread_group);

		//Mutex to protect IO operations
		boost::mutex io_mutex;

		//Initialize thread parameters
		std::vector<int> boundaries;
		if (localfittingobject_sp->m_multicore)
		{
			int thread_step = static_cast<int>(localfittingobject_sp->m_data_q.size() / (localfittingobject_sp->m_num_cores-1));
			int index_delta = localfittingobject_sp->m_data_q.size() - thread_step*(localfittingobject_sp->m_num_cores - 1);
			//Taking into account index offset due to core number
			//Putting indexboundaries into boundaries vector of size cores+1	
			boundaries.push_back(0);
			int offset = 1;
			for (int i = 0; i < (localfittingobject_sp->m_num_cores - 1); i++)
			{
				if (i < index_delta) offset = 1;
				else offset = 0;
				boundaries.push_back(boundaries[i] + thread_step + offset);
			}
		}
		//Define helper vector for temporary storage of results
		std::vector<double> old_model_i;
		old_model_i.resize(localfittingobject_sp->m_model_I.size());

		//Define temporary QVector for data transfer to GUI
		QVector<double> send_fitted_I;
		send_fitted_I.resize(localfittingobject_sp->m_fitted_I.size());

		//Define temporary QVectors for model transfer to GUI
		QVector<double> qvec_model_x;
		qvec_model_x.resize(localfittingobject_sp->m_model.size());
		QVector<double> qvec_model_y;
		qvec_model_y.resize(localfittingobject_sp->m_model.size());
		QVector<double> qvec_model_z;
		qvec_model_z.resize(localfittingobject_sp->m_model.size());

		//Define helper vector for temporary contact storage
		std::vector<int> nr_contacts_vector;
		std::vector<int> nr_old_contacts_vector;
		for (int i = 0; i < localfittingobject_sp->m_model.size(); i++)
		{
			nr_contacts_vector.push_back(localfittingobject_sp->m_model[i]->m_nr_contacts);
			nr_old_contacts_vector.push_back(localfittingobject_sp->m_model[i]->m_nr_contacts);
		}



		//Initialize futher variables
		saxs::coordinate move = (0, 0, 0);
		boost::shared_ptr<saxs::coordinate> movement = boost::make_shared<saxs::coordinate>(move);
		double x, y,z;
		double chi_start = localfittingobject_sp->m_chi;
		double old_chi;
		double old_connectivity;
		double old_compactness;
		double old_coordination;
		double temp_z_coord;
		double temp;
		double old_target_f;
		double new_target_f;
		double connectivityweight;
		double compactnessweight;

		//Potential Field helpers
		double x_pot, y_pot;

		//Output helpers
		std::string str;
		std::stringstream ss;

		for (int i = 0; i < n_runs; i++)
		{
			//Calculate current Temperature
			temp = endtemp + (starttemp - endtemp)*std::pow(deltatemp, i);

			//heating up randomization
			saxs::randseed = (unsigned)time(NULL);
			double rand_scalar_r = temp*localfittingobject_sp->m_diameter;
			double rand_scalar_h = temp*localfittingobject_sp->m_stack_spacing / 2.;
			if (localfittingobject_sp->m_num_stacks>1)rand_scalar_h /= 10.;

			//Every 5th iteration, all atoms outside a critical diameter are kicked pack to the center by 2/3
			if (i % 2 == 0 && i != 0 || (i + 1) == n_runs)
			{
				//Adding random movement every 10th step
				if (i % 10 == 0)
				{
					for (int j = 0; j < localfittingobject_sp->m_model.size(); j++)
					{
						localfittingobject_sp->m_model[j]->m_x += 0.4*rand_scalar_r*0.5*(0.5 - saxs::xor128());
						localfittingobject_sp->m_model[j]->m_y += 0.4*rand_scalar_r*0.5*(0.5 - saxs::xor128());
						if (localfittingobject_sp->m_num_stacks > 1)
							localfittingobject_sp->m_model[j]->m_z += rand_scalar_h*(0.5 - saxs::xor128());
						else
							localfittingobject_sp->m_model[j]->m_z += 0.4*rand_scalar_h*(0.5 - saxs::xor128());
					}
				}

				//Moving atoms with no neighbours
				if ((i + 5) % 10 == 0)
				{
					for (int j = 0; j < localfittingobject_sp->m_model.size(); j++)
					{
						if (1 > localfittingobject_sp->m_model[j]->m_nr_contacts)
						{
							movement->m_x = localfittingobject_sp->m_model[j]->m_x / 1.5 - localfittingobject_sp->m_model[j]->m_x;
							movement->m_y = localfittingobject_sp->m_model[j]->m_y / 1.5 - localfittingobject_sp->m_model[j]->m_y;
							if (localfittingobject_sp->m_num_stacks > 1)
								movement->m_z =0;
							else
								movement->m_z = localfittingobject_sp->m_model[j]->m_z / 1.5 - localfittingobject_sp->m_model[j]->m_z;
							localfittingobject_sp->m_model[j]->m_x += movement->m_x;
							localfittingobject_sp->m_model[j]->m_y += movement->m_y;
							localfittingobject_sp->m_model[j]->m_z += movement->m_z;
						}
					}
				}

				//Critical diameter check
				if (localfittingobject_sp->m_num_stacks > 1)
					for (int j = 0; j < localfittingobject_sp->m_model.size(); j++)
					{
						x = localfittingobject_sp->m_model[j]->m_x;
						y = localfittingobject_sp->m_model[j]->m_y;
						if (x*x + y*y > localfittingobject_sp->m_diameter*localfittingobject_sp->m_diameter / 2.)
						{
							movement->m_x = localfittingobject_sp->m_model[j]->m_x / 3. - localfittingobject_sp->m_model[j]->m_x;
							movement->m_y = localfittingobject_sp->m_model[j]->m_y / 3. - localfittingobject_sp->m_model[j]->m_y;					
							localfittingobject_sp->m_model[j]->m_x += movement->m_x;
							localfittingobject_sp->m_model[j]->m_y += movement->m_y;
						}
					}
				else
				{
					for (int j = 0; j < localfittingobject_sp->m_model.size(); j++)
					{
						x = localfittingobject_sp->m_model[j]->m_x;
						y = localfittingobject_sp->m_model[j]->m_y;
						z = localfittingobject_sp->m_model[j]->m_z;
						if (x*x + y*y + z*z> localfittingobject_sp->m_diameter*localfittingobject_sp->m_diameter)
						{
							movement->m_x = localfittingobject_sp->m_model[j]->m_x / 3. - localfittingobject_sp->m_model[j]->m_x;
							movement->m_y = localfittingobject_sp->m_model[j]->m_y / 3. - localfittingobject_sp->m_model[j]->m_y;
							movement->m_z = localfittingobject_sp->m_model[j]->m_z / 3. - localfittingobject_sp->m_model[j]->m_z;
							localfittingobject_sp->m_model[j]->m_x += movement->m_x;
							localfittingobject_sp->m_model[j]->m_y += movement->m_y;
							localfittingobject_sp->m_model[j]->m_z += movement->m_z;
						}
					}
				}

				//Recalculation of scattering curve
				if (localfittingobject_sp->m_multicore)
				{
					//Instantiate and start threads
					for (unsigned int k = 0; k<(localfittingobject_sp->m_num_cores - 1); k++)
					{
						threadGroup_sp->add_thread(new boost::thread(saxs::calculate_debye, boost::ref(io_mutex), boundaries[k], boundaries[k + 1], localfittingobject_sp));
					}

					//Wait until the end of their jobs
					if (threadGroup_sp)
						threadGroup_sp->join_all();
				}
				else
					saxs::calculate_debye(io_mutex, 0, localfittingobject_sp->m_model_I.size(), localfittingobject_sp);

				saxs::fit_modeltoexp(localfittingobject_sp);
				//This causes the system to be seen as overcomapcted => new calculation of critical cylinder diameter
				localfittingobject_sp->m_compactness = 0.001;
				localfittingobject_sp->m_diameter = localfittingobject_sp->m_diameter / 0.95;
			}

			// Dynamic adjustment of critical cylinder diameter
			if (localfittingobject_sp->m_compactness < 0.3)
			{
				localfittingobject_sp->m_diameter = localfittingobject_sp->m_diameter*0.95;
				localfittingobject_sp->m_contact_d_sq = std::pow(2 * saxs::get_critradius_from_coordinates(localfittingobject_sp->m_diameter,
					localfittingobject_sp->m_stack_spacing, localfittingobject_sp->m_model.size()), 2);
				saxs::calc_contacts_of_model(localfittingobject_sp);
				localfittingobject_sp->m_mean_nr_contacts = saxs::return_mean_contacts_of_model(localfittingobject_sp->m_model);
				localfittingobject_sp->m_mean_connectivity = saxs::return_mean_connectivity_of_fittingobject(localfittingobject_sp);
				localfittingobject_sp->m_compactness = saxs::return_compactness_of_fittingobject(localfittingobject_sp);
				for (int i = 0; i < localfittingobject_sp->m_model.size(); i++)
				{
					nr_contacts_vector[i] = (localfittingobject_sp->m_model[i]->m_nr_contacts);
				}
			}
			if (localfittingobject_sp->m_compactness > 0.8)
			{
				localfittingobject_sp->m_diameter = localfittingobject_sp->m_diameter* 1.05;
				localfittingobject_sp->m_contact_d_sq = std::pow(2 * saxs::get_critradius_from_coordinates(localfittingobject_sp->m_diameter,
					localfittingobject_sp->m_stack_spacing, localfittingobject_sp->m_model.size()), 2);
				saxs::calc_contacts_of_model(localfittingobject_sp);
				localfittingobject_sp->m_mean_nr_contacts = saxs::return_mean_contacts_of_model(localfittingobject_sp->m_model);
				localfittingobject_sp->m_mean_connectivity = saxs::return_mean_connectivity_of_fittingobject(localfittingobject_sp);
				localfittingobject_sp->m_compactness = saxs::return_compactness_of_fittingobject(localfittingobject_sp);
				for (int i = 0; i < localfittingobject_sp->m_model.size(); i++)
				{
					nr_contacts_vector[i] = (localfittingobject_sp->m_model[i]->m_nr_contacts);
				}
			}

			//Adjust NextNeighbour distance such that not more than 30% have more than 12 neighbours
			while (saxs::return_fraction_contactsaturated_das(nr_contacts_vector) > 0.3)
			{
				localfittingobject_sp->m_contact_d_sq = localfittingobject_sp->m_contact_d_sq*0.95;
				saxs::calc_contacts_of_model(localfittingobject_sp);
				localfittingobject_sp->m_mean_nr_contacts = saxs::return_mean_contacts_of_model(localfittingobject_sp->m_model);
				localfittingobject_sp->m_mean_connectivity = saxs::return_mean_connectivity_of_fittingobject(localfittingobject_sp);
				localfittingobject_sp->m_compactness = saxs::return_compactness_of_fittingobject(localfittingobject_sp);
				for (int i = 0; i < localfittingobject_sp->m_model.size(); i++)
				{
					nr_contacts_vector[i] = (localfittingobject_sp->m_model[i]->m_nr_contacts);
				}
			}

			//Adjust Connectivity and COmpactness Weights
			compactnessweight = localfittingobject_sp->m_betaComp * localfittingobject_sp->m_chi;
			connectivityweight = localfittingobject_sp->m_chi * localfittingobject_sp->m_alphaConn;

			//Now start loop over all DAs
			for (int j = 0; j < reftofittingobject->m_model.size(); j++)
			{
				//Make copy of current status
				old_model_i = localfittingobject_sp->m_model_I;
				old_connectivity = localfittingobject_sp->m_mean_connectivity;
				old_compactness = localfittingobject_sp->m_compactness;
				nr_old_contacts_vector = nr_contacts_vector;
				old_chi = localfittingobject_sp->m_chi;
				old_target_f = saxs::tartget_function(old_chi, old_connectivity, connectivityweight,
					old_compactness, compactnessweight, localfittingobject_sp);

				//Random movement using potential field
				x_pot = std::cos(2.*saxs::pi / localfittingobject_sp->m_stack_spacing*localfittingobject_sp->m_model[j]->m_z);
				y_pot = std::sin(2.*saxs::pi / localfittingobject_sp->m_stack_spacing*localfittingobject_sp->m_model[j]->m_z);

				//Generate new random movement
				//Was 0.3
				movement->m_x = rand_scalar_r*(0.2*temp * x_pot + 1*(0.5 - saxs::xor128()));
				movement->m_y = rand_scalar_r*(0.2*temp * y_pot + 1*(0.5 - saxs::xor128()));
				movement->m_z = rand_scalar_h*(0.5 - saxs::xor128());

				//Recalc changes and contacts

				if (localfittingobject_sp->m_multicore)
				{
					//Instantiate and start threads
					threadGroup_sp->add_thread(new boost::thread(saxs::mc_recalc_change_regparams, localfittingobject_sp, j, boost::ref(movement),
						nr_contacts_vector));
					for (unsigned int k = 0; k<(localfittingobject_sp->m_num_cores - 1); k++)
					{
						threadGroup_sp->add_thread(new boost::thread(saxs::mc_recalc_change_I, localfittingobject_sp, j, boost::ref(movement),
							boundaries[k], boundaries[k + 1]));
					}

					//Wait until the end of their jobs
					if (threadGroup_sp)
						threadGroup_sp->join_all();
				}
				else
					//Recalc changes and contacts
					saxs::thread_recalc_change(localfittingobject_sp, j, movement, nr_contacts_vector);

				//Recalc fitted intesity using old lin-reg params
				saxs::recalc_fitted_i(localfittingobject_sp);

				//Return Connectivity
				localfittingobject_sp->m_mean_connectivity = saxs::return_mean_connectivity_of_contactvector(nr_contacts_vector, localfittingobject_sp);
				new_target_f = saxs::tartget_function(localfittingobject_sp->m_chi,
					localfittingobject_sp->m_mean_connectivity, connectivityweight,
					localfittingobject_sp->m_compactness, compactnessweight, localfittingobject_sp);

				//Make final comparison if movement is decreasing target function
				if (new_target_f < old_target_f)
				{
					//YES!! Now accept changes
					localfittingobject_sp->m_model[j]->m_x += movement->m_x;
					localfittingobject_sp->m_model[j]->m_y += movement->m_y;
					if (localfittingobject_sp->m_num_stacks>1)
					{
						temp_z_coord = localfittingobject_sp->m_model[j]->m_z + movement->m_z;
						if (temp_z_coord < 0)  movement->m_z += localfittingobject_sp->m_stack_spacing;
						if (temp_z_coord > localfittingobject_sp->m_stack_spacing)  movement->m_z -= localfittingobject_sp->m_stack_spacing;
					}
					localfittingobject_sp->m_model[j]->m_z += movement->m_z;
					old_target_f = new_target_f;
				}
				else
				{
					//NO!!! Restore old status
					localfittingobject_sp->m_model_I = old_model_i;
					localfittingobject_sp->m_chi = old_chi;
					localfittingobject_sp->m_mean_connectivity = old_connectivity;
					nr_contacts_vector = nr_old_contacts_vector;
					localfittingobject_sp->m_compactness = old_compactness;
				}

			}//endfor j coords

			 //Now all DAs have been moved

			 //Redo linear regression
			fit_modeltoexp(localfittingobject_sp);
			localfittingobject_sp->m_mean_nr_contacts = saxs::return_mean_contacts_of_contacvector(nr_contacts_vector);

			//Write status in tracer
			result_tracer[coreid]->m_chi_tracer.push_back(localfittingobject_sp->m_chi);
			result_tracer[coreid]->m_contact_tracer.push_back(localfittingobject_sp->m_mean_nr_contacts);
			result_tracer[coreid]->m_connectivity_tracer.push_back(localfittingobject_sp->m_mean_connectivity);
			result_tracer[coreid]->m_compactness_tracer.push_back(localfittingobject_sp->m_compactness);
			result_tracer[coreid]->m_target_f_tracer.push_back(old_target_f);

			//Report Progress to GUI
			if (coreid == 0)
			{
				//Generate output for log file in GUI
				ss.str(std::string());
				ss << "Step " << i << ": ";
				ss << std::scientific;
				ss << "chi = ";
				ss << localfittingobject_sp->m_chi;
				ss << ", connect. = ";
				ss << localfittingobject_sp->m_mean_connectivity;
				ss << ", comp. = ";
				ss << localfittingobject_sp->m_compactness;
				str = ss.str();
				emit sig_writelogext(str);


				//Check if ui is ready for processing
				if (ui_ready||i==0)
				{
					//Send 1D data to GUI
					send_fitted_I = QVector<double>::fromStdVector(localfittingobject_sp->m_fitted_I);
					

					for (int k = 0; k < localfittingobject_sp->m_model.size(); k++)
					{
						qvec_model_x[k] = localfittingobject_sp->m_model[k]->m_x;
						qvec_model_y[k] = localfittingobject_sp->m_model[k]->m_y;
						qvec_model_z[k] = localfittingobject_sp->m_model[k]->m_z;
					}
					emit sig_plot_current_thread_Data(send_fitted_I);
					emit sig_live_update_model(qvec_model_x, qvec_model_y, qvec_model_z);
					//Update Statusbar
					emit sig_update_statusbar(double(i) / double(n_runs));
					ui_ready = false;
				}	
			}

			//Now check if the Fitting process was aborted by user
			if (stop_thread)
			{
				if (coreid == 0)
				{
					str = "MANUAL ABORTION!!!";
					emit sig_writelogext(str);
				}
				emit sig_fittingthread_done(coreid);
				break;
			}
		}//endfor i runs

		 //Now the fitting ist done

		 //Now recenter the model such that com_xy = (0,0)
		saxs::recenter_fittingobject(localfittingobject_sp);
		localfittingobject_sp->m_compactness = saxs::return_compactness_of_fittingobject(localfittingobject_sp);

		//Copy data from local temporary vector to results
		for (unsigned int k = 0; k < reftofittingobject->m_model.size(); k++)
		{
			result_tracer[coreid]->m_coordinate[k]->m_x = localfittingobject_sp->m_model[k]->m_x;
			result_tracer[coreid]->m_coordinate[k]->m_y = localfittingobject_sp->m_model[k]->m_y;
			result_tracer[coreid]->m_coordinate[k]->m_z = localfittingobject_sp->m_model[k]->m_z;
			result_tracer[coreid]->m_coordinate[k]->m_nr_contacts = nr_contacts_vector[k];
		}
		result_tracer[coreid]->m_chi = localfittingobject_sp->m_chi;
		result_tracer[coreid]->m_mean_connectivity = localfittingobject_sp->m_mean_connectivity;
		result_tracer[coreid]->m_mean_nr_contacts = localfittingobject_sp->m_mean_nr_contacts;
		result_tracer[coreid]->m_compactness = localfittingobject_sp->m_compactness;
		result_tracer[coreid]->m_model_I = localfittingobject_sp->m_model_I;
		result_tracer[coreid]->m_fitted_I = localfittingobject_sp->m_fitted_I;

		//Send finish signal
		emit sig_fittingthread_done(coreid);
	}

#endif /*!THREAD_UTILS_H*/