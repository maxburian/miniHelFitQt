#ifndef MERGER_UTILS_H
#define MERGER_UTILS_H
//Includes
//STD
#include <vector>
//BOOST
#include <boost\shared_ptr.hpp>
#include <boost/thread.hpp>
//Self-Written
#include "data_def.h"
#include "data_utils.h"

namespace saxs
{

	/***************************************************************************************************************/
	// Merger dedicated helper functions
	/***************************************************************************************************************/

	//returns the mean distance between all points
	double get_mean_distance(const std::vector<saxs::coordinate_sp>& data)
	{
		int counter = 0;
		double r;
		double mean_dist = 0;

		for (int i = 0; i < data.size() - 1; i++) 
		{
			for (int j = 1; j < data.size(); j++)
			{
				r = std::sqrt(
					std::pow(data[i]->m_x - data[j]->m_x, 2) + 
					std::pow(data[i]->m_y - data[j]->m_y, 2) +
					std::pow(data[i]->m_z - data[j]->m_z, 2));
				mean_dist += r;
				counter++;
			}
		}
		return mean_dist / double(counter);
	}

	//returns NSD of two models, one specified by coordiante and one given by vectors
	double get_NSDofmodels(std::vector<saxs::coordinate_sp>& basemodel, double basemodel_fineness,
		std::vector<double>& model2_x, std::vector<double>& model2_y, std::vector<double>& model2_z, double model2_fineness)
	{
		//Helper vars
		double temp_min;
		double temp_dist;
		double nsd_1=0;
		double nsd_2=0;
		int size_1 = basemodel.size();
		int size_2 = model2_x.size();

		//calculate NSD1
		for (int i = 0; i < size_1; i++)
		{
			temp_min = 999999999;
			//Find minimum distance between i in basemodel and any point in model2
			for (int j = 0; j < size_2; j++)
			{
				temp_dist = std::sqrt(
					std::pow(model2_x[j] - basemodel[i]->m_x, 2) + 
					std::pow(model2_y[j] - basemodel[i]->m_y, 2) +
					std::pow(model2_z[j] - basemodel[i]->m_z, 2));
				if (temp_dist < temp_min) temp_min = temp_dist;
			}
			nsd_1 += temp_min*temp_min;
		}

		//calculate NSD2
		for (int i = 0; i < size_2; i++)
		{
			temp_min = 999999999;
			//Find minimum distance between i in basemodel and any point in model2
			for (int j = 0; j < size_1; j++)
			{
				temp_dist = std::sqrt(
					std::pow(basemodel[j]->m_x - model2_x[i], 2) +
					std::pow(basemodel[j]->m_y - model2_y[i], 2) +
					std::pow(basemodel[j]->m_z - model2_z[i], 2));
				if (temp_dist < temp_min) temp_min = temp_dist;
			}
			nsd_2 += temp_min*temp_min;
		}

		//Return full NSD
		return std::sqrt(0.5*(nsd_1/(double(size_1)*model2_fineness*model2_fineness)+
			nsd_2 / (double(size_2)*basemodel_fineness*basemodel_fineness)));
	}

	//recenter model to xyz=(0,0,0)
	void recenter_model(std::vector<saxs::coordinate_sp>& data)
	{
		int nr_da_inside = 0;
		double x_com = 0;
		double y_com = 0;
		double z_com = 0;
		for (int i = 0; i < data.size(); i++)
		{
			x_com += data[i]->m_x;
			y_com += data[i]->m_y;
			z_com += data[i]->m_z;
		}

		x_com = x_com / double(data.size());
		y_com = y_com / double(data.size());
		z_com = z_com / double(data.size());

		for (int i = 0; i < data.size(); i++)
		{
			data[i]->m_x -= x_com;
			data[i]->m_y -= y_com;
			data[i]->m_z -= z_com;
		}
	}

	//Calculates mean NSD of all permutations and the deviation
	std::pair<double, double> get_mean_NSD(std::vector<saxs::mergeobject_sp>& globalemergeobject)
	{
		double ave=0;
		double dev=0;
		int n = globalemergeobject.size();

		for (int i = 0; i < n; i++)
		{
			for (int j = 0; j < n; j++)
			{
				ave += globalemergeobject[i]->m_NSD[j];
			}
		}
		ave = ave / double(n*n-n);

		for (int i = 0; i < n; i++)
		{
			for (int j = 0; j < n; j++)
			{
				if (j == i)continue;
				dev += std::pow(ave-globalemergeobject[i]->m_NSD[j],2);
			}
		}
		dev = std::sqrt(dev/(double(n*n-n)-1));

		std::pair<double, double> temp;
		temp.first = ave;
		temp.second = dev;

		for (int i = 0; i < n; i++)
		{	
			ave = 0;
			for (int j = 0; j < n; j++)
			{
				ave += globalemergeobject[i]->m_NSD[j];
			}
			globalemergeobject[i]->m_mean_NSD = ave / (double(n - 1));
		}

		for (int i = 0; i < n; i++)
		{
			dev = 0;
			for (int j = 0; j < n; j++)
			{
				if (j == i)continue;
				dev += std::pow(globalemergeobject[i]->m_mean_NSD - globalemergeobject[i]->m_NSD[j], 2);
			}
			globalemergeobject[i]->m_sdev_NSD = std::sqrt(dev / (double(n-2)));
		}

		return temp;
	}

	//Merger function
	void mergemodels(std::vector<saxs::mergeobject_sp>& globalemergeobject, 
			std::vector<saxs::coordinate_sp>& sum_model, const int reference, std::vector<int>& includevector)
	{
		//Make temporary objects for local copy
		std::vector<double> x, y, z;
		for (int k = 0; k < includevector.size(); k++)
		{
			if (k == reference)continue;
			for (int i = 0; i < globalemergeobject[includevector[k]]->m_model.size(); i++)
			{
				x.push_back(globalemergeobject[includevector[k]]->m_model[i]->m_x);
				y.push_back(globalemergeobject[includevector[k]]->m_model[i]->m_y);
				z.push_back(globalemergeobject[includevector[k]]->m_model[i]->m_z);
			}
			convertCylinderToStdVectors(globalemergeobject[includevector[k]]->m_model_cyl, x, y, globalemergeobject[reference]->m_NSD_phi[k]);
			if (globalemergeobject[reference]->m_inverted[includevector[k]])for (int j = 0; j < x.size(); j++) x[j] = -1.0*x[j];

			for (int i = 0; i < globalemergeobject[includevector[k]]->m_model.size(); i++)
			{
				sum_model.push_back(coordinate_sp(new coordinate(x[i],y[i],z[i],0)));
			}
			x.clear();
			y.clear();
			z.clear();
		}

	}


	/***************************************************************************************************************/
	// Alignment thread functions
	/***************************************************************************************************************/
	//Main trhead function 

	void returnMinNSD(boost::mutex& io_mutex, std::vector<saxs::mergeobject_sp>& globalemergeobject, const int base_id, const int adjust_id)
	{
		//Make temporary objects for local copy
		std::vector<double> x, y, z;
		for (int i = 0; i < globalemergeobject[adjust_id]->m_model.size(); i++)
		{
			x.push_back(globalemergeobject[adjust_id]->m_model[i]->m_x);
			y.push_back(globalemergeobject[adjust_id]->m_model[i]->m_y);
			z.push_back(globalemergeobject[adjust_id]->m_model[i]->m_z);
		}
		//Now reinitiliaze using cylindrical coords
		convertCylinderToStdVectors(globalemergeobject[adjust_id]->m_model_cyl, x, y, 0);

		//Init minimization variables
		int num_steps = 360;
		double phi_min;
		double phi_current;
		double phi_step = saxs::pi * 2.0 / double(num_steps);
		double NSD_min;
		double NSD_current;
		bool inverted_min = false;

		//Get starting parameters
		NSD_current = get_NSDofmodels(globalemergeobject[base_id]->m_model, globalemergeobject[base_id]->m_mean_dist,
			x, y, z, globalemergeobject[adjust_id]->m_mean_dist);
		NSD_min = NSD_current;
		phi_min = 0;

		//Find alignment minimum
		for (int i = 1; i < num_steps; i++)
		{
			phi_current = phi_step*double(i);
			convertCylinderToStdVectors(globalemergeobject[adjust_id]->m_model_cyl, x, y, phi_current);
			NSD_current = get_NSDofmodels(globalemergeobject[base_id]->m_model, globalemergeobject[base_id]->m_mean_dist,
				x,y,z, globalemergeobject[adjust_id]->m_mean_dist);
			if (NSD_current < NSD_min)
			{
				NSD_min = NSD_current;
				phi_min = phi_current;
			}
		}

		//Now see if it is better to mirror the structure
		for (int i = 0; i < num_steps; i++)
		{
			phi_current = phi_step*double(i);
			convertCylinderToStdVectors(globalemergeobject[adjust_id]->m_model_cyl, x, y, phi_current);
			for (int k = 0; k < x.size(); k++) x[k] = -1.0*x[k];
			NSD_current = get_NSDofmodels(globalemergeobject[base_id]->m_model, globalemergeobject[base_id]->m_mean_dist,
				x, y, z, globalemergeobject[adjust_id]->m_mean_dist);
			if (NSD_current < NSD_min)
			{
				NSD_min = NSD_current;
				phi_min = phi_current;
				inverted_min = true;
			}
		}

		//Lock Mutex
		boost::mutex::scoped_lock lock(io_mutex);

		//Write everything into global object
		globalemergeobject[base_id]->m_NSD[adjust_id]= NSD_min;
		globalemergeobject[adjust_id]->m_NSD[base_id] = NSD_min;
		globalemergeobject[base_id]->m_NSD_phi[adjust_id] = phi_min;
		globalemergeobject[adjust_id]->m_NSD_phi[base_id] = -phi_min;
		globalemergeobject[base_id]->m_inverted[adjust_id] = inverted_min;
		globalemergeobject[adjust_id]->m_inverted[base_id] = inverted_min;

		//Unlock mutex
		lock.unlock();

	}


	//Main node function
	void alignLoadedModels(std::vector<saxs::mergeobject_sp>& globalemergeobject)
	{
		//Thread container
		boost::scoped_ptr<boost::thread_group> threadGroup_sp(new boost::thread_group);

		//Mutex to protect IO operations
		boost::mutex io_mutex;

		//Instantiate and start threads
		for (int i = 0; i<globalemergeobject.size()-1; i++)
		{
			for (int j = i+1; j < globalemergeobject.size(); j++)
			{
				threadGroup_sp->add_thread(new boost::thread(returnMinNSD, boost::ref(io_mutex), globalemergeobject, i,j));
			}
		}

		//Wait until the end of their jobs
		if (threadGroup_sp)
			threadGroup_sp->join_all();

	}


}//End of namespace saxs

#endif /*!MERGER_UTILS_H*/