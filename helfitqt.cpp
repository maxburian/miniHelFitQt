//Include database dependencies
//QT
#include <QtWidgets/QMainWindow>
#include <QMouseEvent>
#include <QLabel>
#include <QProgressBar>
#include <qapplication.h>
#include <qfileinfo.h>
#include <qmessagebox.h>
#include <QMouseEvent>
#include <qevent.h>
#include <qobject.h>
#include <QDoubleValidator>
#include <qstring.h>
//STD
#include <qvector.h>
#include <string.h>
#include <vector>
#include <algorithm>
//Boost
#include <boost/thread.hpp>
#include <boost/chrono.hpp>
#include <boost/scoped_ptr.hpp>
#include <boost/lexical_cast.hpp>
//Self-Written
#include "helfitqt.h"
#include "data_def.h"
#include "data_utils.h"
#include "thread_utils.h"
#include "merger_utils.h"


/***************************************************************************************************************/
// Initialization...
/***************************************************************************************************************/

HelFitQt::HelFitQt(QWidget *parent)
	: QMainWindow(parent)
{
	ui.setupUi(this);
}

HelFitQt::~HelFitQt()
{

}

void HelFitQt::init()
{
	initialize_graphs();
	//Initialize boxvalidator
	ui.lineStackSpacing->setValidator(new QDoubleValidator(0, 1000, 4, this));
	ui.lineCalcQmin->setValidator(new QDoubleValidator(0, 20, 4, this));
	ui.lineCalcQmax->setValidator(new QDoubleValidator(0, 20, 4, this));
	ui.lineStartTemp->setValidator(new QDoubleValidator(0.01, 20, 4, this));
	ui.lineDeltaTemp->setValidator(new QDoubleValidator(0, 2, 4, this));
	ui.lineDAsize->setValidator(new QDoubleValidator(0.001, 100, 4, this));
	//Initialize Sinc Lookup
	saxs::sinc_lookup[0] = 1;
	double sigma = 2.*saxs::pi / 0.000001;
	for (std::size_t i = 1; i < saxs::sinc_lookup.size(); ++i)
	{
		//saxs::sinc_lookup[i] = boost::math::sinc_pi((double)i / 1000);
		saxs::sinc_lookup[i] = boost::math::sinc_pi(double(i) / 1000.)*
			std::exp(-((double(i)/1000.)*(double(i) / 1000.)/(sigma*sigma)));
	}
	//Get number of cores for debyecalc
	core_number = boost::thread::hardware_concurrency();
	ui.spbNrCores->setMaximum(core_number);
	
	//Install event filters
	ui.spbDataMax->installEventFilter(this);
	ui.spbDataMin->installEventFilter(this);
	ui.spbFitMax->installEventFilter(this);
	ui.spbFitMin->installEventFilter(this);

	//Log output
	writetolog("Start of Helfitqt....");
	writetolog("------------------------");
	writetolog("Welcome!");
	writetolog("------------------------");
	writetolog(" ");

	//Status bar
	ui.statusBar->setMaximumHeight(20);
	ui.statusBar->setStyleSheet("font: 12px black;");
	statusLabel = new QLabel(this);
	statusProgressBar = new QProgressBar(this);
	statusLabel->setFixedWidth(700);
	statusProgressBar->setFixedWidth(300);
	statusProgressBar->setMinimum(0);
	statusProgressBar->setMaximum(100);
	statusLabel->setText("Welcome to HelFitQt! Waiting for commands...");
	statusProgressBar->setTextVisible(false);
	ui.statusBar->addPermanentWidget(statusLabel);
	ui.statusBar->addPermanentWidget(statusProgressBar, 0);

	//Initialize recommended Fittin parameters
	globalfittingobject->m_alphaConn = 0;
	globalfittingobject->m_tauConn = 2.0;
	globalfittingobject->m_betaComp = 0.01;
	globalfittingobject->m_sigmaComp = 0.5;
	globalfittingobject->m_gammaHelBias = 0.3;
	globalfittingobject->m_rand_seed_scalar = 0.2;
}


/***************************************************************************************************************/
// Object control - Clear/MoveTo/etc...
/***************************************************************************************************************/

//Clears scattering data
void HelFitQt::clear_data()
{
	imported_data.clear();
	x_expdata.clear();
	y_expdata.clear();
	x_fitdata.clear();
	y_fitdata.clear();
}

//Clears current model
void HelFitQt::clear_model()
{
	imported_model_data.clear();
	imported_model.clear();
	plot_model_cyl.clear();
	model_x.clear();
	model_y.clear();
	model_z.clear();
	phi_moodel_rotation = 0;
}

//Move global variables into fittingobject
void HelFitQt::movedatatofittingobject()
{

	//First move data to global object
	globalfittingobject->m_data_q.clear();
	globalfittingobject->m_data_I.clear();
	globalfittingobject->m_data_e.clear();
	globalfittingobject->m_model_I.clear();
	globalfittingobject->m_fitted_I.clear();
	if (data_loaded)
	{
		globalfittingobject->m_data_q.resize(x_fitdata.size());
		globalfittingobject->m_data_q = x_fitdata.toStdVector();
		globalfittingobject->m_data_I.resize(x_fitdata.size());
		globalfittingobject->m_data_I = y_fitdata.toStdVector();
		globalfittingobject->m_data_e.resize(x_fitdata.size());
		globalfittingobject->m_data_e = e_fitdata.toStdVector();
		globalfittingobject->m_model_I.resize(x_fitdata.size());
		globalfittingobject->m_fitted_I.resize(x_fitdata.size());
	}
	else
	{
		//Resize data intensity to 0 (used as boolean later on)
		globalfittingobject->m_data_I.resize(0);
		//Generate linsapce data
		double step = double(calc_qmax - calc_qmin) / double(num_calcpoints - 1);
		double temp = calc_qmin;
		for (int i = 0; i < num_calcpoints; i++)
		{
			globalfittingobject->m_data_q.push_back(temp);
			temp += step;
		}
		globalfittingobject->m_model_I.resize(num_calcpoints);
		globalfittingobject->m_fitted_I.resize(num_calcpoints);

	}
	std::fill(globalfittingobject->m_model_I.begin(), globalfittingobject->m_model_I.end(), 0);
	std::fill(globalfittingobject->m_fitted_I.begin(), globalfittingobject->m_fitted_I.end(), 0);

	//Now move params to global object
	globalfittingobject->m_num_stacks = int(ui.spbCalcNrStacks->value());
	globalfittingobject->m_stack_spacing = double(ui.lineStackSpacing->text().toDouble());
	globalfittingobject->m_num_cores = int(ui.spbNrCores->value());
	globalfittingobject->m_weighing = curvefit_weight;

}


/***************************************************************************************************************/
// Event Handling - Buttons / MenuItems
/***************************************************************************************************************/

//Button - Load Scattering data (auto signal-slot connection)
void HelFitQt::on_btnLoadFile_clicked()
{
	if (data_loaded) clear_data();
	scatterfilepath = saxs::filedialog("Select data file:", "Data (*.qI);;Chi-file (*.chi);;All (*.*)");
	if (scatterfilepath != "none")
	{
		//scatterfilepath = QFileDialog::getOpenFileName(this, tr("Select File"), "/path/to/file/", tr("Data (*.qI);;All (*.*)"),, QFileDialog::DontUseNativeDialog);
		QFileInfo fi(scatterfilepath);
		scatterfilename = fi.fileName();
		ui.txtFilePath->setText(scatterfilename);

		//Import data from file inte scatteringdata object
		if(saxs::import_scatteringdata(scatterfilepath, imported_data)==false) return;
		//check number of points
		saxs::reduce_scatteringdata_size(imported_data);
		//convert object data into globals x_expdata y_expdata
		saxs::convertScatteringToVector(imported_data, x_expdata, y_expdata, e_expdata);
		//rescale data so no values are <1
		saxs::rescaleScatteringData(imported_data, y_expdata, e_expdata);

		//determine datarange globlas
		min_datarange = 0;
		min_fitrange = 0;
		num_datapoints = x_expdata.size();
		max_datarange = x_expdata.size() - 1;
		max_fitrange = x_expdata.size() - 1;
		num_calcpoints = max_fitrange + 1 - min_fitrange;

		///Altering userinterface for plotting
		ui.spbCalcNrPoints->setValue(num_calcpoints);
		ui.spbDataMin->setEnabled(true);
		ui.spbDataMax->setEnabled(true);

		ui.spbDataMin->setMinimum(1);
		ui.spbDataMin->setMaximum(num_datapoints);
		ui.spbDataMin->setValue(1);
		ui.spbDataMax->setMinimum(1);
		ui.spbDataMax->setMaximum(num_datapoints);
		ui.spbDataMax->setValue(num_datapoints);

		ui.spbFitMin->setEnabled(true);
		ui.spbFitMax->setEnabled(true);

		ui.spbFitMin->setMinimum(1);
		ui.spbFitMin->setMaximum(num_datapoints);
		ui.spbFitMin->setValue(1);
		ui.spbFitMax->setMinimum(1);
		ui.spbFitMax->setMaximum(num_datapoints);
		ui.spbFitMax->setValue(num_datapoints);
		ui.lineQmin->setText(QString::number(x_expdata[min_datarange]));
		ui.lineQmax->setText(QString::number(x_expdata[max_datarange]));
		ui.chkbLogLogPlot->setEnabled(true);
		//Altering userinterface for calc
		ui.lineCalcQmin->setEnabled(false);
		ui.lineCalcQmax->setEnabled(false);
		ui.spbCalcNrPoints->setEnabled(false);
		ui.lineCalcQmax->setText(QString::number(x_expdata[max_datarange]));
		ui.lineCalcQmin->setText(QString::number(x_expdata[min_datarange]));
		ui.spbCalcNrPoints->setValue(max_datarange - min_datarange + 1);

		//set global dataload
		data_loaded = true;

		//plot data
		plot_data();
		ui.tabWidget->setCurrentIndex(0);

		//Write info in ouput field
		writetolog(" ");
		std::string str = " ";
		writetolog(("Succesfully loaded data from: " + scatterfilepath.toStdString()));
		str = str+ boost::lexical_cast<std::string>(num_datapoints);
		str = str+ " points\t\t";
		str = str + "Qmin: " + ui.lineQmin->text().toStdString() + "\t\tQmax: " + ui.lineQmax->text().toStdString();
		writetolog(str);

		str = "Succesfully loaded scattering data from: " + scatterfilepath.toStdString();
		statusLabel->setText(QString::fromStdString(str));
	}
	else
	{
		ui.txtFilePath->setText("none selected");
		statusLabel->setText("No file selected...");
	}
}

//Menu - Load Model from file
void HelFitQt::loadModelFromData()
{
	//Load data from file
	QString modelfilepath;
	modelfilepath = saxs::filedialog("Select model file:", "XYZ (*.xyz);;PDB (*.pdb);;All (*.*)");
	if (modelfilepath != "none")
	{
		if (model_loaded) clear_model();

		//Get file extension
		std::string file_ext = modelfilepath.toStdString().substr(modelfilepath.toStdString().find_last_of(".") + 1);

		if (file_ext =="xyz")	if (saxs::import_xyz_model(modelfilepath, imported_model) == false)return;
		if (file_ext == "pdb")	if (saxs::import_pdb_model(modelfilepath, imported_model) == false)return;
		//Copy original data into plotting qvectors
		saxs::convertCoordinateToVector(imported_model, model_x, model_y, model_z);
		std::vector<double> model_boundaries(3);
		model_boundaries = saxs::convertCoordinateToCylinder(imported_model, plot_model_cyl);
		plot_model(model_boundaries);
		phi_moodel_rotation = saxs::pi;
		ui.hsliderPhiRot->setValue(500);
		model_loaded = true;
		ui.menuDebyeCurrModel->setEnabled(true);
		ui.menuExpandCurrModel->setEnabled(true);
		ui.menuSaveModel->setEnabled(true);
		ui.tabWidget->setCurrentIndex(1);

		//Copying data into globalfittingobject
		globalfittingobject->m_model.clear();
		globalfittingobject->m_model.resize(imported_model.size());
		globalfittingobject->m_model = imported_model;
		globalfittingobject->m_contact_d_sq = -1;

		//Write info in ouput field
		writetolog(" ");
		std::string str = " ";
		writetolog(("Succesfully loaded model from: " + modelfilepath.toStdString()));
		str = str + "The model includes ";
		str = str + boost::lexical_cast<std::string>(imported_model.size());
		str = str + " DAs";
		writetolog(str);

		str = "Succesfully loaded model from: " + modelfilepath.toStdString();
		statusLabel->setText(QString::fromStdString(str));
	}
}

//Menu - Run merger
void HelFitQt::actionRunMerger()
{
	//Set active tab to output log
	ui.tabWidget->setCurrentIndex(2);

	//Opening FileDialog
	QString dialogTitle = "Select model file:";
	QString FileFilter = "PDB-files (*.pdb);;All (*.*)";
	QStringList fileNames = saxs::multifiledialog(dialogTitle,FileFilter);

	//Check if more than 1 model is selected
	if (fileNames.size() < 2)
	{
		saxs::error_message("Only one file selected! This does not make any sense to merge...");
		return;
	}

	//Select target filename
	QString defaultFilter("PDB-files (*.pdb)");
	/* Static method approach */
	QString filename = QFileDialog::getSaveFileName(0, "Save files to..", QDir::currentPath(),
		FileFilter, &defaultFilter, QFileDialog::DontUseNativeDialog);
	if (filename == NULL) return;
	if (filename.toStdString().rfind(".pdb")!=-1) savefilename = filename.toStdString().substr(0, filename.size() - 4);
	else savefilename = filename.toStdString();

	//Disabling ui
	disable_ui();
	statusLabel->setText("Merging datafiles...");

	//Create local mergerobject on stack
	std::vector<saxs::mergeobject_sp> loadedmergerobject;

	//Helper
	std::vector<double> model_boundaries(3);

	//Write data in mergeobject
	for (int i = 0; i < fileNames.size(); i++)
	{
		loadedmergerobject.push_back(saxs::mergeobject_sp(new saxs::mergeobject));

		//Load model from file into object
		if (saxs::import_pdb_model(fileNames[i], loadedmergerobject[i]->m_model)==false) break;

		//Recenter model to xy com
		saxs::recenter_model(loadedmergerobject[i]->m_model);

		//Convert coordinates into cylindrical ones
		model_boundaries = saxs::convertCoordinateToCylinder(loadedmergerobject[i]->m_model, loadedmergerobject[i]->m_model_cyl);

		//Calc mean distance
		loadedmergerobject[i]->m_mean_dist = saxs::get_mean_distance(loadedmergerobject[i]->m_model);
	}

	//Initialize result vector
	for (int i = 0; i < fileNames.size(); i++)
	{
		for (int j = 0; j < fileNames.size(); j++)
		{
			loadedmergerobject[i]->m_inverted.push_back(false);
			loadedmergerobject[i]->m_NSD.push_back(0);
			loadedmergerobject[i]->m_NSD_phi.push_back(0);
		}
	}


	//Write info in ouput field
	writetolog(" ");
	writetolog("******* Cylindrical Model Merger *******");
	std::string str;
	str = "Merging a total of " + boost::lexical_cast<std::string>(loadedmergerobject.size()) + " models.";
	writetolog(str);
	


	//--------------Model Alginment-----------------------------

	//Now find min NSD and align all models
	saxs::alignLoadedModels(loadedmergerobject);

	//get mean nsd and deviation
	std::pair<double, double> mean_NSD = saxs::get_mean_NSD(loadedmergerobject);

	//Write results to log
	std::ostringstream ss;
	writetolog(" ");
	writetolog("NSD results:");
	writetolog(str);
	ss.str("");
	ss.clear();
	ss << "Average NSD = " << mean_NSD.first << " +- " << mean_NSD.second;
	writetolog(ss.str());

	//Find best reference
	double min = loadedmergerobject[0]->m_mean_NSD;
	int min_pntr = 0;
	for (int j = 1; j < loadedmergerobject.size(); j++)
	{
		if (loadedmergerobject[j]->m_mean_NSD < min)
		{
			min = loadedmergerobject[j]->m_mean_NSD;
			min_pntr = j;
		}
	}
	ss.str("");
	ss.clear();
	ss << "The reference model is model " << std::setprecision(0) << min_pntr;
	ss << " with mean NSD " << std::setprecision(3) << min;
	writetolog(ss.str());

	ss.str("");
	ss.clear();
	ss << "Thus, models with average NSD > " << min << " +- " << loadedmergerobject[min_pntr]->m_sdev_NSD << " will be dropped";
	writetolog(ss.str());

	//Include vector - saves which models should be merged
	std::vector<int> includevector;

	//Output
	writetolog("model \tmean NSD");
	int k = 0;
	bool dropitem = false;
	while (k < loadedmergerobject.size())
	{
		if (loadedmergerobject[min_pntr]->m_NSD[k]>(min + loadedmergerobject[min_pntr]->m_sdev_NSD)) dropitem = true;
		ss.str("");
		ss.clear();
		ss << std::fixed << std::setprecision(0);
		ss << k;
		ss << std::fixed << std::setprecision(3);
		ss << "\t" << loadedmergerobject[min_pntr]->m_NSD[k];
		if (dropitem)
		{
			ss << "\tdropped!";
			k++;
			dropitem = false;
		}
		else
		{
			includevector.push_back(k);
			ss << "\tincluded!";
			k++;
		}
		writetolog(ss.str());
	}




	//--------------Model Merging ----------------------------
	
	//Merge all models into one
	std::vector<saxs::coordinate_sp> sum_model;
	saxs::mergemodels(loadedmergerobject, sum_model, min_pntr, includevector);

	//Now save summed model to file
	saxs::writemodeltopdb(savefilename+"_sum", sum_model);

	//Output 
	str = "Full model saved to " + savefilename + "_sum.pdb";
	writetolog(str);




	//--------------Occupancy map calculation-------------------
	
	//Get dimensions of new grid
	double rmax = saxs::get_rmax_from_coordinates(sum_model)+1.;
	std::pair<double,double> z_boundaries = saxs::get_zboundaries_from_coordinates(sum_model);
	double height = z_boundaries.second - z_boundaries.first +2.;
	double V_new = (2*rmax)*(2*rmax) *height;

	//Get dimensions of old grid
	double cyl_diam = saxs::get_cyldiameter_from_coordinates(sum_model);
	double V_old = cyl_diam*cyl_diam / 4.*saxs::pi * height;

	//Number of voxels in Grid
	int n_newgrid = int(double(sum_model.size()) * V_new / V_old*2.);
	double calc_resolution = std::pow(V_new/double(n_newgrid),0.333);

	//Check resolution with user
	str = "Specify resolution of occupancy model [nm]:";
	double resolution = saxs::input_message_doub(QString::fromStdString(str), calc_resolution, 0, 100);
	writetolog("");
	ss.str("");
	ss.clear();
	ss << "Resolution of the occupancy map: " << std::setprecision(2) << resolution << " nm";
	writetolog(ss.str());

	//Generate artifical grid
	std::vector<saxs::coordinate_sp> occ_model;
	saxs::generate_art_grid(occ_model, rmax, z_boundaries, resolution);

	//Impose model and grid
	saxs::impose_model_on_grid(sum_model, occ_model);

	//Save occupancy model
	saxs::writeoccmodeltopdb(savefilename + "_occ", occ_model);
	str = "Occupancy model saved to " + savefilename + "_occ.pdb";
	writetolog(str);

	//Final Message
	writetolog(" ");
	writetolog("Model-merge done!");

	//Enabling Interface
	enable_ui();
	statusLabel->setText("Merging done!");

}

//Menu - Change linreg weight
void HelFitQt::actionWeight0()
{
	curvefit_weight = 0;
	ui.menuWeight0->setChecked(true);
	ui.menuWeight1->setChecked(false);
	ui.menuWeight2->setChecked(false);
}

//Menu - Change linreg weight
void HelFitQt::actionWeight1()
{
	curvefit_weight = 1;
	ui.menuWeight0->setChecked(false);
	ui.menuWeight1->setChecked(true);
	ui.menuWeight2->setChecked(false);
}

//Menu - Change linreg weight
void HelFitQt::actionWeight2()
{
	curvefit_weight = 2;
	ui.menuWeight0->setChecked(false);
	ui.menuWeight1->setChecked(false);
	ui.menuWeight2->setChecked(true);
}

//Menu - Calculate scatteringcurve of current model
void HelFitQt::calcDebyeStackCurrent()
{
	//Check if model is loaded
	if (!model_loaded && !model_generated)
	{
		saxs::error_message("No model loaded!");
		return;
	}

	//Make sure current ui values are in globals
	num_calcpoints = int(ui.spbCalcNrPoints->value());
	calc_qmin = double(ui.lineCalcQmin->text().toDouble());
	calc_qmax = double(ui.lineCalcQmax->text().toDouble());

	//get calculation parameters from gui
	//build fittingobject
	movedatatofittingobject();

	writetolog("Number of cores\t = \t" + boost::lexical_cast<std::string>(globalfittingobject->m_num_cores));

	//Write info in ouput field
	writetolog(" ");
	writetolog("******* Debye Stack Calculation *******");
	std::string str = " ";
	str = "Calculating model data between q = " + ui.lineCalcQmin->text().toStdString() + " - " + ui.lineCalcQmax->text().toStdString();
	str = str +  " using " + boost::lexical_cast<std::string>(num_calcpoints) +" points.";
	writetolog(str);
	writetolog("Parameters:");
	writetolog("Stacking distance\t = \t" + boost::lexical_cast<std::string>(globalfittingobject->m_stack_spacing)+" nm");
	writetolog("Number of stacks\t = \t" + boost::lexical_cast<std::string>(globalfittingobject->m_num_stacks));
	writetolog("Number of cores\t = \t" + boost::lexical_cast<std::string>(globalfittingobject->m_num_cores));

	writetolog("START!");
	qApp->processEvents();
	//Get threads start time
	auto start = boost::chrono::system_clock::now();
		
	globalfittingobject->m_diameter = saxs::get_cyldiameter_from_coordinates(globalfittingobject->m_model);
	globalfittingobject->m_contact_d_sq = std::pow(2 * saxs::get_critradius_from_coordinates(globalfittingobject->m_diameter,
		globalfittingobject->m_stack_spacing, globalfittingobject->m_model.size()), 2);
	saxs::calc_contacts_of_model(globalfittingobject);
	globalfittingobject->m_mean_nr_contacts = saxs::return_mean_contacts_of_model(globalfittingobject->m_model);
	globalfittingobject->m_mean_connectivity = saxs::return_mean_connectivity_of_fittingobject(globalfittingobject);
	globalfittingobject->m_compactness = saxs::return_compactness_of_fittingobject(globalfittingobject);

	globalfittingobject->m_sigma_ff = 2.*saxs::pi/double(ui.lineDAsize->text().toDouble());
	saxs::calc_ff(globalfittingobject);

	saxs::calcDebyeStackCurrentModel(globalfittingobject);

	//Get threads end time
	auto end = boost::chrono::system_clock::now();

	//Calculate elapsed time
	boost::uint64_t elapsed_seconds =boost::chrono::duration_cast<boost::chrono::seconds>(end - start).count();
	writetolog("FINISHED!");
	writetolog("Duration\t\t = \t " + boost::lexical_cast<std::string>(elapsed_seconds) + " seconds.");
	writetolog("Results:");
	writetolog("Estimated diam.\t = \t" + boost::lexical_cast<std::string>(globalfittingobject->m_diameter) + " nm");
	writetolog("Critical NN dist.\t = \t" + boost::lexical_cast<std::string>(std::sqrt(globalfittingobject->m_contact_d_sq)) + " nm");
	writetolog("Ave. number of contacts\t = \t" + boost::lexical_cast<std::string>(globalfittingobject->m_mean_nr_contacts));
	writetolog("Ave. coordination\t = \t" + boost::lexical_cast<std::string>(globalfittingobject->m_mean_nr_coordination));
	writetolog("Global connectivity\t = \t" + boost::lexical_cast<std::string>(globalfittingobject->m_mean_connectivity));
	writetolog("Global compactness\t = \t" + boost::lexical_cast<std::string>(globalfittingobject->m_compactness));

	//Clear previous calcdata
	//Old Datastructure!!!!
	imported_model_data.clear();
	x_imported_model_data.clear();
	y_imported_model_data.clear();

	//Move data to plotable qvector and normalize to I[0]
	for (int i = 0; i < globalfittingobject->m_data_q.size(); i++)
	{
		x_imported_model_data.push_back(globalfittingobject->m_data_q[i]);
		y_imported_model_data.push_back(globalfittingobject->m_fitted_I[i]);
		imported_model_data.push_back(saxs::scatteringdata_sp(new saxs::scatteringdata(
			globalfittingobject->m_data_q[i], globalfittingobject->m_fitted_I[i])));
	}

	plotCalcModel();
	
	//UI chagnes
	statusLabel->setText("Calculation done! Please check results...");
	ui.menuSaveData->setEnabled(true);

}

//Menu - Generate custom model
void HelFitQt::actionGenerateRandModel()
{
	dialog_modelvars *mydialog_modelvars;
	mydialog_modelvars = new dialog_modelvars;
	mydialog_modelvars->exec();
	bool dialogOk = mydialog_modelvars->result();
	
	//Check if parameters are accepted
	if (!dialogOk)
	{
		saxs::error_message("No parameters specified..."); 
		//Destroy object
		delete mydialog_modelvars;
		return;
	}
	//Transfer values into HelFitQt class
	int nr_atoms = mydialog_modelvars->m_nr_atoms;
	double height = mydialog_modelvars->m_stackspacing;
	ui.lineStackSpacing->setText(QString::number(height));
	double diameter = mydialog_modelvars->m_diameter;
	
	//Reinit structures
	generated_model.clear();
	plot_model_cyl.clear();
	model_x.clear();
	model_y.clear();
	model_z.clear();

	//Generate ranodm cylindrical coords:
	saxs::gen_cyl_randomseed(plot_model_cyl, height, diameter, nr_atoms);
	//initialize plotting vectors and load z coords
	saxs::convertCoordinateToVector(plot_model_cyl, model_x, model_y, model_z);
	//Now convert variables into cartesion vectors
	saxs::convertCylinderToQVectors(plot_model_cyl, model_x, model_y, 0);
	//Now move everything into generated_model object
	saxs::convertQvectorsToCoordinate(generated_model, model_x, model_y, model_z);

	
	//Plotting
	std::vector<double> model_boundaries(3);
	model_boundaries = saxs::convertCoordinateToCylinder(generated_model, plot_model_cyl);
	plot_model(model_boundaries);
	phi_moodel_rotation = saxs::pi;
	ui.hsliderPhiRot->setValue(500);
	model_generated = true;
	ui.menuDebyeCurrModel->setEnabled(true);
	ui.menuExpandCurrModel->setEnabled(true);
	ui.tabWidget->setCurrentIndex(1);
	ui.menuSaveModel->setEnabled(true);

	//Copying data into globalfittingobject
	globalfittingobject->m_model.clear();
	globalfittingobject->m_model.resize(generated_model.size());
	globalfittingobject->m_model = generated_model;
	globalfittingobject->m_contact_d_sq = -1;

	//Write info in ouput field
	writetolog(" ");
	std::string str = " ";
	writetolog("Succesfully generated random model!");
	str = str + "The model includes ";
	str = str + boost::lexical_cast<std::string>(generated_model.size());
	str = str + " DAs";
	writetolog(str);

}

//Menu - Expands the current model
void HelFitQt::actionExpandCurrentModel()
{
	int new_points = saxs::input_message("...number of points to add:", 0, 0, 2000);
	if (new_points == 0)
	{
		saxs::error_message("No points added...");
		return;
	}
	if (new_points + globalfittingobject->m_model.size() > 4000)
	{
		saxs::error_message("Your model is becoming to large. Reconsider...");
		return;
	}

	//Make sure diameter is known
	globalfittingobject->m_diameter = saxs::get_cyldiameter_from_coordinates(globalfittingobject->m_model);
	//Now expand model
	saxs::expand_fittingobject(globalfittingobject, new_points);
	//Convert new model to cylindrical coords
	std::vector<double> model_boundaries(3);
	model_boundaries = saxs::convertCoordinateToCylinder(globalfittingobject->m_model, plot_model_cyl);
	//initialize plotting vectors and load z coords
	saxs::convertCoordinateToVector(plot_model_cyl, model_x, model_y, model_z);
	//convert variables into cartesion vectors
	saxs::convertCylinderToQVectors(plot_model_cyl, model_x, model_y, 0);

	//Plotting
	plot_model(model_boundaries);
	phi_moodel_rotation = saxs::pi;
	ui.hsliderPhiRot->setValue(500);
	ui.tabWidget->setCurrentIndex(1);

	//Write info in ouput field
	writetolog(" ");
	std::string str = " ";
	writetolog("Succesfully expanded model!");
	str = str + "The model now includes ";
	str = str + boost::lexical_cast<std::string>(globalfittingobject->m_model.size());
	str = str + " DAs";
	writetolog(str);
}

//Menu - Change fitting parameters
void HelFitQt::actionChangeFittingParams()
{
	dialog_fittingvars *mydialog_fittingvars;
	mydialog_fittingvars = new dialog_fittingvars;

	//Load current variable values into UI
	mydialog_fittingvars->alpha_conn = globalfittingobject->m_alphaConn;
	mydialog_fittingvars->tau_conn = globalfittingobject->m_tauConn;
	mydialog_fittingvars->beta_comp = globalfittingobject->m_betaComp;
	mydialog_fittingvars->sigma_comp = globalfittingobject->m_sigmaComp;
	mydialog_fittingvars->gamma_helbias = globalfittingobject->m_gammaHelBias;
	mydialog_fittingvars->random_seed_scalar = globalfittingobject->m_rand_seed_scalar;
	mydialog_fittingvars->write_current_vars();

	//Run Dialog
	mydialog_fittingvars->exec();
	bool dialogOk = mydialog_fittingvars->result();

	//Check if parameters are accepted
	if (!dialogOk)
	{
		saxs::error_message("No parameters specified...");
		//Destroy object
		delete mydialog_fittingvars;
		return;
	}

	//Transfer values into HelFitQt class
	globalfittingobject->m_alphaConn = mydialog_fittingvars->alpha_conn;
	globalfittingobject->m_tauConn = mydialog_fittingvars->tau_conn;
	globalfittingobject->m_betaComp = mydialog_fittingvars->beta_comp;
	globalfittingobject->m_sigmaComp = mydialog_fittingvars->sigma_comp;
	globalfittingobject->m_gammaHelBias = mydialog_fittingvars->gamma_helbias;
	globalfittingobject->m_rand_seed_scalar = mydialog_fittingvars->random_seed_scalar;

	//Write info in ouput field
	writetolog(" ");
	std::string str = " ";
	writetolog("New fitting parameters specified: ");
	str = "Conn. weight alpha\t = \t";
	writetolog(str + boost::lexical_cast<std::string>(globalfittingobject->m_alphaConn));
	str = "Conn. potential tau\t = \t";
	writetolog(str + boost::lexical_cast<std::string>(globalfittingobject->m_tauConn));
	str = "Comp. weight alpha\t = \t";
	writetolog(str + boost::lexical_cast<std::string>(globalfittingobject->m_betaComp));
	str = "Comp. potential sig.\t = \t";
	writetolog(str + boost::lexical_cast<std::string>(globalfittingobject->m_sigmaComp));
	str = "Helical bias weight \t = \t";
	writetolog(str + boost::lexical_cast<std::string>(globalfittingobject->m_gammaHelBias));
	str = "Random seed scalar \t = \t";
	writetolog(str + boost::lexical_cast<std::string>(globalfittingobject->m_rand_seed_scalar));
}

//Menu - Saves the currently displayed model to pdb files
void HelFitQt::actionSaveModel()
{
	//Select target filename
	QString filters("PDB-files (*.pdb);;All files (*.*)");
	QString defaultFilter("PDB-files (*.pdb)");
	/* Static method approach */
	QString filename = QFileDialog::getSaveFileName(0, "Save files to..", QDir::currentPath(),
		filters, &defaultFilter);
	if (filename == NULL) return;
	savefilename = filename.toStdString().substr(0, filename.size() - 4);

	//Put data into fittingobject
	saxs::convertQvectorsToCoordinate(globalfittingobject->m_model, model_x, model_y, model_z);
	globalfittingobject->m_num_stacks = ui.spbCalcNrStacks->value();

	//Now write to file
	saxs::writemodeltopdb(savefilename, globalfittingobject->m_model);
	saxs::writestackedmodeltopdb(savefilename, globalfittingobject->m_model, globalfittingobject);

	writetolog("Models saved to " + savefilename + ".pdb and *_stck.pdb");
}

//Menu - Saves the currently display data to chi
void HelFitQt::actionSaveData()
{
	//Select target filename
	QString filters("Chi-files (*.chi);;All files (*.*)");
	QString defaultFilter("Chi-files (*.chi)");
	/* Static method approach */
	QString filename = QFileDialog::getSaveFileName(0, "Save files to..", QDir::currentPath(),
		filters, &defaultFilter);
	if (filename == NULL) return;
	savefilename = filename.toStdString().substr(0, filename.size() - 4);

	//Now write to file
	if (data_loaded) 
		saxs::writedatatochi(savefilename, x_fitdata.toStdVector(), y_fitdata.toStdVector(), y_imported_model_data.toStdVector());
	else
	{
		std::vector<double> helper(x_imported_model_data.size(), 0);
		saxs::writedatatochi(savefilename, x_imported_model_data.toStdVector(), helper, y_imported_model_data.toStdVector());
	}

	writetolog("Data saved to " + savefilename + ".chi");
}

//Button - Start Fitting Procedure (auto signal-slot connection)
void HelFitQt::on_btnFit_clicked()
{
	//Check if model is loaded
	if (!model_loaded && !model_generated)
	{
		saxs::error_message("No model loaded!");
		return;
	}

	if (!data_loaded)
	{
		saxs::error_message("No data loaded!");
		return;
	}

	//Select target filename
	QString filters("Chi-files (*.chi);;PDB-files (*.pdb);;All files (*.*)");
	QString defaultFilter("Chi-files (*.chi)");
	/* Static method approach */
	QString filename = QFileDialog::getSaveFileName(0, "Save files to..", QDir::currentPath(),
		filters, &defaultFilter, QFileDialog::DontUseNativeDialog);
	if (filename == NULL) return;
	savefilename = filename.toStdString().substr(0, filename.size() - 4);
	if (filename.toStdString().rfind(".pdb") != -1 || filename.toStdString().rfind(".chi")!=-1)
									savefilename = filename.toStdString().substr(0, filename.size() - 4);
	else savefilename = filename.toStdString();
	if (savefilename.rfind("_") >(savefilename.size()-4)) savefilename = savefilename.substr(0, savefilename.rfind("_"));

	//Make sure current ui values are in globals
	num_calcpoints = int(ui.spbCalcNrPoints->value());
	calc_qmin = double(ui.lineCalcQmin->text().toDouble());
	calc_qmax = double(ui.lineCalcQmax->text().toDouble());
	globalfittingobject->m_start_temp = double(ui.lineStartTemp->text().toDouble());
	globalfittingobject->m_delta_temp = double(ui.lineDeltaTemp->text().toDouble());
	globalfittingobject->m_num_runs = int(ui.spbNrRuns->value());
	globalfittingobject->m_multicore = ui.menuMultiCore->isChecked();


	//build fittingobject
	movedatatofittingobject();

	//first recenter model if only bb is fitted
	saxs::recenter_fittingobject(globalfittingobject);
	
	globalfittingobject->m_diameter = saxs::get_cyldiameter_from_coordinates(globalfittingobject->m_model);
	globalfittingobject->m_contact_d_sq = std::pow(2 * saxs::get_critradius_from_coordinates(globalfittingobject->m_diameter, 
					globalfittingobject->m_stack_spacing, globalfittingobject->m_model.size()), 2);
	saxs::calc_contacts_of_model(globalfittingobject);
	globalfittingobject->m_mean_nr_contacts = saxs::return_mean_contacts_of_model(globalfittingobject->m_model);
	globalfittingobject->m_mean_connectivity = saxs::return_mean_connectivity_of_fittingobject(globalfittingobject);
	globalfittingobject->m_compactness = saxs::return_compactness_of_fittingobject(globalfittingobject);

	globalfittingobject->m_sigma_ff = 2.*saxs::pi / double(ui.lineDAsize->text().toDouble());
	saxs::calc_ff(globalfittingobject);

	saxs::calcDebyeStackCurrentModel(globalfittingobject);

	//Plot to Ui
	//Clear previous calcdata
	imported_model_data.clear();
	x_imported_model_data.clear();
	y_imported_model_data.clear();

	//Move data to plotable qvector and normalize to I[0]
	for (int i = 0; i < globalfittingobject->m_data_q.size(); i++)
	{
		x_imported_model_data.push_back(globalfittingobject->m_data_q[i]);
		y_imported_model_data.push_back(globalfittingobject->m_fitted_I[i]);
		imported_model_data.push_back(saxs::scatteringdata_sp(new saxs::scatteringdata(
			globalfittingobject->m_data_q[i], globalfittingobject->m_fitted_I[i])));
	}
	plotCalcModel();

	writetolog("Number of cores\t = \t" + boost::lexical_cast<std::string>(globalfittingobject->m_num_cores));

	//Write info in ouput field
	writetolog(" ");
	writetolog("******* Debye Fit *******");
	std::string str = " ";
	str = "Calculating model data between q = " + ui.lineCalcQmin->text().toStdString() + " - " + ui.lineCalcQmax->text().toStdString();
	str = str + " using " + boost::lexical_cast<std::string>(num_calcpoints) + " points.";
	writetolog(str);
	writetolog("Parameters:");
	writetolog("Stacking distance\t = \t" + boost::lexical_cast<std::string>(globalfittingobject->m_stack_spacing) + " nm");
	writetolog("Number of stacks\t = \t" + boost::lexical_cast<std::string>(globalfittingobject->m_num_stacks));
	writetolog("Number of cores\t = \t" + boost::lexical_cast<std::string>(globalfittingobject->m_num_cores));
	writetolog("Model attributes before fit:");
	writetolog("Estimated diam.\t = \t" + boost::lexical_cast<std::string>(globalfittingobject->m_diameter) + " nm");
	writetolog("Critical NN dist.\t = \t" + boost::lexical_cast<std::string>(std::sqrt(globalfittingobject->m_contact_d_sq)) + " nm");
	writetolog("Ave. number of contacts\t = \t" + boost::lexical_cast<std::string>(globalfittingobject->m_mean_nr_contacts));
	writetolog("Ave. coordination\t = \t" + boost::lexical_cast<std::string>(globalfittingobject->m_mean_nr_coordination));
	writetolog("Global connectivity\t = \t" + boost::lexical_cast<std::string>(globalfittingobject->m_mean_connectivity));
	writetolog("Global compactness\t = \t" + boost::lexical_cast<std::string>(globalfittingobject->m_compactness));

	switch (globalfittingobject->m_weighing) {
	case (0) : str = "Chisquare (I(q)*q)\t = \t"; break;
	case (1) : str = "Chisquare (I(q)*q^2)\t = \t"; break;
	case (2) : str = "Chisquare (-)\t = \t"; break;
	}
	writetologext("Initial " + str + boost::lexical_cast<std::string>(globalfittingobject->m_chi));

	writetolog("START!");
	qApp->processEvents();
	//Get threads start time
	start = boost::chrono::system_clock::now();

	//Disabling interface
	disable_ui();
	ui.btnStop->setEnabled(true);
	statusLabel->setText("Fitting data....");
	qApp->processEvents();

	//Get approximate model dimensions for contiuous update
	stat_model_boundaries.resize(3);
	stat_model_boundaries = saxs::convertCoordinateToCylinder(globalfittingobject->m_model, plot_model_cyl);
	stat_model_boundaries[0] = stat_model_boundaries[0] * 2.;
	if (globalfittingobject->m_num_stacks < 2)
	{
		stat_model_boundaries[1] = stat_model_boundaries[1] * 2.;
		stat_model_boundaries[2] = stat_model_boundaries[2] * 2.;
	}

	//Failsafe for multicore mode
	if (globalfittingobject->m_multicore && globalfittingobject->m_num_cores == 1)globalfittingobject->m_num_cores = 2;

	//Send to signaling function
	startDebyeFitCurrentModel();
}

//NO ACTION - Function called when fitting is done - merges results
void HelFitQt::finishedDebyeFitCurrentModel()
{
	//Save results to files
	std::string tempstring;

	int num_threads;
	if (globalfittingobject->m_multicore) num_threads = 1;
	else num_threads = globalfittingobject->m_num_cores;

	for (unsigned int i = 0; i < num_threads; i++)
	{
		tempstring = savefilename + "_" + std::to_string(i);
		saxs::writemodeltopdb(tempstring, result_tracer[i]->m_coordinate);
		saxs::writestackedmodeltopdb(tempstring, result_tracer[i]->m_coordinate, globalfittingobject);
		saxs::writedatatochi(tempstring, globalfittingobject->m_data_q, globalfittingobject->m_data_I, result_tracer[i]->m_fitted_I);
		saxs::writechisquretofile(tempstring, i, result_tracer);
	}

	//Check which had the best result
	int chipntr = 0;
	if (num_threads > 1)
	{
		double chimin = globalfittingobject->m_chi;

		//Check which solution was best
		for (unsigned int i = 0; i < globalfittingobject->m_num_cores; i++)
		{
			if (result_tracer[i]->m_chi <= chimin)
			{
				chimin = result_tracer[i]->m_chi;
				chipntr = i;
			}
		}
	}
	//If there is a better solution write best one into globalfittingobject
	globalfittingobject->m_chi = result_tracer[chipntr]->m_chi;
	globalfittingobject->m_model_I = result_tracer[chipntr]->m_model_I;
	globalfittingobject->m_mean_connectivity = result_tracer[chipntr]->m_mean_connectivity;
	globalfittingobject->m_mean_nr_contacts = result_tracer[chipntr]->m_mean_nr_contacts;
	globalfittingobject->m_compactness = result_tracer[chipntr]->m_compactness;
	for (int j = 0; j < globalfittingobject->m_model.size(); j++)
	{
		globalfittingobject->m_model[j]->m_x = result_tracer[chipntr]->m_coordinate[j]->m_x;
		globalfittingobject->m_model[j]->m_y = result_tracer[chipntr]->m_coordinate[j]->m_y;
		globalfittingobject->m_model[j]->m_z = result_tracer[chipntr]->m_coordinate[j]->m_z;
		globalfittingobject->m_model[j]->m_nr_contacts = result_tracer[chipntr]->m_coordinate[j]->m_nr_contacts;
	}

	//Output to extfield
	std::stringstream ss;
	ss << "The best result has core  " << chipntr << " with chi=";
	ss << std::scientific;
	ss << result_tracer[chipntr]->m_chi;
	writetolog(ss.str());

	writetolog("Results saved in " + savefilename + "_*.chi");

	//Get threads end time
	boost::chrono::time_point<boost::chrono::system_clock> end = boost::chrono::system_clock::now();

	//Calculate elapsed time
	boost::uint64_t elapsed_seconds = boost::chrono::duration_cast<boost::chrono::seconds>(end - start).count();
	writetolog("FINISHED!");
	writetolog("Duration\t\t = \t " + boost::lexical_cast<std::string>(elapsed_seconds) + " seconds.");

	//To make sure, recalculate DebyeStack
	globalfittingobject->m_diameter = saxs::get_cyldiameter_from_coordinates(globalfittingobject->m_model);
	globalfittingobject->m_contact_d_sq = std::pow(2 * saxs::get_critradius_from_coordinates(globalfittingobject->m_diameter,
		globalfittingobject->m_stack_spacing, globalfittingobject->m_model.size()), 2);
	
	saxs::calc_contacts_of_model(globalfittingobject);
	globalfittingobject->m_mean_nr_contacts = saxs::return_mean_contacts_of_model(globalfittingobject->m_model);
	globalfittingobject->m_mean_connectivity = saxs::return_mean_connectivity_of_fittingobject(globalfittingobject);
	globalfittingobject->m_compactness = saxs::return_compactness_of_fittingobject(globalfittingobject);
	writetolog("Ave. number of contacts\t = \t" + boost::lexical_cast<std::string>(globalfittingobject->m_mean_nr_contacts));
	writetolog("Global connectivity\t = \t" + boost::lexical_cast<std::string>(globalfittingobject->m_mean_connectivity));
	writetolog("Global compactness\t = \t" + boost::lexical_cast<std::string>(globalfittingobject->m_compactness));

	saxs::calcDebyeStackCurrentModel(globalfittingobject);

	//Clear previous calcdata
	//Old Datastructure!!!!
	imported_model_data.clear();
	x_imported_model_data.clear();
	y_imported_model_data.clear();

	//Move data to plotable qvector and normalize to I[0]
	for (int i = 0; i < globalfittingobject->m_data_q.size(); i++)
	{
		x_imported_model_data.push_back(globalfittingobject->m_data_q[i]);
		y_imported_model_data.push_back(globalfittingobject->m_fitted_I[i]);
		imported_model_data.push_back(saxs::scatteringdata_sp(new saxs::scatteringdata(
			globalfittingobject->m_data_q[i], globalfittingobject->m_fitted_I[i])));
	}

	//GUI changes
	plotCalcModel();
	UpdateFittedModel(globalfittingobject);

	//Re-enabking interface
	enable_ui();
	ui.btnStop->setEnabled(false);
	statusLabel->setText("Fitting done! Please check results...");
	statusProgressBar->setValue(0);

}


/***************************************************************************************************************/
// Event Handling - Sliders / Spinboxes
/***************************************************************************************************************/

//Disabling right clicks
void HelFitQt::mousePressEvent(QMouseEvent *e)
{
	if (e->button() == Qt::RightButton)
	{
		saxs::error_message("No right clicks allowed!");
	}
}

//Eventfilter that opens PopUp window if data-range box is pressed
bool HelFitQt::eventFilter(QObject *obj, QEvent *ev)
{
	//Filtering out keyboard input and opening popup
	if (obj == ui.spbDataMax || obj == ui.spbDataMin || obj == ui.spbFitMin || obj == ui.spbFitMax)
	{
		if (ev->type() == QEvent::KeyPress || ev->type() == QEvent::MouseButtonDblClick) {

			//Check if Tab is clicked
			QKeyEvent *keyEvent = static_cast<QKeyEvent*>(ev);
			if (keyEvent->key() == Qt::Key_Tab) return false;

			QSpinBox * spb = qobject_cast<QSpinBox * > (obj);
			int point;
			int min = spb->minimum();
			int max = spb->maximum();
			int curr = spb->value();
			QString text = " ";
			if (obj == ui.spbDataMax) text = "...maximum plotting range:";
			if (obj == ui.spbDataMin) text = "...minimum plotting range:";
			if (obj == ui.spbFitMax) text = "...maximum fitting range:";
			if (obj == ui.spbFitMin) text = "...minimum fitting range:";
			point = saxs::input_message(text, curr, min, max);
			spb->setValue(point);
			return true;
		}
		else return false;
	}
	else {
		// standard event processing
		return HelFitQt::eventFilter(obj, ev);
	}
}

//Spinbox - Change data/fitting range
void HelFitQt::changedMinDataRange()
{
	if (int(ui.spbDataMin->value() - 1)<max_datarange)
	{
		min_datarange = int(ui.spbDataMin->value()-1);
		ui.spbFitMin->setMinimum(min_datarange + 1);
	}
	else
	{
		ui.spbDataMin->setValue(min_datarange + 1);
	}
	ui.wdgtDataPlot->xAxis->setRange(x_expdata[min_datarange], x_expdata[max_datarange]);
	double ymin = *std::min_element(y_expdata.constBegin() + min_datarange, y_expdata.constBegin() + max_datarange);
	double ymax = *std::max_element(y_expdata.constBegin() + min_datarange, y_expdata.constBegin() + max_datarange);
	ui.wdgtDataPlot->yAxis->setRange(ymin*0.5, 1.5*ymax);
	ui.wdgtDataPlot->replot();
}

//Spinbox - Change data/fitting range
void HelFitQt::changedMaxDataRange()
{
	if (int(ui.spbDataMax->value() - 1) > min_datarange)
	{
		max_datarange = int(ui.spbDataMax->value() - 1);
		ui.spbFitMax->setMaximum(max_datarange + 1);
	}
	else
	{
		ui.spbDataMax->setValue(max_datarange + 1);
	}
	ui.wdgtDataPlot->xAxis->setRange(x_expdata[min_datarange], x_expdata[max_datarange]);
	double ymin = *std::min_element(y_expdata.constBegin() + min_datarange, y_expdata.constBegin() + max_datarange);
	double ymax = *std::max_element(y_expdata.constBegin() + min_datarange, y_expdata.constBegin() + max_datarange);
	ui.wdgtDataPlot->yAxis->setRange(ymin*0.5, 1.5*ymax);
	ui.wdgtDataPlot->replot();
}

//Spinbox - Change data/fitting range
void HelFitQt::changedMinFitRange()
{
	if (int(ui.spbFitMin->value() - 1) < max_fitrange)
	{
		min_fitrange = int(ui.spbFitMin->value() - 1);
		num_calcpoints = max_fitrange  + 1 - min_fitrange;
		replotData();
		ui.lineQmin->setText(QString::number(x_expdata[min_fitrange]));
		ui.lineCalcQmin->setText(QString::number(x_expdata[min_fitrange]));
		ui.spbCalcNrPoints->setValue(num_calcpoints);
	}
	else
	{
		ui.spbFitMin->setValue(min_fitrange + 1);
	}
}

//Spinbox - Change data/fitting range
void HelFitQt::changedMaxFitRange()
{
	if (int(ui.spbFitMax ->value()-1)>min_fitrange)
	{
		max_fitrange = int(ui.spbFitMax->value() - 1);
		num_calcpoints = max_fitrange + 1 - min_fitrange;
		replotData();
		ui.lineQmax->setText(QString::number(x_expdata[max_fitrange]));
		ui.lineCalcQmax->setText(QString::number(x_expdata[max_fitrange]));
		ui.spbCalcNrPoints->setValue(num_calcpoints);
	}
	else
	{
		ui.spbFitMax->setValue(max_fitrange + 2);
	}
}

//Spinbox - Change data/fitting range
void HelFitQt::changedModelRotationPhi(int sliderpos)
{
	phi_moodel_rotation = double(sliderpos) / 500 * saxs::pi;
	if (model_loaded || model_generated) replotModel();
}



/***************************************************************************************************************/
// Graphing/Output functions
/***************************************************************************************************************/

//Initialization
void HelFitQt::initialize_graphs()
{
	//data plotting graphs
	ui.wdgtDataPlot->addGraph();
	// give the axes some labels:
	ui.wdgtDataPlot->xAxis->setLabel("q (1/nm)");
	ui.wdgtDataPlot->yAxis->setLabel("I(q) (-)");
	ui.wdgtDataPlot->yAxis->setRange(0, 1);
	ui.wdgtDataPlot->yAxis->setScaleType(QCPAxis::stLogarithmic);
	ui.wdgtDataPlot->yAxis->setNumberFormat("eb"); // e = exponential, b = beautiful decimal powers
	ui.wdgtDataPlot->yAxis->setNumberPrecision(0); // makes sure "1*10^4" is displayed only as "10^4"
	ui.wdgtDataPlot->yAxis->setScaleLogBase(10);
	ui.wdgtDataPlot->replot();

	//model plotting graphs
	ui.wdgtXZplot->addGraph();
	ui.wdgtXZplot->xAxis->setRange(-20, 20);
	ui.wdgtXZplot->xAxis->grid()->setVisible(false);
	ui.wdgtXZplot->xAxis2->setVisible(true);
	ui.wdgtXZplot->xAxis2->setTickLabels(false);
	ui.wdgtXZplot->yAxis->setRange(0, 50);
	ui.wdgtXZplot->yAxis2->setVisible(true);
	ui.wdgtXZplot->yAxis2->setTickLabels(false);
	ui.wdgtXZplot->yAxis->grid()->setVisible(false);
	ui.wdgtXZplot->replot();

	ui.wdgtYZplot->addGraph();
	ui.wdgtYZplot->xAxis->setRange(-20, 20);
	ui.wdgtYZplot->xAxis->grid()->setVisible(false);
	ui.wdgtYZplot->xAxis2->setVisible(true);
	ui.wdgtYZplot->xAxis2->setTickLabels(false);
	ui.wdgtYZplot->yAxis->setRange(0, 50);
	ui.wdgtYZplot->yAxis2->setVisible(true);
	ui.wdgtYZplot->yAxis2->setTickLabels(false);
	ui.wdgtYZplot->yAxis->grid()->setVisible(false);
	ui.wdgtYZplot->replot();

	ui.wdgtXYplot->addGraph();
	ui.wdgtXYplot->xAxis->setRange(-20, 20);
	ui.wdgtXYplot->xAxis->grid()->setVisible(false);
	ui.wdgtXYplot->xAxis2->setVisible(true);
	ui.wdgtXYplot->xAxis2->setTickLabels(false);
	ui.wdgtXYplot->yAxis->setRange(-20, 20);
	ui.wdgtXYplot->yAxis2->setVisible(true);
	ui.wdgtXYplot->yAxis2->setTickLabels(false);
	ui.wdgtXYplot->yAxis->grid()->setVisible(false);
	ui.wdgtXYplot->replot();
}

//Plotting scattering data after dataload
void HelFitQt::plot_data()
{
	x_fitdata.resize(num_calcpoints);
	y_fitdata.resize(num_calcpoints);
	e_fitdata.resize(num_calcpoints);
	x_fitdata = x_expdata;
	y_fitdata = y_expdata;
	e_fitdata = e_expdata;

	//add experimental data
	ui.wdgtDataPlot->addGraph();
	ui.wdgtDataPlot->graph(0)->setPen(QPen(Qt::blue));
	ui.wdgtDataPlot->graph(0)->setData(x_expdata, y_expdata);

	ui.wdgtDataPlot->addGraph();
	QPen RedPen;
	RedPen.setColor(Qt::red);
	RedPen.setWidthF(2);
	ui.wdgtDataPlot->graph(1)->setPen(RedPen);
	ui.wdgtDataPlot->graph(1)->setScatterStyle(QCPScatterStyle(QCPScatterStyle::ssCircle, Qt::red, Qt::red, 1));
	ui.wdgtDataPlot->graph(1)->setErrorType(QCPGraph::etValue);
	ui.wdgtDataPlot->graph(1)->setErrorPen(QPen(Qt::red));
	ui.wdgtDataPlot->graph(1)->setDataValueError(x_fitdata, y_fitdata, e_fitdata);

	// give the axes some labels:
	ui.wdgtDataPlot->xAxis->setLabel("q (1/nm)");
	ui.wdgtDataPlot->yAxis->setLabel("I(q) (-)");


	// set axes ranges, so we see all data:
	double xmin = x_expdata[min_datarange];
	double xmax = x_expdata[max_datarange];
	ui.wdgtDataPlot->xAxis->setRange(xmin, xmax);
	if (ui.chkbLogLogPlot->checkState() == Qt::Checked)
	{
		ui.wdgtDataPlot->xAxis->setScaleLogBase(10);
		ui.wdgtDataPlot->xAxis->setScaleType(QCPAxis::stLogarithmic);
	}
	else
	{
		ui.wdgtDataPlot->xAxis->setScaleType(QCPAxis::stLinear);
	}

	double ymin = *std::min_element(y_expdata.constBegin() + min_datarange, y_expdata.constBegin() + max_datarange);
	double ymax = *std::max_element(y_expdata.constBegin() + min_datarange, y_expdata.constBegin() + max_datarange);
	ui.wdgtDataPlot->yAxis->setRange(ymin*0.5, 1.5*ymax);
	ui.wdgtDataPlot->yAxis->setScaleType(QCPAxis::stLogarithmic);
	ui.wdgtDataPlot->yAxis->setNumberFormat("eb"); // e = exponential, b = beautiful decimal powers
	ui.wdgtDataPlot->yAxis->setNumberPrecision(0); // makes sure "1*10^4" is displayed only as "10^4"
	ui.wdgtDataPlot->yAxis->setScaleLogBase(10);
	ui.wdgtDataPlot->replot();
}

//Plot scattering data of DebyeCalc model
void HelFitQt::plotCalcModel()
{
	if (data_loaded)
	{
		double chi = globalfittingobject->m_chi;
		std::string str;
		switch (curvefit_weight) {
		case (0) : str = "Chisquare (I(q)*q)\t = \t"; break;
		case (1) : str = "Chisquare (I(q)*q^2)\t = \t"; break;
		case (2) : str = "Chisquare (-)\t = \t"; break;
		}
		writetolog(str + boost::lexical_cast<std::string>(chi));
	}
	int k = 0;
	ui.wdgtDataPlot->addGraph();
	if (data_loaded) { k = 2; }
	QPen BlackPen;
	BlackPen.setColor(Qt::black);
	BlackPen.setWidthF(2);
	ui.wdgtDataPlot->graph(k)->setPen(BlackPen);
	ui.wdgtDataPlot->graph(k)->setData(x_imported_model_data, y_imported_model_data);
	if (!data_loaded)
	{
		double xmin = x_imported_model_data[0];
		double xmax = x_imported_model_data[x_imported_model_data.size() - 1];
		ui.wdgtDataPlot->xAxis->setRange(xmin, xmax);
		double ymin = *std::min_element(y_imported_model_data.constBegin(), y_imported_model_data.constEnd() - 1);
		double ymax = *std::max_element(y_imported_model_data.constBegin(), y_imported_model_data.constEnd() - 1);
		//ui.wdgtDataPlot->yAxis->setRange(0.00001,1.5);
		ui.wdgtDataPlot->yAxis->setRange(0.5*ymin, 1.5*ymax);
		ui.wdgtDataPlot->yAxis->setScaleType(QCPAxis::stLogarithmic);
		ui.wdgtDataPlot->yAxis->setNumberFormat("eb"); // e = exponential, b = beautiful decimal powers
		ui.wdgtDataPlot->yAxis->setNumberPrecision(0); // makes sure "1*10^4" is displayed only as "10^4"
		ui.wdgtDataPlot->yAxis->setScaleLogBase(10);
	}

	ui.wdgtDataPlot->replot();
	ui.tabWidget->setCurrentIndex(0);
}

//Update scattering data of fitted model (Live update during Fitting)
void HelFitQt::UpdateFittedModel(saxs::fittingobject_sp& fitobject)
{

	clear_model();
	imported_model = fitobject->m_model;
	//Copy original data into plotting qvectors
	saxs::convertCoordinateToVector(imported_model, model_x, model_y, model_z);
	std::vector<double> model_boundaries(3);
	model_boundaries = saxs::convertCoordinateToCylinder(imported_model, plot_model_cyl);
	phi_moodel_rotation = saxs::pi;
	ui.hsliderPhiRot->setValue(500);
	model_loaded = true;
	ui.menuDebyeCurrModel->setEnabled(true);
	plot_model(model_boundaries);
	ui.tabWidget->setCurrentIndex(1);
}

//Replot scattering data if range is changed
void HelFitQt::replotData()
{
	x_fitdata.clear();
	y_fitdata.clear();
	e_fitdata.clear();
	x_fitdata.resize(num_calcpoints);
	y_fitdata.resize(num_calcpoints);
	e_fitdata.resize(num_calcpoints);

	for (int i = min_fitrange; i < max_fitrange+1; i++)
	{
		x_fitdata[i- min_fitrange] = x_expdata[i];
		y_fitdata[i- min_fitrange] = y_expdata[i];
		e_fitdata[i - min_fitrange] = e_expdata[i];
	}
	if (ui.chkbLogLogPlot->checkState() == Qt::Checked)
	{
		ui.wdgtDataPlot->xAxis->setScaleLogBase(10);
		ui.wdgtDataPlot->xAxis->setScaleType(QCPAxis::stLogarithmic);
	}
	else
	{
		ui.wdgtDataPlot->xAxis->setScaleType(QCPAxis::stLinear);
	}
	ui.wdgtDataPlot->graph(1)->setDataValueError(x_fitdata, y_fitdata, e_fitdata);
	ui.wdgtDataPlot->replot();
}

//Plotting point model after import
void HelFitQt::plot_model(std::vector<double>& boundaries)
{
	double rmax = boundaries[0] * 1.1;
	double zmax = boundaries[1];
	double zmin = boundaries[2];
	QPen RedPen;
	RedPen.setWidth(2);
	RedPen.setColor(Qt::green);
	ui.wdgtXZplot->addGraph();
	ui.wdgtXZplot->graph(1)->setPen(RedPen);
	ui.wdgtXZplot->graph(1)->setLineStyle(QCPGraph::lsNone);
	ui.wdgtXZplot->graph(1)->setScatterStyle(QCPScatterStyle(QCPScatterStyle::ssDisc, 1));
	ui.wdgtXZplot->graph(1)->setData(model_x, model_z);
	ui.wdgtXZplot->xAxis->setRange(-rmax, rmax);
	ui.wdgtXZplot->yAxis->setRange(zmin, zmax);
	ui.wdgtXZplot->replot();

	ui.wdgtYZplot->addGraph();
	ui.wdgtYZplot->graph(1)->setPen(RedPen);
	ui.wdgtYZplot->graph(1)->setLineStyle(QCPGraph::lsNone);
	ui.wdgtYZplot->graph(1)->setScatterStyle(QCPScatterStyle(QCPScatterStyle::ssDisc, 1));
	ui.wdgtYZplot->graph(1)->setData(model_y, model_z);
	ui.wdgtYZplot->xAxis->setRange(-rmax, rmax);
	ui.wdgtYZplot->yAxis->setRange(zmin, zmax);
	ui.wdgtYZplot->replot();

	ui.wdgtXYplot->addGraph();
	ui.wdgtXYplot->graph(1)->setPen(RedPen);
	ui.wdgtXYplot->graph(1)->setLineStyle(QCPGraph::lsNone);
	ui.wdgtXYplot->graph(1)->setScatterStyle(QCPScatterStyle(QCPScatterStyle::ssDisc, 1));
	ui.wdgtXYplot->graph(1)->setData(model_x, model_y);
	ui.wdgtXYplot->xAxis->setRange(-rmax, rmax);
	ui.wdgtXYplot->yAxis->setRange(-rmax, rmax);
	ui.wdgtXYplot->replot();
}

//Replot Model if rotation slider is moved
void HelFitQt::replotModel() 
{
	saxs::convertCylinderToQVectors(plot_model_cyl,model_x, model_y,phi_moodel_rotation);
	ui.wdgtXZplot->graph(1)->setData(model_x, model_z);
	ui.wdgtXZplot->replot();
	ui.wdgtYZplot->graph(1)->setData(model_y, model_z);
	ui.wdgtYZplot->replot();
	ui.wdgtXYplot->graph(1)->setData(model_x, model_y);
	ui.wdgtXYplot->replot();
}

//Write string into logging box
void HelFitQt::writetolog(std::string text)
{
	std::string time = saxs::now_str();
	QString qstr = QString::fromStdString(time + ": " + text);
	ui.textOutput->append(qstr);
	qApp->processEvents();
}

//External signal handler to write string into logging box (Live update during Fitting)
void HelFitQt::writetologext(std::string text)
{
	writetolog(text);
}

//Disabling interface
void HelFitQt::disable_ui()
{
	ui.menuBar->setEnabled(false);
	ui.groupBox_1->setEnabled(false);
	ui.groupBox_2->setEnabled(false);
	ui.groupBox_3->setEnabled(false);
	ui.btnFit->setEnabled(false);
	ui.hsliderPhiRot->setEnabled(false);
}

//Enabling Interface
void HelFitQt::enable_ui()
{
	ui.menuBar->setEnabled(true);
	ui.groupBox_1->setEnabled(true);
	ui.groupBox_2->setEnabled(true);
	ui.groupBox_3->setEnabled(true);
	ui.btnFit->setEnabled(true);
	ui.hsliderPhiRot->setEnabled(true);
}


/***************************************************************************************************************/
// Signal Handling
/***************************************************************************************************************/

//Receives "Done" signals from threads - if all done run "finishedDebyeFitCurrentModel"
void HelFitQt::fittingthread_done(int finished_thread)
{
	result_tracer[finished_thread]->m_thread_status = true;
	bool helper = true;
	for (int i = 0; i < result_tracer.size(); i++)
	{
		if (!result_tracer[i]->m_thread_status) helper = false;
	}
	if (helper) finishedDebyeFitCurrentModel();
}

//Receives live scattering data from fitting threads (of core 0)
void HelFitQt::plot_current_thread_Data(QVector<double> received_fitted_I)
{
	y_imported_model_data = received_fitted_I;
	QPen BlackPen;
	BlackPen.setColor(Qt::black);
	BlackPen.setWidthF(2);
	ui.wdgtDataPlot->graph(2)->setPen(BlackPen);
	ui.wdgtDataPlot->graph(2)->setData(x_imported_model_data, y_imported_model_data);
	ui.wdgtDataPlot->replot();
}

//Sends stop signal to fitting threads
void HelFitQt::on_btnStop_clicked()
{
	emit sig_stop_thread();
	ui.btnStop->setEnabled(false);
	statusProgressBar->setValue(100);
}

//Receives progress to update status bar
void HelFitQt::update_statusbar(double progress)
{
	QString str;
	int progress_int = int(progress * 100);

	str = "Fitting... " + QString::number(progress_int) + "% done...";
	statusLabel->setText(str);
	statusProgressBar->setValue(progress_int);
	
}

//Receives live model data as x,y,z qvectors
void HelFitQt::live_update_model(QVector<double> qvec_x, QVector<double>qvec_y, QVector<double>qvec_z)
{
	model_x = qvec_x;
	model_y = qvec_y;
	model_z = qvec_z;
	plot_model(stat_model_boundaries);
	emit sig_ui_ready();
}