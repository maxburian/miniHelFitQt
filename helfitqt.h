#ifndef HELFITQT_H
#define HELFITQT_H

//Include database dependencies
//QT
//IF STATIC BUILD:

#include <QtPlugin>

Q_IMPORT_PLUGIN(QWindowsIntegrationPlugin);

//ENDIF STATIC BUILD:
#include <QtWidgets/QMainWindow>
#include <QLabel>
#include <QProgressBar>
#include <qapplication.h>
#include <qobject.h>
#include <qstring.h>
#include <qvector.h>
#include <QMouseEvent>
#include <qevent.h>
//STD
#include <string.h>
#include <vector>
//Boost
#include <boost/chrono.hpp>
#include <boost/scoped_ptr.hpp>
//Self-Written
#include "data_def.h"
#include "ui_helfitqt.h"
#include "dialog_modelvars.h"
#include "dialog_fittingvars.h"

class HelFitQt : public QMainWindow
{
	Q_OBJECT

public:
	HelFitQt(QWidget *parent = 0);
	~HelFitQt();
	void initialize_graphs();
	void init();

	//The main Fittingobject!!
	saxs::fittingobject_sp globalfittingobject{new saxs::fittingobject};
	std::vector<saxs::model_sp> result_tracer;

	//Fileloadparams
	QString scatterfilepath;
	QString scatterfilename;
	std::string savefilename;

	//RandomModelGeneration
	std::vector<saxs::coordinate_sp> generated_model;
	bool model_generated = false;

	//Scatteringdata objets
	std::vector<saxs::scatteringdata_sp> imported_data;
	QVector<double> x_expdata;
	QVector<double> y_expdata;
	QVector<double> e_expdata;
	bool data_loaded = false;
	QVector<double> x_fitdata;
	QVector<double> y_fitdata;
	QVector<double> e_fitdata;
	int min_datarange;
	int max_datarange;
	int num_datapoints;
	int min_fitrange;
	int max_fitrange;

	//Imported Model objects
	std::vector<saxs::scatteringdata_sp> imported_model_data;
	std::vector<saxs::coordinate_sp> imported_model;
	std::vector<saxs::coordinate_sp> plot_model_cyl;
	QVector<double> model_x;
	QVector<double> model_y;
	QVector<double> model_z;
	double phi_moodel_rotation;
	bool model_loaded = false;
	bool new_model=true;

	//Helperobjects
	boost::chrono::time_point<boost::chrono::system_clock> start;
	std::vector<double> stat_model_boundaries;

	//Imported model calculation objects
	//ranges from x_imported_model_data min_calcrange to max_calcrange
	QVector<double> x_imported_model_data;
	QVector<double> y_imported_model_data;
	std::pair<double, double> linreg_results;
	double calc_qmin;
	double calc_qmax;
	int num_calcpoints;
	int num_stacks;
	int core_number;

	//Weight for curvefitting: 0=Is, 1=Iss, 2=I
	int curvefit_weight=0;

public slots:
	void writetologext(std::string text);
	void fittingthread_done(int);
	void plot_current_thread_Data(QVector<double>);
	void update_statusbar(double);
	void live_update_model(QVector<double>, QVector<double>, QVector<double>);

private slots:
	//data load gui
	void clear_data();
	void on_btnLoadFile_clicked();
	void plot_data();
	void changedMinDataRange();
	void changedMaxDataRange();
	void changedMinFitRange();
	void changedMaxFitRange();
	void mousePressEvent(QMouseEvent *e);
	void replotData();

	//Fittingweight
	void actionWeight0();
	void actionWeight1();
	void actionWeight2();

	//model load gui
	void clear_model();
	void plot_model(std::vector<double>& boundaries);
	void loadModelFromData();
	void changedModelRotationPhi(int sliderpos);
	void replotModel();
	void plotCalcModel();
	void actionGenerateRandModel();
	void actionExpandCurrentModel();

	//calc gui
	void calcDebyeStackCurrent();
	void actionChangeFittingParams();

	//Interface handler
	void disable_ui();
	void enable_ui();

	//Merger
	void actionRunMerger();

	//Fitting
	void on_btnFit_clicked();
	void UpdateFittedModel(saxs::fittingobject_sp& fitobject);
	void startDebyeFitCurrentModel();
	void finishedDebyeFitCurrentModel();
	void on_btnStop_clicked();

	//Others
	void movedatatofittingobject();
	void writetolog(std::string text);
	void actionSaveModel();
	void actionSaveData();
	

private:
	Ui::HelFitQtClass ui;
	bool eventFilter(QObject *obj, QEvent *ev);

	//Status Bar
	QLabel *statusLabel;
	QProgressBar *statusProgressBar;

signals:
	void sig_stop_thread();
	void sig_ui_ready();
	


};
#endif // HELFITQT_H
