#ifndef DATA_DEF_H
#define DATA_DEF_H

#include <boost\thread.hpp>
#include <boost\shared_ptr.hpp>
#include <QtCore>
#include <qthread.h>

Q_DECLARE_METATYPE(std::string);

namespace saxs
{
	//Initiliazation of the sinc lookup table
	static std::vector<double> sinc_lookup(100001, 1.0f);
	//Pi lookup
	const double pi = 3.14159265359;


//------------------------------------------------------------------------------
//	Coordinate struct definition and typedef
//------------------------------------------------------------------------------
struct coordinate
{
	coordinate(double x = 0, double y = 0, double z = 0, int nr_contacts=0) :
		m_x(x), m_y(y), m_z(z),m_nr_contacts(nr_contacts), m_nr_coorditations(0) { }

	double m_x;
	double m_y;
	double m_z;
	int m_nr_contacts;
	int m_nr_coorditations;
};

typedef boost::shared_ptr<coordinate> coordinate_sp;

//------------------------------------------------------------------------------
//	Model struct definition and typedef
//------------------------------------------------------------------------------
struct model
{
	std::vector<coordinate_sp> m_coordinate;
	double m_chi;
	double m_mean_nr_contacts;
	double m_mean_connectivity;
	double m_compactness;
	double m_targetf;
	std::vector<double> m_model_I;
	std::vector<double> m_fitted_I;
	std::vector<double> m_chi_tracer;
	std::vector<double> m_contact_tracer;
	std::vector<double> m_connectivity_tracer;
	std::vector<double> m_compactness_tracer;
	std::vector<double> m_target_f_tracer;
	bool m_thread_status; //false - running / true - done
};

typedef boost::shared_ptr<model> model_sp;

//------------------------------------------------------------------------------
//	Scatteringdata struct definition and typedef
//------------------------------------------------------------------------------
struct scatteringdata
{
	scatteringdata(double q = 0, double I = 0, double e=0) :
		m_q(q), m_I(I), m_e(e) { }

	double m_q;
	double m_I;
	double m_e;
};

typedef boost::shared_ptr<scatteringdata> scatteringdata_sp;

//------------------------------------------------------------------------------
//	Fitting object struct definition and typedef
//------------------------------------------------------------------------------
struct fittingobject
{
	/*fittingobject(std::vector<saxs::coordinate_sp> model, std::vector<double> data_q, \
		std::vector<double> data_I, std::vector<double> model_I, double chi=0) : 
		m_model(model), m_data_q (data_q), m_data_I(data_I), m_model_I(model_I), m_chi(chi) {}
		*/
	std::vector<coordinate_sp> m_model;
	std::vector<double> m_data_q;
	std::vector<double> m_data_I;
	std::vector<double> m_ff_I;
	std::vector<double> m_data_e;
	std::vector<double> m_model_I;
	std::vector<double> m_fitted_I;
	double m_chi;
	double m_compactness;
	double m_mean_nr_contacts;
	double m_mean_connectivity;
	double m_mean_nr_coordination;
	std::pair<double, double> m_linreg_results;
	double m_norm;
	double m_norm_offset;
	//Parameters
	int m_num_stacks;
	double m_stack_spacing;
	double m_diameter;
	double m_contact_d_sq;
	double m_sigma_ff;
	//SA Parameters
	int m_num_cores;
	bool m_multicore;
	double m_start_temp;
	double m_delta_temp;
	int m_num_runs;
	int m_weighing;
	//PotentialFields
	double m_alphaConn;
	double m_tauConn;
	double m_betaComp;
	double m_sigmaComp;
};
typedef boost::shared_ptr<fittingobject> fittingobject_sp;

//------------------------------------------------------------------------------
//	Merge object struct definition and typedef
//------------------------------------------------------------------------------
struct mergeobject 
{
	std::vector<saxs::coordinate_sp> m_model;
	std::vector<saxs::coordinate_sp> m_model_cyl;
	std::vector<bool> m_inverted;
	double m_mean_dist;
	std::vector<double> m_NSD;
	double m_mean_NSD;
	double m_sdev_NSD;
	std::vector<double> m_NSD_phi;
};
typedef boost::shared_ptr<mergeobject> mergeobject_sp;

}	//End of namespace


/***************************************************************************************************************/
// QT Class:DebyeFitThread constructor
/***************************************************************************************************************/
class DebyeFitThread : public QThread
{
	Q_OBJECT

	public:
		saxs::fittingobject_sp reftofittingobject;
		std::vector<saxs::model_sp> result_tracer;
		int coreid;
		int n_runs;
		double starttemp;
		double endtemp;
		double deltatemp;
		bool stop_thread = false;
		bool ui_ready = false;

		void run();
		void run2();

	public slots:
		void stop();
		void receive_ui_ready();
		

	signals:
		void sig_writelogext(std::string);
		void sig_fittingthread_done(int);
		void sig_plot_current_thread_Data(QVector<double>);
		void sig_update_statusbar(double);
		void sig_live_update_model(QVector<double>, QVector<double>, QVector<double>);

};

#endif /*!DATA_DEF_H*/