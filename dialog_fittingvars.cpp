#include "dialog_fittingvars.h"
#include "dialog_fittingvars.h"
#include <QDoubleValidator>

dialog_fittingvars::dialog_fittingvars(QWidget *parent)
	: QDialog(parent),
	ui(new Ui::dialog_fittingvars_class)
{
	ui->setupUi(this);
	//Initialize boxvalidator
	ui->lineAlphaConn->setValidator(new QDoubleValidator(-1,10000,4,this));
	ui->lineTauConn->setValidator(new QDoubleValidator(0.0001, 1000, 4, this));
	ui->lineBetaComp->setValidator(new QDoubleValidator(-1, 1000, 4, this));
	ui->lineSigmaComp->setValidator(new QDoubleValidator(0.0001, 1000, 4, this));
	ui->lineGammaHelBias->setValidator(new QDoubleValidator(0, 2, 4, this));
	ui->lineRandomSeedScalar->setValidator(new QDoubleValidator(0, 2, 4, this));
}

dialog_fittingvars::~dialog_fittingvars()
{
	delete ui;
}

void dialog_fittingvars::write_current_vars()
{
	//Write current values into ui
	ui->lineAlphaConn->setText(QString::number(alpha_conn));
	ui->lineTauConn->setText(QString::number(tau_conn));
	ui->lineBetaComp->setText(QString::number(beta_comp));
	ui->lineSigmaComp->setText(QString::number(sigma_comp));
	ui->lineGammaHelBias->setText(QString::number(gamma_helbias));
	ui->lineRandomSeedScalar->setText(QString::number(random_seed_scalar));
}

//Custom functions
void dialog_fittingvars::on_okButton_clicked()
{
	dialogOK = true;
	alpha_conn = double(ui->lineAlphaConn->text().toDouble());
	tau_conn = double(ui->lineTauConn->text().toDouble());
	beta_comp = double(ui->lineBetaComp->text().toDouble());
	sigma_comp = double(ui->lineSigmaComp->text().toDouble());
	gamma_helbias = double(ui->lineGammaHelBias->text().toDouble());
	random_seed_scalar = double(ui->lineRandomSeedScalar->text().toDouble());

	dialog_fittingvars::accept();
}

void dialog_fittingvars::on_cancelButton_clicked()
{
	dialogOK = false;
	dialog_fittingvars::reject();
}
