#ifndef DIALOG_FITTINGVARS_H
#define DIALOG_FITTINGVARS_H

#include <helfitqt.h>
#include "ui_dialog_fittingvars.h"

class dialog_fittingvars : public QDialog
{
	Q_OBJECT

public:
	dialog_fittingvars(QWidget *parent = 0);
	~dialog_fittingvars();
	bool dialogOK = false;
	double alpha_conn;
	double tau_conn;
	double beta_comp;
	double sigma_comp;

	void write_current_vars();

private slots:
	void on_okButton_clicked();
	void on_cancelButton_clicked();

private:
	Ui::dialog_fittingvars_class *ui;
};

#endif // DIALOG_MODELVARS_H#pragma once
