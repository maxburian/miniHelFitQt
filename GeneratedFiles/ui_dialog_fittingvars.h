/********************************************************************************
** Form generated from reading UI file 'dialog_fittingvars.ui'
**
** Created by: Qt User Interface Compiler version 5.6.0
**
** WARNING! All changes made in this file will be lost when recompiling UI file!
********************************************************************************/

#ifndef UI_DIALOG_FITTINGVARS_H
#define UI_DIALOG_FITTINGVARS_H

#include <QtCore/QVariant>
#include <QtWidgets/QAction>
#include <QtWidgets/QApplication>
#include <QtWidgets/QButtonGroup>
#include <QtWidgets/QDialog>
#include <QtWidgets/QFrame>
#include <QtWidgets/QHBoxLayout>
#include <QtWidgets/QHeaderView>
#include <QtWidgets/QLabel>
#include <QtWidgets/QLineEdit>
#include <QtWidgets/QPushButton>
#include <QtWidgets/QSpacerItem>
#include <QtWidgets/QWidget>

QT_BEGIN_NAMESPACE

class Ui_dialog_fittingvars_class
{
public:
    QWidget *layoutWidget;
    QHBoxLayout *hboxLayout;
    QSpacerItem *spacerItem;
    QPushButton *okButton;
    QPushButton *cancelButton;
    QSpacerItem *spacer;
    QLabel *label_15;
    QLineEdit *lineAlphaConn;
    QLabel *label;
    QLabel *label_2;
    QLabel *label_3;
    QLineEdit *lineTauConn;
    QLabel *label_4;
    QFrame *line_2;
    QLabel *label_5;
    QLineEdit *lineBetaComp;
    QLabel *label_6;
    QLabel *label_7;
    QLineEdit *lineSigmaComp;
    QLabel *label_8;
    QLineEdit *lineGammaHelBias;
    QLabel *label_9;
    QFrame *line_3;
    QLineEdit *lineRandomSeedScalar;
    QLabel *label_11;

    void setupUi(QDialog *dialog_fittingvars_class)
    {
        if (dialog_fittingvars_class->objectName().isEmpty())
            dialog_fittingvars_class->setObjectName(QStringLiteral("dialog_fittingvars_class"));
        dialog_fittingvars_class->resize(432, 381);
        dialog_fittingvars_class->setContextMenuPolicy(Qt::NoContextMenu);
        layoutWidget = new QWidget(dialog_fittingvars_class);
        layoutWidget->setObjectName(QStringLiteral("layoutWidget"));
        layoutWidget->setGeometry(QRect(0, 340, 431, 33));
        hboxLayout = new QHBoxLayout(layoutWidget);
        hboxLayout->setSpacing(6);
        hboxLayout->setObjectName(QStringLiteral("hboxLayout"));
        hboxLayout->setContentsMargins(0, 0, 0, 0);
        spacerItem = new QSpacerItem(131, 31, QSizePolicy::Expanding, QSizePolicy::Minimum);

        hboxLayout->addItem(spacerItem);

        okButton = new QPushButton(layoutWidget);
        okButton->setObjectName(QStringLiteral("okButton"));

        hboxLayout->addWidget(okButton);

        cancelButton = new QPushButton(layoutWidget);
        cancelButton->setObjectName(QStringLiteral("cancelButton"));

        hboxLayout->addWidget(cancelButton);

        spacer = new QSpacerItem(131, 31, QSizePolicy::Expanding, QSizePolicy::Minimum);

        hboxLayout->addItem(spacer);

        label_15 = new QLabel(dialog_fittingvars_class);
        label_15->setObjectName(QStringLiteral("label_15"));
        label_15->setGeometry(QRect(30, 20, 311, 20));
        QFont font;
        font.setPointSize(10);
        font.setBold(true);
        font.setWeight(75);
        label_15->setFont(font);
        lineAlphaConn = new QLineEdit(dialog_fittingvars_class);
        lineAlphaConn->setObjectName(QStringLiteral("lineAlphaConn"));
        lineAlphaConn->setGeometry(QRect(90, 110, 61, 22));
        label = new QLabel(dialog_fittingvars_class);
        label->setObjectName(QStringLiteral("label"));
        label->setGeometry(QRect(20, 60, 391, 16));
        label_2 = new QLabel(dialog_fittingvars_class);
        label_2->setObjectName(QStringLiteral("label_2"));
        label_2->setGeometry(QRect(40, 110, 40, 20));
        label_2->setLayoutDirection(Qt::LeftToRight);
        label_2->setAlignment(Qt::AlignRight|Qt::AlignTrailing|Qt::AlignVCenter);
        label_3 = new QLabel(dialog_fittingvars_class);
        label_3->setObjectName(QStringLiteral("label_3"));
        label_3->setGeometry(QRect(20, 90, 161, 16));
        lineTauConn = new QLineEdit(dialog_fittingvars_class);
        lineTauConn->setObjectName(QStringLiteral("lineTauConn"));
        lineTauConn->setGeometry(QRect(290, 110, 61, 22));
        label_4 = new QLabel(dialog_fittingvars_class);
        label_4->setObjectName(QStringLiteral("label_4"));
        label_4->setGeometry(QRect(240, 110, 40, 20));
        label_4->setLayoutDirection(Qt::LeftToRight);
        label_4->setAlignment(Qt::AlignRight|Qt::AlignTrailing|Qt::AlignVCenter);
        line_2 = new QFrame(dialog_fittingvars_class);
        line_2->setObjectName(QStringLiteral("line_2"));
        line_2->setGeometry(QRect(10, 190, 411, 20));
        line_2->setFrameShape(QFrame::HLine);
        line_2->setFrameShadow(QFrame::Sunken);
        label_5 = new QLabel(dialog_fittingvars_class);
        label_5->setObjectName(QStringLiteral("label_5"));
        label_5->setGeometry(QRect(240, 160, 40, 20));
        label_5->setLayoutDirection(Qt::LeftToRight);
        label_5->setAlignment(Qt::AlignRight|Qt::AlignTrailing|Qt::AlignVCenter);
        lineBetaComp = new QLineEdit(dialog_fittingvars_class);
        lineBetaComp->setObjectName(QStringLiteral("lineBetaComp"));
        lineBetaComp->setGeometry(QRect(90, 160, 61, 22));
        label_6 = new QLabel(dialog_fittingvars_class);
        label_6->setObjectName(QStringLiteral("label_6"));
        label_6->setGeometry(QRect(40, 160, 40, 20));
        label_6->setLayoutDirection(Qt::LeftToRight);
        label_6->setAlignment(Qt::AlignRight|Qt::AlignTrailing|Qt::AlignVCenter);
        label_7 = new QLabel(dialog_fittingvars_class);
        label_7->setObjectName(QStringLiteral("label_7"));
        label_7->setGeometry(QRect(20, 140, 161, 16));
        lineSigmaComp = new QLineEdit(dialog_fittingvars_class);
        lineSigmaComp->setObjectName(QStringLiteral("lineSigmaComp"));
        lineSigmaComp->setGeometry(QRect(290, 160, 61, 22));
        label_8 = new QLabel(dialog_fittingvars_class);
        label_8->setObjectName(QStringLiteral("label_8"));
        label_8->setGeometry(QRect(20, 210, 401, 16));
        lineGammaHelBias = new QLineEdit(dialog_fittingvars_class);
        lineGammaHelBias->setObjectName(QStringLiteral("lineGammaHelBias"));
        lineGammaHelBias->setGeometry(QRect(90, 230, 61, 22));
        label_9 = new QLabel(dialog_fittingvars_class);
        label_9->setObjectName(QStringLiteral("label_9"));
        label_9->setGeometry(QRect(40, 230, 40, 20));
        label_9->setLayoutDirection(Qt::LeftToRight);
        label_9->setAlignment(Qt::AlignRight|Qt::AlignTrailing|Qt::AlignVCenter);
        line_3 = new QFrame(dialog_fittingvars_class);
        line_3->setObjectName(QStringLiteral("line_3"));
        line_3->setGeometry(QRect(10, 260, 411, 20));
        line_3->setFrameShape(QFrame::HLine);
        line_3->setFrameShadow(QFrame::Sunken);
        lineRandomSeedScalar = new QLineEdit(dialog_fittingvars_class);
        lineRandomSeedScalar->setObjectName(QStringLiteral("lineRandomSeedScalar"));
        lineRandomSeedScalar->setGeometry(QRect(90, 300, 61, 22));
        label_11 = new QLabel(dialog_fittingvars_class);
        label_11->setObjectName(QStringLiteral("label_11"));
        label_11->setGeometry(QRect(20, 280, 401, 16));

        retranslateUi(dialog_fittingvars_class);
        QObject::connect(okButton, SIGNAL(clicked()), dialog_fittingvars_class, SLOT(accept()));
        QObject::connect(cancelButton, SIGNAL(clicked()), dialog_fittingvars_class, SLOT(reject()));

        QMetaObject::connectSlotsByName(dialog_fittingvars_class);
    } // setupUi

    void retranslateUi(QDialog *dialog_fittingvars_class)
    {
        dialog_fittingvars_class->setWindowTitle(QApplication::translate("dialog_fittingvars_class", "fitting parameters...", 0));
        okButton->setText(QApplication::translate("dialog_fittingvars_class", "OK", 0));
        cancelButton->setText(QApplication::translate("dialog_fittingvars_class", "Cancel", 0));
        label_15->setText(QApplication::translate("dialog_fittingvars_class", "Adjust fitting parameters:", 0));
        lineAlphaConn->setText(QApplication::translate("dialog_fittingvars_class", "-1", 0));
        label->setText(QApplication::translate("dialog_fittingvars_class", "Target function:       f = chi\302\262 \302\267 (1 +  \316\261 \302\267 conn. + \316\262 \302\267 comp.)", 0));
        label_2->setText(QApplication::translate("dialog_fittingvars_class", "\316\261:", 0));
        label_3->setText(QApplication::translate("dialog_fittingvars_class", "connectivity parameters:", 0));
        lineTauConn->setText(QApplication::translate("dialog_fittingvars_class", "2.0", 0));
        label_4->setText(QApplication::translate("dialog_fittingvars_class", "tau:", 0));
        label_5->setText(QApplication::translate("dialog_fittingvars_class", "sigma:", 0));
        lineBetaComp->setText(QApplication::translate("dialog_fittingvars_class", "-1", 0));
        label_6->setText(QApplication::translate("dialog_fittingvars_class", "\316\262:", 0));
        label_7->setText(QApplication::translate("dialog_fittingvars_class", "compactness parameters:", 0));
        lineSigmaComp->setText(QApplication::translate("dialog_fittingvars_class", "2.0", 0));
        label_8->setText(QApplication::translate("dialog_fittingvars_class", "helical bias parameter:         (only active in stacking mode)", 0));
        lineGammaHelBias->setText(QApplication::translate("dialog_fittingvars_class", "0.3", 0));
        label_9->setText(QApplication::translate("dialog_fittingvars_class", "\316\263:", 0));
        lineRandomSeedScalar->setText(QApplication::translate("dialog_fittingvars_class", "0.2", 0));
        label_11->setText(QApplication::translate("dialog_fittingvars_class", "random seed scalar:            (effective every 10th iteration)", 0));
    } // retranslateUi

};

namespace Ui {
    class dialog_fittingvars_class: public Ui_dialog_fittingvars_class {};
} // namespace Ui

QT_END_NAMESPACE

#endif // UI_DIALOG_FITTINGVARS_H
