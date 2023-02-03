/*Main author: Marco Montecchi
             ENEA (Italy)
             email: marco.montecchi@enea.it
Porting to Windows and advanced oscillators by
             Alberto Mittiga
             ENEA (Italy)
             email: alberto.mittiga@enea.it


kSEMAW is a workspace for the analysis of
Spectrophotometric (SP), Ellipsometric (ELI) and
Photothermal Deflection Spectroscopy (PDS) measurements


   Copyright (C) 2022  Marco Montecchi

   This program is free software: you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation, either version 3 of the License, or
   (at your option) any later version.

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with this program.  If not, see <https://www.gnu.org/licenses/>.
*/
#ifndef KSEMAWC_H
#define KSEMAWC_H

#include "ui_ksemawc.h"
class QwtPlotZoomer;
class QwtPlotPicker;
class QwtPlotPanner;
class Plot;
class QPolygon;
using namespace std;

namespace Ui{
class ksemawc;
}


class ksemawc : public QWidget//, private Ui::ksemawcDLG
{
    Q_OBJECT

 
private:
    void Setnk(int ifile);
    void Clrnk(int ifile);
    void updateMatenk(int, QString);
    void ReadSetting(QString);
    Ui::ksemawc *ui;
    QwtPlotPicker* m_picker;

public:
    ksemawc(QWidget *parent=0);    
    QMap<QString, QComboBox*> idToComboBox;
    QMap<QString, QCheckBox*> idToCheckBox;
    QMap<QString, QLineEdit*> idToLineEdit;
    QMap<QString, QDoubleSpinBox*> idToDoubleSpinBox;

public Q_SLOTS:
     void LoadProject();
     void SaveProject();
     void LoadFilenk();
     void ClrFnk();
     void SaveFnk();
     void SaveSetting(int iCall);
     void Setnk1();
     void Setnk2();
     void Setnk3();
     void Setnk4();
     void Setnk5();
     void Setnk6();
     void Setnk7();
     void Setnk8();
     void mDw1();
     void mUp1();
     void mDw2();
     void mUp2();
     void mDw3();
     void mUp3();
     void mDw4();
     void mUp4();
     void mDw5();
     void mUp5();
     void mDw6();
     void mUp6();
     void mDw7();
     void mUp7();
     void mDw8();
     void mUp8();
     void mDw9();
     void mUp9();
     void mDwUp(int iLayer, int Dw1UpM1);
     void Clrnk1();
     void Clrnk2();
     void Clrnk3();
     void Clrnk4();
     void Clrnk5();
     void Clrnk6();
     void Clrnk7();
     void Clrnk8();
     void callSetSample();
     void setSample();
     void listMeas(const QString &);
     void pwTn (const int &);
     void pwTp (const int &);
     void pwRn (const int &);
     void pwRp (const int &);
     void pwR1 (const int &);
     void pwApds (const int &);
     void pwE1 (const int &);
     void pwE2 (const int &);
     void pwE3 (const int &);
     void pwE4 (const int &);
     void pwSubE1 (const int &);
     void pwSubE2 (const int &);
     void pwSubE3 (const int &);
     void pwSubE4 (const int &);
     void MCRange ();
     void listOsc();
     void tabChanged();
     void PanFitEnable();
     void PanFitPar();
     void PanFitChoice();
     void AdjTheta();
     void AdjRoughMax();
     void RefreshModel();
     void SetModel(const int &);
     void closeEvent ( QCloseEvent * event );
     void PlotMENK();
     void PlotME();
     void PlotNK(int iRD);
     void Simula();
     void CalcMis(double mc[15][202]);
     void PlotAve();
     void searchNK();
     void RefTrackG();
     void NumericalSearch();
     void RefIbridG();
     void FitN();
     void FitNK();
     void FitE1E2();
     void FitSelExpMeas();
     void IbridPlotFit();
     void IbridPlotIbrid();
     void IbridFit();
     void IbridOne();
     void IbridOneStore();
     void IbridKernel(QString r);
     void ClearTempIbri();
     void GoBest();
     void GoPrevious();
     void GoNext();
     void rangeWL();
     void refreshFitPar();
     void saveNKsim();
     void saveSim();
     void PlotAbsEL();
     void setOsc1();
     void setOsc2();
     void setOsc3();
     void setOsc4();
     void setOsc5();
     void setOsc6();
     void setOsc7();
     void setOsc8();
     void setOsc9();
     void setOsc10();
     void setOsc11();
     void setOsc12();
     void setOsc13();
     void setOsc14();
     void setOsc15();
     void setOsc16();
     void setOsc17();
     void setOsc18();
     void setOsc19();
     void setOsc20();
     void setOscN(int k);
     void manageLEwl();
     void setMat1();
     void setMat2();
     void setMat3();
     void setMat4();
     void setMat5();
     void setMat6();
     void setMat7();
     void setMat8();
     void setMat9();
     void setMat10();
     void setMat11();
     void setMat12();
     void setMat13();
     void setMat14();
     void setMat15();
     void setMat(int nM);
     void setEMA1();
     void setEMA2();
     void setEMA3();
     void setEMA4();
     void setEMA5();
     void setEMA6();
     void setEMA7();
     void setEMA8();
     void setEMA9();
     void setEMA10();
     void setEMA11();
     void setEMA12();
     void setEMA13();
     void setEMA14();
     void setEMA15();
     void setEMA(int nM);
     void setRifMir();
     void SPADA();
     void selectNsol();
     void Selected( const QPolygon polygon);
     void setKindOsc1();
     void setKindOsc2();
     void setKindOsc3();
     void setKindOsc4();
     void setKindOsc5();
     void setKindOsc6();
     void setKindOsc7();
     void setKindOsc8();
     void setKindOsc9();
     void setKindOsc10();
     void setKindOsc11();
     void setKindOsc12();
     void setKindOsc13();
     void setKindOsc14();
     void setKindOsc15();
     void setKindOsc16();
     void setKindOsc17();
     void setKindOsc18();
     void setKindOsc19();
     void setKindOsc20();
     void setFontDia();
     void msgErrLoad(QString where, QString fnERR);
};

#endif
