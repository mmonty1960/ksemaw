/*Author: Marco Montecchi
Department of Energy Technologies
ENEA C.R. Casaccia
Roma - Italy

kSEMAW is a workspace for the analysis of
Spectrophotometric (SP), Ellipsometric (ELI) and
Photothermal Deflection Spectroscopy (PDS) measurements


   Copyright (C) 2020  Marco Montecchi

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
#ifndef KSEMAW_H
#define KSEMAW_H

#include "ui_ksemaw.h"
using namespace std;

namespace Ui{
class ksemaw;
}


class ksemaw : public QWidget//, private Ui::ksemawDLG
{
    Q_OBJECT

 
private:
    void Setnk(int ifile);
    void Clrnk(int ifile);
    void updateMatenk(int, QString);
    void ReadSetting(QString);
    Ui::ksemaw *ui;

public:
    ksemaw(QWidget *parent=0);
    QString fnproject,pathroot,info,fnk[9],fnSample,fileStore,fStdSpect,
    fnTn,fnTp,fnRn,fnRp,fnR1,fnApds,fnE1,fnE2,fnE3,fnE4,fnFnk,ParFitLab[12],fctrl;
    int ink,lastIndex,ifn,npp,ppm[18],n,nlayer,lastTab,occupyPF,iwspj,lastTabB5;
    double pf[8][22],pm[201][6],par[61][6],cnk[17][4],rxy[31][5];
    QMap<QString, QComboBox*> idToComboBox;
    QMap<QString, QCheckBox*> idToCheckBox;
    QMap<QString, QLineEdit*> idToLineEdit;
    QMap<QString, QDoubleSpinBox*> idToDoubleSpinBox;

public slots:
     void LoadProject();
     void SaveProject();
     void LoadFilenk();
     void ClrFnk();
     void SaveFnk();
     void SaveSetting(int iCall);
     void ClearTemp();
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
     //void saveOsc(int);
     void tabChanged();
     void PanFitEnable();
     void PanFitPar();
     void PanFitChoice();
     void AdjTheta();
     void RefreshModel();
     void SetModel(const int &);
     void closeEvent ( QCloseEvent * event );
     void fileChanged(const QString &);
     void PlotME();
     void Simula();
     void PlotAve();
     void readCursor();
     void searchNK();
     void save1NK();
     void RefTrackG();
     void TrackFromWLmin();
     void TrackFromWLmax();
     void saveNKtmp();
     void RefIbridG();
     void FitN();
     void FitNK();
     void FitE1E2();
     void IbridPlotFit();
     void IbridPlotIbrid();
     void IbridFit();
     void IbridOne();
     void IbridOneStore();
     void ClearTempIbri();
     void AdjSave();
     void EnableAllWl();
     void EnableAllNK();
     void EditWl();
     void EditNK();
     void GoBest();
     void GoPrevious();
     void GoNext();
     //void savePanFit();
     void rangeWL();
     void refreshFitPar();
     void saveSim();
     void PlotAbsEL();
     void Ready(bool slwa);
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
};

#endif
