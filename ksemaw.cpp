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
#include <QtGui>
#include "ksemaw.h"
#include<fstream>
#include <stdio.h>
#include <stdlib.h>
#include <qfile.h>
#include <qtextstream.h>
#include <math.h>
#include<iostream>
#include <QFileDialog>
#include <QMessageBox>
#include <unistd.h>

// if we include <QtGui> there is no need to include every class used: <QString>, <QFileDialog>,...


//funzioni invocate
void previewFile(QString filename, QString lab,QString& info,double& wmax,double& wmin);

ksemaw::ksemaw(QWidget *parent) :
    QWidget(parent),
    ui(new Ui::ksemaw)
{
 //setupUi(this); // this sets up GUI
    ui->setupUi(this);

    const QByteArray value = qgetenv("USER");
    QString uName=QString::fromLocal8Bit(value);
    cout << "current user = " << uName.toStdString() <<endl;
    QString f2ck="/home/"+uName+"/Workspace/temp/defau.1.Spj";
    cout << "filetocheck = " << f2ck.toStdString() << endl;
    QStringList filetocheck(f2ck);
 
 QFileSystemWatcher *qfsw = new QFileSystemWatcher(filetocheck);
 
 // signals/slots mechanism in action
 connect( ui->pBloadPro,  SIGNAL( clicked() ), this, SLOT(LoadProject()));
 connect( ui->pBsavePro,  SIGNAL( clicked() ), this, SLOT(SaveProject()));
 connect( ui->pBloadFnk,  SIGNAL( clicked() ), this, SLOT(LoadFilenk()));
 connect( ui->pBclearFnk, SIGNAL( clicked() ), this, SLOT(ClrFnk()));
 connect( ui->pBsaveFnk,  SIGNAL( clicked() ), this, SLOT(SaveFnk()));
 connect( ui->pushButton_clearNKtemp, SIGNAL( clicked() ), this, SLOT(ClearTemp()));
 connect( ui->pBnk1, SIGNAL( clicked() ), this, SLOT(Setnk1()));
 connect( ui->pBnk2, SIGNAL( clicked() ), this, SLOT(Setnk2()));
 connect( ui->pBnk3, SIGNAL( clicked() ), this, SLOT(Setnk3()));
 connect( ui->pBnk4, SIGNAL( clicked() ), this, SLOT(Setnk4()));
 connect( ui->pBnk5, SIGNAL( clicked() ), this, SLOT(Setnk5()));
 connect( ui->pBnk6, SIGNAL( clicked() ), this, SLOT(Setnk6()));
 connect( ui->pBnk7, SIGNAL( clicked() ), this, SLOT(Setnk7()));
 connect( ui->pBnk8, SIGNAL( clicked() ), this, SLOT(Setnk8()));
 connect( ui->pBclearnk1,SIGNAL( clicked() ), this, SLOT(Clrnk1()));
 connect( ui->pBclearnk2,SIGNAL( clicked() ), this, SLOT(Clrnk2()));
 connect( ui->pBclearnk3,SIGNAL( clicked() ), this, SLOT(Clrnk3()));
 connect( ui->pBclearnk4,SIGNAL( clicked() ), this, SLOT(Clrnk4()));
 connect( ui->pBclearnk5,SIGNAL( clicked() ), this, SLOT(Clrnk5()));
 connect( ui->pBclearnk6,SIGNAL( clicked() ), this, SLOT(Clrnk6()));
 connect( ui->pBclearnk7,SIGNAL( clicked() ), this, SLOT(Clrnk7()));
 connect( ui->pBclearnk8,SIGNAL( clicked() ), this, SLOT(Clrnk8()));
 connect( ui->pBsetSample,SIGNAL( clicked() ), this, SLOT(setSample()));
 connect( ui->lineEdit_sample,SIGNAL (textChanged(const QString &)),
		       this, SLOT(listMeas(const QString &)));
 connect( ui->cBmis1,SIGNAL (currentIndexChanged (const int &) ),
		       this, SLOT(pwTn(const int &)));
 connect( ui->cBmis2,SIGNAL (currentIndexChanged (const int &) ),
		       this, SLOT(pwTp(const int &)));
 connect( ui->cBmis3,SIGNAL (currentIndexChanged (const int &) ),
		       this, SLOT(pwRn(const int &)));
 connect( ui->cBmis4,SIGNAL (currentIndexChanged (const int &) ),
		       this, SLOT(pwRp(const int &)));
 connect( ui->cBmis5,SIGNAL (currentIndexChanged (const int &) ),
		       this, SLOT(pwR1(const int &)));
 connect( ui->cBmis6,SIGNAL (currentIndexChanged (const int &) ),
               this, SLOT(pwApds(const int &)));
 connect( ui->cBmis7,SIGNAL (currentIndexChanged (const int &) ),
		       this, SLOT(pwE1(const int &)));
 connect( ui->cBmis9,SIGNAL (currentIndexChanged (const int &) ),
		       this, SLOT(pwE2(const int &)));
 connect( ui->cBmis11,SIGNAL (currentIndexChanged (const int &) ),
		       this, SLOT(pwE3(const int &)));
 connect( ui->cBmis13,SIGNAL (currentIndexChanged (const int &) ),
		       this, SLOT(pwE4(const int &)));
 connect( ui->cBteE1,SIGNAL (currentIndexChanged (const int &) ),
		       this, SLOT(pwSubE1(const int &)));
 connect( ui->cBteE2,SIGNAL (currentIndexChanged (const int &) ),
		       this, SLOT(pwSubE2(const int &)));
 connect( ui->cBteE3,SIGNAL (currentIndexChanged (const int &) ),
		       this, SLOT(pwSubE3(const int &)));
 connect( ui->cBteE4,SIGNAL (currentIndexChanged (const int &) ),
		       this, SLOT(pwSubE4(const int &)));
 connect(ui->checkB_PAR_22_1,SIGNAL(stateChanged (const int &) ),this, SLOT(setRifMir()));
 connect(ui->checkB_mis1_1,SIGNAL(stateChanged (const int &) ),this, SLOT(MCRange()));
 connect(ui->checkB_mis2_1,SIGNAL(stateChanged (const int &) ),this, SLOT(MCRange()));
 connect(ui->checkB_mis3_1,SIGNAL(stateChanged (const int &) ),this, SLOT(MCRange()));
 connect(ui->checkB_mis4_1,SIGNAL(stateChanged (const int &) ),this, SLOT(MCRange()));
 connect(ui->checkB_mis5_1,SIGNAL(stateChanged (const int &) ),this, SLOT(MCRange()));
 connect(ui->checkB_mis6_1,SIGNAL(stateChanged (const int &) ),this, SLOT(MCRange()));
 connect(ui->checkB_mis7_1,SIGNAL(stateChanged (const int &) ),this, SLOT(MCRange()));
 connect(ui->checkB_mis9_1,SIGNAL(stateChanged (const int &) ),this, SLOT(MCRange()));
 connect(ui->checkB_mis11_1,SIGNAL(stateChanged (const int &) ),this, SLOT(MCRange()));
 connect(ui->checkB_mis13_1,SIGNAL(stateChanged (const int &) ),this, SLOT(MCRange()));
 connect(ui->cB_cnk1a,SIGNAL (currentIndexChanged (const int &) ),this, SLOT(setMat1()));
 connect(ui->cB_cnk2a,SIGNAL (currentIndexChanged (const int &) ),this, SLOT(setMat2()));
 connect(ui->cB_cnk3a,SIGNAL (currentIndexChanged (const int &) ),this, SLOT(setMat3()));
 connect(ui->cB_cnk4a,SIGNAL (currentIndexChanged (const int &) ),this, SLOT(setMat4()));
 connect(ui->cB_cnk5a,SIGNAL (currentIndexChanged (const int &) ),this, SLOT(setMat5()));
 connect(ui->cB_cnk6a,SIGNAL (currentIndexChanged (const int &) ),this, SLOT(setMat6()));
 connect(ui->cB_cnk7a,SIGNAL (currentIndexChanged (const int &) ),this, SLOT(setMat7()));
 connect(ui->cB_cnk8a,SIGNAL (currentIndexChanged (const int &) ),this, SLOT(setMat8()));
 connect(ui->cB_cnk9a,SIGNAL (currentIndexChanged (const int &) ),this, SLOT(setMat9()));
 connect(ui->cB_cnk10a,SIGNAL (currentIndexChanged (const int &) ),this, SLOT(setMat10()));
 connect(ui->cB_cnk11a,SIGNAL (currentIndexChanged (const int &) ),this, SLOT(setMat11()));
 connect(ui->cB_cnk12a,SIGNAL (currentIndexChanged (const int &) ),this, SLOT(setMat12()));
 connect(ui->cB_cnk13a,SIGNAL (currentIndexChanged (const int &) ),this, SLOT(setMat13()));
 connect(ui->cB_cnk14a,SIGNAL (currentIndexChanged (const int &) ),this, SLOT(setMat14()));
 connect(ui->cB_cnk15a,SIGNAL (currentIndexChanged (const int &) ),this, SLOT(setMat15()));
 connect(ui->cB_EMA_1,SIGNAL(stateChanged (const int &) ),this, SLOT(setEMA1()));
 connect(ui->cB_EMA_2,SIGNAL(stateChanged (const int &) ),this, SLOT(setEMA2()));
 connect(ui->cB_EMA_3,SIGNAL(stateChanged (const int &) ),this, SLOT(setEMA3()));
 connect(ui->cB_EMA_4,SIGNAL(stateChanged (const int &) ),this, SLOT(setEMA4()));
 connect(ui->cB_EMA_5,SIGNAL(stateChanged (const int &) ),this, SLOT(setEMA5()));
 connect(ui->cB_EMA_6,SIGNAL(stateChanged (const int &) ),this, SLOT(setEMA6()));
 connect(ui->cB_EMA_7,SIGNAL(stateChanged (const int &) ),this, SLOT(setEMA7()));
 connect(ui->cB_EMA_8,SIGNAL(stateChanged (const int &) ),this, SLOT(setEMA8()));
 connect(ui->cB_EMA_9,SIGNAL(stateChanged (const int &) ),this, SLOT(setEMA9()));
 connect(ui->cB_EMA_10,SIGNAL(stateChanged (const int &) ),this, SLOT(setEMA10()));
 connect(ui->cB_EMA_11,SIGNAL(stateChanged (const int &) ),this, SLOT(setEMA11()));
 connect(ui->cB_EMA_12,SIGNAL(stateChanged (const int &) ),this, SLOT(setEMA12()));
 connect(ui->cB_EMA_13,SIGNAL(stateChanged (const int &) ),this, SLOT(setEMA13()));
 connect(ui->cB_EMA_14,SIGNAL(stateChanged (const int &) ),this, SLOT(setEMA14()));
 connect(ui->cB_EMA_15,SIGNAL(stateChanged (const int &) ),this, SLOT(setEMA15()));
 connect(ui->cBosc_1,SIGNAL(stateChanged (const int &) ),this,SLOT(setOsc1()));
 connect(ui->cBosc_2,SIGNAL(stateChanged (const int &) ),this,SLOT(setOsc2()));
 connect(ui->cBosc_3,SIGNAL(stateChanged (const int &) ),this,SLOT(setOsc3()));
 connect(ui->cBosc_4,SIGNAL(stateChanged (const int &) ),this,SLOT(setOsc4()));
 connect(ui->cBosc_5,SIGNAL(stateChanged (const int &) ),this,SLOT(setOsc5()));
 connect(ui->cBosc_6,SIGNAL(stateChanged (const int &) ),this,SLOT(setOsc6()));
 connect(ui->cBosc_7,SIGNAL(stateChanged (const int &) ),this,SLOT(setOsc7()));
 connect(ui->cBosc_8,SIGNAL(stateChanged (const int &) ),this,SLOT(setOsc8()));
 connect(ui->cBosc_9,SIGNAL(stateChanged (const int &) ),this,SLOT(setOsc9()));
 connect(ui->cBosc_10,SIGNAL(stateChanged (const int &) ),this,SLOT(setOsc10()));
 connect(ui->cBosc_11,SIGNAL(stateChanged (const int &) ),this,SLOT(setOsc11()));
 connect(ui->cBosc_12,SIGNAL(stateChanged (const int &) ),this,SLOT(setOsc12()));
 connect(ui->cBosc_13,SIGNAL(stateChanged (const int &) ),this,SLOT(setOsc13()));
 connect(ui->cBosc_14,SIGNAL(stateChanged (const int &) ),this,SLOT(setOsc14()));
 connect(ui->cBosc_15,SIGNAL(stateChanged (const int &) ),this,SLOT(setOsc15()));
 connect(ui->cBosc_16,SIGNAL(stateChanged (const int &) ),this,SLOT(setOsc16()));
 connect(ui->cBosc_17,SIGNAL(stateChanged (const int &) ),this,SLOT(setOsc17()));
 connect(ui->cBosc_18,SIGNAL(stateChanged (const int &) ),this,SLOT(setOsc18()));
 connect(ui->cBosc_19,SIGNAL(stateChanged (const int &) ),this,SLOT(setOsc19()));
 connect(ui->cBosc_20,SIGNAL(stateChanged (const int &) ),this,SLOT(setOsc20()));
 connect(ui->tabWidget,SIGNAL (currentChanged(const int &)),this, SLOT(tabChanged()));
 connect(ui->dSB_PAR_4_1,SIGNAL(valueChanged(const double&)),this,SLOT(rangeWL()));
 connect(ui->dSB_PAR_4_2,SIGNAL(valueChanged(const double&)),this,SLOT(rangeWL()));
 connect(ui->dSB_PAR_6_1,SIGNAL(valueChanged(const double &)),this,SLOT(AdjTheta()));
 connect(ui->dSB_PAR_14_1,SIGNAL(valueChanged(const double &)),this,SLOT(AdjTheta()));
 connect(ui->dSB_PAR_15_1,SIGNAL(valueChanged(const double &)),this,SLOT(AdjTheta()));
 connect(ui->dSB_PAR_16_1,SIGNAL(valueChanged(const double &)),this,SLOT(AdjTheta()));
 connect(ui->dSB_PAR_17_1,SIGNAL(valueChanged(const double &)),this,SLOT(AdjTheta()));
 connect(ui->comB_PAR_51_3,SIGNAL(currentIndexChanged (const int &) ),this, SLOT(RefreshModel()));
 connect(ui->comB_PAR_52_3,SIGNAL(currentIndexChanged (const int &) ),this, SLOT(RefreshModel()));
 connect(ui->comB_PAR_53_3,SIGNAL(currentIndexChanged (const int &) ),this, SLOT(RefreshModel()));
 connect(ui->comB_PAR_54_3,SIGNAL(currentIndexChanged (const int &) ),this, SLOT(RefreshModel()));
 connect(ui->comB_PAR_55_3,SIGNAL(currentIndexChanged (const int &) ),this, SLOT(RefreshModel()));
 connect(ui->comB_PAR_56_3,SIGNAL(currentIndexChanged (const int &) ),this, SLOT(RefreshModel()));
 connect(ui->comB_PAR_57_3,SIGNAL(currentIndexChanged (const int &) ),this, SLOT(RefreshModel()));
 connect(ui->comB_PAR_58_3,SIGNAL(currentIndexChanged (const int &) ),this, SLOT(RefreshModel()));
 connect(ui->comB_PAR_59_3,SIGNAL(currentIndexChanged (const int &) ),this, SLOT(RefreshModel()));
 connect(ui->sB_PAR_51_2,SIGNAL(valueChanged(const int &)),this,SLOT(SetModel(const int &)));
 connect(ui->sB_PAR_34_5,SIGNAL(valueChanged(const int &)),this,SLOT(PanFitPar()));
 connect(ui->chBeParFit_1,SIGNAL(stateChanged (const int &) ),this, SLOT(PanFitEnable()));
 connect(ui->chBeParFit_2,SIGNAL(stateChanged (const int &) ),this, SLOT(PanFitEnable()));
 connect(ui->chBeParFit_3,SIGNAL(stateChanged (const int &) ),this, SLOT(PanFitEnable()));
 connect(ui->chBeParFit_4,SIGNAL(stateChanged (const int &) ),this, SLOT(PanFitEnable()));
 connect(ui->chBeParFit_5,SIGNAL(stateChanged (const int &) ),this, SLOT(PanFitEnable()));
 connect(ui->chBeParFit_6,SIGNAL(stateChanged (const int &) ),this, SLOT(PanFitEnable()));
 connect(ui->chBeParFit_7,SIGNAL(stateChanged (const int &) ),this, SLOT(PanFitEnable()));
 connect(ui->chBeParFit_8,SIGNAL(stateChanged (const int &) ),this, SLOT(PanFitEnable()));
 connect(ui->chBeParFit_9,SIGNAL(stateChanged (const int &) ),this, SLOT(PanFitEnable()));
 connect(ui->chBeParFit_10,SIGNAL(stateChanged (const int &) ),this, SLOT(PanFitEnable()));
 connect(ui->chBeParFit_11,SIGNAL(stateChanged (const int &) ),this, SLOT(PanFitEnable()));
 connect(ui->chBeParFit_12,SIGNAL(stateChanged (const int &) ),this, SLOT(PanFitEnable()));
 connect(ui->chBeParFit_13,SIGNAL(stateChanged (const int &) ),this, SLOT(PanFitEnable()));
 connect(ui->chBeParFit_14,SIGNAL(stateChanged (const int &) ),this, SLOT(PanFitEnable()));
 connect(ui->chBeParFit_15,SIGNAL(stateChanged (const int &) ),this, SLOT(PanFitEnable()));
 connect(ui->chBeParFit_16,SIGNAL(stateChanged (const int &) ),this, SLOT(PanFitEnable()));
 connect(ui->chBeParFit_17,SIGNAL(stateChanged (const int &) ),this, SLOT(PanFitEnable()));
 connect(ui->cBParFit_1,SIGNAL(currentIndexChanged (const int &) ), this, SLOT(PanFitChoice()));
 connect(ui->cBParFit_2,SIGNAL(currentIndexChanged (const int &) ), this, SLOT(PanFitChoice()));
 connect(ui->cBParFit_3,SIGNAL(currentIndexChanged (const int &) ), this, SLOT(PanFitChoice()));
 connect(ui->cBParFit_4,SIGNAL(currentIndexChanged (const int &) ), this, SLOT(PanFitChoice()));
 connect(ui->cBParFit_5,SIGNAL(currentIndexChanged (const int &) ), this, SLOT(PanFitChoice()));
 connect(ui->cBParFit_6,SIGNAL(currentIndexChanged (const int &) ), this, SLOT(PanFitChoice()));
 connect(ui->cBParFit_7,SIGNAL(currentIndexChanged (const int &) ), this, SLOT(PanFitChoice()));
 connect(ui->cBParFit_8,SIGNAL(currentIndexChanged (const int &) ), this, SLOT(PanFitChoice()));
 connect(ui->cBParFit_9,SIGNAL(currentIndexChanged (const int &) ), this, SLOT(PanFitChoice()));
 connect(ui->cBParFit_10,SIGNAL(currentIndexChanged (const int &) ), this, SLOT(PanFitChoice()));
 connect(ui->cBParFit_11,SIGNAL(currentIndexChanged (const int &) ), this, SLOT(PanFitChoice()));
 connect(ui->cBParFit_12,SIGNAL(currentIndexChanged (const int &) ), this, SLOT(PanFitChoice()));
 connect(ui->cBParFit_13,SIGNAL(currentIndexChanged (const int &) ), this, SLOT(PanFitChoice()));
 connect(ui->cBParFit_14,SIGNAL(currentIndexChanged (const int &) ), this, SLOT(PanFitChoice()));
 connect(ui->cBParFit_15,SIGNAL(currentIndexChanged (const int &) ), this, SLOT(PanFitChoice()));
 connect(ui->cBParFit_16,SIGNAL(currentIndexChanged (const int &) ), this, SLOT(PanFitChoice()));
 connect(ui->cBParFit_17,SIGNAL(currentIndexChanged (const int &) ), this, SLOT(PanFitChoice()));
 connect(qfsw, SIGNAL(fileChanged(const QString &)),this, SLOT( fileChanged(const QString &)));
 connect(ui->pushButton_PlotExpMis,SIGNAL(clicked()), this, SLOT(PlotME()));
 connect(ui->pushButton_Simulate,SIGNAL(clicked()), this, SLOT(Simula()));
 connect(ui->pushButton_PlotAve,SIGNAL(clicked()), this, SLOT(PlotAve()));
 connect(ui->pushButton_cursor,SIGNAL(clicked()), this, SLOT(readCursor()));
 connect(ui->pushButton_SearchNK,SIGNAL(clicked()), this, SLOT(searchNK()));
 connect(ui->pushButton_PickSaveNK,SIGNAL(clicked()), this, SLOT(save1NK()));
 connect(ui->pushButton_RefreshGraph,SIGNAL(clicked()), this, SLOT(RefTrackG()));
 connect(ui->pushButton_AutosFromWLmin,SIGNAL(clicked()), this, SLOT(TrackFromWLmin()));
 connect(ui->pushButton_AutosFromWLmax,SIGNAL(clicked()), this, SLOT(TrackFromWLmax()));
 connect(ui->pushButton_StoreNK,SIGNAL(clicked()), this, SLOT(saveNKtmp()));
 connect(ui->pushButton_RefreshGraph_2,SIGNAL(clicked()), this, SLOT(RefIbridG()));
 connect(ui->pushButton_PlotCurrentFit,SIGNAL(clicked()), this, SLOT(IbridPlotFit()));
 connect(ui->pushButton_FitN,SIGNAL(clicked()), this, SLOT(FitN()));
 connect(ui->pushButton_FitNK,SIGNAL(clicked()), this, SLOT(FitNK()));
 connect(ui->pushButton_FitE1E2,SIGNAL(clicked()), this, SLOT(FitE1E2()));
 connect(ui->pushButton_IbridOne,SIGNAL(clicked()), this, SLOT(IbridOne()));
 connect(ui->pushButton_IbridOneErrStore,SIGNAL(clicked()), this, SLOT(IbridOneStore()));
 connect(ui->pushButton_clearNKtempIbri, SIGNAL( clicked() ), this, SLOT(ClearTempIbri()));
 connect(ui->comboBox_StoreNK,SIGNAL(currentIndexChanged (const int &) ), this, SLOT(AdjSave()));
 connect(ui->comboBox_store,SIGNAL(currentIndexChanged (const int &) ), this, SLOT(AdjSave()));
 connect(ui->pushButton_EnableAllWL,SIGNAL(clicked()), this, SLOT(EnableAllWl()));
 connect(ui->pushButton_enableAllNK,SIGNAL(clicked()), this, SLOT(EnableAllNK()));
 connect(ui->pushButton_EditWL,SIGNAL(clicked()), this, SLOT(EditWl()));
 connect(ui->pushButton_EditNK,SIGNAL(clicked()), this, SLOT(EditNK()));
 connect(ui->pushButton_BackToBestSituation,SIGNAL(clicked()), this,SLOT(GoBest()));
 connect(ui->pushButton_previous,SIGNAL(clicked()), this,SLOT(GoPrevious()));
 connect(ui->pushButton_next,SIGNAL(clicked()), this,SLOT(GoNext()));
 connect(ui->pBmDw1, SIGNAL( clicked() ), this, SLOT(mDw1()));
 connect(ui->pBmUp1, SIGNAL( clicked() ), this, SLOT(mUp1()));
 connect(ui->pBmDw2, SIGNAL( clicked() ), this, SLOT(mDw2()));
 connect(ui->pBmUp2, SIGNAL( clicked() ), this, SLOT(mUp2()));
 connect(ui->pBmDw3, SIGNAL( clicked() ), this, SLOT(mDw3()));
 connect(ui->pBmUp3, SIGNAL( clicked() ), this, SLOT(mUp3()));
 connect(ui->pBmDw4, SIGNAL( clicked() ), this, SLOT(mDw4()));
 connect(ui->pBmUp4, SIGNAL( clicked() ), this, SLOT(mUp4()));
 connect(ui->pBmDw5, SIGNAL( clicked() ), this, SLOT(mDw5()));
 connect(ui->pBmUp5, SIGNAL( clicked() ), this, SLOT(mUp5()));
 connect(ui->pBmDw6, SIGNAL( clicked() ), this, SLOT(mDw6()));
 connect(ui->pBmUp6, SIGNAL( clicked() ), this, SLOT(mUp6()));
 connect(ui->pBmDw7, SIGNAL( clicked() ), this, SLOT(mDw7()));
 connect(ui->pBmUp7, SIGNAL( clicked() ), this, SLOT(mUp7()));
 connect(ui->pBmDw8, SIGNAL( clicked() ), this, SLOT(mDw8()));
 connect(ui->pBmUp8, SIGNAL( clicked() ), this, SLOT(mUp8()));
 connect(ui->pBmDw9, SIGNAL( clicked() ), this, SLOT(mDw9()));
 connect(ui->pBmUp9, SIGNAL( clicked() ), this, SLOT(mUp9()));
 connect(ui->pushButton_saveSim, SIGNAL( clicked() ), this, SLOT(saveSim()));
 connect(ui->pushButton_PlotAbsEL,SIGNAL( clicked() ),this, SLOT(PlotAbsEL()));
 connect(ui->comboBox_searchNK,SIGNAL(currentIndexChanged (const int &) ), this, SLOT(manageLEwl()));
  
 //inizializzazione parametri  
 pathroot="/home/"+uName+"/Workspace/";
 fileStore=pathroot+"temp/defau.1.Spj";
 //fctrl=pathroot+"Qt4-prove/ksemaw/ksemaw.ctrl";
 //QString fRefMir=pathroot+"Qt4-prove/ksemaw/referenceMirrors.txt";
 //fStdSpect=pathroot+"Qt4-prove/ksemaw/standardSpectra.txt";
 fctrl=pathroot+"qtSource/ksemaw/ksemaw.ctrl";
 QString fRefMir=pathroot+"qtSource/ksemaw/referenceMirrors.txt";
 fStdSpect=pathroot+"qtSource/ksemaw/standardSpectra.txt";

 lastIndex=0;
 lastTab=0;
 ifn=0;
 occupyPF=0;
 nlayer=0;
 npp=0;
 iwspj=0;
 for(int i=1;i<=5;i++){
  for(int j=1;j<=60;j++){
   par[j][i]=0.;
  }
 } 
 par[1][3]= 1.;    //efficienza raccolta 1a faccia
 par[2][3]= 1.;    //efficienza raccolta 2a faccia
 par[3][3]=.0022;  //fluttuazione back correction SF
 par[4][1]=2000.;  //WLmin
 par[4][2]=30000.; //WLmax
 par[4][3]=.005;   //DRcal/Rcal=errore su R specchio di calibrazione
 par[5][1]=1.0;    //soglia accettazione soluzione (AUTOS)
 par[5][3]=.0005;  //Errore di lettura
 par[6][3]=1.0;    // ... dipendente da wl
 par[7][2]=1;       //icolor
 par[7][3]=1.2505E-7; //cte per calcolo contributo SUB al segnale PDS
 par[8][1]=0.;     //jobtot
 par[8][2]=0.;     //jobview
 par[10][3]=3;     //RIF06 dopo  del 5/dic/94
 par[22][1]=0;     //mr=0 => R assolute
 par[22][2]=1;     //DF(1)
 par[23][1]=-1;    //DF(2)
 par[23][2]=-1;    //DF(3)
 par[21][3]=.001;  //incremento n
 par[22][3]=1.e-4; //incremento k
 par[23][3]=.001;  //incremento n_2
 par[24][3]=.0;    //incremento k_2
 par[25][3]=.0;    //errore n_sub
 par[26][3]=.25;   //distinguibilit soluzioni (autos)
 par[26][1]=1.;    //iwl_MIN (autos)
 par[26][2]=201.;  //iwl_MAX (autos)
 par[27][2]=1.;    //S1P2=1
 par[28][1]=3.;    //N per il calcolo di <f(x)> 
 par[29][1]=23.;   //N. strati discretizzazione film inomogeneo
 par[33][1]=-10;   //procedura ricerca 3a incognita
 par[34][1]=-10;   //procedura ricerca 4a incognita
 par[35][1]=1.;    //imed     (cnfr)
 par[35][2]=3.;    //s1p2u3   (cnfr)
 par[35][4]=0.;    //in SIMULA cnk(1][1]=0 -> nk=cte
// impostazione puntatore di PM: PPM(17)
 par[34][5]=7.;    //pannello di fit a 7 parametri
 par[35][5]=1.;    //solo 1 parametro fittato
 par[36][5]=102;  //C1
 par[37][5]=103;  //E1
 par[38][5]=104;   //D1
 par[39][5]=107;   //C2
 par[40][5]=108;   //E2
 par[41][5]=109;   //D2
 par[42][5]=112;   //C3
// impostazione film
 for(int i=51;i<=60;i++){
  par[i][3]=1.; //strati "bulk"
 }
// matrice costruzione VNK
 for(int i=1;i<=15;i++){
  cnk[i][1]=.0;
  cnk[i][2]=1.;
  cnk[i][3]=.0;
 }
 cnk[10][1]=0.;  //mezzo di ingresso = aria
 cnk[10][2]=1.;
 cnk[10][3]=0.;
 cnk[16][1]=0.; //mezzo di uscita = aria
 cnk[16][2]=1.;
 cnk[16][3]=0.;
 // parametri modello
 for(int i=1;i<=200;i++){
  for(int j=1;j<=5;j++){
   pm[i][j]=0.;
  }
 }
 pm[1][3]=112;    //fit di C3
 pm[112][2]=1.;   //abilito C3 al fit
 pm[100][1]=1;    //FT# 1 di lavoro in ibridone
 pm[101][1]=2.;   //oscillatore quantistico omogeneo in UV 
 pm[102][1]=0.;   //C1    "          "
 pm[103][1]=5.0;  //E1    "          "
 pm[104][1]=0.04; //D1    "          "
 pm[105][1]=1.;   //K1    "          "
 pm[106][1]=2.;   //oscillatore quantistico omogeneo in IR
 pm[107][1]=0.;   //C2
 pm[108][1]=0.3;  //E2    "          "
 pm[109][1]=0.005;//D2    "          "
 pm[110][1]=1.;   //K2    "          "
 pm[111][1]=4.;   //oscillatore FLAT 
 pm[112][1]=1.5;  //C3    "          "
 pm[113][1]=0.;   //E3    "          "
 pm[114][1]=0.;   //D3    "          "
 pm[115][1]=0.;   //K3    "          "      
// matrice opzioni fit
 for(int i=1;i<=21;i++){
     for(int j=1;j<=7;j++){
         pf[j][i]=0.;
     }
 }
 pf[1][1]=3.;// 3 oscillatori
 pf[1][2]=1.;// osc#1
 pf[1][3]=2.;// osc#2
 pf[1][4]=3.;// osc#3
//grafica
 for(int i=1;i<=6;i++){
  rxy[i][1]=.0;
  rxy[i][2]=100.;
  rxy[i][3]=.0;
  rxy[i][4]=.0;
 }
 for(int i=7;i<=13;i+=2 ){
  rxy[i][1]=.0;
  rxy[i][2]=180.;
  rxy[i][3]=.0;
  rxy[i][4]=.0;
 }
 for(int i=8;i<=14;i+=2 ){
  rxy[i][1]=-180;
  rxy[i][2]=180.;
  rxy[i][3]=.0;
  rxy[i][4]=.0;
 }
 rxy[15][1]=0.;
 rxy[15][2]=10.;
 rxy[15][3]=0.;
 rxy[15][4]=10.;
 for(int i=16;i<=19;i++){
  rxy[i][1]=.0;  
  rxy[i][2]=1.0;
  rxy[i][3]=.0;
  rxy[i][4]=1.0;
 }
 rxy[18][2]=1.;
 rxy[20][1]=2000.;
 rxy[20][2]=30000.;
 rxy[20][3]=2000.;
 rxy[20][4]=30000.;
 rxy[21][1]=0.;
 rxy[21][2]=90.;
 rxy[21][3]=.0;
 rxy[21][4]=90.;
 for(int i=22;i<=30;i++){
  rxy[i][1]=.0;  
  rxy[i][2]=.0;
  rxy[i][3]=.0;
  rxy[i][4]=.0;
 }
  rxy[24][1]=9.;//graph width
  rxy[24][2]=0.5;//graph ratio
  rxy[24][3]=3.;//line width
  rxy[25][3]=2.;//plot in eV
  rxy[25][4]=2.;//spline in eV 
  for(int i=26;i<=27;i++){
   rxy[i][1]=.0;
   rxy[i][2]=10.;
   rxy[i][3]=0.;
   rxy[i][4]=10.;
  }
  for(int i=1;i<=25;i++){//setting values in GUI
    for(int j=1;j<=4;j++){
    if(idToLineEdit.contains("DP_RXY_"+QString::number(i)+"_"+QString::number(j)))
        idToLineEdit["DP_RXY_"+QString::number(i)+"_"+QString::number(j)]
        -> setText(QString::number(rxy[i][j]));
    }
 }
// fine inizializzazione
 
 
 idToComboBox["cBParFit_1"] = ui->cBParFit_1;
 idToComboBox["cBParFit_2"] = ui->cBParFit_2;
 idToComboBox["cBParFit_3"] = ui->cBParFit_3;
 idToComboBox["cBParFit_4"] = ui->cBParFit_4;
 idToComboBox["cBParFit_5"] = ui->cBParFit_5;
 idToComboBox["cBParFit_6"] = ui->cBParFit_6;
 idToComboBox["cBParFit_7"] = ui->cBParFit_7;
 idToComboBox["cBParFit_8"] = ui->cBParFit_8;
 idToComboBox["cBParFit_9"] = ui->cBParFit_9;
 idToComboBox["cBParFit_10"] = ui->cBParFit_10;
 idToComboBox["cBParFit_11"] = ui->cBParFit_11;
 idToComboBox["cBParFit_12"] = ui->cBParFit_12;
 idToComboBox["cBParFit_13"] = ui->cBParFit_13;
 idToComboBox["cBParFit_14"] = ui->cBParFit_14;
 idToComboBox["cBParFit_15"] = ui->cBParFit_15;
 idToComboBox["cBParFit_16"] = ui->cBParFit_16;
 idToComboBox["cBParFit_17"] = ui->cBParFit_17;
 idToComboBox["cB_cnk1a"] = ui->cB_cnk1a;
 idToComboBox["cB_cnk1b"] = ui->cB_cnk1b;
 idToComboBox["cB_cnk2a"] = ui->cB_cnk2a;
 idToComboBox["cB_cnk2b"] = ui->cB_cnk2b;
 idToComboBox["cB_cnk3a"] = ui->cB_cnk3a;
 idToComboBox["cB_cnk3b"] = ui->cB_cnk3b;
 idToComboBox["cB_cnk4a"] = ui->cB_cnk4a;
 idToComboBox["cB_cnk4b"] = ui->cB_cnk4b;
 idToComboBox["cB_cnk5a"] = ui->cB_cnk5a;
 idToComboBox["cB_cnk5b"] = ui->cB_cnk5b;
 idToComboBox["cB_cnk6a"] = ui->cB_cnk6a;
 idToComboBox["cB_cnk6b"] = ui->cB_cnk6b;
 idToComboBox["cB_cnk7a"] = ui->cB_cnk7a;
 idToComboBox["cB_cnk7b"] = ui->cB_cnk7b;
 idToComboBox["cB_cnk8a"] = ui->cB_cnk8a;
 idToComboBox["cB_cnk8b"] = ui->cB_cnk8b;
 idToComboBox["cB_cnk9a"] = ui->cB_cnk9a;
 idToComboBox["cB_cnk9b"] = ui->cB_cnk9b;
 idToComboBox["cB_cnk10a"] = ui->cB_cnk10a;
 idToComboBox["cB_cnk10b"] = ui->cB_cnk10b;
 idToComboBox["cB_cnk11a"] = ui->cB_cnk11a;
 idToComboBox["cB_cnk11b"] = ui->cB_cnk11b;
 idToComboBox["cB_cnk12a"] = ui->cB_cnk12a;
 idToComboBox["cB_cnk12b"] = ui->cB_cnk12b;
 idToComboBox["cB_cnk13a"] = ui->cB_cnk13a;
 idToComboBox["cB_cnk13b"] = ui->cB_cnk13b;
 idToComboBox["cB_cnk14a"] = ui->cB_cnk14a;
 idToComboBox["cB_cnk14b"] = ui->cB_cnk14b;
 idToComboBox["cB_cnk15a"] = ui->cB_cnk15a;
 idToComboBox["cB_cnk15b"] = ui->cB_cnk15b;
 idToComboBox["cBmis1"] = ui->cBmis1;
 idToComboBox["cBmis2"] = ui->cBmis2;
 idToComboBox["cBmis3"] = ui->cBmis3;
 idToComboBox["cBmis4"] = ui->cBmis4;
 idToComboBox["cBmis5"] = ui->cBmis5;
 idToComboBox["cBmis6"] = ui->cBmis6;
 idToComboBox["cBmis7"] = ui->cBmis7;
 idToComboBox["cBmis9"] = ui->cBmis9;
 idToComboBox["cBmis11"] = ui->cBmis11;
 idToComboBox["cBmis13"] = ui->cBmis13;
 idToComboBox["cBteE1"] = ui->cBteE1;
 idToComboBox["cBteE2"] = ui->cBteE2;
 idToComboBox["cBteE3"] = ui->cBteE3;
 idToComboBox["cBteE4"] = ui->cBteE4;
 idToComboBox["comB_PAR_51_1"] = ui->comB_PAR_51_1;
 idToComboBox["comB_PAR_52_1"] = ui->comB_PAR_52_1;
 idToComboBox["comB_PAR_53_1"] = ui->comB_PAR_53_1;
 idToComboBox["comB_PAR_54_1"] = ui->comB_PAR_54_1;
 idToComboBox["comB_PAR_55_1"] = ui->comB_PAR_55_1;
 idToComboBox["comB_PAR_56_1"] = ui->comB_PAR_56_1;
 idToComboBox["comB_PAR_57_1"] = ui->comB_PAR_57_1;
 idToComboBox["comB_PAR_58_1"] = ui->comB_PAR_58_1;
 idToComboBox["comB_PAR_59_1"] = ui->comB_PAR_59_1;
 idToComboBox["comB_PAR_51_3"] = ui->comB_PAR_51_3;
 idToComboBox["comB_PAR_52_3"] = ui->comB_PAR_52_3;
 idToComboBox["comB_PAR_53_3"] = ui->comB_PAR_53_3;
 idToComboBox["comB_PAR_54_3"] = ui->comB_PAR_54_3;
 idToComboBox["comB_PAR_55_3"] = ui->comB_PAR_55_3;
 idToComboBox["comB_PAR_56_3"] = ui->comB_PAR_56_3;
 idToComboBox["comB_PAR_57_3"] = ui->comB_PAR_57_3;
 idToComboBox["comB_PAR_58_3"] = ui->comB_PAR_58_3;
 idToComboBox["comB_PAR_59_3"] = ui->comB_PAR_59_3;
 idToComboBox["cBpm_101_1"] = ui->cBpm_101_1;
 idToComboBox["cBpm_106_1"] = ui->cBpm_106_1;
 idToComboBox["cBpm_111_1"] = ui->cBpm_111_1;
 idToComboBox["cBpm_116_1"] = ui->cBpm_116_1;
 idToComboBox["cBpm_121_1"] = ui->cBpm_121_1;
 idToComboBox["cBpm_126_1"] = ui->cBpm_126_1;
 idToComboBox["cBpm_131_1"] = ui->cBpm_131_1;
 idToComboBox["cBpm_136_1"] = ui->cBpm_136_1;
 idToComboBox["cBpm_141_1"] = ui->cBpm_141_1;
 idToComboBox["cBpm_146_1"] = ui->cBpm_146_1;
 idToComboBox["cBpm_151_1"] = ui->cBpm_151_1;
 idToComboBox["cBpm_156_1"] = ui->cBpm_156_1;
 idToComboBox["cBpm_161_1"] = ui->cBpm_161_1;
 idToComboBox["cBpm_166_1"] = ui->cBpm_166_1;
 idToComboBox["cBpm_171_1"] = ui->cBpm_171_1;
 idToComboBox["cBpm_176_1"] = ui->cBpm_176_1;
 idToComboBox["cBpm_181_1"] = ui->cBpm_181_1;
 idToComboBox["cBpm_186_1"] = ui->cBpm_186_1;
 idToComboBox["cBpm_191_1"] = ui->cBpm_191_1;
 idToComboBox["cBpm_196_1"] = ui->cBpm_196_1;
 
 idToCheckBox["cBosc_1"] = ui->cBosc_1;
 idToCheckBox["cBosc_2"] = ui->cBosc_2;
 idToCheckBox["cBosc_3"] = ui->cBosc_3;
 idToCheckBox["cBosc_4"] = ui->cBosc_4;
 idToCheckBox["cBosc_5"] = ui->cBosc_5;
 idToCheckBox["cBosc_6"] = ui->cBosc_6;
 idToCheckBox["cBosc_7"] = ui->cBosc_7;
 idToCheckBox["cBosc_8"] = ui->cBosc_8;
 idToCheckBox["cBosc_9"] = ui->cBosc_9;
 idToCheckBox["cBosc_10"] = ui->cBosc_10;
 idToCheckBox["cBosc_11"] = ui->cBosc_11;
 idToCheckBox["cBosc_12"] = ui->cBosc_12;
 idToCheckBox["cBosc_13"] = ui->cBosc_13;
 idToCheckBox["cBosc_14"] = ui->cBosc_14;
 idToCheckBox["cBosc_15"] = ui->cBosc_15;
 idToCheckBox["cBosc_16"] = ui->cBosc_16;
 idToCheckBox["cBosc_17"] = ui->cBosc_17;
 idToCheckBox["cBosc_18"] = ui->cBosc_18;
 idToCheckBox["cBosc_19"] = ui->cBosc_19;
 idToCheckBox["cBosc_20"] = ui->cBosc_20;
 idToCheckBox["cB_EMA_1"] = ui->cB_EMA_1;
 idToCheckBox["cB_EMA_2"] = ui->cB_EMA_2;
 idToCheckBox["cB_EMA_3"] = ui->cB_EMA_3;
 idToCheckBox["cB_EMA_4"] = ui->cB_EMA_4;
 idToCheckBox["cB_EMA_5"] = ui->cB_EMA_5;
 idToCheckBox["cB_EMA_6"] = ui->cB_EMA_6;
 idToCheckBox["cB_EMA_7"] = ui->cB_EMA_7;
 idToCheckBox["cB_EMA_8"] = ui->cB_EMA_8;
 idToCheckBox["cB_EMA_9"] = ui->cB_EMA_9;
 idToCheckBox["cB_EMA_10"] = ui->cB_EMA_10;
 idToCheckBox["cB_EMA_11"] = ui->cB_EMA_11;
 idToCheckBox["cB_EMA_12"] = ui->cB_EMA_12;
 idToCheckBox["cB_EMA_13"] = ui->cB_EMA_13;
 idToCheckBox["cB_EMA_14"] = ui->cB_EMA_14;
 idToCheckBox["cB_EMA_15"] = ui->cB_EMA_15;
 idToCheckBox["chBeParFit_1"] = ui->chBeParFit_1;
 idToCheckBox["chBeParFit_2"] = ui->chBeParFit_2;
 idToCheckBox["chBeParFit_3"] = ui->chBeParFit_3;
 idToCheckBox["chBeParFit_4"] = ui->chBeParFit_4;
 idToCheckBox["chBeParFit_5"] = ui->chBeParFit_5;
 idToCheckBox["chBeParFit_6"] = ui->chBeParFit_6;
 idToCheckBox["chBeParFit_7"] = ui->chBeParFit_7;
 idToCheckBox["chBeParFit_8"] = ui->chBeParFit_8;
 idToCheckBox["chBeParFit_9"] = ui->chBeParFit_9;
 idToCheckBox["chBeParFit_10"] = ui->chBeParFit_10;
 idToCheckBox["chBeParFit_11"] = ui->chBeParFit_11;
 idToCheckBox["chBeParFit_12"] = ui->chBeParFit_12;
 idToCheckBox["chBeParFit_13"] = ui->chBeParFit_13;
 idToCheckBox["chBeParFit_14"] = ui->chBeParFit_14;
 idToCheckBox["chBeParFit_15"] = ui->chBeParFit_15;
 idToCheckBox["chBeParFit_16"] = ui->chBeParFit_16;
 idToCheckBox["chBeParFit_17"] = ui->chBeParFit_17;
 idToCheckBox["checkB_mis1_1"] = ui->checkB_mis1_1;
 idToCheckBox["checkB_mis2_1"] = ui->checkB_mis2_1;
 idToCheckBox["checkB_mis3_1"] = ui->checkB_mis3_1;
 idToCheckBox["checkB_mis4_1"] = ui->checkB_mis4_1;
 idToCheckBox["checkB_mis5_1"] = ui->checkB_mis5_1;
 idToCheckBox["checkB_mis6_1"] = ui->checkB_mis6_1;
 idToCheckBox["checkB_mis7_1"] = ui->checkB_mis7_1;
 idToCheckBox["checkB_mis9_1"] = ui->checkB_mis9_1;
 idToCheckBox["checkB_mis11_1"] = ui->checkB_mis11_1;
 idToCheckBox["checkB_mis13_1"] = ui->checkB_mis13_1;
 idToCheckBox["checkB_mis1_2"] = ui->checkB_mis1_2;
 idToCheckBox["checkB_mis2_2"] = ui->checkB_mis2_2;
 idToCheckBox["checkB_mis3_2"] = ui->checkB_mis3_2;
 idToCheckBox["checkB_mis4_2"] = ui->checkB_mis4_2;
 idToCheckBox["checkB_mis5_2"] = ui->checkB_mis5_2;
 idToCheckBox["checkB_mis6_2"] = ui->checkB_mis6_2;
 idToCheckBox["checkB_mis7_2"] = ui->checkB_mis7_2;
 idToCheckBox["checkB_mis8_2"] = ui->checkB_mis8_2;
 idToCheckBox["checkB_mis9_2"] = ui->checkB_mis9_2;
 idToCheckBox["checkB_mis10_2"] = ui->checkB_mis10_2;
 idToCheckBox["checkB_mis11_2"] = ui->checkB_mis11_2;
 idToCheckBox["checkB_mis12_2"] = ui->checkB_mis12_2;
 idToCheckBox["checkB_mis13_2"] = ui->checkB_mis13_2;
 idToCheckBox["checkB_mis14_2"] = ui->checkB_mis14_2;
 idToCheckBox["checkB_mis1_3"] = ui->checkB_mis1_3;
 idToCheckBox["checkB_mis2_3"] = ui->checkB_mis2_3;
 idToCheckBox["checkB_mis3_3"] = ui->checkB_mis3_3;
 idToCheckBox["checkB_mis4_3"] = ui->checkB_mis4_3;
 idToCheckBox["checkB_mis5_3"] = ui->checkB_mis5_3;
 idToCheckBox["checkB_mis6_3"] = ui->checkB_mis6_3;
 idToCheckBox["checkB_mis7_3"] = ui->checkB_mis7_3;
 idToCheckBox["checkB_mis8_3"] = ui->checkB_mis8_3;
 idToCheckBox["checkB_mis9_3"] = ui->checkB_mis9_3;
 idToCheckBox["checkB_mis10_3"] = ui->checkB_mis10_3;
 idToCheckBox["checkB_mis11_3"] = ui->checkB_mis11_3;
 idToCheckBox["checkB_mis12_3"] = ui->checkB_mis12_3;
 idToCheckBox["checkB_mis13_3"] = ui->checkB_mis13_3;
 idToCheckBox["checkB_mis14_3"] = ui->checkB_mis14_3;
 
 idToLineEdit["lineEdit1"]=ui->lineEdit1;
 idToLineEdit["lineEdit2"]=ui->lineEdit2;
 idToLineEdit["lineEdit3"]=ui->lineEdit3;
 idToLineEdit["lineEdit4"]=ui->lineEdit4;
 idToLineEdit["lineEdit5"]=ui->lineEdit5;
 idToLineEdit["lineEdit6"]=ui->lineEdit6;
 idToLineEdit["lineEdit7"]=ui->lineEdit7;
 idToLineEdit["lineEdit8"]=ui->lineEdit8;
 idToLineEdit["lineEdit_infoNK_1"]=ui->lineEdit_infoNK_1;
 idToLineEdit["lineEdit_infoNK_2"]=ui->lineEdit_infoNK_2;
 idToLineEdit["lineEdit_infoNK_3"]=ui->lineEdit_infoNK_3;
 idToLineEdit["lineEdit_infoNK_4"]=ui->lineEdit_infoNK_4;
 idToLineEdit["lineEdit_infoNK_5"]=ui->lineEdit_infoNK_5;
 idToLineEdit["lineEdit_infoNK_6"]=ui->lineEdit_infoNK_6;
 idToLineEdit["lineEdit_infoNK_7"]=ui->lineEdit_infoNK_7;
 idToLineEdit["lineEdit_infoNK_8"]=ui->lineEdit_infoNK_8;
 idToLineEdit["WLminNK1"]=ui->WLminNK1;
 idToLineEdit["WLmaxNK1"]=ui->WLmaxNK1;
 idToLineEdit["WLminNK2"]=ui->WLminNK2;
 idToLineEdit["WLmaxNK2"]=ui->WLmaxNK2;
 idToLineEdit["WLminNK3"]=ui->WLminNK3;
 idToLineEdit["WLmaxNK3"]=ui->WLmaxNK3;
 idToLineEdit["WLminNK4"]=ui->WLminNK4;
 idToLineEdit["WLmaxNK4"]=ui->WLmaxNK4;
 idToLineEdit["WLminNK5"]=ui->WLminNK5;
 idToLineEdit["WLmaxNK5"]=ui->WLmaxNK5;
 idToLineEdit["WLminNK6"]=ui->WLminNK6;
 idToLineEdit["WLmaxNK6"]=ui->WLmaxNK6;
 idToLineEdit["WLminNK7"]=ui->WLminNK7;
 idToLineEdit["WLmaxNK7"]=ui->WLmaxNK7;
 idToLineEdit["WLminNK8"]=ui->WLminNK8;
 idToLineEdit["WLmaxNK8"]=ui->WLmaxNK8;
 idToLineEdit["WLmin1"]=ui->WLmin1;
 idToLineEdit["WLmax1"]=ui->WLmax1;
 idToLineEdit["WLmin2"]=ui->WLmin2;
 idToLineEdit["WLmax2"]=ui->WLmax2;
 idToLineEdit["WLmin3"]=ui->WLmin3;
 idToLineEdit["WLmax3"]=ui->WLmax3;
 idToLineEdit["WLmin4"]=ui->WLmin4;
 idToLineEdit["WLmax4"]=ui->WLmax4;
 idToLineEdit["WLmin5"]=ui->WLmin5;
 idToLineEdit["WLmax5"]=ui->WLmax5;
 idToLineEdit["WLmin6"]=ui->WLmin6;
 idToLineEdit["WLmax6"]=ui->WLmax6;
 idToLineEdit["WLmin7"]=ui->WLmin7;
 idToLineEdit["WLmax7"]=ui->WLmax7;
 idToLineEdit["WLmin9"]=ui->WLmin9;
 idToLineEdit["WLmax9"]=ui->WLmax9;
 idToLineEdit["WLmin11"]=ui->WLmin11;
 idToLineEdit["WLmax11"]=ui->WLmax11;
 idToLineEdit["WLmin13"]=ui->WLmin13;
 idToLineEdit["WLmax13"]=ui->WLmax13;
 idToLineEdit["DP_RXY_1_1"]=ui->DP_RXY_1_1;
 idToLineEdit["DP_RXY_1_2"]=ui->DP_RXY_1_2;
 idToLineEdit["DP_RXY_1_3"]=ui->DP_RXY_1_3;
 idToLineEdit["DP_RXY_1_4"]=ui->DP_RXY_1_4;
 idToLineEdit["DP_RXY_2_1"]=ui->DP_RXY_2_1;
 idToLineEdit["DP_RXY_2_2"]=ui->DP_RXY_2_2;
 idToLineEdit["DP_RXY_2_3"]=ui->DP_RXY_2_3;
 idToLineEdit["DP_RXY_2_4"]=ui->DP_RXY_2_4;
 idToLineEdit["DP_RXY_3_1"]=ui->DP_RXY_3_1;
 idToLineEdit["DP_RXY_3_2"]=ui->DP_RXY_3_2;
 idToLineEdit["DP_RXY_3_3"]=ui->DP_RXY_3_3;
 idToLineEdit["DP_RXY_3_4"]=ui->DP_RXY_3_4;
 idToLineEdit["DP_RXY_4_1"]=ui->DP_RXY_4_1;
 idToLineEdit["DP_RXY_4_2"]=ui->DP_RXY_4_2;
 idToLineEdit["DP_RXY_4_3"]=ui->DP_RXY_4_3;
 idToLineEdit["DP_RXY_4_4"]=ui->DP_RXY_4_4;
 idToLineEdit["DP_RXY_5_1"]=ui->DP_RXY_5_1;
 idToLineEdit["DP_RXY_5_2"]=ui->DP_RXY_5_2;
 idToLineEdit["DP_RXY_5_3"]=ui->DP_RXY_5_3;
 idToLineEdit["DP_RXY_5_4"]=ui->DP_RXY_5_4;
 idToLineEdit["DP_RXY_6_1"]=ui->DP_RXY_6_1;
 idToLineEdit["DP_RXY_6_2"]=ui->DP_RXY_6_2;
 idToLineEdit["DP_RXY_6_3"]=ui->DP_RXY_6_3;
 idToLineEdit["DP_RXY_6_4"]=ui->DP_RXY_6_4;
 idToLineEdit["DP_RXY_7_1"]=ui->DP_RXY_7_1;
 idToLineEdit["DP_RXY_7_2"]=ui->DP_RXY_7_2;
 idToLineEdit["DP_RXY_7_3"]=ui->DP_RXY_7_3;
 idToLineEdit["DP_RXY_7_4"]=ui->DP_RXY_7_4;
 idToLineEdit["DP_RXY_8_1"]=ui->DP_RXY_8_1;
 idToLineEdit["DP_RXY_8_2"]=ui->DP_RXY_8_2;
 idToLineEdit["DP_RXY_8_3"]=ui->DP_RXY_8_3;
 idToLineEdit["DP_RXY_8_4"]=ui->DP_RXY_8_4;
 idToLineEdit["DP_RXY_9_1"]=ui->DP_RXY_9_1;
 idToLineEdit["DP_RXY_9_2"]=ui->DP_RXY_9_2;
 idToLineEdit["DP_RXY_9_3"]=ui->DP_RXY_9_3;
 idToLineEdit["DP_RXY_9_4"]=ui->DP_RXY_9_4;
 idToLineEdit["DP_RXY_10_1"]=ui->DP_RXY_10_1;
 idToLineEdit["DP_RXY_10_2"]=ui->DP_RXY_10_2;
 idToLineEdit["DP_RXY_10_3"]=ui->DP_RXY_10_3;
 idToLineEdit["DP_RXY_10_4"]=ui->DP_RXY_10_4;
 idToLineEdit["DP_RXY_11_1"]=ui->DP_RXY_11_1;
 idToLineEdit["DP_RXY_11_2"]=ui->DP_RXY_11_2;
 idToLineEdit["DP_RXY_11_3"]=ui->DP_RXY_11_3;
 idToLineEdit["DP_RXY_11_4"]=ui->DP_RXY_11_4;
 idToLineEdit["DP_RXY_12_1"]=ui->DP_RXY_12_1;
 idToLineEdit["DP_RXY_12_2"]=ui->DP_RXY_12_2;
 idToLineEdit["DP_RXY_12_3"]=ui->DP_RXY_12_3;
 idToLineEdit["DP_RXY_12_4"]=ui->DP_RXY_12_4;
 idToLineEdit["DP_RXY_13_1"]=ui->DP_RXY_13_1;
 idToLineEdit["DP_RXY_13_2"]=ui->DP_RXY_13_2;
 idToLineEdit["DP_RXY_13_3"]=ui->DP_RXY_13_3;
 idToLineEdit["DP_RXY_13_4"]=ui->DP_RXY_13_4;
 idToLineEdit["DP_RXY_14_1"]=ui->DP_RXY_14_1;
 idToLineEdit["DP_RXY_14_2"]=ui->DP_RXY_14_2;
 idToLineEdit["DP_RXY_14_3"]=ui->DP_RXY_14_3;
 idToLineEdit["DP_RXY_14_4"]=ui->DP_RXY_14_4;
 idToLineEdit["DP_RXY_15_1"]=ui->DP_RXY_15_1;
 idToLineEdit["DP_RXY_15_2"]=ui->DP_RXY_15_2;
 idToLineEdit["DP_RXY_15_3"]=ui->DP_RXY_15_3;
 idToLineEdit["DP_RXY_15_4"]=ui->DP_RXY_15_4;
 idToLineEdit["DP_RXY_16_1"]=ui->DP_RXY_16_1;
 idToLineEdit["DP_RXY_16_2"]=ui->DP_RXY_16_2;
 idToLineEdit["DP_RXY_16_3"]=ui->DP_RXY_16_3;
 idToLineEdit["DP_RXY_16_4"]=ui->DP_RXY_16_4;
 idToLineEdit["DP_RXY_17_1"]=ui->DP_RXY_17_1;
 idToLineEdit["DP_RXY_17_2"]=ui->DP_RXY_17_2;
 idToLineEdit["DP_RXY_17_3"]=ui->DP_RXY_17_3;
 idToLineEdit["DP_RXY_17_4"]=ui->DP_RXY_17_4;
 idToLineEdit["DP_RXY_18_1"]=ui->DP_RXY_18_1;
 idToLineEdit["DP_RXY_18_2"]=ui->DP_RXY_18_2;
 idToLineEdit["DP_RXY_18_3"]=ui->DP_RXY_18_3;
 idToLineEdit["DP_RXY_18_4"]=ui->DP_RXY_18_4;
//  idToLineEdit["DP_RXY_19_1"]=ui->DP_RXY_19_1;
//  idToLineEdit["DP_RXY_19_2"]=ui->DP_RXY_19_2;
//  idToLineEdit["DP_RXY_19_3"]=ui->DP_RXY_19_3;
//  idToLineEdit["DP_RXY_19_4"]=ui->DP_RXY_19_4;
 idToLineEdit["DP_RXY_20_1"]=ui->DP_RXY_20_1;
 idToLineEdit["DP_RXY_20_2"]=ui->DP_RXY_20_2;
 idToLineEdit["DP_RXY_20_3"]=ui->DP_RXY_20_3;
 idToLineEdit["DP_RXY_20_4"]=ui->DP_RXY_20_4;
 idToLineEdit["DP_RXY_21_1"]=ui->DP_RXY_21_1;
 idToLineEdit["DP_RXY_21_2"]=ui->DP_RXY_21_2;
 idToLineEdit["DP_RXY_21_3"]=ui->DP_RXY_21_3;
 idToLineEdit["DP_RXY_21_4"]=ui->DP_RXY_21_4;
 idToLineEdit["DP_RXY_22_1"]=ui->DP_RXY_22_1;
 idToLineEdit["DP_RXY_22_2"]=ui->DP_RXY_22_2;
 idToLineEdit["DP_RXY_22_3"]=ui->DP_RXY_22_3;
 idToLineEdit["DP_RXY_22_4"]=ui->DP_RXY_22_4;
 idToLineEdit["DP_RXY_23_1"]=ui->DP_RXY_23_1;
 idToLineEdit["DP_RXY_23_2"]=ui->DP_RXY_23_2;
 idToLineEdit["DP_RXY_23_3"]=ui->DP_RXY_23_3;
 idToLineEdit["DP_RXY_23_4"]=ui->DP_RXY_23_4;
 idToLineEdit["DP_RXY_24_1"]=ui->DP_RXY_24_1;
 idToLineEdit["DP_RXY_24_2"]=ui->DP_RXY_24_2;
 idToLineEdit["DP_RXY_24_3"]=ui->DP_RXY_24_3;
 idToLineEdit["DP_RXY_24_4"]=ui->DP_RXY_24_4;
//  idToLineEdit["DP_RXY_25_1"]=ui->DP_RXY_25_1;
//  idToLineEdit["DP_RXY_25_2"]=ui->DP_RXY_25_2;
//  idToLineEdit["DP_RXY_25_3"]=ui->DP_RXY_25_3;
//  idToLineEdit["DP_RXY_25_4"]=ui->DP_RXY_25_4;
 idToLineEdit["DP_RXY_26_1"]=ui->DP_RXY_26_1;
 idToLineEdit["DP_RXY_26_2"]=ui->DP_RXY_26_2;
 idToLineEdit["DP_RXY_26_3"]=ui->DP_RXY_26_3;
 idToLineEdit["DP_RXY_26_4"]=ui->DP_RXY_26_4;
 idToLineEdit["DP_RXY_27_1"]=ui->DP_RXY_27_1;
 idToLineEdit["DP_RXY_27_2"]=ui->DP_RXY_27_2;
 idToLineEdit["DP_RXY_27_3"]=ui->DP_RXY_27_3;
 idToLineEdit["DP_RXY_27_4"]=ui->DP_RXY_27_4;
 idToLineEdit["LEcnk1_2"]=ui->LEcnk1_2;
 idToLineEdit["LEcnk1_3"]=ui->LEcnk1_3;
 idToLineEdit["LEcnk2_2"]=ui->LEcnk2_2;
 idToLineEdit["LEcnk2_3"]=ui->LEcnk2_3;
 idToLineEdit["LEcnk3_2"]=ui->LEcnk3_2;
 idToLineEdit["LEcnk3_3"]=ui->LEcnk3_3;
 idToLineEdit["LEcnk4_2"]=ui->LEcnk4_2;
 idToLineEdit["LEcnk4_3"]=ui->LEcnk4_3;
 idToLineEdit["LEcnk5_2"]=ui->LEcnk5_2;
 idToLineEdit["LEcnk5_3"]=ui->LEcnk5_3;
 idToLineEdit["LEcnk6_2"]=ui->LEcnk6_2;
 idToLineEdit["LEcnk6_3"]=ui->LEcnk6_3;
 idToLineEdit["LEcnk7_2"]=ui->LEcnk7_2;
 idToLineEdit["LEcnk7_3"]=ui->LEcnk7_3;
 idToLineEdit["LEcnk8_2"]=ui->LEcnk8_2;
 idToLineEdit["LEcnk8_3"]=ui->LEcnk8_3;
 idToLineEdit["LEcnk9_2"]=ui->LEcnk9_2;
 idToLineEdit["LEcnk9_3"]=ui->LEcnk9_3;
 idToLineEdit["LEcnk10_2"]=ui->LEcnk10_2;
 idToLineEdit["LEcnk10_3"]=ui->LEcnk10_3;
 idToLineEdit["LEcnk11_2"]=ui->LEcnk11_2;
 idToLineEdit["LEcnk11_3"]=ui->LEcnk11_3;
 idToLineEdit["LEcnk12_2"]=ui->LEcnk12_2;
 idToLineEdit["LEcnk12_3"]=ui->LEcnk12_3;
 idToLineEdit["LEcnk13_2"]=ui->LEcnk13_2;
 idToLineEdit["LEcnk13_3"]=ui->LEcnk13_3;
 idToLineEdit["LEcnk14_2"]=ui->LEcnk14_2;
 idToLineEdit["LEcnk14_3"]=ui->LEcnk14_3;
 idToLineEdit["LEcnk15_2"]=ui->LEcnk15_2;
 idToLineEdit["LEcnk15_3"]=ui->LEcnk15_3;
 idToLineEdit["LEpm_102_1"]=ui->LEpm_102_1;
 idToLineEdit["LEpm_103_1"]=ui->LEpm_103_1;
 idToLineEdit["LEpm_104_1"]=ui->LEpm_104_1;
 idToLineEdit["LEpm_105_1"]=ui->LEpm_105_1;
 idToLineEdit["LEpm_107_1"]=ui->LEpm_107_1;
 idToLineEdit["LEpm_108_1"]=ui->LEpm_108_1;
 idToLineEdit["LEpm_109_1"]=ui->LEpm_109_1;
 idToLineEdit["LEpm_110_1"]=ui->LEpm_110_1;
 idToLineEdit["LEpm_112_1"]=ui->LEpm_112_1;
 idToLineEdit["LEpm_113_1"]=ui->LEpm_113_1;
 idToLineEdit["LEpm_114_1"]=ui->LEpm_114_1;
 idToLineEdit["LEpm_115_1"]=ui->LEpm_115_1;
 idToLineEdit["LEpm_117_1"]=ui->LEpm_117_1;
 idToLineEdit["LEpm_118_1"]=ui->LEpm_118_1;
 idToLineEdit["LEpm_119_1"]=ui->LEpm_119_1;
 idToLineEdit["LEpm_120_1"]=ui->LEpm_120_1;
 idToLineEdit["LEpm_122_1"]=ui->LEpm_122_1;
 idToLineEdit["LEpm_123_1"]=ui->LEpm_123_1;
 idToLineEdit["LEpm_124_1"]=ui->LEpm_124_1;
 idToLineEdit["LEpm_125_1"]=ui->LEpm_125_1;
 idToLineEdit["LEpm_127_1"]=ui->LEpm_127_1;
 idToLineEdit["LEpm_128_1"]=ui->LEpm_128_1;
 idToLineEdit["LEpm_129_1"]=ui->LEpm_129_1;
 idToLineEdit["LEpm_130_1"]=ui->LEpm_130_1;
 idToLineEdit["LEpm_132_1"]=ui->LEpm_132_1;
 idToLineEdit["LEpm_133_1"]=ui->LEpm_133_1;
 idToLineEdit["LEpm_134_1"]=ui->LEpm_134_1;
 idToLineEdit["LEpm_135_1"]=ui->LEpm_135_1;
 idToLineEdit["LEpm_137_1"]=ui->LEpm_137_1;
 idToLineEdit["LEpm_138_1"]=ui->LEpm_138_1;
 idToLineEdit["LEpm_139_1"]=ui->LEpm_139_1;
 idToLineEdit["LEpm_140_1"]=ui->LEpm_140_1;
 idToLineEdit["LEpm_142_1"]=ui->LEpm_142_1;
 idToLineEdit["LEpm_143_1"]=ui->LEpm_143_1;
 idToLineEdit["LEpm_144_1"]=ui->LEpm_144_1;
 idToLineEdit["LEpm_145_1"]=ui->LEpm_145_1;
 idToLineEdit["LEpm_147_1"]=ui->LEpm_147_1;
 idToLineEdit["LEpm_148_1"]=ui->LEpm_148_1;
 idToLineEdit["LEpm_149_1"]=ui->LEpm_149_1;
 idToLineEdit["LEpm_150_1"]=ui->LEpm_150_1;
 idToLineEdit["LEpm_152_1"]=ui->LEpm_152_1;
 idToLineEdit["LEpm_153_1"]=ui->LEpm_153_1;
 idToLineEdit["LEpm_154_1"]=ui->LEpm_154_1;
 idToLineEdit["LEpm_155_1"]=ui->LEpm_155_1;
 idToLineEdit["LEpm_157_1"]=ui->LEpm_157_1;
 idToLineEdit["LEpm_158_1"]=ui->LEpm_158_1;
 idToLineEdit["LEpm_159_1"]=ui->LEpm_159_1;
 idToLineEdit["LEpm_160_1"]=ui->LEpm_160_1;
 idToLineEdit["LEpm_162_1"]=ui->LEpm_162_1;
 idToLineEdit["LEpm_163_1"]=ui->LEpm_163_1;
 idToLineEdit["LEpm_164_1"]=ui->LEpm_164_1;
 idToLineEdit["LEpm_165_1"]=ui->LEpm_165_1;
 idToLineEdit["LEpm_167_1"]=ui->LEpm_167_1;
 idToLineEdit["LEpm_168_1"]=ui->LEpm_168_1;
 idToLineEdit["LEpm_169_1"]=ui->LEpm_169_1;
 idToLineEdit["LEpm_170_1"]=ui->LEpm_170_1;
 idToLineEdit["LEpm_172_1"]=ui->LEpm_172_1;
 idToLineEdit["LEpm_173_1"]=ui->LEpm_173_1;
 idToLineEdit["LEpm_174_1"]=ui->LEpm_174_1;
 idToLineEdit["LEpm_175_1"]=ui->LEpm_175_1;
 idToLineEdit["LEpm_177_1"]=ui->LEpm_177_1;
 idToLineEdit["LEpm_178_1"]=ui->LEpm_178_1;
 idToLineEdit["LEpm_179_1"]=ui->LEpm_179_1;
 idToLineEdit["LEpm_180_1"]=ui->LEpm_180_1;
 idToLineEdit["LEpm_182_1"]=ui->LEpm_182_1;
 idToLineEdit["LEpm_183_1"]=ui->LEpm_183_1;
 idToLineEdit["LEpm_184_1"]=ui->LEpm_184_1;
 idToLineEdit["LEpm_185_1"]=ui->LEpm_185_1;
 idToLineEdit["LEpm_187_1"]=ui->LEpm_187_1;
 idToLineEdit["LEpm_188_1"]=ui->LEpm_188_1;
 idToLineEdit["LEpm_189_1"]=ui->LEpm_189_1;
 idToLineEdit["LEpm_190_1"]=ui->LEpm_190_1;
 idToLineEdit["LEpm_192_1"]=ui->LEpm_192_1;
 idToLineEdit["LEpm_193_1"]=ui->LEpm_193_1;
 idToLineEdit["LEpm_194_1"]=ui->LEpm_194_1;
 idToLineEdit["LEpm_195_1"]=ui->LEpm_195_1;
 idToLineEdit["LEpm_197_1"]=ui->LEpm_197_1;
 idToLineEdit["LEpm_198_1"]=ui->LEpm_198_1;
 idToLineEdit["LEpm_199_1"]=ui->LEpm_199_1;
 idToLineEdit["LEpm_200_1"]=ui->LEpm_200_1;
 idToLineEdit["DPparFitV_1"]=ui->DPparFitV_1;
 idToLineEdit["DPparFitV_2"]=ui->DPparFitV_2;
 idToLineEdit["DPparFitV_3"]=ui->DPparFitV_3;
 idToLineEdit["DPparFitV_4"]=ui->DPparFitV_4;
 idToLineEdit["DPparFitV_5"]=ui->DPparFitV_5;
 idToLineEdit["DPparFitV_6"]=ui->DPparFitV_6;
 idToLineEdit["DPparFitV_7"]=ui->DPparFitV_7;
 idToLineEdit["DPparFitV_8"]=ui->DPparFitV_8;
 idToLineEdit["DPparFitV_9"]=ui->DPparFitV_9;
 idToLineEdit["DPparFitV_10"]=ui->DPparFitV_10;
 idToLineEdit["DPparFitV_11"]=ui->DPparFitV_11;
 idToLineEdit["DPparFitV_12"]=ui->DPparFitV_12;
 idToLineEdit["DPparFitV_13"]=ui->DPparFitV_13;
 idToLineEdit["DPparFitV_14"]=ui->DPparFitV_14;
 idToLineEdit["DPparFitV_15"]=ui->DPparFitV_15;
 idToLineEdit["DPparFitV_16"]=ui->DPparFitV_16;
 idToLineEdit["DPparFitV_17"]=ui->DPparFitV_17;
 idToLineEdit["DPparFitErr_1"]=ui->DPparFitErr_1;
 idToLineEdit["DPparFitErr_2"]=ui->DPparFitErr_2;
 idToLineEdit["DPparFitErr_3"]=ui->DPparFitErr_3;
 idToLineEdit["DPparFitErr_4"]=ui->DPparFitErr_4;
 idToLineEdit["DPparFitErr_5"]=ui->DPparFitErr_5;
 idToLineEdit["DPparFitErr_6"]=ui->DPparFitErr_6;
 idToLineEdit["DPparFitErr_7"]=ui->DPparFitErr_7;
 idToLineEdit["DPparFitErr_8"]=ui->DPparFitErr_8;
 idToLineEdit["DPparFitErr_9"]=ui->DPparFitErr_9;
 idToLineEdit["DPparFitErr_10"]=ui->DPparFitErr_10;
 idToLineEdit["DPparFitErr_11"]=ui->DPparFitErr_11;
 idToLineEdit["DPparFitErr_12"]=ui->DPparFitErr_12;
 idToLineEdit["DPparFitErr_13"]=ui->DPparFitErr_13;
 idToLineEdit["DPparFitErr_14"]=ui->DPparFitErr_14;
 idToLineEdit["DPparFitErr_15"]=ui->DPparFitErr_15;
 idToLineEdit["DPparFitErr_16"]=ui->DPparFitErr_16;
 idToLineEdit["DPparFitErr_17"]=ui->DPparFitErr_17;
 idToLineEdit["DPparFitGC_1"]=ui->DPparFitGC_1;
 idToLineEdit["DPparFitGC_2"]=ui->DPparFitGC_2;
 idToLineEdit["DPparFitGC_3"]=ui->DPparFitGC_3;
 idToLineEdit["DPparFitGC_4"]=ui->DPparFitGC_4;
 idToLineEdit["DPparFitGC_5"]=ui->DPparFitGC_5;
 idToLineEdit["DPparFitGC_6"]=ui->DPparFitGC_6;
 idToLineEdit["DPparFitGC_7"]=ui->DPparFitGC_7;
 idToLineEdit["DPparFitGC_8"]=ui->DPparFitGC_8;
 idToLineEdit["DPparFitGC_9"]=ui->DPparFitGC_9;
 idToLineEdit["DPparFitGC_10"]=ui->DPparFitGC_10;
 idToLineEdit["DPparFitGC_11"]=ui->DPparFitGC_11;
 idToLineEdit["DPparFitGC_12"]=ui->DPparFitGC_12;
 idToLineEdit["DPparFitGC_13"]=ui->DPparFitGC_13;
 idToLineEdit["DPparFitGC_14"]=ui->DPparFitGC_14;
 idToLineEdit["DPparFitGC_15"]=ui->DPparFitGC_15;
 idToLineEdit["DPparFitGC_16"]=ui->DPparFitGC_16;
 idToLineEdit["DPparFitGC_17"]=ui->DPparFitGC_17;
 idToLineEdit["LEpar_30_1"]=ui->LEpar_30_1;
 idToLineEdit["LEpar_30_2"]=ui->LEpar_30_2;
 idToLineEdit["LEpar_30_3"]=ui->LEpar_30_3;
 idToLineEdit["LEpar_30_4"]=ui->LEpar_30_4;
 idToLineEdit["LEpar_30_5"]=ui->LEpar_30_5;
 idToLineEdit["LEpar_36_1"]=ui->LEpar_36_1;
 idToLineEdit["LEpar_37_1"]=ui->LEpar_37_1;
 idToLineEdit["LEpar_38_1"]=ui->LEpar_38_1;
 idToLineEdit["LEpar_39_1"]=ui->LEpar_39_1;
 idToLineEdit["LEpar_40_1"]=ui->LEpar_40_1;
 idToLineEdit["LEpar_41_1"]=ui->LEpar_41_1;
 idToLineEdit["LEpar_42_1"]=ui->LEpar_42_1;
 idToLineEdit["LEpar_43_1"]=ui->LEpar_43_1;
 idToLineEdit["LEpar_44_1"]=ui->LEpar_44_1;
 idToLineEdit["LEpar_45_1"]=ui->LEpar_45_1;
 idToLineEdit["LEpar_46_1"]=ui->LEpar_46_1;
 idToLineEdit["LEpar_47_1"]=ui->LEpar_47_1;
 idToLineEdit["LEpar_48_1"]=ui->LEpar_48_1;
 idToLineEdit["LEpar_49_1"]=ui->LEpar_49_1;
 idToLineEdit["LEpar_36_2"]=ui->LEpar_36_2;
 idToLineEdit["LEpar_37_2"]=ui->LEpar_37_2;
 idToLineEdit["LEpar_38_2"]=ui->LEpar_38_2;
 idToLineEdit["LEpar_39_2"]=ui->LEpar_39_2;
 idToLineEdit["LEpar_40_2"]=ui->LEpar_40_2;
 idToLineEdit["LEpar_41_2"]=ui->LEpar_41_2;
 idToLineEdit["LEpar_42_2"]=ui->LEpar_42_2;
 idToLineEdit["LEpar_43_2"]=ui->LEpar_43_2;
 idToLineEdit["LEpar_44_2"]=ui->LEpar_44_2;
 idToLineEdit["LEpar_45_2"]=ui->LEpar_45_2;
 idToLineEdit["LEpar_46_2"]=ui->LEpar_46_2;
 idToLineEdit["LEpar_47_2"]=ui->LEpar_47_2;
 idToLineEdit["LEpar_48_2"]=ui->LEpar_48_2;
 idToLineEdit["LEpar_49_2"]=ui->LEpar_49_2;
 idToLineEdit["LEpar_36_3"]=ui->LEpar_36_3;
 idToLineEdit["LEpar_37_3"]=ui->LEpar_37_3;
 idToLineEdit["LEpar_38_3"]=ui->LEpar_38_3;
 idToLineEdit["LEpar_39_3"]=ui->LEpar_39_3;
 idToLineEdit["LEpar_40_3"]=ui->LEpar_40_3;
 idToLineEdit["LEpar_41_3"]=ui->LEpar_41_3;
 idToLineEdit["LEpar_42_3"]=ui->LEpar_42_3;
 idToLineEdit["LEpar_43_3"]=ui->LEpar_43_3;
 idToLineEdit["LEpar_44_3"]=ui->LEpar_44_3;
 idToLineEdit["LEpar_45_3"]=ui->LEpar_45_3;
 idToLineEdit["LEpar_46_3"]=ui->LEpar_46_3;
 idToLineEdit["LEpar_47_3"]=ui->LEpar_47_3;
 idToLineEdit["LEpar_48_3"]=ui->LEpar_48_3;
 idToLineEdit["LEpar_49_3"]=ui->LEpar_49_3;
 
 idToDoubleSpinBox["dSB_cnk1"]=ui->dSB_cnk1;
 idToDoubleSpinBox["dSB_cnk2"]=ui->dSB_cnk2;
 idToDoubleSpinBox["dSB_cnk3"]=ui->dSB_cnk3;
 idToDoubleSpinBox["dSB_cnk4"]=ui->dSB_cnk4;
 idToDoubleSpinBox["dSB_cnk5"]=ui->dSB_cnk5;
 idToDoubleSpinBox["dSB_cnk6"]=ui->dSB_cnk6;
 idToDoubleSpinBox["dSB_cnk7"]=ui->dSB_cnk7;
 idToDoubleSpinBox["dSB_cnk8"]=ui->dSB_cnk8;
 idToDoubleSpinBox["dSB_cnk9"]=ui->dSB_cnk9;
 idToDoubleSpinBox["dSB_cnk10"]=ui->dSB_cnk10;
 idToDoubleSpinBox["dSB_cnk11"]=ui->dSB_cnk11;
 idToDoubleSpinBox["dSB_cnk12"]=ui->dSB_cnk12;
 idToDoubleSpinBox["dSB_cnk13"]=ui->dSB_cnk13;
 idToDoubleSpinBox["dSB_cnk14"]=ui->dSB_cnk14;
 idToDoubleSpinBox["dSB_cnk15"]=ui->dSB_cnk15;
 idToDoubleSpinBox["dSB_PM_1_1"]=ui->dSB_PM_1_1;
 idToDoubleSpinBox["dSB_PM_2_1"]=ui->dSB_PM_2_1;
 idToDoubleSpinBox["dSB_PM_3_1"]=ui->dSB_PM_3_1;
 idToDoubleSpinBox["dSB_PM_4_1"]=ui->dSB_PM_4_1;
 idToDoubleSpinBox["dSB_PM_5_1"]=ui->dSB_PM_5_1;
 idToDoubleSpinBox["dSB_PM_6_1"]=ui->dSB_PM_6_1;
 idToDoubleSpinBox["dSB_PM_7_1"]=ui->dSB_PM_7_1;
 idToDoubleSpinBox["dSB_PM_8_1"]=ui->dSB_PM_8_1;
 idToDoubleSpinBox["dSB_PM_9_1"]=ui->dSB_PM_9_1;
 idToDoubleSpinBox["dSB_PM_11_1"]=ui->dSB_PM_11_1;
 idToDoubleSpinBox["dSB_PM_12_1"]=ui->dSB_PM_12_1;
 idToDoubleSpinBox["dSB_PM_13_1"]=ui->dSB_PM_13_1;
 idToDoubleSpinBox["dSB_PM_14_1"]=ui->dSB_PM_14_1;
 idToDoubleSpinBox["dSB_PM_15_1"]=ui->dSB_PM_15_1;
 idToDoubleSpinBox["dSB_PM_16_1"]=ui->dSB_PM_16_1;
 idToDoubleSpinBox["dSB_PM_17_1"]=ui->dSB_PM_17_1;
 idToDoubleSpinBox["dSB_PM_18_1"]=ui->dSB_PM_18_1;
 idToDoubleSpinBox["dSB_PM_19_1"]=ui->dSB_PM_19_1;
 idToDoubleSpinBox["dSB_PM_21_1"]=ui->dSB_PM_21_1;
 idToDoubleSpinBox["dSB_PM_22_1"]=ui->dSB_PM_22_1;
 idToDoubleSpinBox["dSB_PM_23_1"]=ui->dSB_PM_23_1;
 idToDoubleSpinBox["dSB_PM_24_1"]=ui->dSB_PM_24_1;
 idToDoubleSpinBox["dSB_PM_25_1"]=ui->dSB_PM_25_1;
 idToDoubleSpinBox["dSB_PM_26_1"]=ui->dSB_PM_26_1;
 idToDoubleSpinBox["dSB_PM_27_1"]=ui->dSB_PM_27_1;
 idToDoubleSpinBox["dSB_PM_28_1"]=ui->dSB_PM_28_1;
 idToDoubleSpinBox["dSB_PM_29_1"]=ui->dSB_PM_29_1;
 idToDoubleSpinBox["dSB_PM_31_1"]=ui->dSB_PM_31_1;
 idToDoubleSpinBox["dSB_PM_32_1"]=ui->dSB_PM_32_1;
 idToDoubleSpinBox["dSB_PM_33_1"]=ui->dSB_PM_33_1;
 idToDoubleSpinBox["dSB_PM_34_1"]=ui->dSB_PM_34_1;
 idToDoubleSpinBox["dSB_PM_35_1"]=ui->dSB_PM_35_1;
 idToDoubleSpinBox["dSB_PM_36_1"]=ui->dSB_PM_36_1;
 idToDoubleSpinBox["dSB_PM_37_1"]=ui->dSB_PM_37_1;
 idToDoubleSpinBox["dSB_PM_38_1"]=ui->dSB_PM_38_1;
 idToDoubleSpinBox["dSB_PM_39_1"]=ui->dSB_PM_39_1;
 idToDoubleSpinBox["dSB_PM_41_1"]=ui->dSB_PM_41_1;
 idToDoubleSpinBox["dSB_PM_42_1"]=ui->dSB_PM_42_1;
 idToDoubleSpinBox["dSB_PM_43_1"]=ui->dSB_PM_43_1;
 idToDoubleSpinBox["dSB_PM_44_1"]=ui->dSB_PM_44_1;
 idToDoubleSpinBox["dSB_PM_45_1"]=ui->dSB_PM_45_1;
 idToDoubleSpinBox["dSB_PM_46_1"]=ui->dSB_PM_46_1;
 idToDoubleSpinBox["dSB_PM_47_1"]=ui->dSB_PM_47_1;
 idToDoubleSpinBox["dSB_PM_48_1"]=ui->dSB_PM_48_1;
 idToDoubleSpinBox["dSB_PM_49_1"]=ui->dSB_PM_49_1;
 idToDoubleSpinBox["dSB_PM_51_1"]=ui->dSB_PM_51_1;
 idToDoubleSpinBox["dSB_PM_52_1"]=ui->dSB_PM_52_1;
 idToDoubleSpinBox["dSB_PM_53_1"]=ui->dSB_PM_53_1;
 idToDoubleSpinBox["dSB_PM_54_1"]=ui->dSB_PM_54_1;
 idToDoubleSpinBox["dSB_PM_55_1"]=ui->dSB_PM_55_1;
 idToDoubleSpinBox["dSB_PM_56_1"]=ui->dSB_PM_56_1;
 idToDoubleSpinBox["dSB_PM_57_1"]=ui->dSB_PM_57_1;
 idToDoubleSpinBox["dSB_PM_58_1"]=ui->dSB_PM_58_1;
 idToDoubleSpinBox["dSB_PM_59_1"]=ui->dSB_PM_59_1;
 idToDoubleSpinBox["dSB_PM_61_1"]=ui->dSB_PM_61_1;
 idToDoubleSpinBox["dSB_PM_62_1"]=ui->dSB_PM_62_1;
 idToDoubleSpinBox["dSB_PM_63_1"]=ui->dSB_PM_63_1;
 idToDoubleSpinBox["dSB_PM_64_1"]=ui->dSB_PM_64_1;
 idToDoubleSpinBox["dSB_PM_65_1"]=ui->dSB_PM_65_1;
 idToDoubleSpinBox["dSB_PM_66_1"]=ui->dSB_PM_66_1;
 idToDoubleSpinBox["dSB_PM_67_1"]=ui->dSB_PM_67_1;
 idToDoubleSpinBox["dSB_PM_68_1"]=ui->dSB_PM_68_1;
 idToDoubleSpinBox["dSB_PM_69_1"]=ui->dSB_PM_69_1;

//oscillator list
 for(int k=1;k<=20;k++){
    idToComboBox["cBpm_"+QString::number(100+1+(k-1)*5)+"_1"]-> clear();
    idToComboBox["cBpm_"+QString::number(100+1+(k-1)*5)+"_1"]->addItem("Lorentz");
    idToComboBox["cBpm_"+QString::number(100+1+(k-1)*5)+"_1"]->addItem("Quant-homo");
    idToComboBox["cBpm_"+QString::number(100+1+(k-1)*5)+"_1"]->addItem("Quant-inhomo");
    idToComboBox["cBpm_"+QString::number(100+1+(k-1)*5)+"_1"]->addItem("Flat");
    idToComboBox["cBpm_"+QString::number(100+1+(k-1)*5)+"_1"]->addItem("Drude");
    idToComboBox["cBpm_"+QString::number(100+1+(k-1)*5)+"_1"]->addItem("Indirect-Gap-Cody");
    idToComboBox["cBpm_"+QString::number(100+1+(k-1)*5)+"_1"]->addItem("Indirect-Gap-Tauc");
    idToComboBox["cBpm_"+QString::number(100+1+(k-1)*5)+"_1"]->addItem("Direct-Gap-Cody");
    idToComboBox["cBpm_"+QString::number(100+1+(k-1)*5)+"_1"]->addItem("Direct-Gap-Tauc");
    idToComboBox["cBpm_"+QString::number(100+1+(k-1)*5)+"_1"]->setCurrentIndex(-1);
 }

//Refrence mirror list
 QFile file0(fRefMir);
 if (!file0.open (QIODevice::ReadOnly | QIODevice::Text))
      return;
 QTextStream stream ( &file0 );
 QString line;
 int iB103=0;
 ui -> comboB_PAR_10_3-> clear();
 ui -> comboB_PAR_10_3->addItem("absolute");
 do{
     iB103++;
     line = stream.readLine();
     ui -> comboB_PAR_10_3->addItem(line);
     line = stream.readLine();
     printf("-> RefMir #%d -> %s\n",iB103,line.toStdString().c_str());
 }while(!stream.atEnd());
 file0.close();

//standard spectra for mean value
 QFile file1(fStdSpect);
 if (!file1.open (QIODevice::ReadOnly | QIODevice::Text))
      return;
 QTextStream stream1 ( &file1 );
 ui->comboBox_average -> clear();
 int iStdS=0;
 do{
     line = stream1.readLine();
     ui->comboBox_average -> addItem(line);
     printf("-> StdSpectrum #%d: %s",iStdS,line.toStdString().c_str());
     line = stream1.readLine();
     printf(" -> %s\n",line.toStdString().c_str());
     iStdS++;
 }while(!stream1.atEnd());
 file1.close();

 ParFitLab[1]="_d";
 ParFitLab[2]="_Dn/<n>";
 ParFitLab[3]="_(nav/<n>-1)";
 ParFitLab[4]="_Dk/<k>";
 ParFitLab[5]="_(kav/<k>-1)";
 ParFitLab[6]="_sigma_roughness";
 ParFitLab[7]="_slope_Dn";
 ParFitLab[8]="_C";
 ParFitLab[9]="_E";
 ParFitLab[10]="_D";
 ParFitLab[11]="_K";
 
 QFile file(fctrl);
  if (!file.open(QIODevice::WriteOnly | QIODevice::Text))
  return;
 QTextStream out(&file);
 out << "pause" << "\n";
 file.close();

 SaveSetting(-1);//save
}


void ksemaw::fileChanged (const QString & ){
 QString scmd;
 do{
     QFile filectrl(fctrl);
     if (!filectrl.open(QIODevice::ReadOnly | QIODevice::Text))
         return;
     QTextStream streamctrl ( &filectrl );
     streamctrl >> scmd;
     filectrl.close();
 } while(scmd.contains("proce"));
 //printf("filechanged scmd=%s occupyPF=%d\n",(scmd.toStdString()).c_str(),occupyPF);
 if(scmd.contains("done!") && occupyPF!=1){
  ReadSetting(fileStore);
  QFile filectrl(fctrl);
  if (!filectrl.open(QIODevice::WriteOnly | QIODevice::Text))
  return;
  QTextStream out(&filectrl);
  out << "pause" << "\n";
  filectrl.close();
 }
}



void ksemaw::closeEvent ( QCloseEvent * event )
{
 int ierr=system("pkill ksemawf");
 if(ierr!=0) printf("Error killing ksemawf!\n");
}




void ksemaw::LoadProject(){
 //SaveSetting(-1);
 printf("-> LoadProject\n");
 for(int i=1;i<=8;i++) Clrnk(i);
 ClrFnk();
 fnproject = QFileDialog::getOpenFileName(
        this,
        "Choose a SEMAW project", //titolo della finestra
        pathroot, //directory iniziale
        "Semaw Project (*.Spj)"); //tipi di file da cercare
// printf("fnprojet=%s\n",fnproject.toStdString().c_str());
 ReadSetting(fnproject);
}
 
 
 
void ksemaw::ReadSetting(QString filename){
 if(occupyPF!=0) return;
 occupyPF=1;
 QString scmd;
 QFile filectrl(fctrl);
 if (!filectrl.open(QIODevice::ReadOnly | QIODevice::Text)){
    occupyPF=0;
    return;
 }
 QTextStream streamctrl ( &filectrl );
 streamctrl >> scmd;
 filectrl.close();
//  printf("scmd=%s\n",scmd.toStdString().c_str());
//  printf("filename=%s\n",filename.toStdString().c_str());
 if(filename==fileStore && !scmd.contains("done!")){
     occupyPF=0;
     return;
 }
 printf("-> ReadSetting from %s\n",(filename.toStdString()).c_str());
 Qt::CheckState state;
 QFile file(filename);
 if (!file.open (QIODevice::ReadOnly | QIODevice::Text)){
      occupyPF=0;
      return;
 }
 QTextStream stream ( &file );
 QString line,line2;
 int irx=30;
 line = stream.readLine();
 if(line.contains("iVspj=1"))
     line = stream.readLine();
 else{
     irx=25;
     printf("-> Spj old version\n");
 }
 ui->lineEdit_infoP -> setText(line.simplified());
 line2=fnproject.section(pathroot, 1, 1);
 ui->lineEdit_P -> setText(line2);
 
 for(int i=1;i<=8;i++){
  line = stream.readLine();
  fnk[i]=pathroot+line.simplified()+".nk";
 }
 line = stream.readLine();
 fnFnk=pathroot+line.simplified()+".nk";
 line = stream.readLine();//path-file std spectrum
// if(line.contains("mate/aa999.9", Qt::CaseInsensitive))
//  ui->comboBox_average ->setCurrentIndex(0);
// else if(line.contains("ufo_/media.1", Qt::CaseInsensitive))
//  ui->comboBox_average ->setCurrentIndex(1);
// else if(line.contains("ufo_/media.2", Qt::CaseInsensitive))
//  ui->comboBox_average ->setCurrentIndex(2);
// else if(line.contains("ufo_/media.3", Qt::CaseInsensitive))
//  ui->comboBox_average ->setCurrentIndex(3);
// else if(line.contains("ufo_/media.4", Qt::CaseInsensitive))
//  ui->comboBox_average ->setCurrentIndex(4);
// else if(line.contains("ufo_/media.5", Qt::CaseInsensitive))
//  ui->comboBox_average ->setCurrentIndex(5);
// else if(line.contains("ufo_/media.6", Qt::CaseInsensitive))
//  ui->comboBox_average ->setCurrentIndex(6);
// else if(line.contains("ufo_/pwtdr.1", Qt::CaseInsensitive))
//  ui->comboBox_average ->setCurrentIndex(7);
// else if(line.contains("ufo_/pwiqe.1", Qt::CaseInsensitive))
//  ui->comboBox_average ->setCurrentIndex(8);
// else if(line.contains("ufo_/media.7", Qt::CaseInsensitive))
//  ui->comboBox_average ->setCurrentIndex(9);
// else if(line.contains("ufo_/media.9", Qt::CaseInsensitive))
//  ui->comboBox_average ->setCurrentIndex(10);
 line = stream.readLine();
 fnSample=pathroot+line.simplified();
 line = stream.readLine();
 fnE1=pathroot+line.simplified()+".el";
 line = stream.readLine();
 fnE2=pathroot+line.simplified()+".el";
 line = stream.readLine();
 fnE3=pathroot+line.simplified()+".el";
 line = stream.readLine();
 fnE4=pathroot+line.simplified()+".el";
 ifn=1;
 Setnk1();
 Setnk2();
 Setnk3();
 Setnk4();
 Setnk5();
 Setnk6();
 Setnk7();
 Setnk8();
 LoadFilenk();
 setSample();
 //RXY
 ifn=0;
 for(int i=1;i<=irx;i++){
  for(int j=1;j<=4;j++){
   stream>> rxy[i][j];
   if(i==24 && j==1 && rxy[i][j]< 1.)
       rxy[i][j]=6.;
   if(i==24 && j==2 && rxy[i][j]< 0.5)
       rxy[i][j]=0.5;
   if(i==24 && j==3 && (rxy[i][j]< 1 || rxy[i][j]> 201))
       rxy[i][j]=3.;
   if(idToLineEdit.contains("DP_RXY_"+QString::number(i)+"_"+QString::number(j)))
       idToLineEdit["DP_RXY_"+QString::number(i)+"_"+QString::number(j)]
       -> setText(QString::number(rxy[i][j]));
  }
 }
 if(int(rxy[25][3]+0.5)==1)
   ui->checkB_RXY_25_3 -> setCheckState ( Qt::Unchecked );
 else
   ui->checkB_RXY_25_3 -> setCheckState ( Qt::Checked );
 if(int(rxy[25][4]+0.5)==1)
   ui->checkB_RXY_25_4 -> setCheckState ( Qt::Unchecked );
 else
   ui->checkB_RXY_25_4 -> setCheckState ( Qt::Checked );
 //CNK
 int i0=0,i1=0,i2=0;
 double f2=0.;
 for(int i=1;i<=15;i++){
  for(int j=1;j<=3;j++){
   stream>> cnk[i][j];
  }
 }
 //PAR
 for(int i=1;i<=60;i++){
  for(int j=1;j<=5;j++){
   stream>> par[i][j];
  }
 }
 if(filename!=fileStore){
  par[8][1]=0.;
  par[8][2]=0.;
 }
 int JJ=1,Jitem=-1,nMis;
 for(int i=1;i<=14;i++){
  if(int(par[i][4]+0.5)==-1){
   idToCheckBox["checkB_mis"+QString::number(i)+"_1"] -> setCheckState ( Qt::Unchecked );
   idToCheckBox["checkB_mis"+QString::number(i)+"_2"] -> setCheckState ( Qt::Checked );
   idToCheckBox["checkB_mis"+QString::number(i)+"_3"] -> setCheckState ( Qt::Unchecked );
  } 
  else if(int(par[i][4]+0.5)==0){
   idToCheckBox["checkB_mis"+QString::number(i)+"_1"] -> setCheckState ( Qt::Unchecked );
   idToCheckBox["checkB_mis"+QString::number(i)+"_2"] -> setCheckState ( Qt::Unchecked );
   idToCheckBox["checkB_mis"+QString::number(i)+"_3"] -> setCheckState ( Qt::Unchecked );
  } 
  else if(int(par[i][4]+0.5)==1){
   idToCheckBox["checkB_mis"+QString::number(i)+"_1"] -> setCheckState ( Qt::Checked );
   idToCheckBox["checkB_mis"+QString::number(i)+"_2"] -> setCheckState ( Qt::Checked );
   idToCheckBox["checkB_mis"+QString::number(i)+"_3"] -> setCheckState ( Qt::Unchecked );
  } 
  else if(int(par[i][4]+0.5)==2){
   idToCheckBox["checkB_mis"+QString::number(i)+"_1"] -> setCheckState ( Qt::Checked );
   idToCheckBox["checkB_mis"+QString::number(i)+"_2"] -> setCheckState ( Qt::Checked );
   idToCheckBox["checkB_mis"+QString::number(i)+"_3"] -> setCheckState ( Qt::Checked );
  }
  nMis=int(par[20+JJ][4]+0.5)-1;
  if(nMis>=0) idToComboBox["cBmis"+QString::number(i)] -> setCurrentIndex(nMis);
  if(JJ>=7){
   Jitem=idToComboBox["cBteE"+QString::number(JJ-6)] -> findText(QString::number(par[13+JJ-6][1]),Qt::MatchExactly);
   if(Jitem>=0) idToComboBox["cBteE"+QString::number(JJ-6)] ->setCurrentIndex(Jitem);
   i++;
   if(int(par[i][4]+0.5)==-1){
    idToCheckBox["checkB_mis"+QString::number(i)+"_2"] -> setCheckState ( Qt::Checked );
    idToCheckBox["checkB_mis"+QString::number(i)+"_3"] -> setCheckState ( Qt::Unchecked );
   } 
   else if(int(par[i][4]+0.5)==0){
    idToCheckBox["checkB_mis"+QString::number(i)+"_2"] -> setCheckState ( Qt::Unchecked );
    idToCheckBox["checkB_mis"+QString::number(i)+"_3"] -> setCheckState ( Qt::Unchecked );
   } 
   else if(int(par[i][4]+0.5)==1){
    idToCheckBox["checkB_mis"+QString::number(i)+"_2"] -> setCheckState ( Qt::Checked );
    idToCheckBox["checkB_mis"+QString::number(i)+"_3"] -> setCheckState ( Qt::Unchecked );
   } 
   else if(int(par[i][4]+0.5)==2){
    idToCheckBox["checkB_mis"+QString::number(i)+"_2"] -> setCheckState ( Qt::Checked );
    idToCheckBox["checkB_mis"+QString::number(i)+"_3"] -> setCheckState ( Qt::Checked );
   }
  }
  JJ++;
 }
 ui->dSB_PAR_3_3 -> setValue(par[3][3]);
 ui->dSB_PAR_4_1 -> setValue(par[4][1]);
 ui->dSB_PAR_4_2 -> setValue(par[4][2]);
 ui->dSB_PAR_4_3 -> setValue(par[4][3]);
 ui->dSB_PAR_5_3 -> setValue(par[5][3]);
 ui->dSB_PAR_6_1 -> setValue(par[6][1]);
 ui->sB_PAR_8_1 -> setValue(int(par[8][1]+.5));
 if(par[8][2]>par[8][1])
     par[8][2]=par[8][1];
 ui->sB_PAR_8_2 -> setValue(int(par[8][2]+.5));
 if(int(par[9][1]+0.5)==1)
     ui->checkBox_setPsoK -> setCheckState ( Qt::Checked );
 else
     ui->checkBox_setPsoK -> setCheckState ( Qt::Unchecked );
 if(int(par[10][1]+0.5)==0)
     ui->cBox_PAR_10_1 -> setCheckState (Qt::Unchecked);
 else
     ui->cBox_PAR_10_1 -> setCheckState (Qt::Checked);
 ui->LEpar_12_2 -> setText(QString::number(par[12][2]));
 ui->LEpar_13_1 -> setText(QString::number(par[13][1]));
 ui->LEpar_13_2 -> setText(QString::number(par[13][2]));
 int irif=int(par[10][3]+.5);
   ui->comboB_PAR_10_3 -> setCurrentIndex(irif);
 int iverboso=int(par[18][1]+.5);
   if(iverboso==1)
    ui->checkB_par_18_1 -> setCheckState ( Qt::Checked );
   else
    ui->checkB_par_18_1 -> setCheckState ( Qt::Unchecked );
 int isv=int(par[18][2]+0.5);
 ui->comboBox_StoreNK -> setCurrentIndex(isv);
 ui->comboBox_store -> setCurrentIndex(isv);
 ui->sBenabledWL -> setValue(int(par[21][1]+0.5));
 ui->sBenabledNK -> setValue(int(par[21][2]+0.5));
 int imultiply=int(par[22][1]+.5);
   if(imultiply==0)
    ui->checkB_PAR_22_1 -> setCheckState ( Qt::Unchecked );
   else
    ui->checkB_PAR_22_1 -> setCheckState ( Qt::Checked );
 ui->DP_PAR_27_1  -> setText(QString::number(par[27][1]));
 int s1p2=int(par[27][2]+.5);
 ui->cB_PAR_27_2 -> setCurrentIndex(s1p2-1);
 ui->sB_PAR_28_1 -> setValue(int(par[28][1]+0.5));
 ui->sB_PAR_29_1 -> setValue(int(par[29][1]+0.5));
 ui->sBtotWL -> setValue(int(par[28][2]+0.5));
 ui->sBtotNK -> setValue(int(par[29][2]+0.5));
 for(int i=1;i<=5;i++){
  idToLineEdit["LEpar_30_"+QString::number(i)] -> setText(QString::number(par[30][i]));
 }
 ui->cB_PAR_31_1 -> setCurrentIndex(int(par[31][1]+0.5));
 ui->cB_PAR_31_2 -> setCurrentIndex(int(par[31][2]+0.5));
 ui->DP_PAR_31_3 -> setText(QString::number(par[31][3]));
 ui->DP_PAR_31_4 -> setText(QString::number(par[31][4]));
 if(int(par[31][5]+0.5)==0)
     ui->checkBox_logScale->setCheckState(Qt::Unchecked);
 else
     ui->checkBox_logScale->setCheckState(Qt::Checked);
 ui->cB_PAR_32_1 -> setCurrentIndex(int(par[32][1]+0.5));
 ui->cB_PAR_32_2 -> setCurrentIndex(int(par[32][2]+0.5));
 ui->DP_PAR_32_3 -> setText(QString::number(par[32][3]));
 ui->DP_PAR_32_4 -> setText(QString::number(par[32][4]));
 ui->DP_PAR_21_3 -> setText(QString::number(par[21][3]));
 ui->DP_PAR_22_3 -> setText(QString::number(par[22][3]));
 nlayer=int(par[51][2]+0.5);//N, layer
 int isimmetry=int(par[52][2]+.5);
   if(isimmetry==0)
    ui->checkB_PAR_52_2 -> setCheckState ( Qt::Unchecked );
   else
    ui->checkB_PAR_52_2 -> setCheckState ( Qt::Checked );
 int ihemi=int(par[54][2]+.5);
   if(ihemi==0)
    ui->checkB_PAR_54_2 -> setCheckState ( Qt::Unchecked );
   else
    ui->checkB_PAR_54_2 -> setCheckState ( Qt::Checked );
 ui->DP_PAR_55_2 -> setText(QString::number(par[55][2]));
 npp=int(par[34][5]+0.5);//N par fit
 ui->cB_PAR_35_2 -> setCurrentIndex(int(par[35][2]-1.+0.5));
 n=int(par[35][5]+0.5);//N par enabled for fit
 ui->sB_PAR_34_5 -> setValue(npp);
 ui->sB_PAR_35_5 -> setValue(n);
 ui->sB_free_deg -> setValue(npp-n);
 for(int i=1;i<=14;i++){
  state=idToCheckBox["checkB_mis"+QString::number(i)+"_2"] -> checkState();
  if(state==Qt::Checked){
   idToLineEdit["LEpar_"+QString::number(35+i)+"_1"] -> setText(QString::number(par[35+i][1]));
   idToLineEdit["LEpar_"+QString::number(35+i)+"_2"] -> setText(QString::number(par[35+i][2]));
   idToLineEdit["LEpar_"+QString::number(35+i)+"_3"] -> setText(QString::number(par[35+i][3]));
  }
  else{
   idToLineEdit["LEpar_"+QString::number(35+i)+"_1"] -> clear();
   idToLineEdit["LEpar_"+QString::number(35+i)+"_2"] -> clear();
   idToLineEdit["LEpar_"+QString::number(35+i)+"_3"] -> clear();
  }
 }
 ui->DPenabledNK_WLmin -> setText(QString::number(par[38][1]));
 ui->DPenabledNK_WLmax -> setText(QString::number(par[38][2]));
 ui->DPenabledWL_WLmin -> setText(QString::number(par[38][3]));
 ui->DPenabledWL_WLmax -> setText(QString::number(par[38][4]));

 ui->comboBox_average->setCurrentIndex(int(par[35][1]+0.5));//set std spectrum

 for(int i=1;i<=npp;i++)
   ppm[i]=int(par[35+i][5]+0.5);
 //PM 
 for(int i=1;i<=200;i++){
  for(int j=1;j<=5;j++){
   stream>> pm[i][j];
  }
  //printf("pm[%d][1-5]: %f %f %f %f %f\n",i,pm[i][1],pm[i][2],pm[i][3],pm[i][4],pm[i][5]);
 }
 //printf("\tReadSetting: pm[100][1]= %f\n",pm[100][1]);
 //PF
 for(int i=1;i<=21;i++){
  for(int j=1;j<=7;j++){
   stream>> pf[j][i];
  }
 }
 file.close();
 //inserimento valori cnk
 for(int i=1;i<=15;i++){
  i0=int(cnk[i][1]+0.5);
  if(i0 >=0 && i0<=17){
   idToComboBox["cB_cnk"+QString::number(i)+"a"] -> setCurrentIndex(i0);
   idToCheckBox["cB_EMA_"+QString::number(i)] -> setCheckState ( Qt::Unchecked );
  }
  else{
   i1=int(cnk[i][1]/1000.+0.5);
   i2=int((cnk[i][1]-i1*1000.)/10.+0.5);
   f2=cnk[i][1]-i1*1000.-i2*10.;
   idToComboBox["cB_cnk"+QString::number(i)+"a"] -> setCurrentIndex(i1);
   idToComboBox["cB_cnk"+QString::number(i)+"b"] -> setCurrentIndex(i2);
   idToDoubleSpinBox["dSB_cnk"+QString::number(i)] ->setValue(f2);
   idToCheckBox["cB_EMA_"+QString::number(i)] -> setCheckState ( Qt::Checked );
  }
  idToLineEdit["LEcnk"+QString::number(i)+"_2"] -> setText(QString::number(cnk[i][2]));
  idToLineEdit["LEcnk"+QString::number(i)+"_3"] -> setText(QString::number(cnk[i][3]));
 }
 //aggiornamento valori PanFit
 for(int i=1;i<=npp;i++){
  state=idToCheckBox["chBeParFit_"+QString::number(i)]-> checkState();
  int ip=int(ppm[i]+0.5);
  idToLineEdit["DPparFitV_"+QString::number(i)] -> setText(QString::number(pm[ip][1]));
  if( state == Qt::Checked ){
   idToLineEdit["DPparFitErr_"+QString::number(i)] -> setText(QString::number(pm[ip][4]));
   idToLineEdit["DPparFitGC_"+QString::number(i)] -> setText(QString::number(pm[ip][5]));
  } 
  else{
   idToLineEdit["DPparFitErr_"+QString::number(i)] -> setText(QString::number(0));
   idToLineEdit["DPparFitGC_"+QString::number(i)] -> setText(QString::number(0));
  }
 }
 
 i0=int(pm[100][1]+0.5);
 //printf("i0= %d\n",i0);
 ui->cB_cnk1a -> setCurrentIndex(i0);//fit option
 ui->sB_PAR_51_2 -> setValue(nlayer);
 occupyPF=0;
 SetModel(nlayer);
 listOsc();
 PanFitPar();
}


void ksemaw::setRifMir(){
    Qt::CheckState state;
    state=ui->checkB_PAR_22_1 -> checkState();
    ui->comboB_PAR_10_3-> setEnabled(state==Qt::Checked);
}



void ksemaw::AdjTheta(){
 double theta=ui->dSB_PAR_6_1 -> value();
 ui->dSB_PAR_6_1bis -> setValue(theta);
 ui->dSB_PAR_6_1tris -> setValue(theta);
 ui->dSB_PAR_6_1quater -> setValue(theta);
 ui->dSB_PAR_6_1quinto -> setValue(theta);
 theta=ui->dSB_PAR_14_1 -> value();
 ui->dSB_PAR_14_1bis -> setValue(theta);
 ui->dSB_PAR_14_1tris -> setValue(theta);
 theta=ui->dSB_PAR_15_1 -> value();
 ui->dSB_PAR_15_1bis -> setValue(theta);
 ui->dSB_PAR_15_1tris -> setValue(theta);
 theta=ui->dSB_PAR_16_1 -> value();
 ui->dSB_PAR_16_1bis -> setValue(theta);
 ui->dSB_PAR_16_1tris -> setValue(theta);
 theta=ui->dSB_PAR_17_1 -> value();
 ui->dSB_PAR_17_1bis -> setValue(theta);
 ui->dSB_PAR_17_1tris -> setValue(theta);
}



void ksemaw::SaveProject(){
 SaveSetting(-1);
 int ierr;
 occupyPF=1;
 QString command, subfnproject;
 subfnproject=ui->lineEdit_P -> text();
 if(subfnproject=="mate/aa999.9"){
  fnproject = QFileDialog::getSaveFileName(
          this,
          "Filename to save",
          pathroot,
          "Semaw Project (*.Spj)");
  if(!fnproject.contains(".Spj"))
   fnproject=fnproject+".Spj";
  subfnproject=fnproject.section(pathroot, 1, 1);
  ui->lineEdit_P -> setText(subfnproject);
//  ui->lineEdit_P -> setText(subfnproject.section('.', 0, 1));
 }
 else{
   if(!subfnproject.contains(".Spj")){
    fnproject=pathroot+subfnproject+".Spj";
    ui->lineEdit_P -> setText(subfnproject+".Spj");
   }
   else
    fnproject=pathroot+subfnproject;
 }
 QFile file(fnproject);
 if( file.exists() ){
   QMessageBox msgBox;
   msgBox.setText("The file already exist!");
   msgBox.setInformativeText("Do you want to save anyway?");
   msgBox.setStandardButtons(QMessageBox::Save | QMessageBox::Cancel);
   msgBox.setDefaultButton(QMessageBox::Save);
   int ret = msgBox.exec();
   switch (ret) {
    case QMessageBox::Save:
       // Save was clicked
       command="cp "+fileStore+" "+fnproject;
       ierr=system((command.toStdString()).c_str());
       if(ierr != 0)
           printf("Error copying project!!!\n");
       else
           printf("Project saved as %s\n",(fnproject.toStdString()).c_str());
       break;
    case QMessageBox::Cancel:
       // Cancel was clicked
       printf("as you requested file was not saved!\n");
       break;
   }
 }
 else{
   command="cp "+fileStore+" "+fnproject;
   ierr= system((command.toStdString()).c_str());
   if(ierr != 0) printf("Error copying project!!!\n");
   printf("Project saved as %s\n",(fnproject.toStdString()).c_str());
 }  
 occupyPF=0;
}


void ksemaw::SaveSetting(int iCall){
    if(occupyPF!=0) return;
    occupyPF=1;
    printf("-> SaveSetting su %s\n",(fileStore.toStdString()).c_str());
    Qt::CheckState state,state1,state2;
    QString stringa,lab,svalue;
    int itab;
    if(iCall<0)
        itab=ui->tabWidget -> currentIndex();
    else
        itab=iCall;
    if(itab==1){
        state=ui->checkBox_setPsoK -> checkState ();
        if( state == Qt::Checked )
            par[9][1]=1.;
        else
            par[9][1]=0.;
        for(int i=1;i<=9;i++){
            if(i <= nlayer){
                par[50+i][1]=idToComboBox["comB_PAR_5"+QString::number(i)+"_1"] -> currentIndex();
                par[50+i][1]++;
                par[50+i][3]=idToComboBox["comB_PAR_5"+QString::number(i)+"_3"] -> currentIndex();
                par[50+i][3]++;
                pm[50+i][1]=idToDoubleSpinBox["dSB_PM_"+QString::number(50+i)+"_1"] -> value();
                if(int(par[50+i][3]+0.5)==1){
                    pm[i][1]=idToDoubleSpinBox["dSB_PM_"+QString::number(i)+"_1"] -> value();
                    pm[i][1]=pm[i][1]/1.E-7;
                } else if(int(par[50+i][3]+0.5)==2){
                    pm[i][1]=idToDoubleSpinBox["dSB_PM_"+QString::number(i)+"_1"] -> value();
                } else if(int(par[50+i][3]+0.5)==3){
                    pm[i][1]=idToDoubleSpinBox["dSB_PM_"+QString::number(i)+"_1"] -> value();
                    pm[10+i][1]=idToDoubleSpinBox["dSB_PM_"+QString::number(10+i)+"_1"] -> value();
                    pm[20+i][1]=idToDoubleSpinBox["dSB_PM_"+QString::number(20+i)+"_1"] -> value();
                    pm[30+i][1]=idToDoubleSpinBox["dSB_PM_"+QString::number(30+i)+"_1"] -> value();
                    pm[40+i][1]=idToDoubleSpinBox["dSB_PM_"+QString::number(40+i)+"_1"] -> value();
                    pm[60+i][1]=idToDoubleSpinBox["dSB_PM_"+QString::number(60+i)+"_1"] -> value();
                } else if(int(par[50+i][3]+0.5)==4 || int(par[50+i][3]+0.5)==5){
                    pm[i][1]=idToDoubleSpinBox["dSB_PM_"+QString::number(i)+"_1"] -> value();
                    pm[10+i][1]=idToDoubleSpinBox["dSB_PM_"+QString::number(10+i)+"_1"] -> value();
                }
            }
        }
        //aggiornamento valori PanFit
        for(int i=1;i<=npp;i++){
         state=idToCheckBox["chBeParFit_"+QString::number(i)]-> checkState();
         int ip=int(ppm[i]+0.5);
         idToLineEdit["DPparFitV_"+QString::number(i)] -> setText(QString::number(pm[ip][1]));
         if( state == Qt::Checked ){
          idToLineEdit["DPparFitErr_"+QString::number(i)] -> setText(QString::number(pm[ip][4]));
          idToLineEdit["DPparFitGC_"+QString::number(i)] -> setText(QString::number(pm[ip][5]));
         }
         else{
          idToLineEdit["DPparFitErr_"+QString::number(i)] -> setText(QString::number(0));
          idToLineEdit["DPparFitGC_"+QString::number(i)] -> setText(QString::number(0));
         }
        }
    }
    else if(itab==2){
        int ioptFit=ui->cB_cnk1a -> currentIndex();
        if(ioptFit>0 && ioptFit<8){
            printf("-> SaveOsc pf[%d,i]\n",ioptFit);
            for(int i=1;i<=21;i++) pf[ioptFit][i]=0.;
            int j=2;
            Qt::CheckState state;
            for(int i=1;i<=20;i++){
                state=idToCheckBox["cBosc_"+QString::number(i)]-> checkState();
                if( state == Qt::Checked ) {
                    pf[ioptFit][j]=i;
                    pm[100+1+(i-1)*5][1]=idToComboBox["cBpm_"+QString::number(100+1+(i-1)*5)+"_1"] -> currentIndex();
                    pm[100+1+(i-1)*5][1]++;
                    pm[100+2+(i-1)*5][1]=idToLineEdit["LEpm_"+QString::number(100+2+(i-1)*5)+"_1"] -> text().toDouble();
                    pm[100+3+(i-1)*5][1]=idToLineEdit["LEpm_"+QString::number(100+3+(i-1)*5)+"_1"] -> text().toDouble();
                    pm[100+4+(i-1)*5][1]=idToLineEdit["LEpm_"+QString::number(100+4+(i-1)*5)+"_1"] -> text().toDouble();
                    pm[100+5+(i-1)*5][1]=idToLineEdit["LEpm_"+QString::number(100+5+(i-1)*5)+"_1"] -> text().toDouble();
                    j++;
                    printf("pm[%d][1]: %f %f %f %f %f\n",100+2+(i-1)*5,pm[100+1+(i-1)*5][1],
                            pm[100+2+(i-1)*5][1],pm[100+3+(i-1)*5][1],
                            pm[100+4+(i-1)*5][1],pm[100+5+(i-1)*5][1]);
                }
            }
            pf[ioptFit][1]=j-2;
        }
    }
    else if(itab==4){
        printf("-> savePanFit\n");
        int n,n1,n2,ip,jpf;
        QString Lj,Valore;
        Qt::CheckState state;
        n=17;
        npp=17;
        jpf=0;
        for(int j=1;j<=17;j++){
          Lj=idToComboBox["cBParFit_"+QString::number(j)] -> currentText();
          state=idToCheckBox["chBeParFit_"+QString::number(j)]-> checkState();
          if(Lj=="none"){
            n--;
            npp--;
          }
          else{
            n1=(Lj.at(1)).digitValue();
            if(Lj.at(2).isNumber()) {
              n1=n1*10;
              n1=n1+((Lj.at(2)).digitValue());
            }
            n2=0;
            ip=0;
            do{
              n2++;
            } while(!Lj.contains(ParFitLab[n2],Qt::CaseSensitive));
            if(1<=n2 && n2<=7)
              ip=10*(n2-1)+n1;
            else if(8<=n2 && n2<=11)
              ip=100+(n1-1)*5+(n2-7)+1;
               printf("SavePanFit:  n1=%d n2=%d ip=%d\n",n1,n2,ip);
            ppm[j]=ip;
            Valore=idToLineEdit["DPparFitV_"+QString::number(j)] -> text();
            pm[ip][1]=Valore.toDouble();
            if( state == Qt::Unchecked ){
              pm[ip][2]=0;
              n--;
            }
            else{
              jpf++;
              pm[ip][2]=jpf;
              pm[jpf][3]=ip;
            }
          }
        }
        int nCeck=ui->sB_PAR_34_5 -> value();
        if(nCeck!=npp){
            QMessageBox msgBox;
            msgBox.setText("Please select each individual fit-parameter before fitting!!!");
            msgBox.setStandardButtons(QMessageBox::Ok);
            msgBox.exec();
            occupyPF=0;
            return;
        }
        par[34][5]=npp;//N par fit
        par[35][5]=n;//N par enabled
        for(int i=1;i<=npp;i++){
            par[35+i][5]=ppm[i];
        }
        ui->sB_PAR_35_5 -> setValue(n);
        ui->sB_free_deg -> setValue(npp-n);
        ui->sB_PAR_34_5 -> setValue(npp);
        SetModel(nlayer);
        listOsc();
    }
    int ci,i0,i1,i2;
    double f2;
    QFile file(fileStore);
    if (!file.open(QIODevice::WriteOnly | QIODevice::Text)){
        occupyPF=0;
        return;
    }
    QTextStream out(&file);
    out << "iVspj=1" << "\n";
    stringa=ui->lineEdit_infoP -> text();
    out << stringa.simplified() << "\n";
    //file nk
    for(int i=1;i<=8;i++){
        stringa=idToLineEdit["lineEdit"+QString::number(i)] -> text();
        out << stringa << "\n";
    }
    stringa=ui->lineEdit_Fnk -> text();
    stringa=stringa.section('.', 0, 1);
    out << stringa << "\n";
    //media
    ci=ui->comboBox_average -> currentIndex();
    par[35][1]=ci;//N# std spectrum
//standard spectra for mean value
    QFile file1(fStdSpect);
    if (!file1.open (QIODevice::ReadOnly | QIODevice::Text)){
        occupyPF=0;
        return;
    }
    QTextStream stream1 ( &file1 );
    QString line;
    for(int i=0;i<=ci;i++){
        line = stream1.readLine();
        line = stream1.readLine();
    }
    out << line <<Qt::endl;
    //sampleName
    stringa=ui->lineEdit_sample-> text();
    out << stringa << "\n";
    //ellipsometric measurements
    lab=ui->cBmis7 -> currentText();
    state=ui->checkB_mis7_1->checkState();
    if(lab.isEmpty() || state==Qt::Unchecked)
        out << "mate/aa999.9" << "\n";
    else{
        out << stringa+"."+lab << "\n";
        ui->dSB_PAR_14_1-> setValue(ui->cBteE1 -> currentText().toDouble());
    }
    lab=ui->cBmis9 -> currentText();
    state=ui->checkB_mis9_1->checkState();
    if(lab.isEmpty() || state==Qt::Unchecked)
        out << "mate/aa999.9" << "\n";
    else{
        out << stringa+"."+lab << "\n";
        ui->dSB_PAR_15_1-> setValue(ui->cBteE2 -> currentText().toDouble());
    }
    lab=ui->cBmis11 -> currentText();
    state=ui->checkB_mis11_1->checkState();
    if(lab.isEmpty() || state==Qt::Unchecked)
        out << "mate/aa999.9" << "\n";
    else{
        out << stringa+"."+lab << "\n";
        ui->dSB_PAR_16_1-> setValue(ui->cBteE3 -> currentText().toDouble());
    }
    lab=ui->cBmis13 -> currentText();
    state=ui->checkB_mis13_1->checkState();
    if(lab.isEmpty() || state==Qt::Unchecked)
        out << "mate/aa999.9" << "\n";
    else{
        out << stringa+"."+lab << "\n";
        ui->dSB_PAR_17_1-> setValue(ui->cBteE4 -> currentText().toDouble());
    }
    //RXY
    state=ui->checkB_RXY_25_3 -> checkState();
    if(state==Qt::Unchecked)
        rxy[25][3]=1;
    else
        rxy[25][3]=2;
    state=ui->checkB_RXY_25_4 -> checkState();
    if(state==Qt::Unchecked)
        rxy[25][4]=1;
    else
        rxy[25][4]=2;
    for(int i=1;i<=30;i++){
        for(int j=1;j<=4;j++){
            if(idToLineEdit.contains("DP_RXY_"+QString::number(i)+"_"+QString::number(j))){
                svalue=idToLineEdit["DP_RXY_"+QString::number(i)+"_"+QString::number(j)] -> text();
                rxy[i][j]=svalue.toDouble();
            }
        }
    }
    //CNK
    for(int i=1;i<=15;i++){
        state=idToCheckBox["cB_EMA_"+QString::number(i)] -> checkState();
        if(state==Qt::Unchecked){
            i0=idToComboBox["cB_cnk"+QString::number(i)+"a"] -> currentIndex();
            cnk[i][1]=i0;
        }
        else{
            i1=idToComboBox["cB_cnk"+QString::number(i)+"a"] -> currentIndex();
            i2=idToComboBox["cB_cnk"+QString::number(i)+"b"] -> currentIndex();
            f2=idToDoubleSpinBox["dSB_cnk"+QString::number(i)] -> value();
            cnk[i][1]=i1*1000.+i2*10.+f2;
        }
        svalue=idToLineEdit["LEcnk"+QString::number(i)+"_2"] -> text();
        cnk[i][2]=svalue.toDouble();
        svalue=idToLineEdit["LEcnk"+QString::number(i)+"_3"] -> text();
        cnk[i][3]=svalue.toDouble();
    }
    //PAR
    //  setting measurements
    par[23][1]=-1.0; //no SF
    par[23][2]=-1.0; //no ELI
    par[24][1]=1.0;  //IUVIR=1
    int JJ=1;
    for(int i=1;i<=14;i++){
        state =idToCheckBox["checkB_mis"+QString::number(i)+"_1"] -> checkState ();
        state1=idToCheckBox["checkB_mis"+QString::number(i)+"_2"] -> checkState ();
        state2=idToCheckBox["checkB_mis"+QString::number(i)+"_3"] -> checkState ();
        if(state==Qt::Unchecked && state1==Qt::Unchecked)
            par[i][4]=0;
        else if(state==Qt::Unchecked && state1==Qt::Checked)
            par[i][4]=-1;
        else if(state==Qt::Checked && state2==Qt::Unchecked )
            par[i][4]=1;
        else if(state==Qt::Checked && state2==Qt::Checked)
            par[i][4]=2;
        if(i<=6 && state==Qt::Checked && par[23][1]<0.)
            par[23][1]=1.0;
        if(i >6 && state==Qt::Checked && par[23][2]<0.)
            par[23][2]=1.0;
        if(i >6 && state==Qt::Checked){
            lab=idToComboBox["cBmis"+QString::number(i)] ->currentText();
            if(lab.at(0).toLatin1() == 'i') par[24][1]=2.0;//SF IR
        }
        par[20+JJ][4]=idToComboBox["cBmis"+QString::number(i)] -> currentIndex();
        par[20+JJ][4]++;
        if(JJ>=7){
            i++;
            state1=idToCheckBox["checkB_mis"+QString::number(i)+"_2"] -> checkState ();
            state2=idToCheckBox["checkB_mis"+QString::number(i)+"_3"] -> checkState ();
            if(state==Qt::Unchecked && state1==Qt::Unchecked)
                par[i][4]=0;
            else if(state==Qt::Unchecked && state1==Qt::Checked)
                par[i][4]=-1;
            else if(state==Qt::Checked && state2==Qt::Unchecked )
                par[i][4]=1;
            else if(state==Qt::Checked && state2==Qt::Checked)
                par[i][4]=2;
        }
        JJ++;
    }
    par[3][3]=ui->dSB_PAR_3_3 -> value();
    par[4][1]=ui->dSB_PAR_4_1 -> value();
    par[4][2]=ui->dSB_PAR_4_2 -> value();
    //if(par[4][1]>rxy[20][1]) rxy[20][1]=par[4][1];
    //if(par[4][2]<rxy[20][2]) rxy[20][2]=par[4][2];
    par[4][3]=ui->dSB_PAR_4_3 -> value();
    par[5][3]=ui->dSB_PAR_5_3 -> value();
    par[6][1]=ui->dSB_PAR_6_1 -> value();
    par[8][1]=ui->sB_PAR_8_1 -> value();
    par[8][2]=ui->sB_PAR_8_2 -> value();
    state=ui->cBox_PAR_10_1 -> checkState ();
    if(state==Qt::Unchecked)
        par[10][1]=0;
    else
        par[10][1]=1;//plot eps1 eps2
    par[10][3]=ui->comboB_PAR_10_3 -> currentIndex();
    par[14][1]=ui->dSB_PAR_14_1 -> value();
    par[15][1]=ui->dSB_PAR_15_1 -> value();
    par[16][1]=ui->dSB_PAR_16_1 -> value();
    par[17][1]=ui->dSB_PAR_17_1 -> value();
    state=ui->checkB_par_18_1 -> checkState ( );//verboso
    if(state==Qt::Unchecked)
        par[18][1]=0;
    else
        par[18][1]=1;
    par[18][2]=double(ui->comboBox_StoreNK -> currentIndex());
    state=ui->checkB_PAR_22_1 -> checkState ();//moltiplicazione per Rrif
    if(state==Qt::Unchecked)
        par[22][1]=0;
    else
        par[22][1]=1;
    par[27][2]=ui->cB_PAR_27_2 -> currentIndex();//S1P2
    par[27][2]++;
    par[28][1]=ui->sB_PAR_28_1 -> value();
    par[29][1]=ui->sB_PAR_29_1 -> value();
    svalue=ui->DP_PAR_31_3 -> text();
    par[31][1]=ui->cB_PAR_31_1 -> currentIndex();
    par[31][2]=ui->cB_PAR_31_2 -> currentIndex();
    par[31][3]=svalue.toDouble();
    svalue=ui->DP_PAR_31_4 -> text();
    par[31][4]=svalue.toDouble();
    state=ui->checkBox_logScale -> checkState();
    if(state==Qt::Unchecked)
        par[31][5]=0.;
    else{
        par[31][5]=1.;
        if(rxy[17][1]<=0.)
            rxy[17][1]=rxy[17][2]/10000.;
    }
    par[32][1]=ui->cB_PAR_32_1 -> currentIndex();
    par[32][2]=ui->cB_PAR_32_2 -> currentIndex();
    svalue=ui->DP_PAR_32_3 -> text();
    par[32][3]=svalue.toDouble();
    svalue=ui->DP_PAR_32_4 -> text();
    par[32][4]=svalue.toDouble();
    svalue=ui->DP_PAR_21_3 -> text();
    par[21][3]=svalue.toDouble();
    if(fabs(par[21][3]) < 0.003*fabs(rxy[16][2]-rxy[16][1]))
        par[21][3]=0.003*fabs(rxy[16][2]-rxy[16][1]);
    svalue=ui->DP_PAR_22_3 -> text();
    par[22][3]=svalue.toDouble();
    if(fabs(par[22][3]) < 0.003*fabs(rxy[17][2]-rxy[17][1]))
        par[22][3]=0.003*fabs(rxy[17][2]-rxy[17][1]);
    par[35][4]=cnk[1][1];//usato in simula
    //  setting the model
    nlayer=ui->sB_PAR_51_2 -> value();
    par[51][2]=nlayer;
    state=ui->checkB_PAR_52_2 -> checkState ();
    if(state==Qt::Unchecked)
        par[52][2]=0;
    else
        par[52][2]=1;
    state=ui->checkB_PAR_54_2 -> checkState ();
    if(state==Qt::Unchecked)
        par[54][2]=0;
    else
        par[54][2]=1;
    for(int i=1;i<=9;i++){
        if(i<=nlayer){
            par[50+i][1]=idToComboBox["comB_PAR_5"+QString::number(i)+"_1"] -> currentIndex();
            par[50+i][1]++;
            par[50+i][3]=idToComboBox["comB_PAR_5"+QString::number(i)+"_3"] -> currentIndex();
            par[50+i][3]++;
            if(int(par[50+i][1]+0.5)==1) par[53][2]=i;//puntatore strato incognito
        }
        else{
            par[50+i][1]=0;
            par[50+i][3]=0;
        }
    }
    svalue=ui->DP_PAR_55_2 -> text();
    par[55][2]=svalue.toDouble();
    par[34][5]=npp;
    par[35][2]=ui->cB_PAR_35_2 ->currentIndex();
    par[35][2]++;
    par[35][5]=n;
    par[51][2]=nlayer;

    pm[100][1]=ui->cB_cnk1a -> currentIndex();//fit option

    //save rxy
    for(int i=1;i<=30;i++){
        for(int j=1;j<=4;j++){
            out << rxy[i][j] << "\t";
        }
        out << "\n";
    }
    //save cnk
    for(int i=1;i<=15;i++){
        for(int j=1;j<=3;j++){
            out << QString::number(cnk[i][j],  'f',  5) << "\t";
        }
        out << "\n";
    }
    //save par
    for(int i=1;i<=60;i++){
        for(int j=1;j<=5;j++){
            out << par[i][j] << "\t";
        }
        out << "\n";
    }
    //save pm
    for(int i=1;i<=200;i++){
        for(int j=1;j<=5;j++){
            out << pm[i][j] << "\t";
        }
        out << "\n";
    }
    //save pf
    for(int i=1;i<=21;i++){
        for(int j=1;j<=7;j++){
            out << pf[j][i] <<"\t";
        }
        out << "\n";
    }
    file.close();
    occupyPF=0;
}

void ksemaw::LoadFilenk(){
 SaveSetting(-1);
// printf("ifn=%d fnFnk=%s\n",ifn,(fnFnk.toStdString()).c_str());
 if(ifn==0){
  fnFnk = QFileDialog::getOpenFileName(
	 this,
	 "Choose a file-nk", //titolo della finestra
	 pathroot, //directory iniziale
	 "file-nk (*.nk)"); //tipi di file da cercare
 }
 QFile file(fnFnk);
 if (!file.open (QIODevice::ReadOnly | QIODevice::Text))
      return;
 QTextStream stream ( &file );
 QString line,line2;
 line = stream.readLine();
 ui->lineEdit_infoFnk -> setText(line.simplified());
 file.close();
 line2=fnFnk.section(pathroot, 1, 1);
 ui->lineEdit_Fnk -> setText(line2);
// ui->lineEdit_Fnk -> setText(line2.section('.', 0, 1));
 //printf("LoadFilenk: occupyPF = %d\n",occupyPF);
 if(occupyPF!=1) {
  SaveSetting(-1);
  QFile file2(fctrl);
  if (!file2.open(QIODevice::WriteOnly | QIODevice::Text))
   return;
  QTextStream out(&file2);
  out << "iofnk" << "\n";
  out << "2" << "\n";
  file2.close();
 }
}

void ksemaw::ClrFnk(){
 ui->lineEdit_Fnk -> setText("mate/aa999.9");
 ui->lineEdit_infoFnk ->setText("");
// inviare comando "0"
}


void ksemaw::SaveFnk(){
 SaveSetting(-1);
 QString line;
 line=ui->lineEdit_Fnk -> text();
 if(line=="mate/aa999.9"){
   fnFnk = QFileDialog::getSaveFileName(
          this,
          "Filename to save",
          pathroot,
          "nk file (*.nk)");
   if(!fnFnk.contains(".nk"))
    fnFnk=fnFnk+".nk";
  line=fnFnk.section(pathroot, 1, 1);
  ui->lineEdit_Fnk -> setText(line);
//  ui->lineEdit_Fnk -> setText(line.section('.', 0, 1));
 }
 else{
  if(!line.contains(".nk")){
      occupyPF=1;
      line=line+".nk";
      ui->lineEdit_Fnk -> setText(line);
      usleep(500);
      occupyPF=0;
  }
  fnFnk=pathroot+line;
 }
 QFile file(fnFnk);
 int iok=1;
 if( file.exists() ){
   QMessageBox msgBox;
   msgBox.setText("The file already exist!");
   msgBox.setInformativeText("Do you want to save anyway?");
   msgBox.setStandardButtons(QMessageBox::Save | QMessageBox::Cancel);
   msgBox.setDefaultButton(QMessageBox::Save);
   int ret = msgBox.exec();
   switch (ret) {
    case QMessageBox::Save:
       // Save was clicked
       iok=1;
       break;
    case QMessageBox::Cancel:
       // Cancel was clicked
       iok=0;
       break;
   }
 }
 if(iok==1){
  SaveSetting(-1);
  QFile file(fctrl);
  if (!file.open(QIODevice::WriteOnly | QIODevice::Text))
  return;
  QTextStream out(&file);
  out << "iofnk" << "\n";
  out << "1" << "\n";
  line=ui->lineEdit_infoFnk -> text();
  out << line << "\n";
  file.close();
 }
}


void ksemaw::Setnk(int ifile){
 if(ifn==0){
  fnk[ifile]=QFileDialog::getOpenFileName(
	 this,
	 "Choose a file-nk", //titolo della finestra
	 pathroot, //directory iniziale
	 "file-nk (*.nk)"); //tipi di file da cercare
 }
 QFile file(fnk[ifile]);
 if (!file.open (QIODevice::ReadOnly | QIODevice::Text))
      return;
 double val;
 QTextStream stream ( &file );
 QString line,line2;
 line = stream.readLine();
 idToLineEdit["lineEdit_infoNK_"+QString::number(ifile)] -> setText(line.simplified());
 line2=fnk[ifile].section(pathroot, 1, 1);
 idToLineEdit["lineEdit"+QString::number(ifile)] -> setText(line2.section('.', 0, 1));
 line = stream.readLine();
 line = stream.readLine();
 line=line.simplified();
 val=line.section(' ', 0, 0).toDouble();
 idToLineEdit["WLminNK"+QString::number(ifile)] -> setText(QString::number(int(val+0.5)));
 line2=line;
 do {
  line=line2;
  line2 = stream.readLine();
 } while (!line2.isNull());
 line=line.simplified();
 val=line.section(' ', 0, 0).toDouble();
 idToLineEdit["WLmaxNK"+QString::number(ifile)] -> setText(QString::number(int(val+0.5)));
 file.close();
 MCRange();
 line2=idToLineEdit["lineEdit"+QString::number(ifile)] -> text();
 updateMatenk(7+ifile,line2);
}

void ksemaw::Setnk1(){
 Setnk(1);
}

void ksemaw::Setnk2(){
 Setnk(2);
}

void ksemaw::Setnk3(){
 Setnk(3);
}

void ksemaw::Setnk4(){
 Setnk(4);
}

void ksemaw::Setnk5(){
 Setnk(5);
}

void ksemaw::Setnk6(){
 Setnk(6);
}

void ksemaw::Setnk7(){
 Setnk(7);
}

void ksemaw::Setnk8(){
 Setnk(8);
}

void ksemaw::mDw1(){
    mDwUp(1,1);
}

void ksemaw::mUp1(){
    mDwUp(1,-1);
}

void ksemaw::mDw2(){
    mDwUp(2,1);
}

void ksemaw::mUp2(){
    mDwUp(2,-1);
}

void ksemaw::mDw3(){
    mDwUp(3,1);
}

void ksemaw::mUp3(){
    mDwUp(3,-1);
}

void ksemaw::mDw4(){
    mDwUp(4,1);
}

void ksemaw::mUp4(){
    mDwUp(4,-1);
}

void ksemaw::mDw5(){
    mDwUp(5,1);
}

void ksemaw::mUp5(){
    mDwUp(5,-1);
}

void ksemaw::mDw6(){
    mDwUp(6,1);
}

void ksemaw::mUp6(){
    mDwUp(6,-1);
}

void ksemaw::mDw7(){
    mDwUp(7,1);
}

void ksemaw::mUp7(){
    mDwUp(7,-1);
}

void ksemaw::mDw8(){
    mDwUp(8,1);
}

void ksemaw::mUp8(){
    mDwUp(8,-1);
}

void ksemaw::mDw9(){
    mDwUp(9,1);
}

void ksemaw::mUp9(){
    mDwUp(9,-1);
}

void ksemaw::mDwUp(int iLayer, int Dw1UpM1){
    double tmp1=par[50+iLayer][1];
    double tmp3=par[50+iLayer][3];
    par[50+iLayer][1]=par[50+iLayer+Dw1UpM1][1];
    par[50+iLayer][3]=par[50+iLayer+Dw1UpM1][3];
    par[50+iLayer+Dw1UpM1][1]=tmp1;
    par[50+iLayer+Dw1UpM1][3]=tmp3;
    for(int j=0;j<=6;j++){
        for(int k=1;k<=5;k++){
            double tmp=pm[10*j+iLayer][k];
            pm[10*j+iLayer][k]=pm[10*j+iLayer+Dw1UpM1][k];
            pm[10*j+iLayer+Dw1UpM1][k]=tmp;
        }
    }
    SetModel(nlayer);
}


void ksemaw::Clrnk(int ifile){
 idToLineEdit["lineEdit"+QString::number(ifile)] -> setText("mate/aa999.9");
 idToLineEdit["lineEdit_infoNK_"+QString::number(ifile)] ->setText("");
 idToLineEdit["WLminNK"+QString::number(ifile)] -> setText("");
 idToLineEdit["WLmaxNK"+QString::number(ifile)] -> setText("");
 MCRange();
 updateMatenk(7+ifile,"nk-"+QString::number(ifile));
}

void ksemaw::Clrnk1(){
 Clrnk(1);
}

void ksemaw::Clrnk2(){
 Clrnk(2);
}

void ksemaw::Clrnk3(){
 Clrnk(3);
}

void ksemaw::Clrnk4(){
 Clrnk(4);
}

void ksemaw::Clrnk5(){
 Clrnk(5);
}

void ksemaw::Clrnk6(){
 Clrnk(6);
}

void ksemaw::Clrnk7(){
 Clrnk(7);
}

void ksemaw::Clrnk8(){
 Clrnk(8);
}

void ksemaw::setSample(){
 QString line2;
 if(ifn==0){
  fnSample=QFileDialog::getOpenFileName(
        this,
        "Choose a sample", //titolo della finestra
        pathroot, //directory iniziale
        "file-exp (*.tn *.tp *.rn *.rp *.r1 *.an *.el)"); //tipi di file da cercare
 }
 fnSample=fnSample.section('.', 0, 0);
 line2=fnSample.section(pathroot, 1, 1);
 ui->lineEdit_sample -> setText(line2);
}

void ksemaw::listMeas(const QString &){
 int i,j;
 QString lab,estens[7];
 estens[1]=".tn";
 estens[2]=".tp";
 estens[3]=".rn";
 estens[4]=".rp";
 estens[5]=".r1";
 estens[6]=".an";
 
 for(int imis=1;imis<=6;imis++){
  idToComboBox["cBmis"+QString::number(imis)] -> clear();
  for(j=0;j<2;j++){
   lab="v";
   if(j==1) lab="i";
   for(i=0;i<10;i++){
    QFile file(fnSample+"."+lab+QString::number(i)+estens[imis]);
    if(file.exists()) {
     idToComboBox["cBmis"+QString::number(imis)] -> addItem(lab+QString::number(i));
     file.close();
    }
   }
  }
 }
 
 ui->cBmis7 -> clear();
 ui->cBmis9 -> clear();
 ui->cBmis11 -> clear();
 ui->cBmis13 -> clear();
 for(i=0;i<10;i++){
  QFile file(fnSample+"."+QString::number(i)+".el");
  if(file.exists()) {
    ui->cBmis7 -> addItem(QString::number(i));
    ui->cBmis9 -> addItem(QString::number(i));
    ui->cBmis11 -> addItem(QString::number(i));
    ui->cBmis13 -> addItem(QString::number(i));
    file.close();
  }
 }
 
}

void ksemaw::pwTn(const int &){
 double wmin=0.,wmax=0.;
 QString lab,info;
 lab=ui->cBmis1 -> currentText();
 fnTn=fnSample+"."+lab+".tn";
 previewFile(fnTn,lab,info,wmin,wmax);
 ui->lineEdit_Tn -> setText(info);
 if(wmin>0. && wmax>0.){
  ui->WLmin1 -> setText(QString::number(wmin));
  ui->WLmax1 -> setText(QString::number(wmax));
  }
 else{
  ui-> WLmin1 -> setText("");
  ui->WLmax1 -> setText("");
 }
 MCRange();
}

void ksemaw::pwTp(const int &){
 double wmin=0.,wmax=0.;
 QString lab,info;
 lab=ui->cBmis2 -> currentText();
 fnTp=fnSample+"."+lab+".tp";
 previewFile(fnTp,lab,info,wmin,wmax);
 ui->lineEdit_Tp -> setText(info);
 if(wmin>0. && wmax>0.){
  ui->WLmin2 -> setText(QString::number(wmin));
  ui->WLmax2 -> setText(QString::number(wmax));
  }
 else{
  ui->WLmin2 -> setText("");
  ui->WLmax2 -> setText("");
 }
 MCRange();
}

void ksemaw::pwRn(const int &){
 double wmin=0.,wmax=0.;
 QString lab,info;
 lab=ui->cBmis3 -> currentText();
 fnRn=fnSample+"."+lab+".rn";
 previewFile(fnRn,lab,info,wmin,wmax);
 ui->lineEdit_Rn -> setText(info);
 if(wmin>0. && wmax>0.){
  ui->WLmin3 -> setText(QString::number(wmin));
  ui->WLmax3 -> setText(QString::number(wmax));
  }
 else{
  ui->WLmin3 -> setText("");
  ui->WLmax3 -> setText("");
 }
 MCRange();
}

void ksemaw::pwRp(const int &){
 double wmin=0.,wmax=0.;
 QString lab,info;
 lab=ui->cBmis4 -> currentText();
 fnRp=fnSample+"."+lab+".rp";
 previewFile(fnRp,lab,info,wmin,wmax);
 ui->lineEdit_Rp -> setText(info);
 if(wmin>0. && wmax>0.){
  ui->WLmin4 -> setText(QString::number(wmin));
  ui->WLmax4 -> setText(QString::number(wmax));
 }
 else{
  ui->WLmin4 -> setText("");
  ui->WLmax4 -> setText("");
 }
 MCRange();
}

void ksemaw::pwR1(const int &){
 double wmin=0.,wmax=0.;
 QString lab,info;
 lab=ui->cBmis5 -> currentText();
 fnR1=fnSample+"."+lab+".r1";
 previewFile(fnR1,lab,info,wmin,wmax);
 ui->lineEdit_R1 -> setText(info);
 if(wmin>0. && wmax>0.){
  ui->WLmin5 -> setText(QString::number(wmin));
  ui->WLmax5 -> setText(QString::number(wmax));
 }
 else{
  ui->WLmin5 -> setText("");
  ui->WLmax5 -> setText("");
 }
 MCRange();
}

void ksemaw::pwApds(const int &){
 double wmin=0.,wmax=0.;
 QString lab,info;
 lab=ui->cBmis6 -> currentText();
 fnApds=fnSample+"."+lab+".an";
 previewFile(fnApds,lab,info,wmin,wmax);
 ui->lineEdit_Apds -> setText(info);
 if(wmin>0. && wmax>0.){
  ui->WLmin6 -> setText(QString::number(wmin));
  ui->WLmax6 -> setText(QString::number(wmax));
 }
 else{
  ui->WLmin6 -> setText("");
  ui->WLmax6 -> setText("");
 }
 MCRange();
}

void ksemaw::pwE1(const int &){
 int ntel=0,ndat=0,i,index;
 double teta=0;
 double wmin=0.,wmax=0.;
 ui->cBteE1 -> clear();
 QString lab,info;
 index=ui->cBmis7 ->currentIndex();
 if(index>=0){
  lab=ui->cBmis7 -> currentText();
  fnE1=fnSample+"."+lab+".el";
  QFile file(fnE1);
  if (!file.open (QIODevice::ReadOnly | QIODevice::Text))
       return;
  else if(file.exists()) {
   QTextStream stream ( &file );
   info = stream.readLine();
   ui->lineEdit_E1 -> setText(info);
   stream >> ntel >> ndat >> info;
   for(i=0;i<ntel;i++){
    stream >> teta >> ndat >> wmin >> wmax;
    ui->cBteE1 -> addItem(QString::number(teta));
    ui->WLmin7 -> setText(QString::number(wmin*10.));
    ui->WLmax7 -> setText(QString::number(wmax*10.));
   }
   file.close();
  }
 }
 else{
  fnE1="";
  ui->lineEdit_E1 -> clear();
 }
}

void ksemaw::pwE2(const int &){
 int ntel=0,ndat=0,i,index;
 double teta=0;
 double wmin=0.,wmax=0.;
 ui->cBteE2 -> clear();
 QString lab,info;
 index=ui->cBmis9 ->currentIndex();
 if(index>=0){
  lab=ui->cBmis9 -> currentText();
  fnE2=fnSample+"."+lab+".el";
  QFile file(fnE2);
  if (!file.open (QIODevice::ReadOnly | QIODevice::Text))
       return;
  else if(file.exists()) {
   QTextStream stream ( &file );
   info = stream.readLine();
   ui->lineEdit_E2 -> setText(info);
   stream >> ntel >> ndat >> info;
   for(i=0;i<ntel;i++){
    stream >> teta >> ndat >> wmin >> wmax;
    ui->cBteE2 -> addItem(QString::number(teta));
    ui->WLmin9 -> setText(QString::number(wmin*10.));
    ui->WLmax9 -> setText(QString::number(wmax*10.));
   }
   file.close();
  }
 }
 else{
  fnE2="";
  ui->lineEdit_E2 -> clear();
 }
}

void ksemaw::pwE3(const int &){
 int ntel=0,ndat=0,i,index;
 double teta=0;
 double wmin=0.,wmax=0.;
 ui->cBteE3 -> clear();
 QString lab,info;
 index=ui->cBmis11 ->currentIndex();
 if(index>=0){
  lab=ui->cBmis11 -> currentText();
  fnE3=fnSample+"."+lab+".el";
  QFile file(fnE3);
  if (!file.open (QIODevice::ReadOnly | QIODevice::Text))
       return;
  else if(file.exists()) {
   QTextStream stream ( &file );
   info = stream.readLine();
   ui->lineEdit_E3 -> setText(info);
   stream >> ntel >> ndat >> info;
   for(i=0;i<ntel;i++){
    stream >> teta >> ndat >> wmin >> wmax;
    ui->cBteE3 -> addItem(QString::number(teta));
    ui->WLmin11 -> setText(QString::number(wmin*10.));
    ui->WLmax11 -> setText(QString::number(wmax*10.));
   }
   file.close();
  }
 }
 else{
  fnE3="";
  ui->lineEdit_E3 -> clear();
 }
}

void ksemaw::pwE4(const int &){
 int ntel=0,ndat=0,i,index;
 double teta=0;
 double wmin=0.,wmax=0.;
 ui->cBteE4 -> clear();
 QString lab,info;
 index=ui->cBmis13 ->currentIndex();
 if(index>=0){
  lab=ui->cBmis13 -> currentText();
  fnE4=fnSample+"."+lab+".el";
  QFile file(fnE4);
  if (!file.open (QIODevice::ReadOnly | QIODevice::Text))
       return;
  else if(file.exists()) {
   QTextStream stream ( &file );
   info = stream.readLine();
   ui->lineEdit_E4 -> setText(info);
   stream >> ntel >> ndat >> info;
   for(i=0;i<ntel;i++){
    stream >> teta >> ndat >> wmin >> wmax;
    ui->cBteE4 -> addItem(QString::number(teta));
    ui->WLmin13 -> setText(QString::number(wmin*10.));
    ui->WLmax13 -> setText(QString::number(wmax*10.));
   }
   file.close();
  }
 }
 else{
  fnE4="";
  ui->lineEdit_E4 -> clear();
 }
}

void ksemaw::pwSubE1(const int &){
 int index,i,ndat;
 double teta=0;
 double wmin=0.,wmax=0.;
 QString line;
 index=ui->cBteE1 -> currentIndex();
 if(index>=0){
  QFile file(fnE1);
  if (!file.open (QIODevice::ReadOnly | QIODevice::Text))
       return;
  else if(file.exists()) {
   QTextStream stream ( &file );
   for(i=0;i<index+2;i++) line = stream.readLine();
   stream >> teta >> ndat >> wmin >> wmax;
   wmin=wmin*10.;
   wmax=wmax*10.;
  }
 }
 if(wmin>0. && wmax>0.){
  ui->WLmin7 -> setText(QString::number(wmin));
  ui->WLmax7 -> setText(QString::number(wmax));
 }
 else{
  ui->WLmin7 -> setText("");
  ui->WLmax7 -> setText("");
 }
 MCRange();
}

void ksemaw::pwSubE2(const int &){
 int index,i,ndat;
 float teta=0;
 double wmin=0.,wmax=0.;
 QString line;
 index=ui->cBteE2 -> currentIndex();
 if(index>=0){
  QFile file(fnE2);
  if (!file.open (QIODevice::ReadOnly | QIODevice::Text))
       return;
  else if(file.exists()) {
   QTextStream stream ( &file );
   for(i=0;i<index+2;i++) line = stream.readLine();
   stream >> teta >> ndat >> wmin >> wmax;
   wmin=wmin*10.;
   wmax=wmax*10.;
  }
 }
 if(wmin>0. && wmax>0.){
  ui->WLmin9 -> setText(QString::number(wmin));
  ui->WLmax9 -> setText(QString::number(wmax));
 }
 else{
  ui->WLmin9 -> setText("");
  ui->WLmax9 -> setText("");
 }
 MCRange();
}

void ksemaw::pwSubE3(const int &){
 int index,i,ndat;
 float teta=0;
 double wmin=0.,wmax=0.;
 QString line;
 index=ui->cBteE3 -> currentIndex();
 if(index>=0){
  QFile file(fnE3);
  if (!file.open (QIODevice::ReadOnly | QIODevice::Text))
       return;
  else if(file.exists()) {
   QTextStream stream ( &file );
   for(i=0;i<index+2;i++) line = stream.readLine();
   stream >> teta >> ndat >> wmin >> wmax;
   wmin=wmin*10.;
   wmax=wmax*10.;
  }
 }
 if(wmin>0. && wmax>0.){
  ui->WLmin11 -> setText(QString::number(wmin));
  ui->WLmax11 -> setText(QString::number(wmax));
 }
 else{
  ui->WLmin11 -> setText("");
  ui->WLmax11 -> setText("");
 }
 MCRange();
}

void ksemaw::pwSubE4(const int &){
 int index,i,ndat;
 float teta=0;
 double wmin=0.,wmax=0.;
 QString line;
 index=ui->cBteE4 -> currentIndex();
 if(index>=0){
  QFile file(fnE4);
  if (!file.open (QIODevice::ReadOnly | QIODevice::Text))
       return;
  else if(file.exists()) {
   QTextStream stream ( &file );
   for(i=0;i<index+2;i++) line = stream.readLine();
   stream >> teta >> ndat >> wmin >> wmax;
   wmin=wmin*10.;
   wmax=wmax*10.;
  }
 }
 if(wmin>0. && wmax>0.){
  ui->WLmin13 -> setText(QString::number(wmin));
  ui->WLmax13 -> setText(QString::number(wmax));
 }
 else{
  ui->WLmin13 -> setText("");
  ui->WLmax13 -> setText("");
 }
 MCRange();
}

void ksemaw::MCRange(){
 double Lmin=0.,Lmax=1.E+20,vmin,vmax;
 QString str;
 Qt::CheckState state1;
 for(int ink=1;ink<=8;ink++){
  str=idToLineEdit["WLminNK"+QString::number(ink)] -> text();
  vmin=str.toDouble();
  str=idToLineEdit["WLmaxNK"+QString::number(ink)] -> text();
  vmax=str.toDouble();
  if(vmin != 0.) Lmin=max(Lmin,vmin);
  if(vmax != 0.) Lmax=min(Lmax,vmax);
 }
 for(int imis=1;imis<=14;imis++){
  state1 = idToCheckBox["checkB_mis"+QString::number(imis)+"_1"] -> checkState();
  if( state1 == Qt::Checked ) {
   str=idToLineEdit["WLmin"+QString::number(imis)] -> text();
   vmin=str.toDouble();
   str=idToLineEdit["WLmax"+QString::number(imis)] -> text();
   vmax=str.toDouble();
   if(vmin != 0.) Lmin=max(Lmin,vmin);
   if(vmax != 0.) Lmax=min(Lmax,vmax);
  }
  if(imis>=7) imis++;
 }
 if(Lmin>0. && Lmax<1.E+20){
  ui->dSB_PAR_4_1 -> setValue(Lmin);
  ui->dSB_PAR_4_2 -> setValue(Lmax);
  //ui->DP_RXY_20_1 -> setText(QString::number(Lmin));
  //ui->DP_RXY_20_2 -> setText(QString::number(Lmax));
  ui->DP_RXY_20_3 -> setText(QString::number(Lmin));
  ui->DP_RXY_20_4 -> setText(QString::number(Lmax));
 }
 else{
  ui->dSB_PAR_4_1 -> cleanText();
  ui->dSB_PAR_4_2 -> cleanText();
  ui->DP_RXY_20_1 -> setText(" ");
  ui->DP_RXY_20_2 -> setText(" ");
  ui->DP_RXY_20_3 -> setText(" ");
  ui->DP_RXY_20_4 -> setText(" ");
 }
}

void ksemaw::updateMatenk(int index, QString newName){
 for(int i=1;i<=15;i++){
  idToComboBox["cB_cnk"+QString::number(i)+"a"]-> setItemText(index,newName);
  idToComboBox["cB_cnk"+QString::number(i)+"b"]-> setItemText(index,newName);
 }
}

void previewFile(QString filename, QString lab,QString& info,double& wmin,double& wmax){
 int i,idelta,ndati,ilinrim;
 double ridelta,dini,dfin,dmin,dmax,val;
 QString line,line0;
 QFile file(filename);
 if (!file.open (QIODevice::ReadOnly | QIODevice::Text))
      return;
 else if(file.exists()) {
  QTextStream stream ( &file );
  info = stream.readLine();
  line0 = stream.readLine();
//  printf("line0= %s \n",(line0.toStdString()).c_str());
  if(info.contains("PE UV", Qt::CaseInsensitive)){
// trattasi di file Perkin Elmer L900 - L950
    for(i=0;i<6;i++) line = stream.readLine();
    info = stream.readLine();
    for(i=0;i<5;i++) line = stream.readLine();
    if(line.contains("PerkinElmer UV WinLab 5", Qt::CaseInsensitive))
     ilinrim=65;
    else
     ilinrim=69;
    for(i=0;i<ilinrim;i++) line = stream.readLine();
    stream >> wmax;
    stream >> ridelta;
    if(ridelta > 0.)
     idelta=int(ridelta+0.5);
    else
     idelta=int(ridelta-0.5);
    stream >> ndati;
    wmin=wmax+double(idelta*(ndati-1));
    wmax=wmax*10.;
    wmin=wmin*10.;
//     printf("file Perkin Elmer L900 - L950 \n");
//     printf(" ridelta=%f ndati=%d wmin= %f wmax= %f \n",ridelta,ndati,wmin,wmax);
  }
  else if(line0.contains("#####SCALED") || line0.contains("#####SCALEDA") ||
          line0.contains("#####SCALED%")){
// trattasi di file scalato
    stream >> dini;
    stream >> dfin;
    if(dini < dfin){
      wmin=dini;
      wmax=dfin;
    }
    else{
      wmin=dfin;
      wmax=dini;
    }
    wmax=wmax*10.;
    wmin=wmin*10.;
  }
  else if(line0.contains("##XYUNITS= W")){
// trattasi di file l/9 o ftir vecchio tipo
    stream >> line >> dmax >> line >> val;
    stream >> line >> dmin >> line >> val;
//    printf("file vecchio tipo dmin= %f dmax= %f \n",dmin,dmax);
    if(lab.contains("v")){
     wmin=dmin*10.;
     wmax=dmax*10.;
    }
    else{
     wmin=(1.E+8)/dmax;
     wmax=(1.E+8)/dmin;
    }
  }
  else{
   line0 = stream.readLine();
   if(line0.contains("##DATA TYPE= UL")){
//          trattasi file l/9 nuovo tipo
    for(i=0;i<8;i++) line = stream.readLine();
    stream >> line >> dmin;
    stream >> line >> dmax;
    wmin=dmin*10.;
    wmax=dmax*10.;
   }
   else if(line0.contains("##DATA TYPE= UV")){
// trattasi file l/19
     for(i=0;i<11;i++) line = stream.readLine();
     stream >> line >> dmax;
     stream >> line >> dmin;
     wmin=dmin*10.;
     wmax=dmax*10.;
   }
   else if(line0.contains("##DATA TYPE= IN")){
    // trattasi file ftir nuovo tipo
     for(i=0;i<14;i++) line = stream.readLine();
     stream >> line >> dmin;
     stream >> line >> dmax;
     wmin=(1.E+8)/dmin;
     wmax=(1.E+8)/dmax;          
   }
  }
 }
// printf("wmin= %f wmax= %f \n",wmin,wmax);
}

void ksemaw::SetModel(const int &){
    int itab=ui->tabWidget -> currentIndex();
    if(occupyPF!=0 && itab!=4) return;
    printf("-> SetModel\n");
    occupyPF=1;
    int nlayerOld=nlayer;
    nlayer=ui->sB_PAR_51_2 -> value();
    int Dnl=nlayer-nlayerOld;
    if(Dnl>0){
        for(int i=0;i<nlayerOld;i++){
            par[50+nlayer-i][1]=par[50+nlayerOld-i][1];
            par[50+nlayer-i][3]=par[50+nlayerOld-i][3];
            for(int j=0;j<=6;j++){
                for(int k=1;k<=5;k++)
                    pm[10*j+nlayer-i][k]=pm[10*j+nlayerOld-i][k];
            }
        }
    }
    else{
        for(int i=1;i<=nlayer;i++){
            par[50+i][1]=par[50+i-Dnl][1];
            par[50+i][3]=par[50+i-Dnl][3];
            for(int j=0;j<=6;j++){
                for(int k=1;k<=5;k++)
                    pm[10*j+i][k]=pm[10*j+i-Dnl][k];
            }
        }
    }
    ui->pBmDw1->setEnabled(false);
    ui->pBmDw2->setEnabled(false);
    ui->pBmUp2->setEnabled(false);
    ui->pBmDw3->setEnabled(false);
    ui->pBmUp3->setEnabled(false);
    ui->pBmDw4->setEnabled(false);
    ui->pBmUp4->setEnabled(false);
    ui->pBmDw5->setEnabled(false);
    ui->pBmUp5->setEnabled(false);
    ui->pBmDw6->setEnabled(false);
    ui->pBmUp6->setEnabled(false);
    ui->pBmDw7->setEnabled(false);
    ui->pBmUp7->setEnabled(false);
    ui->pBmDw8->setEnabled(false);
    ui->pBmUp8->setEnabled(false);
    ui->pBmDw9->setEnabled(false);
    ui->pBmUp9->setEnabled(false);
    if(nlayer>0)
        ui->pBmDw1->setEnabled(true);
    if(nlayer>1){
        if(nlayer>2) ui->pBmDw2->setEnabled(true);
        ui->pBmUp2->setEnabled(true);
    }
    if(nlayer>2){
        if(nlayer>3) ui->pBmDw3->setEnabled(true);
        ui->pBmUp3->setEnabled(true);
    }
    if(nlayer>3){
        if(nlayer>4) ui->pBmDw4->setEnabled(true);
        ui->pBmUp4->setEnabled(true);
    }
    if(nlayer>4){
        if(nlayer>5) ui->pBmDw5->setEnabled(true);
        ui->pBmUp5->setEnabled(true);
    }
    if(nlayer>5){
        if(nlayer>6) ui->pBmDw6->setEnabled(true);
        ui->pBmUp6->setEnabled(true);
    }
    if(nlayer>6){
        if(nlayer>7) ui->pBmDw7->setEnabled(true);
        ui->pBmUp7->setEnabled(true);
    }
    if(nlayer>7){
        if(nlayer>8) ui->pBmDw8->setEnabled(true);
        ui->pBmUp8->setEnabled(true);
    }
    if(nlayer>8){
        ui->pBmUp9->setEnabled(true);
    }
    for(int i=1;i<=9;i++){
        if(i <= nlayer){
            if(int(par[50+i][3]+0.5)<1 || int(par[50+i][3]+0.5)>5) par[50+i][3]=1.;
            if(int(par[50+i][1]+0.5)<1 || int(par[50+i][1]+0.5)>9) par[50+i][1]=1.;
            idToComboBox["comB_PAR_5"+QString::number(i)+"_1"] -> setEnabled(true);
            idToComboBox["comB_PAR_5"+QString::number(i)+"_3"] -> setEnabled(true);
            idToComboBox["comB_PAR_5"+QString::number(i)+"_1"] -> setCurrentIndex(int(par[50+i][1]+0.5)-1);
            idToComboBox["comB_PAR_5"+QString::number(i)+"_3"] -> setCurrentIndex(int(par[50+i][3]+0.5)-1);
            if(int(par[50+i][3]+0.5)==1){
                idToDoubleSpinBox["dSB_PM_"+QString::number(i)+"_1"]    -> setEnabled(true);
                idToDoubleSpinBox["dSB_PM_"+QString::number(i)+"_1"]    -> setValue(pm[i][1]*1.E-7);
                idToDoubleSpinBox["dSB_PM_"+QString::number(10+i)+"_1"] -> setEnabled(false);
                idToDoubleSpinBox["dSB_PM_"+QString::number(20+i)+"_1"] -> setEnabled(false);
                idToDoubleSpinBox["dSB_PM_"+QString::number(30+i)+"_1"] -> setEnabled(false);
                idToDoubleSpinBox["dSB_PM_"+QString::number(40+i)+"_1"] -> setEnabled(false);
                idToDoubleSpinBox["dSB_PM_"+QString::number(60+i)+"_1"] -> setEnabled(false);
                idToDoubleSpinBox["dSB_PM_"+QString::number(50+i)+"_1"] -> setEnabled(true);
                idToDoubleSpinBox["dSB_PM_"+QString::number(50+i)+"_1"] -> setValue(pm[50+i][1]);
            } else if(int(par[50+i][3]+0.5)==2){
                idToDoubleSpinBox["dSB_PM_"+QString::number(i)+"_1"]    -> setEnabled(true);
                idToDoubleSpinBox["dSB_PM_"+QString::number(i)+"_1"]    -> setValue(pm[i][1]);
                idToDoubleSpinBox["dSB_PM_"+QString::number(10+i)+"_1"] -> setEnabled(false);
                idToDoubleSpinBox["dSB_PM_"+QString::number(20+i)+"_1"] -> setEnabled(false);
                idToDoubleSpinBox["dSB_PM_"+QString::number(30+i)+"_1"] -> setEnabled(false);
                idToDoubleSpinBox["dSB_PM_"+QString::number(40+i)+"_1"] -> setEnabled(false);
                idToDoubleSpinBox["dSB_PM_"+QString::number(60+i)+"_1"] -> setEnabled(false);
                idToDoubleSpinBox["dSB_PM_"+QString::number(50+i)+"_1"] -> setEnabled(true);
                idToDoubleSpinBox["dSB_PM_"+QString::number(50+i)+"_1"] -> setValue(pm[50+i][1]);
            } else if(int(par[50+i][3]+0.5)==3){
                idToDoubleSpinBox["dSB_PM_"+QString::number(i)+"_1"]    -> setEnabled(true);
                idToDoubleSpinBox["dSB_PM_"+QString::number(10+i)+"_1"] -> setEnabled(true);
                idToDoubleSpinBox["dSB_PM_"+QString::number(20+i)+"_1"] -> setEnabled(true);
                idToDoubleSpinBox["dSB_PM_"+QString::number(30+i)+"_1"] -> setEnabled(true);
                idToDoubleSpinBox["dSB_PM_"+QString::number(40+i)+"_1"] -> setEnabled(true);
                idToDoubleSpinBox["dSB_PM_"+QString::number(50+i)+"_1"] -> setEnabled(true);
                idToDoubleSpinBox["dSB_PM_"+QString::number(60+i)+"_1"] -> setEnabled(true);
                idToDoubleSpinBox["dSB_PM_"+QString::number(i)+"_1"]    -> setValue(pm[i][1]);
                idToDoubleSpinBox["dSB_PM_"+QString::number(10+i)+"_1"] -> setValue(pm[10+i][1]);
                idToDoubleSpinBox["dSB_PM_"+QString::number(20+i)+"_1"] -> setValue(pm[20+i][1]);
                idToDoubleSpinBox["dSB_PM_"+QString::number(30+i)+"_1"] -> setValue(pm[30+i][1]);
                idToDoubleSpinBox["dSB_PM_"+QString::number(40+i)+"_1"] -> setValue(pm[40+i][1]);
                idToDoubleSpinBox["dSB_PM_"+QString::number(50+i)+"_1"] -> setValue(pm[50+i][1]);
                idToDoubleSpinBox["dSB_PM_"+QString::number(60+i)+"_1"] -> setValue(pm[60+i][1]);
            }else if(int(par[50+i][3]+0.5)==4 || int(par[50+i][3]+0.5)==5){
                idToDoubleSpinBox["dSB_PM_"+QString::number(i)+"_1"]    -> setEnabled(true);
                idToDoubleSpinBox["dSB_PM_"+QString::number(10+i)+"_1"] -> setEnabled(true);
                idToDoubleSpinBox["dSB_PM_"+QString::number(20+i)+"_1"] -> setEnabled(false);
                idToDoubleSpinBox["dSB_PM_"+QString::number(30+i)+"_1"] -> setEnabled(false);
                idToDoubleSpinBox["dSB_PM_"+QString::number(40+i)+"_1"] -> setEnabled(false);
                idToDoubleSpinBox["dSB_PM_"+QString::number(50+i)+"_1"] -> setEnabled(true);
                idToDoubleSpinBox["dSB_PM_"+QString::number(60+i)+"_1"] -> setEnabled(false);
                idToDoubleSpinBox["dSB_PM_"+QString::number(i)+"_1"]    -> setValue(pm[i][1]);
                idToDoubleSpinBox["dSB_PM_"+QString::number(10+i)+"_1"] -> setValue(pm[10+i][1]);
                idToDoubleSpinBox["dSB_PM_"+QString::number(50+i)+"_1"] -> setValue(pm[50+i][1]);
            }
        } else{
            idToComboBox["comB_PAR_5"+QString::number(i)+"_1"] -> setEnabled(false);
            idToComboBox["comB_PAR_5"+QString::number(i)+"_3"] -> setEnabled(false);
            idToDoubleSpinBox["dSB_PM_"+QString::number(i)+"_1"] -> setEnabled(false);
            idToDoubleSpinBox["dSB_PM_"+QString::number(10+i)+"_1"] -> setEnabled(false);
            idToDoubleSpinBox["dSB_PM_"+QString::number(20+i)+"_1"] -> setEnabled(false);
            idToDoubleSpinBox["dSB_PM_"+QString::number(30+i)+"_1"] -> setEnabled(false);
            idToDoubleSpinBox["dSB_PM_"+QString::number(40+i)+"_1"] -> setEnabled(false);
            idToDoubleSpinBox["dSB_PM_"+QString::number(50+i)+"_1"] -> setEnabled(false);
            idToDoubleSpinBox["dSB_PM_"+QString::number(60+i)+"_1"] -> setEnabled(false);
        }
    }
    occupyPF=0;
}


void ksemaw::RefreshModel(){
 if(occupyPF!=1){
  printf("-> RefreshModel\n");
  for(int i=1;i<=9;i++){
   if(i <= nlayer){    
    par[50+i][3]=idToComboBox["comB_PAR_5"+QString::number(i)+"_3"] -> currentIndex();
    par[50+i][3]++;
    par[50+i][1]=idToComboBox["comB_PAR_5"+QString::number(i)+"_1"] -> currentIndex();
    par[50+i][1]++;
    if(int(par[50+i][3]+0.5)==1){
     idToDoubleSpinBox["dSB_PM_"+QString::number(i)+"_1"] -> setEnabled(true);
//     idToDoubleSpinBox["dSB_PM_"+QString::number(i)+"_1"] -> setValue(pm[i][1]*1.E-7);
     idToDoubleSpinBox["dSB_PM_"+QString::number(10+i)+"_1"] -> setEnabled(false);
     idToDoubleSpinBox["dSB_PM_"+QString::number(20+i)+"_1"] -> setEnabled(false);
     idToDoubleSpinBox["dSB_PM_"+QString::number(30+i)+"_1"] -> setEnabled(false);
     idToDoubleSpinBox["dSB_PM_"+QString::number(40+i)+"_1"] -> setEnabled(false);
     idToDoubleSpinBox["dSB_PM_"+QString::number(60+i)+"_1"] -> setEnabled(false);
     idToDoubleSpinBox["dSB_PM_"+QString::number(50+i)+"_1"] -> setEnabled(true);
//     idToDoubleSpinBox["dSB_PM_"+QString::number(i)+"_1"] -> setValue(pm[50+i][1]);
    } else if(int(par[50+i][3]+0.5)==2){
     idToDoubleSpinBox["dSB_PM_"+QString::number(i)+"_1"] -> setEnabled(true);
//     idToDoubleSpinBox["dSB_PM_"+QString::number(i)+"_1"] -> setValue(pm[i][1]);
     idToDoubleSpinBox["dSB_PM_"+QString::number(10+i)+"_1"] -> setEnabled(false);
     idToDoubleSpinBox["dSB_PM_"+QString::number(20+i)+"_1"] -> setEnabled(false);
     idToDoubleSpinBox["dSB_PM_"+QString::number(30+i)+"_1"] -> setEnabled(false);
     idToDoubleSpinBox["dSB_PM_"+QString::number(40+i)+"_1"] -> setEnabled(false);
     idToDoubleSpinBox["dSB_PM_"+QString::number(60+i)+"_1"] -> setEnabled(false);
     idToDoubleSpinBox["dSB_PM_"+QString::number(50+i)+"_1"] -> setEnabled(true);
     idToDoubleSpinBox["dSB_PM_"+QString::number(50+i)+"_1"] -> setValue(pm[50+i][1]);
    } else if(int(par[50+i][3]+0.5)==3){
     idToDoubleSpinBox["dSB_PM_"+QString::number(i)+"_1"] -> setEnabled(true);
     idToDoubleSpinBox["dSB_PM_"+QString::number(10+i)+"_1"] -> setEnabled(true);
     idToDoubleSpinBox["dSB_PM_"+QString::number(20+i)+"_1"] -> setEnabled(true);
     idToDoubleSpinBox["dSB_PM_"+QString::number(30+i)+"_1"] -> setEnabled(true);
     idToDoubleSpinBox["dSB_PM_"+QString::number(40+i)+"_1"] -> setEnabled(true);
     idToDoubleSpinBox["dSB_PM_"+QString::number(50+i)+"_1"] -> setEnabled(true);
     idToDoubleSpinBox["dSB_PM_"+QString::number(60+i)+"_1"] -> setEnabled(true);
//     idToDoubleSpinBox["dSB_PM_"+QString::number(i)+"_1"] -> setValue(pm[i][1]);
     idToDoubleSpinBox["dSB_PM_"+QString::number(10+i)+"_1"] -> setValue(pm[10+i][1]);
     idToDoubleSpinBox["dSB_PM_"+QString::number(20+i)+"_1"] -> setValue(pm[20+i][1]);
     idToDoubleSpinBox["dSB_PM_"+QString::number(30+i)+"_1"] -> setValue(pm[30+i][1]);
     idToDoubleSpinBox["dSB_PM_"+QString::number(40+i)+"_1"] -> setValue(pm[40+i][1]);
     idToDoubleSpinBox["dSB_PM_"+QString::number(50+i)+"_1"] -> setValue(pm[50+i][1]);
     idToDoubleSpinBox["dSB_PM_"+QString::number(60+i)+"_1"] -> setValue(pm[60+i][1]);
    }else if(int(par[50+i][3]+0.5)==4 || int(par[50+i][3]+0.5)==5){
        idToDoubleSpinBox["dSB_PM_"+QString::number(i)+"_1"] -> setEnabled(true);
        idToDoubleSpinBox["dSB_PM_"+QString::number(10+i)+"_1"] -> setEnabled(true);
        idToDoubleSpinBox["dSB_PM_"+QString::number(20+i)+"_1"] -> setEnabled(false);
        idToDoubleSpinBox["dSB_PM_"+QString::number(30+i)+"_1"] -> setEnabled(false);
        idToDoubleSpinBox["dSB_PM_"+QString::number(40+i)+"_1"] -> setEnabled(false);
        idToDoubleSpinBox["dSB_PM_"+QString::number(50+i)+"_1"] -> setEnabled(true);
        idToDoubleSpinBox["dSB_PM_"+QString::number(60+i)+"_1"] -> setEnabled(false);
   //     idToDoubleSpinBox["dSB_PM_"+QString::number(i)+"_1"] -> setValue(pm[i][1]);
        idToDoubleSpinBox["dSB_PM_"+QString::number(10+i)+"_1"] -> setValue(pm[10+i][1]);
        idToDoubleSpinBox["dSB_PM_"+QString::number(50+i)+"_1"] -> setValue(pm[50+i][1]);
       }
   } else{
    idToComboBox["comB_PAR_5"+QString::number(i)+"_1"] -> setEnabled(false);
    idToComboBox["comB_PAR_5"+QString::number(i)+"_3"] -> setEnabled(false);
    idToDoubleSpinBox["dSB_PM_"+QString::number(i)+"_1"] -> setEnabled(false);
    idToDoubleSpinBox["dSB_PM_"+QString::number(10+i)+"_1"] -> setEnabled(false);
    idToDoubleSpinBox["dSB_PM_"+QString::number(20+i)+"_1"] -> setEnabled(false);
    idToDoubleSpinBox["dSB_PM_"+QString::number(30+i)+"_1"] -> setEnabled(false);
    idToDoubleSpinBox["dSB_PM_"+QString::number(40+i)+"_1"] -> setEnabled(false);
    idToDoubleSpinBox["dSB_PM_"+QString::number(50+i)+"_1"] -> setEnabled(false);
    idToDoubleSpinBox["dSB_PM_"+QString::number(60+i)+"_1"] -> setEnabled(false);
   }   
  }
 }
}


void ksemaw::setMat1(){
    setMat(1);
    listOsc();
}
void ksemaw::setMat2(){
    setMat(2);
}
void ksemaw::setMat3(){
    setMat(3);
}
void ksemaw::setMat4(){
    setMat(4);
}
void ksemaw::setMat5(){
    setMat(5);
}
void ksemaw::setMat6(){
    setMat(6);
}
void ksemaw::setMat7(){
    setMat(7);
}
void ksemaw::setMat8(){
    setMat(8);
}
void ksemaw::setMat9(){
    setMat(9);
}
void ksemaw::setMat10(){
    setMat(10);
}
void ksemaw::setMat11(){
    setMat(11);
}
void ksemaw::setMat12(){
    setMat(12);
}
void ksemaw::setMat13(){
    setMat(13);
}
void ksemaw::setMat14(){
    setMat(14);
}
void ksemaw::setMat15(){
    setMat(15);
}
void ksemaw::setMat(int nM){
    printf("->setMat called by cnk %d\n",nM);
    int curI=idToComboBox["cB_cnk"+QString::number(nM)+"a"] -> currentIndex();
    idToLineEdit["LEcnk"+QString::number(nM)+"_2"] -> setEnabled(curI==0);
    idToLineEdit["LEcnk"+QString::number(nM)+"_3"] -> setEnabled(curI==0);
}

void ksemaw::setEMA1(){
    setEMA(1);
}
void ksemaw::setEMA2(){
    setEMA(2);
}
void ksemaw::setEMA3(){
    setEMA(3);
}
void ksemaw::setEMA4(){
    setEMA(4);
}
void ksemaw::setEMA5(){
    setEMA(5);
}
void ksemaw::setEMA6(){
    setEMA(6);
}
void ksemaw::setEMA7(){
    setEMA(7);
}
void ksemaw::setEMA8(){
    setEMA(8);
}
void ksemaw::setEMA9(){
    setEMA(9);
}
void ksemaw::setEMA10(){
    setEMA(10);
}
void ksemaw::setEMA11(){
    setEMA(11);
}
void ksemaw::setEMA12(){
    setEMA(12);
}
void ksemaw::setEMA13(){
    setEMA(13);
}
void ksemaw::setEMA14(){
    setEMA(14);
}
void ksemaw::setEMA15(){
    setEMA(15);
}
void ksemaw::setEMA(int nM){
    Qt::CheckState state;
    state=idToCheckBox["cB_EMA_"+QString::number(nM)] -> checkState();
    printf("->setEMA called by checkBox %d\n",nM);
    idToComboBox["cB_cnk"+QString::number(nM)+"b"] -> setEnabled(state==Qt::Checked);
    idToDoubleSpinBox["dSB_cnk"+QString::number(nM)] -> setEnabled(state==Qt::Checked);
}

void ksemaw::setOsc1(){
    setOscN(1);
}

void ksemaw::setOsc2(){
    setOscN(2);
}

void ksemaw::setOsc3(){
    setOscN(3);
}

void ksemaw::setOsc4(){
    setOscN(4);
}

void ksemaw::setOsc5(){
    setOscN(5);
}

void ksemaw::setOsc6(){
    setOscN(6);
}

void ksemaw::setOsc7(){
    setOscN(7);
}

void ksemaw::setOsc8(){
    setOscN(8);
}

void ksemaw::setOsc9(){
    setOscN(9);
}

void ksemaw::setOsc10(){
    setOscN(10);
}

void ksemaw::setOsc11(){
    setOscN(11);
}

void ksemaw::setOsc12(){
    setOscN(12);
}

void ksemaw::setOsc13(){
    setOscN(13);
}

void ksemaw::setOsc14(){
    setOscN(14);
}

void ksemaw::setOsc15(){
    setOscN(15);
}

void ksemaw::setOsc16(){
    setOscN(16);
}

void ksemaw::setOsc17(){
    setOscN(17);
}

void ksemaw::setOsc18(){
    setOscN(18);
}

void ksemaw::setOsc19(){
    setOscN(19);
}

void ksemaw::setOsc20(){
    setOscN(20);
}

void ksemaw::setOscN(int k){
    //int itab=ui->tabWidget -> currentIndex();
    if(occupyPF!=0) return;
    int ioptFit=ui->cB_cnk1a -> currentIndex();//Fit#
    int Nosc=pf[ioptFit][1];//N. Osc
    Qt::CheckState state;
    state=idToCheckBox["cBosc_"+QString::number(k)]->checkState();
    printf("->setOscN ioptFit=%d N.Osc=%d call-by-oscN=%d ",ioptFit,Nosc,k);
    if(state==Qt::Checked){
        printf("to be add\n");
        pf[ioptFit][Nosc+1+1]=k;
        pf[ioptFit][1]=Nosc+1;
    }
    else{
        printf("to be deleted\n");
        int jj=1;
        while(pf[ioptFit][jj+1]!=k && jj<Nosc){
            jj++;
        }
        for (int i=jj;i<=Nosc;i++) {
            pf[ioptFit][i+1]=pf[ioptFit][i+2];
        }
        pf[ioptFit][Nosc+1]=0;
        pf[ioptFit][1]=Nosc-1;
    }
    listOsc();
}

void ksemaw::listOsc(){
 //int itab=ui->tabWidget -> currentIndex();
 if(occupyPF!=0) return;
 occupyPF=1;
 int iok;
 int ioptFit=ui->cB_cnk1a -> currentIndex();
 printf("->ListOsc: ioptFit= %d\n",ioptFit);
 if(ioptFit==0 || ioptFit>7){
     for(int k=1;k<=20;k++){
         idToComboBox["cBpm_"+QString::number(100+1+(k-1)*5)+"_1"] -> setEnabled(false);
         idToLineEdit["LEpm_"+QString::number(100+2+(k-1)*5)+"_1"] -> setEnabled(false);
         idToLineEdit["LEpm_"+QString::number(100+3+(k-1)*5)+"_1"] -> setEnabled(false);
         idToLineEdit["LEpm_"+QString::number(100+4+(k-1)*5)+"_1"] -> setEnabled(false);
         idToLineEdit["LEpm_"+QString::number(100+5+(k-1)*5)+"_1"] -> setEnabled(false);
     }
 }
 if(ioptFit>0 && ioptFit<8){
  for(int k=1;k<=20;k++){
   iok=0;
   for(int i=2;i<=21;i++) {
    if(int(pf[ioptFit][i]+0.5)==k)iok=1;
   }
   if(iok==1){
    idToCheckBox["cBosc_"+QString::number(k)]-> setCheckState ( Qt::Checked );
    idToComboBox["cBpm_"+QString::number(100+1+(k-1)*5)+"_1"] -> setEnabled(true);
    idToLineEdit["LEpm_"+QString::number(100+2+(k-1)*5)+"_1"] -> setEnabled(true);
    idToLineEdit["LEpm_"+QString::number(100+3+(k-1)*5)+"_1"] -> setEnabled(true);
    idToLineEdit["LEpm_"+QString::number(100+4+(k-1)*5)+"_1"] -> setEnabled(true);
    idToLineEdit["LEpm_"+QString::number(100+5+(k-1)*5)+"_1"] -> setEnabled(true);
    idToComboBox["cBpm_"+QString::number(100+1+(k-1)*5)+"_1"] -> setCurrentIndex(int(pm[100+1+(k-1)*5][1]+0.5)-1);
    idToLineEdit["LEpm_"+QString::number(100+2+(k-1)*5)+"_1"] -> setText(QString::number(abs(pm[100+2+(k-1)*5][1])));
    idToLineEdit["LEpm_"+QString::number(100+3+(k-1)*5)+"_1"] -> setText(QString::number(pm[100+3+(k-1)*5][1]));
    idToLineEdit["LEpm_"+QString::number(100+4+(k-1)*5)+"_1"] -> setText(QString::number(pm[100+4+(k-1)*5][1]));
    idToLineEdit["LEpm_"+QString::number(100+5+(k-1)*5)+"_1"] -> setText(QString::number(pm[100+5+(k-1)*5][1]));
   }
   else{
    idToCheckBox["cBosc_"+QString::number(k)]-> setCheckState ( Qt::Unchecked );
    idToComboBox["cBpm_"+QString::number(100+1+(k-1)*5)+"_1"] -> setEnabled(false);
    idToLineEdit["LEpm_"+QString::number(100+2+(k-1)*5)+"_1"] -> setEnabled(false);
    idToLineEdit["LEpm_"+QString::number(100+3+(k-1)*5)+"_1"] -> setEnabled(false);
    idToLineEdit["LEpm_"+QString::number(100+4+(k-1)*5)+"_1"] -> setEnabled(false);
    idToLineEdit["LEpm_"+QString::number(100+5+(k-1)*5)+"_1"] -> setEnabled(false);
   }
  }
 }
 else{
  for(int k=1;k<=20;k++)
    idToCheckBox["cBosc_"+QString::number(k)]-> setCheckState ( Qt::Unchecked );
 }
 lastIndex=ioptFit;
 occupyPF=0;
}


void ksemaw::rangeWL(){
    double Lmin=ui->dSB_PAR_4_1 -> value();
    double Lmax=ui->dSB_PAR_4_2 -> value();
    Qt::CheckState state;
    state = ui->checkBox_msgr -> checkState();
    if(state == Qt::Unchecked){
        rxy[20][1]=Lmin;
        rxy[20][2]=Lmax;
        ui->DP_RXY_20_1 ->setText(QString::number(Lmin));
        ui->DP_RXY_20_2 ->setText(QString::number(Lmax));
    }
}


void ksemaw::saveSim(){
    int ierr,iok=1;
    QString command, fn2s;
    QString fSim=pathroot+"expo/MisSim.dat";
    fn2s = QFileDialog::getSaveFileName(
                this,
                "Filename to save",
                pathroot+"expo/",
                "dat file (*.dat)");
    QFile file(fn2s);
    if( file.exists() ){
        QMessageBox msgBox;
        msgBox.setText("The file already exist!");
        msgBox.setInformativeText("Do you want to save anyway?");
        msgBox.setStandardButtons(QMessageBox::Save | QMessageBox::Cancel);
        msgBox.setDefaultButton(QMessageBox::Save);
        int ret = msgBox.exec();
        switch (ret) {
        case QMessageBox::Save:
            // Save was clicked
            iok=1;
            break;
        case QMessageBox::Cancel:
            // Cancel was clicked
            iok=0;
            printf("as you requested file was not saved!\n");
            break;
        }
    }
    if(iok==1){
        cout<< "fsim= "<<fSim.toStdString()<<endl;
        cout<< "fn2s= "<<fn2s.toStdString()<<endl;
        command="cp "+fSim+" "+fn2s;
        ierr=system((command.toStdString()).c_str());
        if(ierr != 0)
            printf("Error copying MisSim.dat!!!\n");
        else
            printf("MisSim saved as %s\n",(fn2s.toStdString()).c_str());
    }
}

void ksemaw::tabChanged(){
    QMessageBox msgBox;
    int itab=ui->tabWidget -> currentIndex();
    printf("-> tabChanged:  itab=%d lastTab=%d\n",itab,lastTab);
    SaveSetting(lastTab);
    QString filepro=ui->lineEdit_P -> text();
    if(filepro.contains("mate/aa999")){
        ui->tabWidget -> setCurrentIndex(0);
        if(iwspj==0){
            msgBox.setText("First of all you MUST name the Project!");
            msgBox.setStandardButtons(QMessageBox::Ok);
            msgBox.exec();
            iwspj=1;
        }
        return;
    }
    if(lastTab==0 && itab!=0){
        occupyPF=1;
        QFile file(fctrl);
        if (!file.open(QIODevice::WriteOnly | QIODevice::Text)){
            printf("! in tabChanged failure opening fctrl\n");
            occupyPF=0;
            return;
        }
        QTextStream out(&file);
        out << "valin" << "\n";
        file.close();
        occupyPF=0;
    }
    if(itab==4){
        int ioptFit=ui->cB_cnk1a -> currentIndex();
        if(ioptFit==0 || ioptFit>7){
            ui->tabWidget -> setCurrentIndex(2);
            msgBox.setText("Please select a valid FIT# and SAVE!");
            msgBox.setStandardButtons(QMessageBox::Ok);
            msgBox.exec();
            return;
        }
        refreshFitPar();
    }
    if(itab==5)
        lastTabB5=lastTab;
    if(lastTab==5){
        if(lastTabB5==2)
            PlotME();
        else if(lastTabB5==3)
            RefTrackG();
        else if(lastTabB5==4)
            RefIbridG();
    }
    lastTab=itab;
}


void ksemaw::refreshFitPar(){
    int icoherent,klim,ivp,n_fresh,iDecine,ilayer,iosc,Jcombo,ip;
    occupyPF=1;
    printf("-> refreshFitPar\n");

    npp=ui->sB_PAR_34_5 -> value();
    n_fresh=0;
    Qt::CheckState state;
    state = Qt::Checked;
    QString Lj, Lparametro,Valore;
    for(int j=1;j<=17;j++){
        idToComboBox["cBParFit_"+QString::number(j)] -> clear();
        idToComboBox["cBParFit_"+QString::number(j)] -> addItem("none");
    }
    for(int i=1;i<=nlayer;i++){
        icoherent=idToComboBox["comB_PAR_5"+QString::number(i)+"_3"] -> currentIndex();
        if(icoherent>0){
            Lj="L"+QString::number(i);
            klim=1;
            if(icoherent==2)
                klim=7;
            else if(icoherent == 3 || icoherent == 4)
                klim=2;
            for(int j=1;j<=npp;j++){
                for(int k=1;k<=klim;k++){
                    idToComboBox["cBParFit_"+QString::number(j)]-> addItem(Lj+ParFitLab[k]);
                }
            }
        }
    }
    for(int i=1;i<=20;i++){
        Lj="O"+QString::number(i);
        state=idToCheckBox["cBosc_"+QString::number(i)]-> checkState();
        if( state == Qt::Checked ) {
            for(int j=1;j<=npp;j++){
                for(int k=8;k<=11;k++){
                    idToComboBox["cBParFit_"+QString::number(j)]-> addItem(Lj+ParFitLab[k]);
                }
            }
        }
    }
    for(int j=1;j<=17;j++){
        if(j<=npp){
            //if(lastTab != itab){
                ivp=int(pm[ppm[j]][2]+0.5);
                ip=ppm[j];
                if(ip<=59){
                    iDecine=int(double(ip)/10.);
                    ilayer=ip-iDecine*10;
                    Lparametro="L"+QString::number(ilayer)+ParFitLab[iDecine+1];
                    //      printf("j=%d ip=%d iDecine=%d ilayer=%d Lparametro=%s\n",j,ip,iDecine,ilayer,(Lparametro.toStdString()).c_str());
                }
                else{
                    iosc=int(double(ip-101)/5.)+1;
                    Lparametro="O"+QString::number(iosc)+ParFitLab[7+ip-101-5*(iosc-1)];
                    //      printf("j=%d ip=%d iosc=%d Lparametro=%s\n",j,ip,iosc,(Lparametro.toStdString()).c_str());
                }
                Jcombo=idToComboBox["cBParFit_"+QString::number(j)] -> findText(Lparametro);
                if(Jcombo>=0) idToComboBox["cBParFit_"+QString::number(j)] -> setCurrentIndex(Jcombo);
                idToLineEdit["DPparFitV_"+QString::number(j)] -> setText(QString::number(pm[ip][1]));
                if(ivp==0)
                    idToCheckBox["chBeParFit_"+QString::number(j)] -> setCheckState ( Qt::Unchecked );
                else{
                    idToCheckBox["chBeParFit_"+QString::number(j)] -> setCheckState ( Qt::Checked );
                    n_fresh++;
                    idToLineEdit["DPparFitErr_"+QString::number(j)] -> setText(QString::number(pm[ip][4]));
                    idToLineEdit["DPparFitGC_"+QString::number(j)] -> setText(QString::number(pm[ip][5]));
                }
            //}
            //else{
            //    state=idToCheckBox["chBeParFit_"+QString::number(j)]-> checkState();
            //    if( state == Qt::Checked ) n_fresh++;
            //}
        }
        else
            idToCheckBox["chBeParFit_"+QString::number(j)] -> setCheckState ( Qt::Unchecked );
    }
    n=n_fresh;
    ui->sB_PAR_35_5 -> setValue(n);
    ui->sB_free_deg -> setValue(npp-n);
    occupyPF=0;
    SaveSetting(-1);
}


//void ksemaw::savePanFit(){
//  printf("-> savePanFit\n");
//  int n,n1,n2,ip,jpf;
//  QString Lj,Valore;
//  Qt::CheckState state;
//  n=17;
//  npp=17;
//  jpf=0;
//  for(int j=1;j<=17;j++){
//    Lj=idToComboBox["cBParFit_"+QString::number(j)] -> currentText();
//    state=idToCheckBox["chBeParFit_"+QString::number(j)]-> checkState();
//    if(Lj=="none"){
//      n--;
//      npp--;
//    }
//    else{
//      n1=(Lj.at(1)).digitValue();
//      if(Lj.at(2).isNumber()) {
//        n1=n1*10;
//        n1=n1+((Lj.at(2)).digitValue());
//      }
//      n2=0;
//      ip=0;
//      do{
//        n2++;
//      } while(!Lj.contains(ParFitLab[n2],Qt::CaseSensitive));
//      if(1<=n2 && n2<=7)
//        ip=10*(n2-1)+n1;
//      else if(8<=n2 && n2<=11)
//        ip=100+(n1-1)*5+(n2-7)+1;
//         printf("SavePanFit:  n1=%d n2=%d ip=%d\n",n1,n2,ip);
//      ppm[j]=ip;
//      Valore=idToLineEdit["DPparFitV_"+QString::number(j)] -> text();
//      pm[ip][1]=Valore.toDouble();
//      if( state == Qt::Unchecked ){
//        pm[ip][2]=0;
//        n--;
//      }
//      else{
//        jpf++;
//        pm[ip][2]=jpf;
//        pm[jpf][3]=ip;
//      }
//    }
//  }
//  par[34][5]=npp;//N par fit
//  par[35][5]=n;//N par enabled
//  sB_PAR_35_5 -> setValue(n);
//  sB_free_deg -> setValue(npp-n);
//  sB_PAR_34_5 -> setValue(npp);
//}

void ksemaw::PanFitEnable(){
 printf("-> PanFitEnable\n");
 int n_fresh,n1,n2,ip;
 int itab=ui->tabWidget -> currentIndex();
 
 if(occupyPF==0){
  occupyPF=1;
  npp=ui->sB_PAR_34_5 -> value();
  n_fresh=0;
  Qt::CheckState state;
  state = Qt::Checked;
  QString Lj,Valore;
  if(itab==4 && itab==lastTab){
//   printf("Start PanFitEnable\n");
   for(int j=1;j<=17;j++){
    if(j<=npp){
     state=idToCheckBox["chBeParFit_"+QString::number(j)]-> checkState();
     if( state == Qt::Checked ) n_fresh++;
     Lj=idToComboBox["cBParFit_"+QString::number(j)] -> currentText();
     if(Lj.at(0)=='L'){
      n1=(Lj.at(1)).digitValue();
      n2=0;
      ip=0;
      do{
       n2++;
      } while(!Lj.contains(ParFitLab[n2],Qt::CaseSensitive) && n2<7);
      ip=10*(n2-1)+n1;
     }
     else{
      n1=(Lj.at(1)).digitValue();
      if(Lj.at(2).isNumber()) {
       n1=n1*10;
       n1=n1+((Lj.at(2)).digitValue());
      }
      n2=6;
      ip=0;
      do{
       n2++;
      } while(!Lj.contains(ParFitLab[n2],Qt::CaseSensitive) && n2<11);
      ip=100+(n1-1)*5+(n2-7)+1;
     }
     ppm[j]=ip;
     Valore=idToLineEdit["DPparFitV_"+QString::number(j)] -> text();
     pm[ip][1]=Valore.toDouble();
     if( state == Qt::Unchecked )
      pm[ip][2]=0;
     else{
      pm[ip][2]=j;
      pm[n_fresh][3]=ip;
     }
//     printf("j=%d ip=%d pm[%d][1]=%f pm[%d][2]=%f pm[%d][3]=%f\n",j,ip,ip,pm[ip][1],ip,pm[ip][2],n_fresh,pm[n_fresh][3]);
    } 
    else
     idToCheckBox["chBeParFit_"+QString::number(j)] -> setCheckState ( Qt::Unchecked );
   }
   n=n_fresh;
   ui->sB_PAR_35_5 -> setValue(n);
   ui->sB_free_deg -> setValue(npp-n);
//   printf("END PanFitEnable\n");
  }
  occupyPF=0;
 }
}


void ksemaw::PanFitPar(){
    if(occupyPF!=0) return;
    printf("-> PanFitPar\n");
    int icoherent,klim,nppNew;
    int itab=ui->tabWidget -> currentIndex();
    if(occupyPF==0){
        occupyPF=1;
        nppNew=ui->sB_PAR_34_5 -> value();
        Qt::CheckState state;
        state = Qt::Checked;
        QString Lj, Lparametro;
        if(itab==4 && itab==lastTab){
           printf("\tStart PanFitPar: nppNew=%d npp=%d\n",nppNew,npp);
            if(nppNew > npp){
                for(int j=npp+1;j<=nppNew;j++){
                    idToComboBox["cBParFit_"+QString::number(j)] -> clear();
                    idToComboBox["cBParFit_"+QString::number(j)] -> addItem("none");
                    idToLineEdit["DPparFitV_"+QString::number(j)]-> clear();
                    idToLineEdit["DPparFitErr_"+QString::number(j)]-> clear();
                    idToLineEdit["DPparFitGC_"+QString::number(j)]-> clear();
                }
                for(int i=1;i<=nlayer;i++){
                    icoherent=idToComboBox["comB_PAR_5"+QString::number(i)+"_3"] -> currentIndex();
                    if(icoherent>0){
                        Lj="L"+QString::number(i);
                        klim=1;
                        if(icoherent==2)
                            klim=7;
                        else if(icoherent==3 || icoherent==4)
                            klim=2;
                        for(int j=npp+1;j<=nppNew;j++){
                            for(int k=1;k<=klim;k++){
                                idToComboBox["cBParFit_"+QString::number(j)]-> addItem(Lj+ParFitLab[k]);
                            }
                        }
                    }
                }
                for(int i=1;i<=20;i++){
                    Lj="O"+QString::number(i);
                    state=idToCheckBox["cBosc_"+QString::number(i)]-> checkState();
                    if( state == Qt::Checked ) {
                        for(int j=npp+1;j<=nppNew;j++){
                            for(int k=8;k<=11;k++){
                                idToComboBox["cBParFit_"+QString::number(j)]-> addItem(Lj+ParFitLab[k]);
                            }
                        }
                    }
                }
            }
            else{
                for(int j=nppNew+1;j<=npp;j++){
                    idToComboBox["cBParFit_"+QString::number(j)] -> clear();
                    idToComboBox["cBParFit_"+QString::number(j)] -> addItem("none");
                    idToCheckBox["chBeParFit_"+QString::number(j)] -> setCheckState ( Qt::Unchecked );
                    idToLineEdit["DPparFitV_"+QString::number(j)] -> setText("");
                    idToLineEdit["DPparFitErr_"+QString::number(j)] -> clear();
                    idToLineEdit["DPparFitGC_"+QString::number(j)] -> clear();
                }
            }
            //   printf("End PanFitPar\n");
            npp=nppNew;
            PanFitEnable();
        }
        occupyPF=0;
    }
}


void ksemaw::PanFitChoice(){
    int n1,n2,ip;
    int itab=ui->tabWidget -> currentIndex();
    QString Lj;
    Qt::CheckState state;
    if(occupyPF==0){
        printf("-> PanFitChoice\n");
        occupyPF=1;
        if(itab==4 && itab==lastTab){
            npp=ui->sB_PAR_34_5 -> value();
            for(int j=1;j<=npp;j++){
                Lj=idToComboBox["cBParFit_"+QString::number(j)] -> currentText();
                if(Lj=="none"){
                    idToLineEdit["DPparFitV_"+QString::number(j)]   -> clear();
                    idToLineEdit["DPparFitErr_"+QString::number(j)] -> clear();
                    idToLineEdit["DPparFitGC_"+QString::number(j)]  -> clear();
                    idToCheckBox["chBeParFit_"+QString::number(j)] -> setCheckState ( Qt::Unchecked );
                }
                else{
                    n1=(Lj.at(1)).digitValue();
                    if(Lj.at(2).isNumber()) {
                        n1=n1*10;
                        n1=n1+((Lj.at(2)).digitValue());
                    }
                    n2=0;
                    ip=0;
                    do{
                        n2++;
                    } while(!Lj.contains(ParFitLab[n2],Qt::CaseSensitive) && n2<=11);
                    if(1<=n2 && n2<=7)
                        ip=10*(n2-1)+n1;
                    else if(8<=n2 && n2<=11)
                        ip=100+(n1-1)*5+(n2-7)+1;
                    printf("PanFitChoice j=%d Lj=%s n1=%d n2=%d ip=%d\n",j,(Lj.toStdString()).c_str(),n1,n2,ip);
                    if(1<=ip && ip<=200){
                        idToLineEdit["DPparFitV_"+QString::number(j)] -> setText(QString::number(pm[ip][1]));
                        state=idToCheckBox["chBeParFit_"+QString::number(j)]-> checkState();
                        if( state == Qt::Checked ){
                            idToLineEdit["DPparFitErr_"+QString::number(j)] -> setText(QString::number(pm[ip][4]));
                            idToLineEdit["DPparFitGC_"+QString::number(j)] -> setText(QString::number(pm[ip][5]));
                        }
                        else{
                            idToLineEdit["DPparFitErr_"+QString::number(j)] -> clear();
                            idToLineEdit["DPparFitGC_"+QString::number(j)]  -> clear();
                        }
                    }
                }
            }
            //   printf("End PanFitChoice\n");
        }
        occupyPF=0;
    }
}



void ksemaw::PlotME(){
 SaveSetting(-1);
 QFile file(fctrl);
  if (!file.open(QIODevice::WriteOnly | QIODevice::Text))
  return;
 QTextStream out(&file);
 printf("simul & G\n");
 out << "simul" << "\n";
 out << "G" << "\n";
 file.close();
}

void ksemaw::Simula(){
 //SaveSetting(-1);
 Qt::CheckState state;
 //int ioptFit=cB_cnk1a -> currentIndex();
   //saveOsc(ioptFit);
 state=ui->cBox_PAR_19_1 -> checkState();
 int iok=0;
 if(state == Qt::Checked){
  QMessageBox msgBox;
  msgBox.setText("You are imposing Simulation -> Experimental!");
  msgBox.setInformativeText("Are you sure?");
  msgBox.setStandardButtons(QMessageBox::Yes | QMessageBox::Cancel);
  msgBox.setDefaultButton(QMessageBox::Cancel);
  int ret = msgBox.exec();
  switch (ret) {
    case QMessageBox::Yes:
       // Save was clicked
       iok=1;
       break;
    case QMessageBox::Cancel:
       // Cancel was clicked
       iok=0;
       break;
  }
 }
 if(iok==1)
  par[19][1]=1;
 else
  par[19][1]=0;
 state=ui->cBox_PAR_9_2 -> checkState();
 par[9][2]=0;
 if(state == Qt::Checked)
     par[9][2]=1;
 SaveSetting(-1);
 QFile file(fctrl);
  if (!file.open(QIODevice::WriteOnly | QIODevice::Text))
  return;
 QTextStream out(&file);
 out << "simul" << "\n";
 out << "g" << "\n";
 file.close();
 ui->cBox_PAR_19_1 -> setCheckState(Qt::Unchecked);
 ui->cBox_PAR_9_2  -> setCheckState(Qt::Unchecked);
// printf("checkpoint: occupyPF = %d\n",occupyPF);
}

void ksemaw::PlotAve(){
 SaveSetting(-1);
 QFile file(fctrl);
  if (!file.open(QIODevice::WriteOnly | QIODevice::Text))
  return;
 QTextStream out(&file);
 out << "simul" << "\n";
 out << "gt" << "\n";
 file.close();
}

void ksemaw::PlotAbsEL(){
 SaveSetting(-1);
 QFile file(fctrl);
  if (!file.open(QIODevice::WriteOnly | QIODevice::Text))
  return;
 QTextStream out(&file);
 out << "simul" << "\n";
 out << "pA" << "\n";
 file.close();
}

void ksemaw::readCursor(){
 par[12][1]=ui->sB_par_12_1 -> value();
 SaveSetting(-1);
 QFile file(fctrl);
  if (!file.open(QIODevice::WriteOnly | QIODevice::Text))
  return;
 QTextStream out(&file);
 out << "curso" << "\n";
 file.close();
}

void ksemaw::manageLEwl(){
   int iChoiceWL=ui->comboBox_searchNK -> currentIndex();
   ui->dSB_WLsearchNK->setEnabled(iChoiceWL==0);
}

void ksemaw::searchNK(){
 int iChoiceWL=ui->comboBox_searchNK -> currentIndex();
 if(iChoiceWL==0)
  par[7][1]=ui->dSB_WLsearchNK -> value();
 else if(iChoiceWL==1)
  par[7][1]=ui->dSB_PAR_4_1 -> value();
 else if(iChoiceWL==2)
  par[7][1]=ui->dSB_PAR_4_2 -> value();
 SaveSetting(-1);
 QFile file(fctrl);
  if (!file.open(QIODevice::WriteOnly | QIODevice::Text))
  return;
 QTextStream out(&file);
 out << "autos" << "\n";
 out << "x" << "\n";
 file.close();
}

void ksemaw::save1NK(){
 QFile file(fctrl);
  if (!file.open(QIODevice::WriteOnly | QIODevice::Text))
  return;
 QTextStream out(&file);
 out << "autos" << "\n";
 out << "i" << "\n";
 file.close();
}

void ksemaw::RefTrackG(){
 SaveSetting(-1);
 QFile file(fctrl);
  if (!file.open(QIODevice::WriteOnly | QIODevice::Text))
  return;
 QTextStream out(&file);
 out << "autos" << "\n";
 out << "G" << "\n";
 file.close();
}

void ksemaw::TrackFromWLmin(){
 SaveSetting(-1);
 QFile file(fctrl);
  if (!file.open(QIODevice::WriteOnly | QIODevice::Text))
  return;
 QTextStream out(&file);
 out << "autos" << "\n";
 out << "a" << "\n";
 file.close();
 Ready(false);
 QString scmd;
 do{
     usleep(100);
     QFile filectrl(fctrl);
     if (!filectrl.open(QIODevice::ReadOnly | QIODevice::Text))
         return;
     QTextStream streamctrl ( &filectrl );
     streamctrl >> scmd;
     filectrl.close();
 } while(!scmd.contains("done!"));
 Ready(true);
}

void ksemaw::TrackFromWLmax(){
 SaveSetting(-1);
 QFile file(fctrl);
  if (!file.open(QIODevice::WriteOnly | QIODevice::Text))
  return;
 QTextStream out(&file);
 out << "autos" << "\n";
 out << "A" << "\n";
 file.close();
 Ready(false);
 QString scmd;
 do{
     usleep(100);
     QFile filectrl(fctrl);
     if (!filectrl.open(QIODevice::ReadOnly | QIODevice::Text))
         return;
     QTextStream streamctrl ( &filectrl );
     streamctrl >> scmd;
     filectrl.close();
 } while(!scmd.contains("done!"));
 Ready(true);
}

void ksemaw::Ready(bool slwa){
    cout<<"-> set Ready "<<slwa<<endl;
    ksemaw::setEnabled(slwa);
    ksemaw::update();
}

void ksemaw::saveNKtmp(){
 QFile file(fctrl);
  if (!file.open(QIODevice::WriteOnly | QIODevice::Text))
  return;
 QTextStream out(&file);
 out << "autos" << "\n";
 out << "M" << "\n";
 file.close();
}

void ksemaw::ClearTemp(){
 QFile file(fctrl);
  if (!file.open(QIODevice::WriteOnly | QIODevice::Text))
  return;
 QTextStream out(&file);
 out << "autos" << "\n";
 out << "0" << "\n";
 file.close();
}

void ksemaw::RefIbridG(){
 SaveSetting(-1);
 QFile file(fctrl);
  if (!file.open(QIODevice::WriteOnly | QIODevice::Text))
  return;
 QTextStream out(&file);
 out << "ibrid" << "\n";
 out << "G" << "\n";
 file.close();
}

void ksemaw::IbridPlotFit(){
 SaveSetting(-1);
 QFile file(fctrl);
  if (!file.open(QIODevice::WriteOnly | QIODevice::Text))
  return;
 QTextStream out(&file);
 out << "ibrid" << "\n";
 out << "g" << "\n";
 file.close();
}

void ksemaw::IbridPlotIbrid(){
 SaveSetting(-1);
 QFile file(fctrl);
  if (!file.open(QIODevice::WriteOnly | QIODevice::Text))
  return;
 QTextStream out(&file);
 out << "ibrid" << "\n";
 out << "gi" << "\n";
 file.close();
}

void ksemaw::FitN(){
    par[32][5]=0.;//Fit n&K in IbridOne
    IbridFit();
}

void ksemaw::FitNK(){
    par[32][5]=1.;//Fit n&K in IbridOne
    IbridFit();
}

void ksemaw::FitE1E2(){
    par[10][1]=1.;//Plot epsi1 and epsi2
    par[32][5]=2.;//Fit n&K in IbridOne
    IbridFit();
}

void ksemaw::IbridFit(){
 Qt::CheckState state;
 npp=ui->sB_PAR_34_5 -> value();
 for(int j=1;j<=npp;j++){
   state=idToCheckBox["chBeParFit_"+QString::number(j)]-> checkState();
   if(ppm[j]<100 && state==Qt::Checked)
       idToCheckBox["chBeParFit_"+QString::number(j)] -> setCheckState ( Qt::Unchecked );
 }
 PanFitEnable();
 SaveSetting(-1);
 QFile file(fctrl);
  if (!file.open(QIODevice::WriteOnly | QIODevice::Text))
  return;
 QTextStream out(&file);
 out << "ibrid" << "\n";
 out << "f" << "\n";
 file.close();
}

void ksemaw::IbridOne(){
 SaveSetting(-1);
 QFile file(fctrl);
  if (!file.open(QIODevice::WriteOnly | QIODevice::Text))
  return;
 QTextStream out(&file);
 out << "ibrid" << "\n";
 out << "i" << "\n";
 file.close();
}

void ksemaw::IbridOneStore(){
 SaveSetting(-1);
 QFile file(fctrl);
  if (!file.open(QIODevice::WriteOnly | QIODevice::Text))
  return;
 QTextStream out(&file);
 out << "ibrid" << "\n";
 out << "M" << "\n";
 file.close();
}

void ksemaw::ClearTempIbri(){
 QFile file(fctrl);
  if (!file.open(QIODevice::WriteOnly | QIODevice::Text))
  return;
 QTextStream out(&file);
 out << "ibrid" << "\n";
 out << "0" << "\n";
 file.close();
}

void ksemaw::AdjSave(){
 int itab,isv;
 itab=ui->tabWidget -> currentIndex();
 if(itab==3){
  isv=ui->comboBox_StoreNK -> currentIndex();
  ui->comboBox_store -> setCurrentIndex(isv);
 } else{
  isv=ui->comboBox_store -> currentIndex();
  ui->comboBox_StoreNK -> setCurrentIndex(isv);
 }
 par[18][2]=isv;
}

void ksemaw::EnableAllWl(){
 QFile file(fctrl);
  if (!file.open(QIODevice::WriteOnly | QIODevice::Text))
  return;
 QTextStream out(&file);
 out << "ibrid" << "\n";
 out << "eR" << "\n";
 file.close();
}

void ksemaw::EnableAllNK(){
 QFile file(fctrl);
  if (!file.open(QIODevice::WriteOnly | QIODevice::Text))
  return;
 QTextStream out(&file);
 out << "ibrid" << "\n";
 out << "eN" << "\n";
 file.close();
}

void ksemaw::EditWl(){
 QFile file(fctrl);
  if (!file.open(QIODevice::WriteOnly | QIODevice::Text))
  return;
 QTextStream out(&file);
 out << "ibrid" << "\n";
 out << "em" << "\n";
 file.close();
}

void ksemaw::EditNK(){
 QFile file(fctrl);
  if (!file.open(QIODevice::WriteOnly | QIODevice::Text))
  return;
 QTextStream out(&file);
 out << "ibrid" << "\n";
 out << "e " << "\n";
 file.close();
}

void ksemaw::GoBest(){
 QFile file(fctrl);
  if (!file.open(QIODevice::WriteOnly | QIODevice::Text))
  return;
 QTextStream out(&file);
 out << "ibrid" << "\n";
 out << "R " << "\n";
 file.close();
}

void ksemaw::GoPrevious(){
 int jobview=ui->sB_PAR_8_2->value();
 int jobtot=ui->sB_PAR_8_1->value();
 if(jobview<=1)
     return;
 jobview--;
 QString ftmp=pathroot+"temp/semaw"+QString::number(jobview).rightJustified(4, '0')+".Spj";
 ReadSetting(ftmp);
 ui->sB_PAR_8_1->setValue(jobtot);
 ui->sB_PAR_8_2->setValue(jobview);
}

void ksemaw::GoNext(){
 int jobview=ui->sB_PAR_8_2->value();
 int jobtot=ui->sB_PAR_8_1->value();
 if(jobview>=jobtot)
     return;
 jobview++;
 QString ftmp=pathroot+"temp/semaw"+QString::number(jobview).rightJustified(4, '0')+".Spj";
 ReadSetting(ftmp);
 ui->sB_PAR_8_1->setValue(jobtot);
 ui->sB_PAR_8_2->setValue(jobview);
}
