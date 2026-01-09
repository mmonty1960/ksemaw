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


   Copyright (C) 2024  Marco Montecchi

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


/*
//C   par[60][5] matrix where several different parameters are stored
//C
//C      [1][1],[1][2]........NMIN, NMAX [solution search range ]
//C      [2][1],[2][2]........KMIN, KMAX [  "      "        "   ]
//C      [3][1],[3][2]........ /  ,  /
//C      [4][1],[4][2]........LAMBDAmin,LAMBDAmax of loaded data-file
//C      [5][1],[5][2]........ /  ,  /
//C      [6][1],[6][2]........TETA,  /
//C      [7][1],[7][2]........ /  ,  /
//C      [8][1],[8][2]........ jobtot ,  jobcurrent
//C      [9][1],[9][2]........if 1 k>=0, if 1 nk_sim->temp from scratch
//C     [10][1],[10][2].......if 1 plot eps1 eps2 , if 1 no 3-points limit in ELI resampling
//C     [11][1],[11][2].......attenuation coeff. of k by Fit#N , /
//C     [12][1],[12][2]....... / , /
//C     [13][1],[13][2]....... /, iDelCon (0->as computed, 1->keep connected)
//C     [14][1],[14][2]....... ThetaEli1 , /
//C     [15][1],[15][2]....... ThetaEli2 , /
//C     [16][1],[16][2]....... ThetaEli3 , /
//C     [17][1],[17][2]....... ThetaEli4 , /
//C     [18][1],[18][2].......verbose [0=no, 1=yes], /
//C     [19][1],[19][2]....... / , /
//C     [20][1],[20][2]....... / , /
//C     [21][1],[21][2].......mwla[N. wl enabled], mda[N. nk enabled]
//C     [22][1],[22][2].......mr[0=>R absolute, 1=> R*Rrif] , Nmis enabled
//C     [23][1],[23][2].......1 to load SF files, 1 to load ELI files  (0 no files to load)
//C     [24][1],[24][2].......IUVIR,index-lambda
//C     [25][1],[25][2]....... / , /
//C     [26][1],[26][2]....... / , /
//C     [27][1],[27][2].......chi2_initial ,S1P2
//C     [28][1] .............N discretization for computing <f[x]>
//C            [28][2].......mwl [Ntotal wl]
//C     [29][1],.............N sublayers modeling a single inhomogeneous film
//C            [29][2]....... /
//C     [30][1-5] ...........{/}
//C     [31][1-    4]........{/}
//C     [32][1-    4]........{/}
//C     [33][1-    4]........{/}
//C     [34][1-    4]........{/}
//C     [35][1-    4]... standard-spectrum_for_weighting , S1P2U3 , / , CNK[1][1] in simulation
//C
//C     [51-60][1]      pointer to VNK -> nk layer #j=1,9
//C         [51][2]     Nlayer [max 9]
//C         [52][2]     symmetric multilayer on the back face of the last layer of Nlayer:
//C                     0 = no, 1 = yes
//C         [53][2]     pointer to the unknown layer
//C         [54][2]     type of SF spectra: specular/direct 0, hemispherical 1
//C         [55][2]     chi2fin [last]
//C         [56][2]     fredeg [degree of freedom]
//C       [51-60][3]    type of layer: 1 = bulk [not coherent]
//C                                    2 = homogeneous film
//C                                    3 = inhomogeneous film
//C       [ 1-20][3]    INSTR[20] [see below]
//C       [21-22][3]    Delta_n , Delta_k [numerical increment]
//C        [ 1-14][4]   DATO[1-14]
//C        [21-35][4]   selected number of a given measurement; as an example in baseName.v2.tn is 2
//C            [31][5]  log scale of k-plot [0->linear 1->logarithmic]
//C            [32][5]  0->fit_n   1->fit_nk   2->fit epsi1 and epsi2
//C            [33][5]  Delta & Psi units [0-> Delta(deg) Psi(deg) 1-> cos(Delta) tan (Psi)
//C            [34][5]  N parameters displayed in dat-fit tab
//C            [35][5]  N parameters enabled to fit
//C
//C      [36][1-3] Tn:    Delta/error[rms] <Experimental>  <Simulation>
//C      [37][1-3] Tp:    Delta/error[rms] <Experimental>  <Simulation>
//C      [38][1-    4]... WLmin_solnk, WLmax_solnk,WLmis_mis,WLmax_mis enabled
//C        .....
//C      [49][1-3] PSI-4: Delta/error[rms] <Experimental>  <Simulation>
//C
//C         [36-53][5]  PPM[17]
//C
//C
//C        IUVIR : kind of wavelength range = 1 UV-VIS-NIR
//C                                           2 IR
//C
//C
//C   INSTR[20] instrument parameters saved in par:
//C                         INSTR[i]=par[i][3]
//C   INSTR[1] : a[deff]= efficiency of capture of back face reflection [1st order]
//C   INSTR[2] : a[2*deff]= idem but 2nd order
//C   INSTR[3] : DBK/BK= drift of SF baseline
//C   INSTR[4] : DRcal/Rcal=relative error of reference mirror
//C   INSTR[5] : Reading error
//C   INSTR[6] :   "      "     not [0] or depending [1] on wavelength
//C   INSTR[7] : cte substrate contribution to PDS signal
//C   INSTR[8] : /
//C   INSTR[9] : /
//C   INSTR[10]: 0<-> absolute     user is allowed to customzize the refence mirror list
//C              1<-> RIF05
//C              2<-> RIF06 before 5/dic/94
//C              3<-> RIF06 after  5/dic/94
//C              4<-> RIF05 after 18/mar/96
//C              5<-> RIF08
//C              6<-> RIF08 after 17/gen/2018
//C   INSTR[11]: /
//C   INSTR[12]: /
//C   INSTR[13]: /
//C   INSTR[14]: /
//C   INSTR[15]: /
//C   INSTR[16]: /
//C   INSTR[17]: /
//C   INSTR[18]: /
//C   INSTR[19]: /
//C   INSTR[20]: /
//C
//C
//C
//C   DATO[14] loaded experimental measurements saved in par:
//C                         DATO[i]=par[i][4]
//C       i=1: Tnormal
//C       i=2: Tpolarised
//C       i=3: Rnormal
//C       i=4: Rpolarised
//C       i=5: R1_normal
//C       i=6: Apds
//C       i=7: DELTA_1
//C       i=8: PSI_1
//C      ....
//C       i=14: PSI_4
//C
//C       DATO[i] = 0 -> no measure
//C                 1 -> the measure is loaded
//C                 2 -> the measure is loaded and is processed
//C
//C
//C
//C
//C
//c **** pm[201][6] generalized multilayer & oscillators ****
//c        [i][1] = parameter value
//c        [i][2] = managment in data-fit: 0 -> constant
//c                                        j -> P[j] where P is a fit parameter
//c        [j][3] = pointer P: P[j] = PM[int(PM[j][3])][1]
//c        [i][4] = ERROR
//c        [i][5] = global correlation
//c
//c   1- 9 thickness
//c  11-19 Dn/<n> n-gradient
//c  21-29 n-curvature
//c  31-39 Dk/<k> k-gradient
//c  41-49 k-curvature
//c  51-59 roughness
//c  61-69 slope_Dn/<n> in 1/eV
//c  71-85 f-EMA
//c  86-89 ThetaEli_1, ThetaEli_2,ThetaEli_3,ThetaEli_4
//c  91-99 nonunifornity of thickness (%)
//c 100  Fit# used in data-fit [1,2 ...,7]
//c 101  fdisp#1
//c 102  C1
//c 103  E1
//c 104  D1
//c 105  K1
//c etc etc up to the 20th oscillator
//C
//C
//C
//C   ppm[nparmax=17] pointer pm stored in par:
//C                       PPM[i]=par[35+i][5] i=1,17
//C                       n_parameters_fit=par[35][5]
//C
//C
//C
//C   MIS[18][1000][3] matrix SF spectra and nk-data loaded as file-nk
//C      [i][j][k]
//C       i=1: Tnormal     j=1,1001 <-> Lambda   k=1 <-> value
//C       i=2: Tpolarised                         2 <-> error
//C       i=3: Rnormal
//C       i=4: Rpolarised
//C       i=5: R1_normal
//C       i=6: Apds
//C       i=7: Lambda
//C       i=8-15:  material #1-8   k=1 <-> n
//C                                  2 <-> k
//C       i=16: n[k=1],k[k=2] IbridOne
//C       i=17: enabling in data-fit
//C
//C   ELI[8][1000][2] matrix ellipsometric spectra
//C      [i][j][k]
//C       i=1: DELTA_1   j=1,1001 <-> Lambda   k=1 <-> value
//C       i=2: PSI_1                            2 <-> error
//C      .. ..
//C       i=7: DELTA_4
//C       i=8: PSI_4
//C
//C
//C   NANK[1,..,8] : nome "mate/aa999.9" nk-known
//C   NANK[9]      : nome "mate/aa999.9" nk-solution
//C   NANK[10]     : nome "mate/aa999  " custom weights for mean computing
//C   NANK[11]     : nome "mate/aa999" SF spectra
//C   NANK[12,.,15]: nome "mate/aa999" ELI spectra
//C   NANK[16]     : nome "mate/aa999.9" ksemaw project
//C
//C   CNK[17][7]
//C   CNK[J][K]: managment matrix of VNK[J][H]
//C       J=1    <-> VNK[1][1 o 2]    = n,k working
//C       J=2..8 <-> VNK[2..8][1 o 2] = n,k known
//C       J=10   <-> VNK[10][1 o 2]   = n,k input medium Tnormal,..,R1_normal
//C       J=11   <-> VNK[11][1 o 2]   = n,k input medium PDS
//C       J=12   <-> VNK[12][1 o 2]   = n,k input medium PSI,DELTA #1
//C       J=13   <-> VNK[13][1 o 2]   = n,k input medium PSI,DELTA #2
//C       J=14   <-> VNK[14][1 o 2]   = n,k input medium PSI,DELTA #3
//C       J=15   <-> VNK[15][1 o 2]   = n,k input medium PSI,DELTA #
//C       J=16   <-> VNK[9][1 o 2]    = n,k output medium Tnormal,..,R1_normal
//C    CNK[J][1 and 4]=0       =>  n,k cte  n=CNK[J][2 and 5]  k=CNK[J][3 and 6]
//C    CNK[J][1 and 4]=1..5    =>  n,k by  Fit#
//C    CNK[J][1 and 4]=8..15   =>  n,k by file-nk loaded in MIS[K][L][NK]
//C    CNK[J][1 and 4]=16      =>  n,k Material#1
//C
//C
//C    rxy[30][4]: managment matrix of plots
//C       [J ][K]
//C
//C       J= 1 <-> Tnormal
//C       J= 2 <-> Tpolarised
//C       J= 3 <-> Rnormal
//C       J= 4 <-> Rpolarised
//C       J= 5 <-> R1_normal
//C       J= 6 <-> Apds
//C       J= 7 <-> Delta_1
//C       J= 8 <-> Psi_1
//C       J= 9 <-> Delta_2
//C       J=10 <-> Psi_2
//C       J=11 <-> Delta_3
//C       J=12 <-> Psi_3
//C       J=13 <-> Delta_4
//C       J=14 <-> Psi_4
//C       J=15 <-> A_back=1-Tn-R1
//C       J=16 <-> n
//C       J=17 <-> k
//C       J=18 <-> A_front=1-Tn-Rn
//C       J=19 <-> cos(Delta)
//C       J=20 <-> Lambda [nm]
//C       J=21 <-> Teta [deg]
//C       J=22 <-> tan(Psi)
//C       J=23 <-> Tau,Rho,Rho1
//C       J=24 <-> Size of graph [width & height & aspect]
//C       J=25 <->  /, /   ,wl/eV, spline in eV
//C       J=26 <-> eps1
//C       J=27 <-> eps2
//C
//C          for 1=< J <=23 && 26,27
//C              rxy[J][1]  = Xmin in plot
//C              rxy[J][2]  = Xmax in plot
//C              rxy[J][3]  = Vmin, minimum value of J matter
//C              rxy[J][4]  = Vmax, maximum value of J matter
//C         int(rxy[25][1]) = index of matter plotted along X axes
//C         int(rxy[25][2]) = index of matter plotted along Y axes
//C         int(rxy[25][3]) = 1 => nm
//C         int(rxy[25][3]) = 2 => eV
//C         int(rxy[25][4]) = 1 => nm step
//C         int(rxy[25][4]) = 2 => eV step
//C
//C    ARSE[500][2] service array
*/


#include "ksemawc.h"
#include <stdio.h>
#include <stdlib.h>
#include <qfile.h>
#include <qtextstream.h>
#include <math.h>
#include <complex>
#include <iostream>
#include <QFileDialog>
#include <QInputDialog>
#include <QMessageBox>
#include <unistd.h>
#include <cminpack.h>
#include <QtGlobal>
#include <qwt_plot.h>
#include <qwt_plot_curve.h>
#include <qwt_series_data.h>
#include <qwt_plot_grid.h>
#include <qwt_legend.h>
#include <qwt_plot_item.h>
#include <qwt_plot_marker.h>
#include <qwt_symbol.h>
#include <qwt_interval.h>
#include <qwt_plot_intervalcurve.h>
#include <qwt_interval_symbol.h>
#include <QwtPickerMachine>
#include <qwt_picker_machine.h>
#include <qwt_picker.h>
#include <qwt_plot_picker.h>
#include <QwtPlotZoomer>
#include <QDateTime>
#include <QPen>
#include <QGraphicsPolygonItem>
#include <QPolygonF>
#include <QFontDialog>
#include <QApplication>
#include <QFontInfo>
#include <QSettings>
#include <gsl/gsl_sf_result.h>
#include <gsl/gsl_complex_math.h>
#include <gsl/gsl_specfunc.h>



//global variables
QString fnproject,pathroot,fileStore,fileStore0,fStdSpect,fRefMir,fNKsim,fMisSim,filechi2,lastAction,fileStoreSpjName,caller;
QString info,fnk[9],fnSample,fnTn,fnTp,fnRn,fnRp,fnR1,fnApds,fnE1,fnE2,fnE3,fnE4,fnFnk,ParFitLab[14],NANK[17];

QwtPlot *G1_Tn, *G2_Tp, *G3_Rn, *G4_Rp, *G5_R1, *G6_Apds, *G7_D, *G8_P, *G9_Af, *G10_tTR, *G11_nk, *G12_wn, *G13_wk, *G14_we1, *G15_we2, *G16_Ab;
QPolygon polygon;
QPolygonF polygonF;
QColor myColor[7]={Qt::black,Qt::blue,Qt::cyan,Qt::green,Qt::magenta,Qt::red,Qt::yellow};
QString labelQwt;//label to use in absorptance plot
int iColor=0;
int iSelected=0;//used to stop recursive n-selection by polygon

int ink,lastIndex,ifn,npp,ppm[18],nPar,nlayer,lastTab,occupyPF,occupyPEj,DATO[15],IXW[17],iw=0,L1E2=2,iRecChi2=0;
int ifirstcall=0;//used to initialize fit
int ifirstWarning=0;//used to warn about A= 1-T-R<0
int jobtot=0;
int nOpenGraph=3;
int NeV; // number of interpolated points in Ev or lambda
int iDelCon=0;
int iNewSample=0;
int iwl2print=100;//index of WL at which print when it is due

double MIS[18][1000][3],ELI[9][1000][3],pf[8][22],pm[201][6],par[61][6],CNK[17][7],rxy[31][5],ARSE[501][3];
double SOL[4000][7]; // solutions from the search algorithm
double PIG=acos(-1.);
double deg2rad=PIG/180.;
double Dperiod=180.;//period of DELTA
double Nema,Kema,sqn,sqk,EpsiR,EpsiI;
double Vservice[3];//for transporting
double ABSmax=20.;//max value of abs(imag(DELTA)), beyond the layer is considered as exit medium
double cDAW[7];//c[NMAX+1] with NMAX=6 use in DAWS function
double w[6][7]={
    {0.,0.              ,0.             ,0.         ,0.         ,0.         ,0.},
    {0.,9.78264917E-1   ,1.08675416E-2  ,0.         ,0.         ,0.         ,0.},
    {0.,5.9826E-1       ,1.9423E-1      ,6.6460E-3  ,0.         ,0.         ,0.},
    {0.,0.39905         ,0.242036       ,5.40056E-2 ,4.43305E-3 ,0.         ,0.},
    {0.,2.99373E-1      ,2.25978E-1     ,9.7192E-2  ,2.38179E-2 ,3.32573E-3 ,0.},
    {0.,2.39559E-1      ,2.00097E-1     ,1.16606E-1 ,4.74085E-2 ,1.34476E-2 ,2.66126E-3}
};
double dx[6][7]={
    {0.,0. ,0.   ,0. ,0.   ,0.  ,0.},
    {0.,0. ,3.   ,0. ,0.   ,0.  ,0.},
    {0.,0. ,1.5  ,3. ,0.   ,0.  ,0.},
    {0.,0. ,1.   ,2. ,3.   ,0.  ,0.},
    {0.,0. ,0.75 ,1.5,2.25 ,3.  ,0.},
    {0.,0. ,0.6  ,1.2,1.8  ,2.4 ,3.}
};
double EliTab[5][11][4];//summary table of the selected ELI data files
//           [j][k][l]
//            j=1, ., 4 index of loadable ELI slots (0 is not used)
//           [j][0][0] ->nTheta at j
//           [j][k][l] k=1,..,9 is the theta index in the ELI selected file
//                  l=0->theta
//                    1->nDat
//                    2->wlMin
//                    3->wlMax

complex<double> freCoeff[6][2];// Fresnell coefficient
//                     Poynting-s  Poynting-p
//                     rhoS         rhoP
//                     rho1S        rho1P
//                     n-i*k[NFA]   /
//                     Bs           Cs
//                     Bp           Cp

bool iStop=false;

//invoked functions
void previewFile(QString filename, QString lab,QString& info,double& wmax,double& wmin);
int FPAR(void *p, int m, int n, const double *x, double *fvec, int iflag);
void CONVER(double X[12000],double Y[12000],int NDATI,int N1,int STEP,int I,int IUVIR);
void PLOTline1bar2(int iL1B2,int iRD,int iCol,int ic,int Ndata,double *Xp,double *Yp,double *ErrXp,double *ErrYp);
void COSVNK(double VNK[17][3],int L);
void SETVNK(int io,int J,double VNK[17][3],int L,int Kcte);
void EMA(double N,double K,double NA,double KA,double FA);
void FDISP(int iopt,double eV);
long double I2CODYRE(long double E,long double E0,long double Ep,long double D);
long double ICCODYRE(long double E,long double Ep,long double D);
long double I2CODYIM(long double E,long double E0,long double Ep,long double D);
long double ICCODYIM(long double E,long double Ep,long double D);
long double ICTAUCRE(long double E,long double Ep,long double D);
long double I2TAUCRE(long double E,long double E0,long double Ep,long double D);
long double ICTAUCIM(long double E,long double Ep,long double D);
long double I2TAUCIM(long double E,long double E0,long double Ep,long double D);
long double DAWS(long double x);
long double ReFdirGap(long double E,long double E0,long double D,long double x);
long double aImFdirGap(long double E,long double E0,long double D,long double x);
long double ReM0M3(long double E,long double E0,long double E3,long double reD);
long double aImM0M3(long double E,long double E0,long double E3,long double reD);
long double ReCodyM1M2(long double E,long double C,long double E0,long double E3,long double reD);
long double ImCodyM1M2(long double E,long double C,long double E0,long double E3,long double reD);
long double ReTaucM1M2(long double E,long double C,long double E0,long double E3,long double reD);
long double ImTaucM1M2(long double E,long double C,long double E0,long double E3,long double reD);
long double ReTaucL(long double E,long double E0,long double D,long double E3);
long double aImTaucL(long double E,long double E0,long double D,long double E3);
long double ReDirCodyUrb(long double E,long double C,long double E0,long double E3,long double D);
long double ImDirCodyUrb(long double E,long double C,long double E0,long double E3,long double D);
long double ReDirTaucUrb(long double E,long double C,long double E0,long double E3,long double D);
long double ImDirTaucUrb(long double E,long double C,long double E0,long double E3,long double D);
long double ImIndirCodyUrb(long double E,long double C,long double E0,long double E3,long double D);
long double ReIndirCodyUrb(long double E,long double C,long double E0,long double E3,long double D);
long double ImIndirTaucUrb(long double E,long double C,long double E0,long double E3,long double D,long double Ex,long double Ey,long double TUx,long double TUy, long double KTi);
long double ReIndirTaucUrb(long double E,long double C,long double E0,long double E3,long double D,long double Ex,long double Ey,long double TUx,long double TUy, long double KTi);
void ASSEMBLER(int iwl, int ikind,double teta,double vot[6][3]);
void BUILDER(int iwl,int ikind,int ifst,int ncoe, complex<double> pq,double vosi[6][3]);
void CALFRE(int NFA,double wl,complex<double> pq,complex<double> IR[999],double d[999],complex<double> out[9][3]);
int SOLVE(int imis,int iWL,double *Xp,double *Yp,double *ErrYp);
double FMER(double k);
double FindRoot(double (*FT)(double),double t,double dt,double dlim,double tol);
int FSQ(void *p, int m, int n, const double *x, double *fvec, int iflag);
int FSEM(void *p, int m, int n, const double *x, double *fvec, int iflag);
int FRCK(void *p, int m, int n, const double *x, double *fvec, int iflag);
double DELTAT(double k);
int NINT(double x);
void MATINV(int n, int np,double **as,double **b);
void lubksb(int np,double **a,int n,int *indx,double *b);
void ludcmp(int np,double **a,int n,int *indx,double d);
//void AbsLayInco(double atten, double Ra, double R1a, double Rb);
void nextColor();

static struct pointToFit2{
    int Nt;
}pTF2[1];


ksemawc::ksemawc(QWidget *parent) :
    QWidget(parent),
    ui(new Ui::ksemawc)
{
    //setupUi(this); // this sets up GUI
    ui->setupUi(this);
    printf("****************************************************\n");
    printf("              Program C++ kSEMAW\n\n");
    printf("Spectro-Ellipsometric Measurement Analysis Workbench\n");
    printf("  (spectrophotometric, ellipsometric and PDS)\n\n");
    printf("         version 3.4 9 January 2026\n\n");
    printf("       Main author: Marco Montecchi, ENEA (Italy)\n");
    printf("          email: marco.montecchi@enea.it\n");
    printf("          Porting to Windows and advanced oscillators by\n");
    printf("               Alberto Mittiga, ENEA (Italy)\n");
    printf("          email: alberto.mittiga@enea.it\n ");
    printf("****************************************************\n");

//    const QByteArray value = qgetenv("USER");
//    QString uName=QString::fromLocal8Bit(value);
//    cout << "current user = " << uName.toStdString() <<endl;

    // signals/slots mechanism in action
    connect( ui->pBloadPro,  SIGNAL( clicked() ), this, SLOT(LoadProject()));
    connect( ui->pBsavePro,  SIGNAL( clicked() ), this, SLOT(SaveProject()));
    connect( ui->pBloadFnk,  SIGNAL( clicked() ), this, SLOT(LoadFilenk()));
    connect( ui->pBclearFnk, SIGNAL( clicked() ), this, SLOT(ClrFnk()));
    connect( ui->pBnk1, &QPushButton::clicked, this, [this]() {Setnk(1);});
    connect( ui->pBnk2, &QPushButton::clicked, this, [this]() {Setnk(2);});
    connect( ui->pBnk3, &QPushButton::clicked, this, [this]() {Setnk(3);});
    connect( ui->pBnk4, &QPushButton::clicked, this, [this]() {Setnk(4);});
    connect( ui->pBnk5, &QPushButton::clicked, this, [this]() {Setnk(5);});
    connect( ui->pBnk6, &QPushButton::clicked, this, [this]() {Setnk(6);});
    connect( ui->pBnk7, &QPushButton::clicked, this, [this]() {Setnk(7);});
    connect( ui->pBnk8, &QPushButton::clicked, this, [this]() {Setnk(8);});
    connect( ui->pBclearnk1, &QPushButton::clicked, this, [this]() {Clrnk(1);});
    connect( ui->pBclearnk2, &QPushButton::clicked, this, [this]() {Clrnk(2);});
    connect( ui->pBclearnk3, &QPushButton::clicked, this, [this]() {Clrnk(3);});
    connect( ui->pBclearnk4, &QPushButton::clicked, this, [this]() {Clrnk(4);});
    connect( ui->pBclearnk5, &QPushButton::clicked, this, [this]() {Clrnk(5);});
    connect( ui->pBclearnk6, &QPushButton::clicked, this, [this]() {Clrnk(6);});
    connect( ui->pBclearnk7, &QPushButton::clicked, this, [this]() {Clrnk(7);});
    connect( ui->pBclearnk8, &QPushButton::clicked, this, [this]() {Clrnk(8);});
    connect( ui->pBsetSample,SIGNAL( clicked() ), this, SLOT(callSetSample()));
    connect( ui->pBclearFn,SIGNAL( clicked() ), this, SLOT(Clrfn()));
    connect( ui->lineEdit_sample,SIGNAL(textChanged(QString)),this, SLOT(listMeas(QString)));
    connect( ui->cBmis1,SIGNAL(currentIndexChanged(int)),this,SLOT(pwTn(int)));
    connect( ui->cBmis2,SIGNAL(currentIndexChanged(int)),this,SLOT(pwTp(int)));
    connect( ui->cBmis3,SIGNAL(currentIndexChanged(int)),this,SLOT(pwRn(int)));
    connect( ui->cBmis4,SIGNAL(currentIndexChanged(int)),this,SLOT(pwRp(int)));
    connect( ui->cBmis5,SIGNAL(currentIndexChanged(int)),this,SLOT(pwR1(int)));
    connect( ui->cBmis6,SIGNAL(currentIndexChanged(int)),this,SLOT(pwApds(int)));
    connect( ui->cBmis7 ,qOverload<int>(&QComboBox::currentIndexChanged), this, [this]() {pwEj(1);});
    connect( ui->cBmis9 ,qOverload<int>(&QComboBox::currentIndexChanged), this, [this]() {pwEj(2);});
    connect( ui->cBmis11,qOverload<int>(&QComboBox::currentIndexChanged), this, [this]() {pwEj(3);});
    connect( ui->cBmis13,qOverload<int>(&QComboBox::currentIndexChanged), this, [this]() {pwEj(4);});
    connect( ui->cBteE1,qOverload<int>(&QComboBox::currentIndexChanged), this, [this]() {pwSubEj(1);});
    connect( ui->cBteE2,qOverload<int>(&QComboBox::currentIndexChanged), this, [this]() {pwSubEj(2);});
    connect( ui->cBteE3,qOverload<int>(&QComboBox::currentIndexChanged), this, [this]() {pwSubEj(3);});
    connect( ui->cBteE4,qOverload<int>(&QComboBox::currentIndexChanged), this, [this]() {pwSubEj(4);});
    connect(ui->comboB_PAR_10_3,SIGNAL(currentIndexChanged(int)),this,SLOT(setRifMir()));
    connect(ui->checkB_PAR_22_1,SIGNAL(stateChanged(int) ),this, SLOT(setRifMir()));
    connect(ui->checkB_mis1_1,SIGNAL(stateChanged(int)),this,SLOT(MCRange()));
    connect(ui->checkB_mis2_1,SIGNAL(stateChanged(int)),this,SLOT(MCRange()));
    connect(ui->checkB_mis3_1,SIGNAL(stateChanged(int)),this,SLOT(MCRange()));
    connect(ui->checkB_mis4_1,SIGNAL(stateChanged(int)),this,SLOT(MCRange()));
    connect(ui->checkB_mis5_1,SIGNAL(stateChanged(int)),this,SLOT(MCRange()));
    connect(ui->checkB_mis6_1,SIGNAL(stateChanged(int)),this,SLOT(MCRange()));
    connect(ui->checkB_mis7_1,SIGNAL(stateChanged(int)),this,SLOT(MCRange()));
    connect(ui->checkB_mis9_1,SIGNAL(stateChanged(int)),this,SLOT(MCRange()));
    connect(ui->checkB_mis11_1,SIGNAL(stateChanged(int)),this, SLOT(MCRange()));
    connect(ui->checkB_mis13_1,SIGNAL(stateChanged(int)),this, SLOT(MCRange()));
    connect(ui->cB_cnk1a, qOverload<int>(&QComboBox::currentIndexChanged), this, [this]() {
        setMat(1);
        listOsc();
    });
    connect(ui->cB_cnk1b,qOverload<int>(&QComboBox::currentIndexChanged), this, [this]() {setMat(1);});
    connect(ui->cB_cnk2a,qOverload<int>(&QComboBox::currentIndexChanged), this, [this]() {setMat(2);});
    connect(ui->cB_cnk3a,qOverload<int>(&QComboBox::currentIndexChanged), this, [this]() {setMat(3);});
    connect(ui->cB_cnk4a,qOverload<int>(&QComboBox::currentIndexChanged), this, [this]() {setMat(4);});
    connect(ui->cB_cnk5a,qOverload<int>(&QComboBox::currentIndexChanged), this, [this]() {setMat(5);});
    connect(ui->cB_cnk6a,qOverload<int>(&QComboBox::currentIndexChanged), this, [this]() {setMat(6);});
    connect(ui->cB_cnk7a,qOverload<int>(&QComboBox::currentIndexChanged), this, [this]() {setMat(7);});
    connect(ui->cB_cnk8a,qOverload<int>(&QComboBox::currentIndexChanged), this, [this]() {setMat(8);});
    connect(ui->cB_cnk9a,qOverload<int>(&QComboBox::currentIndexChanged), this, [this]() {setMat(9);});
    connect(ui->cB_cnk10a,qOverload<int>(&QComboBox::currentIndexChanged), this, [this]() {setMat(10);});
    connect(ui->cB_cnk11a,qOverload<int>(&QComboBox::currentIndexChanged), this, [this]() {setMat(11);});
    connect(ui->cB_cnk12a,qOverload<int>(&QComboBox::currentIndexChanged), this, [this]() {setMat(12);});
    connect(ui->cB_cnk13a,qOverload<int>(&QComboBox::currentIndexChanged), this, [this]() {setMat(13);});
    connect(ui->cB_cnk14a,qOverload<int>(&QComboBox::currentIndexChanged), this, [this]() {setMat(14);});
    connect(ui->cB_cnk15a,qOverload<int>(&QComboBox::currentIndexChanged), this, [this]() {setMat(15);});
    connect(ui->cB_cnk1b,qOverload<int>(&QComboBox::currentIndexChanged), this, [this]() {setEMA(1);});
    connect(ui->cB_cnk2b,qOverload<int>(&QComboBox::currentIndexChanged), this, [this]() {setEMA(2);});
    connect(ui->cB_cnk3b,qOverload<int>(&QComboBox::currentIndexChanged), this, [this]() {setEMA(3);});
    connect(ui->cB_cnk4b,qOverload<int>(&QComboBox::currentIndexChanged), this, [this]() {setEMA(4);});
    connect(ui->cB_cnk5b,qOverload<int>(&QComboBox::currentIndexChanged), this, [this]() {setEMA(5);});
    connect(ui->cB_cnk6b,qOverload<int>(&QComboBox::currentIndexChanged), this, [this]() {setEMA(6);});
    connect(ui->cB_cnk7b,qOverload<int>(&QComboBox::currentIndexChanged), this, [this]() {setEMA(7);});
    connect(ui->cB_cnk8b,qOverload<int>(&QComboBox::currentIndexChanged), this, [this]() {setEMA(8);});
    connect(ui->cB_cnk9b,qOverload<int>(&QComboBox::currentIndexChanged), this, [this]() {setEMA(9);});
    connect(ui->cB_cnk10b,qOverload<int>(&QComboBox::currentIndexChanged), this, [this]() {setEMA(10);});
    connect(ui->cB_cnk11b,qOverload<int>(&QComboBox::currentIndexChanged), this, [this]() {setEMA(11);});
    connect(ui->cB_cnk12b,qOverload<int>(&QComboBox::currentIndexChanged), this, [this]() {setEMA(12);});
    connect(ui->cB_cnk13b,qOverload<int>(&QComboBox::currentIndexChanged), this, [this]() {setEMA(13);});
    connect(ui->cB_cnk14b,qOverload<int>(&QComboBox::currentIndexChanged), this, [this]() {setEMA(14);});
    connect(ui->cB_cnk15b,qOverload<int>(&QComboBox::currentIndexChanged), this, [this]() {setEMA(15);});
    connect(ui->cB_EMA_1,&QCheckBox::stateChanged, this, [this]() {setEMA(1);});
    connect(ui->cB_EMA_2,&QCheckBox::stateChanged, this, [this]() {setEMA(2);});
    connect(ui->cB_EMA_3,&QCheckBox::stateChanged, this, [this]() {setEMA(3);});
    connect(ui->cB_EMA_4,&QCheckBox::stateChanged, this, [this]() {setEMA(4);});
    connect(ui->cB_EMA_5,&QCheckBox::stateChanged, this, [this]() {setEMA(5);});
    connect(ui->cB_EMA_6,&QCheckBox::stateChanged, this, [this]() {setEMA(6);});
    connect(ui->cB_EMA_7,&QCheckBox::stateChanged, this, [this]() {setEMA(7);});
    connect(ui->cB_EMA_8,&QCheckBox::stateChanged, this, [this]() {setEMA(8);});
    connect(ui->cB_EMA_9,&QCheckBox::stateChanged, this, [this]() {setEMA(9);});
    connect(ui->cB_EMA_10,&QCheckBox::stateChanged, this, [this]() {setEMA(10);});
    connect(ui->cB_EMA_11,&QCheckBox::stateChanged, this, [this]() {setEMA(11);});
    connect(ui->cB_EMA_12,&QCheckBox::stateChanged, this, [this]() {setEMA(12);});
    connect(ui->cB_EMA_13,&QCheckBox::stateChanged, this, [this]() {setEMA(13);});
    connect(ui->cB_EMA_14,&QCheckBox::stateChanged, this, [this]() {setEMA(14);});
    connect(ui->cB_EMA_15,&QCheckBox::stateChanged, this, [this]() {setEMA(15);});
    connect(ui->cBosc_1,&QCheckBox::stateChanged, this, [this]() {setOscN(1);});
    connect(ui->cBosc_2,&QCheckBox::stateChanged, this, [this]() {setOscN(2);});
    connect(ui->cBosc_3,&QCheckBox::stateChanged, this, [this]() {setOscN(3);});
    connect(ui->cBosc_4,&QCheckBox::stateChanged, this, [this]() {setOscN(4);});
    connect(ui->cBosc_5,&QCheckBox::stateChanged, this, [this]() {setOscN(5);});
    connect(ui->cBosc_6,&QCheckBox::stateChanged, this, [this]() {setOscN(6);});
    connect(ui->cBosc_7,&QCheckBox::stateChanged, this, [this]() {setOscN(7);});
    connect(ui->cBosc_8,&QCheckBox::stateChanged, this, [this]() {setOscN(8);});
    connect(ui->cBosc_9,&QCheckBox::stateChanged, this, [this]() {setOscN(9);});
    connect(ui->cBosc_10,&QCheckBox::stateChanged, this, [this]() {setOscN(10);});
    connect(ui->cBosc_11,&QCheckBox::stateChanged, this, [this]() {setOscN(11);});
    connect(ui->cBosc_12,&QCheckBox::stateChanged, this, [this]() {setOscN(12);});
    connect(ui->cBosc_13,&QCheckBox::stateChanged, this, [this]() {setOscN(13);});
    connect(ui->cBosc_14,&QCheckBox::stateChanged, this, [this]() {setOscN(14);});
    connect(ui->cBosc_15,&QCheckBox::stateChanged, this, [this]() {setOscN(15);});
    connect(ui->cBosc_16,&QCheckBox::stateChanged, this, [this]() {setOscN(16);});
    connect(ui->cBosc_17,&QCheckBox::stateChanged, this, [this]() {setOscN(17);});
    connect(ui->cBosc_18,&QCheckBox::stateChanged, this, [this]() {setOscN(18);});
    connect(ui->cBosc_19,&QCheckBox::stateChanged, this, [this]() {setOscN(19);});
    connect(ui->cBosc_20,&QCheckBox::stateChanged, this, [this]() {setOscN(20);});
    connect(ui->tabWidget,SIGNAL(currentChanged(int)),this,SLOT(tabChanged()));
    connect(ui->dSB_PAR_4_1,SIGNAL(valueChanged(double)),this,SLOT(rangeWL()));
    connect(ui->dSB_PAR_4_2,SIGNAL(valueChanged(double)),this,SLOT(rangeWL()));
    connect(ui->dSB_PAR_6_1,SIGNAL(valueChanged(double)),this,SLOT(AdjTheta()));
    connect(ui->dSB_PM_86_1,SIGNAL(valueChanged(double)),this,SLOT(AdjTheta()));
    connect(ui->dSB_PM_87_1,SIGNAL(valueChanged(double)),this,SLOT(AdjTheta()));
    connect(ui->dSB_PM_88_1,SIGNAL(valueChanged(double)),this,SLOT(AdjTheta()));
    connect(ui->dSB_PM_89_1,SIGNAL(valueChanged(double)),this,SLOT(AdjTheta()));
    connect(ui->dSB_PM_51_1,SIGNAL(valueChanged(double)),this,SLOT(AdjRoughMax()));
    connect(ui->dSB_PM_52_1,SIGNAL(valueChanged(double)),this,SLOT(AdjRoughMax()));
    connect(ui->dSB_PM_53_1,SIGNAL(valueChanged(double)),this,SLOT(AdjRoughMax()));
    connect(ui->dSB_PM_54_1,SIGNAL(valueChanged(double)),this,SLOT(AdjRoughMax()));
    connect(ui->dSB_PM_55_1,SIGNAL(valueChanged(double)),this,SLOT(AdjRoughMax()));
    connect(ui->dSB_PM_56_1,SIGNAL(valueChanged(double)),this,SLOT(AdjRoughMax()));
    connect(ui->dSB_PM_57_1,SIGNAL(valueChanged(double)),this,SLOT(AdjRoughMax()));
    connect(ui->dSB_PM_58_1,SIGNAL(valueChanged(double)),this,SLOT(AdjRoughMax()));
    connect(ui->dSB_PM_59_1,SIGNAL(valueChanged(double)),this,SLOT(AdjRoughMax()));
    connect(ui->comB_PAR_51_3,SIGNAL(currentIndexChanged(int)),this,SLOT(RefreshModel()));
    connect(ui->comB_PAR_52_3,SIGNAL(currentIndexChanged(int)),this,SLOT(RefreshModel()));
    connect(ui->comB_PAR_53_3,SIGNAL(currentIndexChanged(int)),this,SLOT(RefreshModel()));
    connect(ui->comB_PAR_54_3,SIGNAL(currentIndexChanged(int)),this,SLOT(RefreshModel()));
    connect(ui->comB_PAR_55_3,SIGNAL(currentIndexChanged(int)),this,SLOT(RefreshModel()));
    connect(ui->comB_PAR_56_3,SIGNAL(currentIndexChanged(int)),this,SLOT(RefreshModel()));
    connect(ui->comB_PAR_57_3,SIGNAL(currentIndexChanged(int)),this,SLOT(RefreshModel()));
    connect(ui->comB_PAR_58_3,SIGNAL(currentIndexChanged(int)),this,SLOT(RefreshModel()));
    connect(ui->comB_PAR_59_3,SIGNAL(currentIndexChanged(int)),this,SLOT(RefreshModel()));
    connect(ui->sB_PAR_51_2,SIGNAL(valueChanged(int)),this,SLOT(SetModel(int)));
    connect(ui->sB_PAR_34_5,SIGNAL(valueChanged(int)),this,SLOT(PanFitPar()));
    connect(ui->chBeParFit_1,SIGNAL(stateChanged(int)),this,SLOT(PanFitEnable()));
    connect(ui->chBeParFit_2,SIGNAL(stateChanged(int)),this,SLOT(PanFitEnable()));
    connect(ui->chBeParFit_3,SIGNAL(stateChanged(int)),this,SLOT(PanFitEnable()));
    connect(ui->chBeParFit_4,SIGNAL(stateChanged(int)),this,SLOT(PanFitEnable()));
    connect(ui->chBeParFit_5,SIGNAL(stateChanged(int)),this,SLOT(PanFitEnable()));
    connect(ui->chBeParFit_6,SIGNAL(stateChanged(int)),this,SLOT(PanFitEnable()));
    connect(ui->chBeParFit_7,SIGNAL(stateChanged(int)),this,SLOT(PanFitEnable()));
    connect(ui->chBeParFit_8,SIGNAL(stateChanged(int)),this,SLOT(PanFitEnable()));
    connect(ui->chBeParFit_9,SIGNAL(stateChanged(int)),this,SLOT(PanFitEnable()));
    connect(ui->chBeParFit_10,SIGNAL(stateChanged(int)),this,SLOT(PanFitEnable()));
    connect(ui->chBeParFit_11,SIGNAL(stateChanged(int)),this,SLOT(PanFitEnable()));
    connect(ui->chBeParFit_12,SIGNAL(stateChanged(int)),this,SLOT(PanFitEnable()));
    connect(ui->chBeParFit_13,SIGNAL(stateChanged(int)),this,SLOT(PanFitEnable()));
    connect(ui->chBeParFit_14,SIGNAL(stateChanged(int)),this,SLOT(PanFitEnable()));
    connect(ui->chBeParFit_15,SIGNAL(stateChanged(int)),this,SLOT(PanFitEnable()));
    connect(ui->chBeParFit_16,SIGNAL(stateChanged(int)),this,SLOT(PanFitEnable()));
    connect(ui->chBeParFit_17,SIGNAL(stateChanged(int)),this,SLOT(PanFitEnable()));
    connect(ui->cBParFit_1,SIGNAL(currentIndexChanged(int)),this,SLOT(PanFitChoice()));
    connect(ui->cBParFit_2,SIGNAL(currentIndexChanged(int)),this,SLOT(PanFitChoice()));
    connect(ui->cBParFit_3,SIGNAL(currentIndexChanged(int)),this,SLOT(PanFitChoice()));
    connect(ui->cBParFit_4,SIGNAL(currentIndexChanged(int)),this,SLOT(PanFitChoice()));
    connect(ui->cBParFit_5,SIGNAL(currentIndexChanged(int)),this,SLOT(PanFitChoice()));
    connect(ui->cBParFit_6,SIGNAL(currentIndexChanged(int)),this,SLOT(PanFitChoice()));
    connect(ui->cBParFit_7,SIGNAL(currentIndexChanged(int)),this,SLOT(PanFitChoice()));
    connect(ui->cBParFit_8,SIGNAL(currentIndexChanged(int)),this,SLOT(PanFitChoice()));
    connect(ui->cBParFit_9,SIGNAL(currentIndexChanged(int)),this,SLOT(PanFitChoice()));
    connect(ui->cBParFit_10,SIGNAL(currentIndexChanged(int)),this,SLOT(PanFitChoice()));
    connect(ui->cBParFit_11,SIGNAL(currentIndexChanged(int)),this,SLOT(PanFitChoice()));
    connect(ui->cBParFit_12,SIGNAL(currentIndexChanged(int)),this,SLOT(PanFitChoice()));
    connect(ui->cBParFit_13,SIGNAL(currentIndexChanged(int)),this,SLOT(PanFitChoice()));
    connect(ui->cBParFit_14,SIGNAL(currentIndexChanged(int)),this,SLOT(PanFitChoice()));
    connect(ui->cBParFit_15,SIGNAL(currentIndexChanged(int)),this,SLOT(PanFitChoice()));
    connect(ui->cBParFit_16,SIGNAL(currentIndexChanged(int)),this,SLOT(PanFitChoice()));
    connect(ui->cBParFit_17,SIGNAL(currentIndexChanged(int)),this,SLOT(PanFitChoice()));
    connect(ui->pushButton_PlotExpMis,SIGNAL(clicked()), this, SLOT(PlotMENK()));
    connect(ui->pushButton_PlotExpMis2,SIGNAL(clicked()), this, SLOT(PlotMENK()));
    connect(ui->pushButton_Simulate,SIGNAL(clicked()), this, SLOT(Simula()));
    connect(ui->pushButton_modelSim,SIGNAL(clicked()), this, SLOT(Simula()));
    connect(ui->pushButton_DataFitSim,SIGNAL(clicked()), this, SLOT(Simula()));
    connect(ui->pushButton_PlotAve,SIGNAL(clicked()), this, SLOT(PlotAve()));
    connect(ui->pushButton_texture,SIGNAL(clicked()), this, SLOT(PlotTexturized()));
    connect(ui->pushButton_SearchNK,SIGNAL(clicked()), this, SLOT(searchNK()));
    connect(ui->pushButton_RefreshGraph,SIGNAL(clicked()), this, SLOT(RefTrackG()));
    connect(ui->pushButton_AutosFromWLmin,SIGNAL(clicked()), this, SLOT(NumericalSearch()));
    connect(ui->pushButton_PlotCurrentFit,SIGNAL(clicked()), this, SLOT(IbridPlotFit()));
    connect(ui->pushButton_FitN,SIGNAL(clicked()), this, SLOT(FitN()));
    connect(ui->pushButton_FitNK,SIGNAL(clicked()), this, SLOT(FitNK()));
    connect(ui->pushButton_FitE1E2,SIGNAL(clicked()), this, SLOT(FitE1E2()));
    connect(ui->pushButton_FitSelExpMeas,SIGNAL(clicked()), this, SLOT(FitSelExpMeas()));
    connect(ui->pushButton_IbridOne,SIGNAL(clicked()), this, SLOT(IbridOne()));
    connect(ui->pushButton_IbridOneErrStore,SIGNAL(clicked()), this, SLOT(IbridOneStore()));
    connect(ui->pushButton_stop,SIGNAL(clicked()), this, SLOT(setStop()));
    connect(ui->pushButton_store,SIGNAL(clicked()), this,SLOT(StoreFitSet()));
    connect(ui->pushButton_BackToBestSituation,SIGNAL(clicked()), this,SLOT(GoBest()));
    connect(ui->pushButton_previous,SIGNAL(clicked()), this,SLOT(GoPrevious()));
    connect(ui->pushButton_next,SIGNAL(clicked()), this,SLOT(GoNext()));
    connect(ui->pushButton_saveNKens,SIGNAL(clicked()), this,SLOT(SaveFnk()));
    connect(ui->pushButton_saveNKbestFit,SIGNAL(clicked()), this,SLOT(SaveFnk()));
    connect(ui->pBmDw1, &QPushButton::clicked, this, [this]() {mDwUp(1, 1);});
    connect(ui->pBmDw2, &QPushButton::clicked, this, [this]() {mDwUp(2, 1);});
    connect(ui->pBmDw3, &QPushButton::clicked, this, [this]() {mDwUp(3, 1);});
    connect(ui->pBmDw4, &QPushButton::clicked, this, [this]() {mDwUp(4, 1);});
    connect(ui->pBmDw5, &QPushButton::clicked, this, [this]() {mDwUp(5, 1);});
    connect(ui->pBmDw6, &QPushButton::clicked, this, [this]() {mDwUp(6, 1);});
    connect(ui->pBmDw7, &QPushButton::clicked, this, [this]() {mDwUp(7, 1);});
    connect(ui->pBmDw8, &QPushButton::clicked, this, [this]() {mDwUp(8, 1);});
    connect(ui->pBmDw9, &QPushButton::clicked, this, [this]() {mDwUp(9, 1);});
    connect(ui->pBmUp1, &QPushButton::clicked, this, [this]() {mDwUp(1,-1);});
    connect(ui->pBmUp2, &QPushButton::clicked, this, [this]() {mDwUp(2,-1);});
    connect(ui->pBmUp3, &QPushButton::clicked, this, [this]() {mDwUp(3,-1);});
    connect(ui->pBmUp4, &QPushButton::clicked, this, [this]() {mDwUp(4,-1);});
    connect(ui->pBmUp5, &QPushButton::clicked, this, [this]() {mDwUp(5,-1);});
    connect(ui->pBmUp6, &QPushButton::clicked, this, [this]() {mDwUp(6,-1);});
    connect(ui->pBmUp7, &QPushButton::clicked, this, [this]() {mDwUp(7,-1);});
    connect(ui->pBmUp8, &QPushButton::clicked, this, [this]() {mDwUp(8,-1);});
    connect(ui->pBmUp9, &QPushButton::clicked, this, [this]() {mDwUp(9,-1);});
    connect(ui->pushButton_saveSim, SIGNAL( clicked() ), this, SLOT(saveSim()));
    connect(ui->pushButton_saveNKsim,SIGNAL( clicked() ), this, SLOT(saveNKsim()));
    connect(ui->pushButton_PlotAbsEL,SIGNAL( clicked() ),this, SLOT(PlotAbsEL()));
    connect(ui->comboBox_searchNK,SIGNAL(currentIndexChanged(int)),this,SLOT(manageLEwl()));
    connect(ui->pushButton_selectN,SIGNAL( clicked() ), this, SLOT(selectNsol()));
    connect(ui->cBpm_101_1,SIGNAL(currentIndexChanged(int)),this,SLOT(setKindOsc1()));
    connect(ui->cBpm_106_1,SIGNAL(currentIndexChanged(int)),this,SLOT(setKindOsc2()));
    connect(ui->cBpm_111_1,SIGNAL(currentIndexChanged(int)),this,SLOT(setKindOsc3()));
    connect(ui->cBpm_116_1,SIGNAL(currentIndexChanged(int)),this,SLOT(setKindOsc4()));
    connect(ui->cBpm_121_1,SIGNAL(currentIndexChanged(int)),this,SLOT(setKindOsc5()));
    connect(ui->cBpm_126_1,SIGNAL(currentIndexChanged(int)),this,SLOT(setKindOsc6()));
    connect(ui->cBpm_131_1,SIGNAL(currentIndexChanged(int)),this,SLOT(setKindOsc7()));
    connect(ui->cBpm_136_1,SIGNAL(currentIndexChanged(int)),this,SLOT(setKindOsc8()));
    connect(ui->cBpm_141_1,SIGNAL(currentIndexChanged(int)),this,SLOT(setKindOsc9()));
    connect(ui->cBpm_146_1,SIGNAL(currentIndexChanged(int)),this,SLOT(setKindOsc10()));
    connect(ui->cBpm_151_1,SIGNAL(currentIndexChanged(int)),this,SLOT(setKindOsc11()));
    connect(ui->cBpm_156_1,SIGNAL(currentIndexChanged(int)),this,SLOT(setKindOsc12()));
    connect(ui->cBpm_161_1,SIGNAL(currentIndexChanged(int)),this,SLOT(setKindOsc13()));
    connect(ui->cBpm_166_1,SIGNAL(currentIndexChanged(int)),this,SLOT(setKindOsc14()));
    connect(ui->cBpm_171_1,SIGNAL(currentIndexChanged(int)),this,SLOT(setKindOsc15()));
    connect(ui->cBpm_176_1,SIGNAL(currentIndexChanged(int)),this,SLOT(setKindOsc16()));
    connect(ui->cBpm_181_1,SIGNAL(currentIndexChanged(int)),this,SLOT(setKindOsc17()));
    connect(ui->cBpm_186_1,SIGNAL(currentIndexChanged(int)),this,SLOT(setKindOsc18()));
    connect(ui->cBpm_191_1,SIGNAL(currentIndexChanged(int)),this,SLOT(setKindOsc19()));
    connect(ui->cBpm_196_1,SIGNAL(currentIndexChanged(int)),this,SLOT(setKindOsc20()));
    connect(ui->comboBox_DeltaPsiScale,SIGNAL(currentIndexChanged(int)),this,SLOT(setRangeEli()));

    connect(ui->pushButton_setFont, SIGNAL( clicked() ), this, SLOT(setFontDia()));

    occupyPF=1;//disable severals call during initial setting
    occupyPEj=1;

    // parameter initialization
#ifdef __unix__
#define IS_POSIX 1
#else
#define IS_POSIX 0
#endif

    QDir dir;  //current directory
    dir.cdUp();//cd ..
    dir.cdUp();//cd ..
    if (IS_POSIX == 1) {
        //Linux path initialization
        //nothing to do
    }
    else {
        //windows path inizialization
        dir.cdUp();//cd ..
    }
    pathroot=dir.absolutePath()+"/";
    fRefMir=pathroot+"qtSource/ksemawc/referenceMirrors.txt";
    fStdSpect=pathroot+"qtSource/ksemawc/standardSpectra.txt";
    fileStore=pathroot+"temp/defau.1.Spj";
    fileStore0=pathroot+"temp/defau.0.Spj";
    fileStoreSpjName=pathroot+"temp/SpjName.txt";
    fNKsim=pathroot+"expo/NKsim.dat";
    fMisSim=pathroot+"expo/MisSim.dat";
    filechi2=pathroot+"temp/ksemawc.log";

    cout<<"pathroot= "<<pathroot.toStdString()<<"\n";

    // stylesheet
    cout << "styleSheet= " << QApplication::style()->metaObject()->className() <<"\n";
    cout << "Current StyleSheet:" << qApp->styleSheet().toStdString().c_str()<<"\n";
    //set specific background color to groupBox of DataAnalysis Tab
    ui->groupBox->setStyleSheet("background-color: #eff0f1;");
    ui->groupBox_2->setStyleSheet("background-color: #ffdede;");
    ui->groupBox_3->setStyleSheet("background-color: #ffffde;");
    //ui->groupBox_4->setStyleSheet("background-color: #d3e4ff;");
    ui->groupBox_5->setStyleSheet("background-color: #deffde;");
    ui->groupBox_6->setStyleSheet("background-color: #ddfdff;");
    ui->groupBox_7->setStyleSheet("background-color: #fce2ce;");
    //ui->tabModel->setStyleSheet("background-color: #deffde;");

    lastAction="";
    lastIndex=0;
    lastTab=0;
    ifn=0;
    nlayer=0;
    npp=0;
    for(int i=1;i<=5;i++){
        for(int j=1;j<=60;j++){
            par[j][i]=0.;
        }
    }
    par[1][3]= 1.;    //capture effiency of reflection from 1st face
    par[2][3]= 1.;    //capture effiency of reflection from 2nd face
    par[3][3]=.0022;  //drift baseline of spectrophotometer
    par[4][1]=200.;  //wavelength MIN (nm)
    par[4][2]=3000.; //wavelength MAX (nm)
    par[4][3]=.005;   //relative error of reference mirror reflectance
    par[5][3]=.0005;  //reading error
    par[6][3]=1.0;    // ... depending on wavelengh
    par[7][3]=1.2505E-7; //cte used for computing the substrate contribution to PDS signal
    par[11][1]=1;     //attenuation coeff. of k by Fit#N
    par[11][2]=0;     //is set to 1 when ludcmp found singular matrix
    par[21][3]=.001;  //n step
    par[22][3]=1.e-4; //k step
    par[28][1]=3.;    //N for computing <f(x)>
    par[28][2]=201;   // number of interpolated points in Ev or lambda
    NeV=par[28][2];
    par[29][1]=23.;   //N. sublayer used to model each inhomogeneous film
    // set the pointer to PM: PPM(17)
    par[34][5]=7.;    //7 parameters
    par[36][5]=102;   //C1
    par[37][5]=103;   //E1
    par[38][5]=104;   //D1
    par[39][5]=107;   //C2
    par[40][5]=108;   //E2
    par[41][5]=109;   //D2
    par[42][5]=112;   //C3
    // initialization of CNK which pilot VNK setting
    for(int i=1;i<=16;i++){
        CNK[i][1]=.0;//option "constant nk"
        CNK[i][2]=1.;//n1 value
        CNK[i][3]=.0;//k1 value
        CNK[i][4]=-1.;//No EMA
        CNK[i][5]=1.;//n2 value
        CNK[i][6]=.0;//k2 value
    }
    // initialization of model parameters pm
    for(int i=1;i<=200;i++){
        for(int j=1;j<=5;j++){
            pm[i][j]=0.;
        }
    }
    pm[86][1]=-1.;
    pm[87][1]=-1.;
    pm[88][1]=-1.;
    pm[89][1]=-1.;
    pm[100][1]=1;    //option FT#1
    pm[101][1]=2.;   //quantistic homogeneous oscillator in UV
    pm[102][1]=0.;   //C1    "          "
    pm[103][1]=5.0;  //E1    "          "
    pm[104][1]=0.04; //D1    "          "
    pm[105][1]=0.;   //K1    "          "
    pm[106][1]=2.;   //quantistic homogeneous oscillator in IR
    pm[107][1]=0.;   //C2
    pm[108][1]=0.3;  //E2    "          "
    pm[109][1]=0.005;//D2    "          "
    pm[110][1]=0.;   //K2    "          "
    pm[111][1]=4.;   //FLAT term
    pm[112][1]=1.5;  //C3    "          "
    //graphic parameters
    for(int i=1;i<=6;i++){//spectrophotometric spectra
        rxy[i][1]=.0;  //Ymin
        rxy[i][2]=100.;//Ymax
        rxy[i][3]=.0;  //Ymin-data
        rxy[i][4]=.0;  //Ymax-data
    }
    for(int i=7;i<=13;i+=2 ){//DELTA ellipsometric spectra
        rxy[i][1]=-180.;  //min
        rxy[i][2]=180.;//max
        rxy[i][3]=.0;  //min-data
        rxy[i][4]=.0;  //max-data
    }
    for(int i=8;i<=14;i+=2 ){//PSI ellipsometric spectra
        rxy[i][1]=0;//min
        rxy[i][2]=90.;//max
        rxy[i][3]=.0;  //min-data
        rxy[i][4]=.0;  //max-data
    }
    for(int i=16;i<=17;i++){//n k
        rxy[i][1]=.0;
        rxy[i][2]=2.0;
        rxy[i][3]=.0;
        rxy[i][4]=.0;
    }
    rxy[15][1]=-10.0;
    rxy[15][2]=100.0;
    rxy[15][3]=.0;
    rxy[15][4]=.0;
    rxy[18][1]=-10.0;
    rxy[18][2]=100.0;
    rxy[18][3]=.0;
    rxy[18][4]=.0;
    rxy[19][1]=-1.0;//cos(Delta)
    rxy[19][2]=1.0;
    rxy[19][3]=.0;
    rxy[19][4]=.0;
    rxy[20][1]=200.; //wlmin (nm)
    rxy[20][2]=3000.;//wlmax
    rxy[20][3]=200.; //wlmin
    rxy[20][4]=3000.;//wlmax
    rxy[21][1]=0.;  //theta_inc min
    rxy[21][2]=90.; //theta_inc max
    rxy[21][3]=.0;  //theta_inc min
    rxy[21][4]=90.; //theta_inc max
    rxy[22][1]=0;//tan(Psi)
    rxy[22][2]=10.0;
    rxy[22][3]=.0;
    rxy[22][4]=.0;
    for(int i=23;i<=30;i++){
        rxy[i][1]=.0;
        rxy[i][2]=.0;
        rxy[i][3]=.0;
        rxy[i][4]=.0;
    }
    rxy[24][1]=500.;//graph-window width
    rxy[24][2]=300.;//graph-window high
    rxy[24][3]=3.;//line width
    rxy[25][3]=2.;//plot in eV
    rxy[25][4]=2.;//spline in eV
    for(int i=26;i<=27;i++){
        rxy[i][1]=.0;
        rxy[i][2]=10.;
        rxy[i][3]=0.;
        rxy[i][4]=10.;
    }
    ui->sB_iwl2print->setValue(iwl2print);

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
    idToLineEdit["lineEdit_E1"]=ui->lineEdit_E1;
    idToLineEdit["lineEdit_E2"]=ui->lineEdit_E2;
    idToLineEdit["lineEdit_E3"]=ui->lineEdit_E3;
    idToLineEdit["lineEdit_E4"]=ui->lineEdit_E4;
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
    // idToLineEdit["DP_RXY_22_1"]=ui->DP_RXY_22_1;
    // idToLineEdit["DP_RXY_22_2"]=ui->DP_RXY_22_2;
    // idToLineEdit["DP_RXY_22_3"]=ui->DP_RXY_22_3;
    // idToLineEdit["DP_RXY_22_4"]=ui->DP_RXY_22_4;
    idToLineEdit["DP_RXY_23_1"]=ui->DP_RXY_23_1;
    idToLineEdit["DP_RXY_23_2"]=ui->DP_RXY_23_2;
    idToLineEdit["DP_RXY_23_3"]=ui->DP_RXY_23_3;
    idToLineEdit["DP_RXY_23_4"]=ui->DP_RXY_23_4;
    idToLineEdit["DP_RXY_24_1"]=ui->DP_RXY_24_1;
    idToLineEdit["DP_RXY_24_2"]=ui->DP_RXY_24_2;
    idToLineEdit["DP_RXY_24_3"]=ui->DP_RXY_24_3;
    //  idToLineEdit["DP_RXY_24_4"]=ui->DP_RXY_24_4;
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
    idToLineEdit["LEcnk1_5"]=ui->LEcnk1_5;
    idToLineEdit["LEcnk1_6"]=ui->LEcnk1_6;
    idToLineEdit["LEcnk2_2"]=ui->LEcnk2_2;
    idToLineEdit["LEcnk2_3"]=ui->LEcnk2_3;
    idToLineEdit["LEcnk2_5"]=ui->LEcnk2_5;
    idToLineEdit["LEcnk2_6"]=ui->LEcnk2_6;
    idToLineEdit["LEcnk3_2"]=ui->LEcnk3_2;
    idToLineEdit["LEcnk3_3"]=ui->LEcnk3_3;
    idToLineEdit["LEcnk3_5"]=ui->LEcnk3_5;
    idToLineEdit["LEcnk3_6"]=ui->LEcnk3_6;
    idToLineEdit["LEcnk4_2"]=ui->LEcnk4_2;
    idToLineEdit["LEcnk4_3"]=ui->LEcnk4_3;
    idToLineEdit["LEcnk4_5"]=ui->LEcnk4_5;
    idToLineEdit["LEcnk4_6"]=ui->LEcnk4_6;
    idToLineEdit["LEcnk5_2"]=ui->LEcnk5_2;
    idToLineEdit["LEcnk5_3"]=ui->LEcnk5_3;
    idToLineEdit["LEcnk5_5"]=ui->LEcnk5_5;
    idToLineEdit["LEcnk5_6"]=ui->LEcnk5_6;
    idToLineEdit["LEcnk6_2"]=ui->LEcnk6_2;
    idToLineEdit["LEcnk6_3"]=ui->LEcnk6_3;
    idToLineEdit["LEcnk6_5"]=ui->LEcnk6_5;
    idToLineEdit["LEcnk6_6"]=ui->LEcnk6_6;
    idToLineEdit["LEcnk7_2"]=ui->LEcnk7_2;
    idToLineEdit["LEcnk7_3"]=ui->LEcnk7_3;
    idToLineEdit["LEcnk7_5"]=ui->LEcnk7_5;
    idToLineEdit["LEcnk7_6"]=ui->LEcnk7_6;
    idToLineEdit["LEcnk8_2"]=ui->LEcnk8_2;
    idToLineEdit["LEcnk8_3"]=ui->LEcnk8_3;
    idToLineEdit["LEcnk8_5"]=ui->LEcnk8_5;
    idToLineEdit["LEcnk8_6"]=ui->LEcnk8_6;
    idToLineEdit["LEcnk9_2"]=ui->LEcnk9_2;
    idToLineEdit["LEcnk9_3"]=ui->LEcnk9_3;
    idToLineEdit["LEcnk9_5"]=ui->LEcnk9_5;
    idToLineEdit["LEcnk9_6"]=ui->LEcnk9_6;
    idToLineEdit["LEcnk10_2"]=ui->LEcnk10_2;
    idToLineEdit["LEcnk10_3"]=ui->LEcnk10_3;
    idToLineEdit["LEcnk10_5"]=ui->LEcnk10_5;
    idToLineEdit["LEcnk10_6"]=ui->LEcnk10_6;
    idToLineEdit["LEcnk11_2"]=ui->LEcnk11_2;
    idToLineEdit["LEcnk11_3"]=ui->LEcnk11_3;
    idToLineEdit["LEcnk11_5"]=ui->LEcnk11_5;
    idToLineEdit["LEcnk11_6"]=ui->LEcnk11_6;
    idToLineEdit["LEcnk12_2"]=ui->LEcnk12_2;
    idToLineEdit["LEcnk12_3"]=ui->LEcnk12_3;
    idToLineEdit["LEcnk12_5"]=ui->LEcnk12_5;
    idToLineEdit["LEcnk12_6"]=ui->LEcnk12_6;
    idToLineEdit["LEcnk13_2"]=ui->LEcnk13_2;
    idToLineEdit["LEcnk13_3"]=ui->LEcnk13_3;
    idToLineEdit["LEcnk13_5"]=ui->LEcnk13_5;
    idToLineEdit["LEcnk13_6"]=ui->LEcnk13_6;
    idToLineEdit["LEcnk14_2"]=ui->LEcnk14_2;
    idToLineEdit["LEcnk14_3"]=ui->LEcnk14_3;
    idToLineEdit["LEcnk14_5"]=ui->LEcnk14_5;
    idToLineEdit["LEcnk14_6"]=ui->LEcnk14_6;
    idToLineEdit["LEcnk15_2"]=ui->LEcnk15_2;
    idToLineEdit["LEcnk15_3"]=ui->LEcnk15_3;
    idToLineEdit["LEcnk15_5"]=ui->LEcnk15_5;
    idToLineEdit["LEcnk15_6"]=ui->LEcnk15_6;
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
    idToDoubleSpinBox["dSB_PM_71_1"]=ui->dSB_PM_71_1;
    idToDoubleSpinBox["dSB_PM_72_1"]=ui->dSB_PM_72_1;
    idToDoubleSpinBox["dSB_PM_73_1"]=ui->dSB_PM_73_1;
    idToDoubleSpinBox["dSB_PM_74_1"]=ui->dSB_PM_74_1;
    idToDoubleSpinBox["dSB_PM_75_1"]=ui->dSB_PM_75_1;
    idToDoubleSpinBox["dSB_PM_76_1"]=ui->dSB_PM_76_1;
    idToDoubleSpinBox["dSB_PM_77_1"]=ui->dSB_PM_77_1;
    idToDoubleSpinBox["dSB_PM_78_1"]=ui->dSB_PM_78_1;
    idToDoubleSpinBox["dSB_PM_79_1"]=ui->dSB_PM_79_1;
    idToDoubleSpinBox["dSB_PM_80_1"]=ui->dSB_PM_80_1;
    idToDoubleSpinBox["dSB_PM_81_1"]=ui->dSB_PM_81_1;
    idToDoubleSpinBox["dSB_PM_82_1"]=ui->dSB_PM_82_1;
    idToDoubleSpinBox["dSB_PM_83_1"]=ui->dSB_PM_83_1;
    idToDoubleSpinBox["dSB_PM_84_1"]=ui->dSB_PM_84_1;
    idToDoubleSpinBox["dSB_PM_85_1"]=ui->dSB_PM_85_1;
    idToDoubleSpinBox["dSB_PM_86_1"]=ui->dSB_PM_86_1;
    idToDoubleSpinBox["dSB_PM_87_1"]=ui->dSB_PM_87_1;
    idToDoubleSpinBox["dSB_PM_88_1"]=ui->dSB_PM_88_1;
    idToDoubleSpinBox["dSB_PM_89_1"]=ui->dSB_PM_89_1;
    idToDoubleSpinBox["dSB_PM_91_1"]=ui->dSB_PM_91_1;
    idToDoubleSpinBox["dSB_PM_92_1"]=ui->dSB_PM_92_1;
    idToDoubleSpinBox["dSB_PM_93_1"]=ui->dSB_PM_93_1;
    idToDoubleSpinBox["dSB_PM_94_1"]=ui->dSB_PM_94_1;
    idToDoubleSpinBox["dSB_PM_95_1"]=ui->dSB_PM_95_1;
    idToDoubleSpinBox["dSB_PM_96_1"]=ui->dSB_PM_96_1;
    idToDoubleSpinBox["dSB_PM_97_1"]=ui->dSB_PM_97_1;
    idToDoubleSpinBox["dSB_PM_98_1"]=ui->dSB_PM_98_1;
    idToDoubleSpinBox["dSB_PM_99_1"]=ui->dSB_PM_99_1;

    for(int i=1;i<=25;i++){//setting values in GUI
        for(int j=1;j<=4;j++){
            if(idToLineEdit.contains("DP_RXY_"+QString::number(i)+"_"+QString::number(j)))
                idToLineEdit["DP_RXY_"+QString::number(i)+"_"+QString::number(j)]
                        -> setText(QString::number(rxy[i][j]));
            if(i==24&&j==3)
                j++;
        }
    }
    for(int i=1;i<=16;i++)
        IXW[i]=-1;//no graph window

    //oscillator list
    for(int k=1;k<=20;k++){
        idToComboBox["cBpm_"+QString::number(100+1+(k-1)*5)+"_1"]-> clear();
        idToComboBox["cBpm_"+QString::number(100+1+(k-1)*5)+"_1"]->addItem("Lorentz");      //#1
        idToComboBox["cBpm_"+QString::number(100+1+(k-1)*5)+"_1"]->addItem("Quant-homo");   //#2
        idToComboBox["cBpm_"+QString::number(100+1+(k-1)*5)+"_1"]->addItem("Quant-inhomo"); //#3
        idToComboBox["cBpm_"+QString::number(100+1+(k-1)*5)+"_1"]->addItem("Flat");         //#4
        idToComboBox["cBpm_"+QString::number(100+1+(k-1)*5)+"_1"]->addItem("Drude");        //#5
        idToComboBox["cBpm_"+QString::number(100+1+(k-1)*5)+"_1"]->addItem("Ind-Gap-Cody"); //#6
        idToComboBox["cBpm_"+QString::number(100+1+(k-1)*5)+"_1"]->addItem("Ind-Gap-Tauc"); //#7
        idToComboBox["cBpm_"+QString::number(100+1+(k-1)*5)+"_1"]->addItem("Dir-Gap-Cody"); //#8
        idToComboBox["cBpm_"+QString::number(100+1+(k-1)*5)+"_1"]->addItem("Dir-Gap-Tauc"); //#9
        idToComboBox["cBpm_"+QString::number(100+1+(k-1)*5)+"_1"]->addItem("Ind-Gap-Cody-Urbach");  //#10
        idToComboBox["cBpm_"+QString::number(100+1+(k-1)*5)+"_1"]->addItem("Ind-Gap-Tauc-Urbach");  //#11
        idToComboBox["cBpm_"+QString::number(100+1+(k-1)*5)+"_1"]->addItem("Dir-Gap-Cody-Urbach");  //#12
        idToComboBox["cBpm_"+QString::number(100+1+(k-1)*5)+"_1"]->addItem("Dir-Gap-Tauc-Urbach");  //#13
        idToComboBox["cBpm_"+QString::number(100+1+(k-1)*5)+"_1"]->addItem("Dir-Gap-Tauc-Exciton"); //#14
        idToComboBox["cBpm_"+QString::number(100+1+(k-1)*5)+"_1"]->addItem("Dir-Gap-Cody-M1M2");    //#15
        idToComboBox["cBpm_"+QString::number(100+1+(k-1)*5)+"_1"]->addItem("Dir-Gap-Tauc-M1M2");    //#16
        idToComboBox["cBpm_"+QString::number(100+1+(k-1)*5)+"_1"]->addItem("Drude-ionized");        //#17
        idToComboBox["cBpm_"+QString::number(100+1+(k-1)*5)+"_1"]->addItem("Lorentz-Dirac");        //#18
        //idToComboBox["cBpm_"+QString::number(100+1+(k-1)*5)+"_1"]->addItem("Volum-Scat");           //#19
        idToComboBox["cBpm_"+QString::number(100+1+(k-1)*5)+"_1"]->setCurrentIndex(-1);
    }


    QString line;
    //Refrence mirror list
    QFile file0(fRefMir);
    if (file0.open (QIODevice::ReadOnly | QIODevice::Text)){
        QTextStream stream ( &file0 );
        int iB103=0;
        ui -> comboB_PAR_10_3-> clear();
        ui -> comboB_PAR_10_3->addItem("none");
        do{
            iB103++;
            line = stream.readLine();
            ui -> comboB_PAR_10_3->addItem(line);
            line = stream.readLine();
            printf("-> RefMir #%d -> %s\n",iB103,line.toStdString().c_str());
        }while(!stream.atEnd());
        file0.close();
    }

    //standard spectra for mean value
    QFile file1(fStdSpect);
    if (file1.open (QIODevice::ReadOnly | QIODevice::Text)){
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
    }

    ParFitLab[0]="_dNonUni";
    ParFitLab[1]="_d";
    ParFitLab[2]="_Dn/<n>";
    ParFitLab[3]="_(nav/<n>-1)";
    ParFitLab[4]="_Dk/<k>";
    ParFitLab[5]="_(kav/<k>-1)";
    ParFitLab[6]="_sigma_rough";
    ParFitLab[7]="_slope_Dn";
    ParFitLab[8]="_C";
    ParFitLab[9]="_E";
    ParFitLab[10]="_D";
    ParFitLab[11]="_W";
    ParFitLab[12]="fEMA_";
    ParFitLab[13]="Theta_";

    //set the latest Size & font if settings exist
    QSettings settings("ksemawc","ksemawc");
    QString usnam=settings.fileName();
    cout<<"setting-file exists! "<<usnam.toStdString()<<"\n";
    if (settings.status() == QSettings::NoError && settings.status() != QSettings::FormatError && settings.status() !=QSettings::AccessError){
        restoreGeometry(settings.value("geometry").toByteArray());
        QVariant fontFamily=settings.value("fontFamily",QString());
        int fontSize = settings.value("pointSize", int()).toInt();
        bool fontIsBold = settings.value("bold", true).toBool();
        bool fontIsItalic = settings.value("italic", false).toBool();
        //cout << "font= "  << (fontFamily.toString()).toStdString()<<"\n";
        //cout << "size = " << fontSize<<"\n";
        //cout << "bold = " << fontIsBold<<"\n";
        //cout << "italic= "<<fontIsItalic<<"\n";
        QFont font;
        font.setFamily(fontFamily.toString());
        font.setPointSize(fontSize);
        font.setBold(fontIsBold);
        font.setItalic(fontIsItalic);
        const QWidgetList allWidgets = QApplication::allWidgets();
        for (QWidget *widget : allWidgets){
            widget->setFont(font);
            widget->update();
        }
        QCoreApplication::processEvents();
    }

    occupyPF=0;
    occupyPEj=0;

    //make your choice
    QMessageBox msgBox;
    msgBox.setWindowTitle("Starting ksemawc project");
    msgBox.setText("Please make a choice:");
    QPushButton* continueButton = msgBox.addButton("Continue the previous session", QMessageBox::AcceptRole);
    QPushButton* newButton = msgBox.addButton("Start a new project", QMessageBox::RejectRole);
    QPushButton* openButton = msgBox.addButton("Open an existing project", QMessageBox::DestructiveRole);
    msgBox.setDefaultButton(newButton);
    int ret=msgBox.exec();
    printf("ret=%d\n",ret);
    QString SpjName;
    QFile file(fileStoreSpjName);
    switch (ret) {

    case QMessageBox::AcceptRole:
        printf("Continue the previous session\n");
        //Load Project name
        if(file.open(QIODevice::ReadOnly | QIODevice::Text)){
            QTextStream stream (&file);
            SpjName=stream.readLine();
            file.close();
        }
        if(SpjName.contains("mate/aa999.9")||SpjName.isEmpty()){
            caller="initialization";
            ReadSetting(fileStore0);
            ui->lineEdit_infoP-> setText("The last session has been loaded. Now please set the project name!");
        }
        else{
            caller="initialization";
            ReadSetting(pathroot+SpjName);
            ui->lineEdit_P -> setText(SpjName);
        }
        break;

    case QMessageBox::RejectRole:
        printf("Start a new project\n");
        break;

    case QMessageBox::DestructiveRole:
        printf("Open an existing project\n");
        LoadProject();
        break;
    }

    SaveSetting(-1);//save
}


void ksemawc::closeEvent ( QCloseEvent * event )
{
    //save Project name
    QString fileSpjName=ui->lineEdit_P -> text();
    QFile file(fileStoreSpjName);
    QTextStream stream (&file);
    if(file.open(QIODevice::WriteOnly | QIODevice::Text)){
        stream << fileSpjName << "\n";
        file.close();
    }

    //save project in default
    QString command;
    int ierr;
    if (IS_POSIX == 1){
        command="cp "+fileStore+" "+fileStore0;
        ierr=system((command.toStdString()).c_str());}
    else{
        command="copy "+fileStore+" "+fileStore0;
        QStringList List;
        List =command.split("/");
        int nV=List.count();
        //printf("\nnV=%d\n",nV);
        QString commandw;
        QString pezzo;
        for(int iv=0;iv<nV;iv++){
            pezzo=List.at(iv).toLocal8Bit().constData();
            //printf("pezzo[%d]= %s\t",iv,pezzo.toStdString().c_str());
            if(iv==0){
                commandw=pezzo;
            }
            else{
                commandw=commandw+"\\"+pezzo;
            }
        }
        //printf("comando=%s\n",(commandw.toStdString()).c_str());
        ierr=system((commandw.toStdString()).c_str());
    }
    if(ierr != 0){
        msgErrLoad("SaveProject",fnproject);
        printf("Error copying project!!!\n");
    }
    else
        printf("Project saved as %s\n",(fnproject.toStdString()).c_str());
    //save widget size
    QSettings settings("ksemawc","ksemawc");
    settings.setValue("geometry", saveGeometry());
    qApp->quit();
}


void ksemawc::setFontDia(){
    bool ok;
    QFont font=QFontDialog::getFont(&ok,this);//,QFont("Noto Sans",8)
    if (ok) {
        // the user clicked OK and font is set to the font the user selected
        const QWidgetList allWidgets = QApplication::allWidgets();
        for (QWidget *widget : allWidgets){
            widget->setFont(font);
            widget->update();
        }
        QCoreApplication::processEvents();
    } else {
        // the user canceled the dialog; font is set to the initial
        // value, in this case Helvetica [Cronyx], 10
    }
    //cout << "fontFamily= "  << font.family().toStdString()<<"\n";
    //cout << "size = " << font.pointSize()<<"\n";
    //cout << "bold = " << font.bold()<<"\n";
    //cout << "italic = " << font.italic()<<"\n";
    //save the choice
    QSettings settings("ksemawc","ksemawc");
    settings.setValue("fontFamily", font.family());
    settings.setValue("pointSize", font.pointSize());
    settings.setValue("bold", font.bold());
    settings.setValue("italic", font.italic());
}


void ksemawc::LoadProject(){
    printf("-> LoadProject\n");
    fnproject = QFileDialog::getOpenFileName(
                this,
                "Choose a SEMAW project", //window title
                pathroot,                 //initial directory
                "Semaw Project (*.Spj)"); //file extension
    // printf("fnprojet=%s\n",fnproject.toStdString().c_str());
    if(fnproject.isEmpty())
        return;

    occupyPF=10;
    ui->tabWidget -> setCurrentIndex(0);
    //thetaEli reset
    pm[86][1]=0.;
    pm[87][1]=0.;
    pm[88][1]=0.;
    pm[89][1]=0.;
    // reset of pm - model parameters
    for(int i=1;i<=200;i++){
        for(int j=1;j<=5;j++){
            pm[i][j]=0.;
        }
    }
    //reset oscillator list
    ui->sB_PAR_34_5 -> setValue(0);
    for(int k=1;k<=20;k++){
        idToComboBox["cBpm_"+QString::number(100+1+(k-1)*5)+"_1"] -> setCurrentIndex(-1);
        idToLineEdit["LEpm_"+QString::number(100+2+(k-1)*5)+"_1"] -> clear();
        idToLineEdit["LEpm_"+QString::number(100+3+(k-1)*5)+"_1"] -> clear();
        idToLineEdit["LEpm_"+QString::number(100+4+(k-1)*5)+"_1"] -> clear();
        idToLineEdit["LEpm_"+QString::number(100+5+(k-1)*5)+"_1"] -> clear();
        idToComboBox["cBpm_"+QString::number(100+1+(k-1)*5)+"_1"] -> setEnabled(false);
        idToLineEdit["LEpm_"+QString::number(100+2+(k-1)*5)+"_1"] -> setEnabled(false);
        idToLineEdit["LEpm_"+QString::number(100+3+(k-1)*5)+"_1"] -> setEnabled(false);
        idToLineEdit["LEpm_"+QString::number(100+4+(k-1)*5)+"_1"] -> setEnabled(false);
        idToLineEdit["LEpm_"+QString::number(100+5+(k-1)*5)+"_1"] -> setEnabled(false);
    }
    //reset panel fit
    for(int j=1;j<=17;j++){
        idToCheckBox["chBeParFit_"+QString::number(j)] -> setEnabled(false);
        idToComboBox["cBParFit_"+QString::number(j)] -> clear();
        idToComboBox["cBParFit_"+QString::number(j)] -> addItem("none");
        idToCheckBox["chBeParFit_"+QString::number(j)] -> setCheckState ( Qt::Unchecked );
        idToLineEdit["DPparFitV_"+QString::number(j)] -> clear();//setText("");
        idToLineEdit["DPparFitErr_"+QString::number(j)] -> clear();
        idToLineEdit["DPparFitGC_"+QString::number(j)] -> clear();
        idToLineEdit["DPparFitV_"+QString::number(j)] -> setEnabled(false);
        idToLineEdit["DPparFitErr_"+QString::number(j)] -> setEnabled(false);
        idToLineEdit["DPparFitGC_"+QString::number(j)] -> setEnabled(false);
    }

    caller="LoadProject";
    ReadSetting(fnproject);
    SPADA();//load file-nk & measurements
    SaveSetting(-1);
    ifirstWarning=0;
}



void ksemawc::ReadSetting(QString filename){
    if(occupyPF!=0 && occupyPF!=10) return;
    printf("-> ReadSetting (occupyPF=%d) from %s\n",occupyPF,(filename.toStdString()).c_str());
    for(int i=1;i<=8;i++)
        Clrnk(i);
    ClrFnk();
    Clrfn();
    Qt::CheckState state;
    QFileInfo info(filename);
    QDateTime dtTm=info.lastModified();
    QString as = dtTm.toString("yyyy-MM-dd HH:mm:ss");
    cout<<"date: "<<as.toStdString()<<"\n";
    QDateTime tLim(QDate(2022, 1, 1), QTime(0, 00, 0));
    QDateTime tLim2(QDate(2023, 4, 11), QTime(0, 00, 0));
    QDateTime tLim3(QDate(2025, 10, 25), QTime(0, 00, 0));
    int old0new1=0;
    if(dtTm<tLim)
        printf("The project is older than January 2022\n");
    else if(dtTm>=tLim && dtTm<tLim2){
        printf("The project is newer than January 2022 and older than 11 April 2023\n");
        old0new1=1;
    }
    else{
        printf("The project version is up-to-date\n");
        old0new1=2;
    }
    double a2nm=1.;
    if(dtTm<tLim3){
        a2nm=0.1;
        printf("The projet is older than 25 October 2025 wl and thickness will be converted in nm!!!\n");
    }
    QFile file(filename);
    if (!file.open (QIODevice::ReadOnly | QIODevice::Text)){
        occupyPF=0;
        return;
    }
    QTextStream stream ( &file );
    QString line,line2,lab;
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
        NANK[i]=line.simplified();
        idToLineEdit["lineEdit"+QString::number(i)] -> setText(NANK[i]);
        fnk[i]=pathroot+line.simplified()+".nk";
        ifn=0;
        if(!NANK[i].contains("mate/aa999.9")){
            ifn=1;
            Setnk(i);
        }
    }

    //nk-solutions
    line = stream.readLine();
    NANK[9]=line.simplified();
    ui->lineEdit_Fnk->setText(NANK[9]+".nk");
    fnFnk=pathroot+line.simplified()+".nk";
    //standard spectrum for mean computing
    line = stream.readLine();
    NANK[10]=line.simplified();
    //base-name of SF measurement
    line = stream.readLine();
    NANK[11]=line.simplified();
    ui->lineEdit_sample-> setText(NANK[11]);
    fnSample=pathroot+line.simplified();
    ifn=0;
    if(!fnSample.contains("mate/aa999")){
        ifn=1;
        listMeas(fnSample);
    }
    line = stream.readLine();
    NANK[12]=line.simplified();
    fnE1=pathroot+line.simplified()+".el";
    line = stream.readLine();
    NANK[13]=line.simplified();
    fnE2=pathroot+line.simplified()+".el";
    line = stream.readLine();
    NANK[14]=line.simplified();
    fnE3=pathroot+line.simplified()+".el";
    line = stream.readLine();
    NANK[15]=line.simplified();
    fnE4=pathroot+line.simplified()+".el";

    QString pezzo;
    //setSample();

    //RXY
    int DP0cosDtanP1=ui->comboBox_DeltaPsiScale->currentIndex();
    for(int i=1;i<=irx;i++){
        line = stream.readLine();
        line=line.simplified();
        QStringList List0;
        List0 =line.split(" ");
        for(int j=1;j<=4;j++){
            pezzo=List0.at(j-1).toLocal8Bit().constData();
            rxy[i][j]=pezzo.toDouble();
            if(i==6 && j==2 && rxy[i][j]<=0.)
                rxy[i][j]=1.;
            if(i==24 && j==1 && rxy[i][j]< 10.)
                rxy[i][j]=500.;
            if(i==24 && j==2 && rxy[i][j]< 10.)
                rxy[i][j]=300;
            if(i==24 && j==3 && (rxy[i][j]< 1 || rxy[i][j]> 20))
                rxy[i][j]=2.;
            if(idToLineEdit.contains("DP_RXY_"+QString::number(i)+"_"+QString::number(j))){
                int iEff=i;
                if(i==7 && DP0cosDtanP1==1)
                    iEff=19;
                else if(i==8 && DP0cosDtanP1==1)
                    iEff=22;
                idToLineEdit["DP_RXY_"+QString::number(i)+"_"+QString::number(j)]
                        -> setText(QString::number(rxy[iEff][j]));
            }
        }
    }
    rxy[20][1]=rxy[20][1]*a2nm;//WLmin in nm
    rxy[20][2]=rxy[20][2]*a2nm;//Wlmax in nm
    if(NINT(rxy[25][3])==1)
        ui->checkB_RXY_25_3 -> setCheckState ( Qt::Unchecked );
    else
        ui->checkB_RXY_25_3 -> setCheckState ( Qt::Checked );
    if(NINT(rxy[25][4])==1)
        ui->checkB_RXY_25_4 -> setCheckState ( Qt::Unchecked );
    else
        ui->checkB_RXY_25_4 -> setCheckState ( Qt::Checked );

    //CNK
    for(int i=1;i<=15;i++){
        line = stream.readLine();
        // printf("line: %s\n",line.toStdString().c_str());
        // fflush(stdout);
        line=line.simplified();
        QStringList List0;
        List0 =line.split(" ");
        int nV=List0.count();
        //printf("nV=%d\n",nV);
        fflush(stdout);
        for(int j=1;j<=nV;j++){
            pezzo=List0.at(j-1).toLocal8Bit().constData();
            CNK[i][j]=pezzo.toDouble();
            //printf("CNK[%d][%d]=%f\n",i,j,CNK[i][j]);
            //fflush(stdout);
        }
        if(nV==3){
            if(CNK[i][1]<=17){
                CNK[i][4]=-1.;
            }
            else{
                int i1=NINT(CNK[i][1]/1000.);
                int i2=NINT((CNK[i][1]-i1*1000.)/10.);
                CNK[i][1]=i1;
                CNK[i][4]=i2;
            }
            //printf("-> CNK[%d][1]=%f CNK[%d][4]=%f\n",i,CNK[i][1],i,CNK[i][4]);
            //fflush(stdout);
        }
    }

    //PAR
    for(int i=1;i<=60;i++){
        for(int j=1;j<=5;j++){
            stream>> par[i][j];
        }
    }
    par[4][1]=par[4][1]*a2nm;//WLmin in nm
    par[4][2]=par[4][2]*a2nm;//WLmax in nm
    int nItem=0;
    if(!fnE1.contains("mate/aa999")){
        line2=fnE1.section('.', 1, 1);
        nItem=ui->cBmis7 -> count();
        int i=-1;
        do{
            i++;
            lab=ui->cBmis7 -> itemText(i);
        }while(!line2.contains(lab) && i<nItem);
        ui->cBmis7 ->setCurrentIndex(i);
    }
    if(!fnE2.contains("mate/aa999")){
        line2=fnE2.section('.', 1, 1);
        nItem=ui->cBmis9 -> count();
        int i=-1;
        do{
            i++;
            lab=ui->cBmis9 -> itemText(i);
        }while(!line2.contains(lab) && i<nItem);
        ui->cBmis9 ->setCurrentIndex(i);
    }
    if(!fnE3.contains("mate/aa999")){
        line2=fnE3.section('.', 1, 1);
        nItem=ui->cBmis11 -> count();
        int i=-1;
        do{
            i++;
            lab=ui->cBmis11 -> itemText(i);
        }while(!line2.contains(lab) && i<nItem);
        ui->cBmis11 ->setCurrentIndex(i);
    }
    if(!fnE4.contains("mate/aa999")){
        line2=fnE4.section('.', 1, 1);
        nItem=ui->cBmis13 -> count();
        int i=-1;
        do{
            i++;
            lab=ui->cBmis13 -> itemText(i);
        }while(!line2.contains(lab) && i<nItem);
        ui->cBmis13 ->setCurrentIndex(i);
    }
    int JJ=1,Jitem=-1,nMis;
    for(int i=1;i<=14;i++){
        if(NINT(par[i][4])==-1){
            idToCheckBox["checkB_mis"+QString::number(i)+"_1"] -> setCheckState ( Qt::Unchecked );
            idToCheckBox["checkB_mis"+QString::number(i)+"_2"] -> setCheckState ( Qt::Checked );
            idToCheckBox["checkB_mis"+QString::number(i)+"_3"] -> setCheckState ( Qt::Unchecked );
            idToCheckBox["checkB_mis"+QString::number(i)+"_3"] -> setEnabled(false);
            DATO[i]=-1;
        }
        else if(NINT(par[i][4])==0){
            idToCheckBox["checkB_mis"+QString::number(i)+"_1"] -> setCheckState ( Qt::Unchecked );
            idToCheckBox["checkB_mis"+QString::number(i)+"_2"] -> setCheckState ( Qt::Unchecked );
            idToCheckBox["checkB_mis"+QString::number(i)+"_3"] -> setCheckState ( Qt::Unchecked );
            idToCheckBox["checkB_mis"+QString::number(i)+"_3"] -> setEnabled(false);
            DATO[i]=0;
        }
        else if(NINT(par[i][4])==1){
            idToCheckBox["checkB_mis"+QString::number(i)+"_1"] -> setCheckState ( Qt::Checked );
            idToCheckBox["checkB_mis"+QString::number(i)+"_2"] -> setCheckState ( Qt::Checked );
            idToCheckBox["checkB_mis"+QString::number(i)+"_3"] -> setCheckState ( Qt::Unchecked );
            idToCheckBox["checkB_mis"+QString::number(i)+"_3"] -> setEnabled(true);
            DATO[i]=1;
        }
        else if(NINT(par[i][4])==2){
            idToCheckBox["checkB_mis"+QString::number(i)+"_1"] -> setCheckState ( Qt::Checked );
            idToCheckBox["checkB_mis"+QString::number(i)+"_2"] -> setCheckState ( Qt::Checked );
            idToCheckBox["checkB_mis"+QString::number(i)+"_3"] -> setCheckState ( Qt::Checked );
            idToCheckBox["checkB_mis"+QString::number(i)+"_3"] -> setEnabled(true);
            DATO[i]=2;
        }
        nMis=NINT(par[20+JJ][4]);
        if(nMis>=0){
            int iCount=idToComboBox["cBmis"+QString::number(i)] -> count();
            for(int ii=0;ii<iCount;ii++){
                QString Lcb=idToComboBox["cBmis"+QString::number(i)] -> itemText(ii);
                if(Lcb.contains(QString::number(nMis)))
                    idToComboBox["cBmis"+QString::number(i)] -> setCurrentIndex(ii);
            }
            //idToComboBox["cBmis"+QString::number(i)] -> setCurrentIndex(nMis);
        }
        if(i==7 || i==9 || i==11 || i==13){//if(JJ>=7)//set angle and Psi checkBox
            Jitem=idToComboBox["cBteE"+QString::number(JJ-6)] -> findText(QString::number(par[13+JJ-6][1]),Qt::MatchExactly);
            //printf("i=%d JJ=%d Jitem=%d\n",i,JJ,Jitem);
            if(Jitem>=0) idToComboBox["cBteE"+QString::number(JJ-6)] ->setCurrentIndex(Jitem);
            i++;
            if(NINT(par[i][4])==-1){
                idToCheckBox["checkB_mis"+QString::number(i)+"_2"] -> setCheckState ( Qt::Checked );
                idToCheckBox["checkB_mis"+QString::number(i)+"_3"] -> setCheckState ( Qt::Unchecked );
                idToCheckBox["checkB_mis"+QString::number(i)+"_3"] -> setEnabled(false);
                DATO[i]=-1;
            }
            else if(NINT(par[i][4])==0){
                idToCheckBox["checkB_mis"+QString::number(i)+"_2"] -> setCheckState ( Qt::Unchecked );
                idToCheckBox["checkB_mis"+QString::number(i)+"_3"] -> setCheckState ( Qt::Unchecked );
                idToCheckBox["checkB_mis"+QString::number(i)+"_3"] -> setEnabled(false);
                DATO[i]=0;
            }
            else if(NINT(par[i][4])==1){
                idToCheckBox["checkB_mis"+QString::number(i)+"_2"] -> setCheckState ( Qt::Checked );
                idToCheckBox["checkB_mis"+QString::number(i)+"_3"] -> setCheckState ( Qt::Unchecked );
                idToCheckBox["checkB_mis"+QString::number(i)+"_3"] -> setEnabled(true);
                DATO[i]=1;
            }
            else if(NINT(par[i][4])==2){
                idToCheckBox["checkB_mis"+QString::number(i)+"_2"] -> setCheckState ( Qt::Checked );
                idToCheckBox["checkB_mis"+QString::number(i)+"_3"] -> setCheckState ( Qt::Checked );
                idToCheckBox["checkB_mis"+QString::number(i)+"_3"] -> setEnabled(true);
                DATO[i]=2;
            }
        }
        JJ++;
    }
    ui->dSB_PAR_3_3 -> setValue(par[3][3]);//Err DBase/Base
    ui->dSB_PAR_4_1 -> setValue(par[4][1]);//Lmin
    ui->dSB_PAR_4_2 -> setValue(par[4][2]);//Lmax
    ui->dSB_PAR_4_3 -> setValue(par[4][3]);//Err DRef/Ref
    ui->dSB_PAR_5_3 -> setValue(par[5][3]);//ErrReading
    ui->dSB_PAR_6_1 -> setValue(par[6][1]);//theta
    if(caller=="LoadProject"){
        ui->sB_PAR_8_1 -> setValue(0);
        ui->sB_PAR_8_2 -> setValue(0);
        jobtot=0;
        ifirstcall=0;
        ui->spinBox_jobBest -> setValue(0);
        ui->lineEdit_chi2Best->setText("1.E+99");
    }
    if(NINT(par[9][1])==1)
        ui->checkBox_setPsoK -> setCheckState ( Qt::Checked );
    else
        ui->checkBox_setPsoK -> setCheckState ( Qt::Unchecked );
    if(NINT(par[10][1])==0)
        ui->cBox_PAR_10_1 -> setCheckState (Qt::Unchecked);
    else
        ui->cBox_PAR_10_1 -> setCheckState (Qt::Checked);
    if(NINT(par[10][2])==0)
        ui->checkBox_3points -> setCheckState (Qt::Checked);
    else
        ui->checkBox_3points -> setCheckState (Qt::Unchecked);
    if(old0new1==0){
        QMessageBox msgBox;
        msgBox.setText("The project kind is old => k=extinction_coefficient=0");
        msgBox.setInformativeText("Would you like to update to k evaluated from epsilon?");
        msgBox.setStandardButtons(QMessageBox::Yes | QMessageBox::Cancel);
        msgBox.setDefaultButton(QMessageBox::Cancel);
        int ret = msgBox.exec();
        switch (ret) {
        case QMessageBox::Yes:
            par[11][1]=1.;
            break;
        case QMessageBox::Cancel:
            break;
        }
    }
    int expA=0;
    if(par[11][1]>0.)
        expA=NINT(log10(par[11][1]));//attenuation coefficient of k by Fit#N
    else
        expA=0;
    double valA=par[11][1]/pow(10.,expA);
    //printf("par[11][1]=%f expA=%d valA=%f\n",par[11][1],expA,valA);
    ui->spinBox_expAENS->setValue(expA);
    ui->doubleSpinBox_valAENS->setValue(valA);
    int irif=NINT(par[10][3]);
    ui->comboB_PAR_10_3 -> setCurrentIndex(irif);
    if(NINT(par[13][2])==0){
        iDelCon=0;
        ui->checkBox_deltaConnect->setCheckState(Qt::Unchecked);
    }
    else{
        iDelCon=1;
        ui->checkBox_deltaConnect->setCheckState(Qt::Checked);
    }
    iw=NINT(par[18][1]);
    if(iw==1)
        ui->checkB_par_18_1 -> setCheckState ( Qt::Checked );
    else
        ui->checkB_par_18_1 -> setCheckState ( Qt::Unchecked );
    int imultiply=NINT(par[22][1]);
    if(imultiply==0)
        ui->checkB_PAR_22_1 -> setCheckState ( Qt::Unchecked );
    else
        ui->checkB_PAR_22_1 -> setCheckState ( Qt::Checked );
    ui->DP_PAR_27_1  -> setText(QString::number(par[27][1]));
    int s1p2=NINT(par[27][2]);
    ui->cB_PAR_27_2 -> setCurrentIndex(s1p2-1);
    ui->sB_PAR_28_1 -> setValue(NINT(par[28][1]));
    ui->sB_PAR_28_2 -> setValue(NINT(par[28][2]));
    ui->sB_PAR_29_1 -> setValue(NINT(par[29][1]));
    if(NINT(par[31][5])==0)
        ui->checkBox_logScale->setCheckState(Qt::Unchecked);
    else
        ui->checkBox_logScale->setCheckState(Qt::Checked);
    ui->comboBox_DeltaPsiScale ->setCurrentIndex(NINT(par[33][5]));
    nlayer=NINT(par[51][2]);//N, layer
    int isimmetry=NINT(par[52][2]);
    if(isimmetry==0)
        ui->checkB_PAR_52_2 -> setCheckState ( Qt::Unchecked );
    else
        ui->checkB_PAR_52_2 -> setCheckState ( Qt::Checked );
    int ihemi=NINT(par[54][2]);
    if(ihemi==0)
        ui->checkB_PAR_54_2 -> setCheckState ( Qt::Unchecked );
    else
        ui->checkB_PAR_54_2 -> setCheckState ( Qt::Checked );
    ui->DP_PAR_55_2 -> setText(QString::number(par[55][2]));
    npp=NINT(par[34][5]);//N par fit
    ui->cB_PAR_35_2 -> setCurrentIndex(NINT(par[35][2]-1.));
    nPar=NINT(par[35][5]);//N par enabled for fit
    ui->sB_PAR_34_5 -> setValue(npp);
    ui->sB_PAR_35_5 -> setValue(nPar);
    for(int i=1;i<=14;i++){
        state=idToCheckBox["checkB_mis"+QString::number(i)+"_2"] -> checkState();
        if(state==Qt::Checked && DATO[i]>0){
            idToLineEdit["LEpar_"+QString::number(35+i)+"_1"] -> setText(QString::number(par[35+i][1]));
            if(i<=6){
                idToLineEdit["LEpar_"+QString::number(35+i)+"_2"] -> setText(QString::number(par[35+i][2]));
                idToLineEdit["LEpar_"+QString::number(35+i)+"_3"] -> setText(QString::number(par[35+i][3]));
            }
            else{
                idToLineEdit["LEpar_"+QString::number(35+i)+"_2"] -> clear();
                idToLineEdit["LEpar_"+QString::number(35+i)+"_3"] -> clear();
            }
        }
        else{
            idToLineEdit["LEpar_"+QString::number(35+i)+"_1"] -> clear();
            idToLineEdit["LEpar_"+QString::number(35+i)+"_2"] -> clear();
            idToLineEdit["LEpar_"+QString::number(35+i)+"_3"] -> clear();
        }
    }

    ui->comboBox_average->setCurrentIndex(NINT(par[35][1]));//set std spectrum

    for(int i=1;i<=npp;i++)
        ppm[i]=NINT(par[35+i][5]);

    //PM
    for(int i=1;i<=200;i++){
        for(int j=1;j<=5;j++){
            stream>> pm[i][j];
        }
        //printf("pm[%d][1-5]: %f %f %f %f %f\n",i,pm[i][1],pm[i][2],pm[i][3],pm[i][4],pm[i][5]);
    }
    if(a2nm<1.){//Angstrom -> nm
        for(int i=1;i<=9;i++){
            pm[i][1]=pm[i][1]*a2nm;//thickness
            pm[50+i][1]=pm[50+i][1]*a2nm;//roughness
        }
    }
    //printf("\tReadSetting: pm[100][1]= %f\n",pm[100][1]);
    if(old0new1<2){//renormalization of oscillator coeff for old project
        for(int i=1;i<=20;i++){
            if(pm[101+(i-1)*5][1]==1){//Lorentz
                pm[102+(i-1)*5][1]=pm[102+(i-1)*5][1]*PIG/2.;
                pm[104+(i-1)*5][1]=pm[104+(i-1)*5][1]/2.;
            }
            else if(pm[101+(i-1)*5][1]==2){//Quantum
                pm[102+(i-1)*5][1]=pm[102+(i-1)*5][1]*PIG*pm[104+(i-1)*5][1]*pm[103+(i-1)*5][1];
            }
            else if(pm[101+(i-1)*5][1]==3){//Inhomo
                double sqrln2=0.832554611;
                long double E0=abs(pm[103+(i-1)*5][1]);
                long double D=abs(pm[104+(i-1)*5][1]);
                long double Kinho=pow((D/sqrln2),2.)/2*exp(-pow(E0*sqrln2/D,2.))+
                                    E0*D/sqrln2/2.*sqrt(PIG)*erfcl(-E0*sqrln2/D);
                pm[102+(i-1)*5][1]=pm[102+(i-1)*5][1]*Kinho;
            }
            else if(pm[101+(i-1)*5][1]==5){//Drude
                double E0=abs(pm[103+(i-1)*5][1]);
                pm[102+(i-1)*5][1]=E0*E0*PIG/2.;
                pm[103+(i-1)*5][1]=0.;
            }
            else if(pm[101+(i-1)*5][1]==6){//Indirect-Cody
                double E0=abs(pm[103+(i-1)*5][1]);
                double De=abs(pm[105+(i-1)*5][1]);
                double KCi=E0*De+De*De/2.;
                pm[102+(i-1)*5][1]=pm[102+(i-1)*5][1]*KCi;
            }
            else if(pm[101+(i-1)*5][1]==7){//Indirect-Tauc
                double E0=abs(pm[103+(i-1)*5][1]);
                double De=abs(pm[105+(i-1)*5][1]);
                double r=E0/De;
                double KTi=2.*(8.*r*r*log((4.*r+1.)/(4.*r)) + 8.*pow((r+1.),2.)*log((4.*r+4.)/(4.*r+3.))
                           - (1.+8.*r*(r+1.))*log((4.*r+3.)/(4.*r+1.)) );
                pm[102+(i-1)*5][1]=pm[102+(i-1)*5][1]*KTi;
            }
            else if(pm[101+(i-1)*5][1]==8){//Direct-Cody
                double E0=abs(pm[103+(i-1)*5][1]);
                double K=abs(pm[105+(i-1)*5][1]);
                double KCd=PIG/16.*K*K*(2.*E0+K);
                pm[102+(i-1)*5][1]=pm[102+(i-1)*5][1]*KCd*PIG;
            }
            else if(pm[101+(i-1)*5][1]==9){//Direct-Tauc
                double E0=abs(pm[103+(i-1)*5][1]);
                double K=abs(pm[105+(i-1)*5][1]);
                double KTd=PIG/2.*(2.*E0+K-2.*sqrt(E0*(E0+K)));
                pm[102+(i-1)*5][1]=pm[102+(i-1)*5][1]*KTd;
            }
        }
    }

    //PF
    for(int i=1;i<=21;i++){
        for(int j=1;j<=7;j++){
            stream>> pf[j][i];
        }
    }
    file.close();
    //inserimento valori cnk
    double f2=0.;
    for(int i=1;i<=15;i++){
        idToComboBox["cB_cnk"+QString::number(i)+"a"] -> setCurrentIndex(NINT(CNK[i][1]));
        if(NINT(CNK[i][1])==0){
            idToLineEdit["LEcnk"+QString::number(i)+"_2"] -> setEnabled(true);
            idToLineEdit["LEcnk"+QString::number(i)+"_3"] -> setEnabled(true);
            idToLineEdit["LEcnk"+QString::number(i)+"_2"] -> setText(QString::number(CNK[i][2]));
            idToLineEdit["LEcnk"+QString::number(i)+"_3"] -> setText(QString::number(CNK[i][3]));
        }
        else{
            idToLineEdit["LEcnk"+QString::number(i)+"_2"] -> setEnabled(false);
            idToLineEdit["LEcnk"+QString::number(i)+"_3"] -> setEnabled(false);
        }
        if(NINT(CNK[i][4])<0.){
            idToCheckBox["cB_EMA_"+QString::number(i)] -> setCheckState ( Qt::Unchecked );
            idToComboBox["cB_cnk"+QString::number(i)+"b"] -> setCurrentIndex(-1);
            idToComboBox["cB_cnk"+QString::number(i)+"b"] -> setEnabled(false);
            idToLineEdit["LEcnk"+QString::number(i)+"_5"] -> setEnabled(false);
            idToLineEdit["LEcnk"+QString::number(i)+"_6"] -> setEnabled(false);
            idToDoubleSpinBox["dSB_PM_"+QString::number(70+i)+"_1"] -> setEnabled(false);
        }
        else{
            idToCheckBox["cB_EMA_"+QString::number(i)] -> setCheckState ( Qt::Checked );
            idToComboBox["cB_cnk"+QString::number(i)+"b"] -> setEnabled(true);
            idToComboBox["cB_cnk"+QString::number(i)+"b"] -> setCurrentIndex(NINT(CNK[i][4]));
            if(NINT(CNK[i][4])==0){
                idToLineEdit["LEcnk"+QString::number(i)+"_5"] -> setEnabled(true);
                idToLineEdit["LEcnk"+QString::number(i)+"_6"] -> setEnabled(true);
                idToLineEdit["LEcnk"+QString::number(i)+"_5"] -> setText(QString::number(CNK[i][5]));
                idToLineEdit["LEcnk"+QString::number(i)+"_6"] -> setText(QString::number(CNK[i][6]));
            }
            else{
                idToLineEdit["LEcnk"+QString::number(i)+"_5"] -> setEnabled(false);
                idToLineEdit["LEcnk"+QString::number(i)+"_6"] -> setEnabled(false);
            }
            f2=pm[70+i][1];
            idToDoubleSpinBox["dSB_PM_"+QString::number(70+i)+"_1"] -> setEnabled(true);
            idToDoubleSpinBox["dSB_PM_"+QString::number(70+i)+"_1"] -> setValue(f2);
        }
    }
    //aggiornamento valori PanFit
    for(int i=1;i<=npp;i++){
        state=idToCheckBox["chBeParFit_"+QString::number(i)]-> checkState();
        int ip=NINT(ppm[i]);
        if(ip>100)
            pm[ip][1]=abs(pm[ip][1]);
        idToLineEdit["DPparFitV_"+QString::number(i)] -> setEnabled(true);
        idToLineEdit["DPparFitErr_"+QString::number(i)] -> setEnabled(true);
        idToLineEdit["DPparFitGC_"+QString::number(i)] -> setEnabled(true);
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
    for(int i=npp+1;i<=17;i++){
        idToLineEdit["DPparFitV_"+QString::number(i)] -> setEnabled(false);
        idToLineEdit["DPparFitErr_"+QString::number(i)] -> setEnabled(false);
        idToLineEdit["DPparFitGC_"+QString::number(i)] -> setEnabled(false);
    }
    //printf("i0= %d\n",i0);
    ui->cB_cnk1a -> setCurrentIndex(NINT(pm[100][1]));//fit option
    ui->sB_PAR_51_2 -> setValue(nlayer);
    occupyPF=0;

    QString sname=ui->lineEdit_sample -> text();
    if(sname.isEmpty() || sname=="mate/aa999")
        Clrfn();
    MCRange();
    SetModel(nlayer);
    listOsc();
    PanFitPar();

    printf("exit from ReadSetting\n");
}


void ksemawc::setRifMir(){
    Qt::CheckState state;
    state=ui->checkB_PAR_22_1 -> checkState();
    ui->comboB_PAR_10_3-> setEnabled(state==Qt::Checked);
    if(state==Qt::Checked){
        int nrif=ui->comboB_PAR_10_3 -> currentIndex();
        printf("->setRifMir nrif=%d\n",nrif);
        if(nrif!=0){
            QString line,line2,fnam;
            QFile file(fRefMir);
            if(file.open(QIODevice::ReadOnly | QIODevice::Text)){
                QTextStream stream (&file);
                for(int irif=1;irif<=nrif;irif++){
                    line = stream.readLine();
                    line2 = stream.readLine();
                }
                file.close();
                fnam=pathroot+line2.simplified();
                double wmin,wmax;
                previewFile(fnam,"v",info,wmin,wmax);
                ui->WLminRef->setText(QString::number(wmin));
                ui->WLmaxRef->setText(QString::number(wmax));
                MCRange();
            }
            else
                msgErrLoad("setRifMir",fRefMir);
        }
    }
    else
        MCRange();
}



void ksemawc::AdjTheta(){
    double theta=ui->dSB_PAR_6_1 -> value();
    ui->dSB_PAR_6_1bis -> setValue(theta);
    ui->dSB_PAR_6_1tris -> setValue(theta);
//    ui->dSB_PAR_6_1quater -> setValue(theta);
//    ui->dSB_PAR_6_1quinto -> setValue(theta);
//    theta=ui->dSB_PM_86_1 -> value();
//    ui->dSB_PAR_14_1bis -> setValue(theta);
//    theta=ui->dSB_PM_87_1 -> value();
//    ui->dSB_PAR_15_1bis -> setValue(theta);
//    theta=ui->dSB_PM_88_1 -> value();
//    ui->dSB_PAR_16_1bis -> setValue(theta);
//    theta=ui->dSB_PM_89_1 -> value();
//    ui->dSB_PAR_17_1bis -> setValue(theta);
}


void ksemawc::AdjRoughMax(){
    double d,dUP,rghMax;
    for(int i=1;i<10;i++){
        int iKind=idToComboBox["comB_PAR_5"+QString::number(i)+"_3"] ->currentIndex();
        if(iKind>0){
            d=idToDoubleSpinBox["dSB_PM_"+QString::number(i)+"_1"] ->value();
            dUP=d;
            if(i>1)
               dUP=idToDoubleSpinBox["dSB_PM_"+QString::number(i-1)+"_1"] ->value();
            rghMax=min(d/3.,dUP/3.);
            idToDoubleSpinBox["dSB_PM_5"+QString::number(i)+"_1"] ->setMaximum(rghMax);
        }
        else
            idToDoubleSpinBox["dSB_PM_5"+QString::number(i)+"_1"] ->setMaximum(100.);
    }
}



void ksemawc::setRangeEli(){
    printf("-> setRangeEli\n");
    int kD,kP;
    int DP0cosDtanP1=ui->comboBox_DeltaPsiScale->currentIndex();
    if(DP0cosDtanP1==0){
        kD=7;
        kP=8;
    }
    else{
        kD=19;
        kP=22;
    }
    for(int j=1;j<=4;j++){//DP_RXY_7_1
        idToLineEdit["DP_RXY_7_"+QString::number(j)]->setText(QString::number(rxy[kD][j]));
        idToLineEdit["DP_RXY_8_"+QString::number(j)]->setText(QString::number(rxy[kP][j]));
    }
}


void ksemawc::readRangeEli(){
    printf("-> readRangeEli\n");
    int kD,kP;
    int DP0cosDtanP1=ui->comboBox_DeltaPsiScale->currentIndex();
    if(DP0cosDtanP1==0){
        kD=7;
        kP=8;
    }
    else{
        kD=19;
        kP=22;
    }
    for(int j=1;j<=2;j++){//DP_RXY_7_1
        rxy[kD][j]=idToLineEdit["DP_RXY_7_"+QString::number(j)]->text().toDouble();
        rxy[kP][j]=idToLineEdit["DP_RXY_8_"+QString::number(j)]->text().toDouble();
    }
}



void ksemawc::SaveProject(){
    SaveSetting(-1);
    int ierr;
    occupyPF=1;
    QString command,subfnproject,sample;
    subfnproject=ui->lineEdit_P -> text();
    if(subfnproject.contains("mate/aa999.9")){
        sample=ui->lineEdit_sample-> text();
        if(sample.isEmpty())
            subfnproject="";
        else
            subfnproject=sample+".1.Spj";
    }
    fnproject = QFileDialog::getSaveFileName(
        this,
        tr("Filename to save"),
        pathroot+"/"+subfnproject,
        "Semaw Project (*.Spj)");
    if(!fnproject.contains(".Spj"))
        fnproject=fnproject+".Spj";
    QFile file(fnproject);
    subfnproject=fnproject.section(pathroot, 1, 1);
    if(fnproject!=".Spj"){
        ui->lineEdit_P -> setText(subfnproject);
        if (IS_POSIX == 1) {
            command="cp "+fileStore+" "+fnproject;
            ierr=system((command.toStdString()).c_str());}
        else
        {
            command="copy "+fileStore+" "+fnproject;
            QStringList List;
            List =command.split("/");
            int nV=List.count();
            //printf("\nnV=%d\n",nV);
            QString commandw;
            QString pezzo;
            for(int iv=0;iv<nV;iv++){
                pezzo=List.at(iv).toLocal8Bit().constData();
                //printf("pezzo[%d]= %s\t",iv,pezzo.toStdString().c_str());
                if(iv==0){
                    commandw=pezzo;
                }
                else{
                    commandw=commandw+"\\"+pezzo;
                }
            }
            //printf("comando=%s\n",(commandw.toStdString()).c_str());
            ierr=system((commandw.toStdString()).c_str());
        }
        if(ierr != 0){
            msgErrLoad("SaveProject",fnproject);
            printf("Error copying project!!!\n");
        }
        else
            printf("Project saved as %s\n",(fnproject.toStdString()).c_str());
    }
    occupyPF=0;
}



void ksemawc::SaveSetting(int iCall){
    if(occupyPF!=0) return;
    printf("-> SaveSetting (occupyPF=%d option_iCall=%d) to %s\n",
           occupyPF,iCall,(fileStore.toStdString()).c_str());
    occupyPF=1;
    Qt::CheckState state,state1,state2;
    QString stringa,lab,svalue;
    int itab;
    if(iCall<0)
        itab=ui->tabWidget -> currentIndex();
    else
        itab=iCall;
    printf("itab=%d\n",itab);
    if(itab==1){//Model TAB
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
                pm[i][1]=idToDoubleSpinBox["dSB_PM_"+QString::number(i)+"_1"] -> value();
                pm[50+i][1]=idToDoubleSpinBox["dSB_PM_"+QString::number(50+i)+"_1"] -> value();
                pm[90+i][1]=idToDoubleSpinBox["dSB_PM_"+QString::number(90+i)+"_1"] -> value();
                if(NINT(par[50+i][3])==1){
                    pm[i][1]=pm[i][1]/1.E-6;//mm->nm
                } else if(NINT(par[50+i][3])==3){
                    pm[10+i][1]=idToDoubleSpinBox["dSB_PM_"+QString::number(10+i)+"_1"] -> value();
                    pm[20+i][1]=idToDoubleSpinBox["dSB_PM_"+QString::number(20+i)+"_1"] -> value();
                    pm[30+i][1]=idToDoubleSpinBox["dSB_PM_"+QString::number(30+i)+"_1"] -> value();
                    pm[40+i][1]=idToDoubleSpinBox["dSB_PM_"+QString::number(40+i)+"_1"] -> value();
                    pm[60+i][1]=idToDoubleSpinBox["dSB_PM_"+QString::number(60+i)+"_1"] -> value();
                }
            }
        }
        for(int i=71;i<=89;i++)
            pm[i][1]=idToDoubleSpinBox["dSB_PM_"+QString::number(i)+"_1"] -> value();//save f2 and ThetaEli values
        //refresh PanFit values
        for(int i=1;i<=npp;i++){
            state=idToCheckBox["chBeParFit_"+QString::number(i)]-> checkState();
            int ip=NINT(ppm[i]);
            //printf("npp=%d i=%d ip=%d\n",npp,i,ip);
            if(ip>100 && ip <=200)
                pm[ip][1]=abs(pm[ip][1]);
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
//        //enabled measures
//        for(int i=1;i<=13;i++){
//            state=idToCheckBox["checkB_mis"+QString::number(i)+"_1"] -> checkState();
//            idToCheckBox["checkB_mis"+QString::number(i)+"_2"]->setCheckState(state);
//            idToCheckBox["checkB_mis"+QString::number(i)+"_3"]->setCheckState(state);
//            if(state==Qt::Unchecked)
//                idToCheckBox["checkB_mis"+QString::number(i)+"_3"] -> setEnabled(false);
//            else
//                idToCheckBox["checkB_mis"+QString::number(i)+"_3"] -> setEnabled(true);
//            if(i==7 || i==9 || i==11 ) i++;
//        }
    }
    else if(itab==2){//Simulation TAB
        setTabSim();
    }
    else if(itab==3){//Data Analysis
        for(int i=1;i<=14;i++){
            state2=idToCheckBox["checkB_mis"+QString::number(i)+"_3"] -> checkState ();
            if(state2==Qt::Checked){
                        DATO[i]=2;
            }
        }
        printf("\t apply savePanFit\n");
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
                n2=-1;
                ip=0;
                do{
                    n2++;
                } while(!Lj.contains(ParFitLab[n2],Qt::CaseSensitive));
                if(n2==0)
                    ip=90+n1;//dnonuni
                else if(1<=n2 && n2<=7)//film parameters
                    ip=10*(n2-1)+n1;
                else if(8<=n2 && n2<=11)//oscillator parameters
                    ip=100+(n1-1)*5+(n2-7)+1;
                else if(n2==12)//fEMA
                    ip=70+Lj.at(5).digitValue();
                else if(n2==13)//Theta
                    ip=85+Lj.at(6).digitValue();

                ppm[j]=ip;
                Valore=idToLineEdit["DPparFitV_"+QString::number(j)] -> text();
                pm[ip][1]=Valore.toDouble();

                //check value
                if((ip>0 && ip<10) || (ip>50 && ip<60) || ip>100)//thickness or roughness or oscillator
                    if(pm[ip][1]<.0)
                        pm[ip][1]=0.;
                if((ip>=71 && ip<=85) || (ip>=91 && ip<=99)){//fEMA or d_nonuni
                    if(pm[ip][1]<.0)
                        pm[ip][1]=0.;
                    if(pm[ip][1]>1.)
                        pm[ip][1]=1.;
                }
                idToLineEdit["DPparFitV_"+QString::number(j)] -> setText(QString::number(pm[ip][1]));//eventually correct the value

                printf("\tSavePanFit:  j=%d n1=%d n2=%d ip=%d pm[%d][1]=%f\n",j,n1,n2,ip,ip,pm[ip][1]);
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
        ui->sB_PAR_34_5 -> setValue(npp);
        SetModel(nlayer);
        listOsc();
    }
    int ci;
    double f2;
    QFile file(fileStore);
    if (!file.open(QIODevice::WriteOnly | QIODevice::Text)){
        msgErrLoad("SaveSetting-fileStore",fileStore);
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
        NANK[i]=stringa;
        out << stringa << "\n";
    }
    stringa=ui->lineEdit_Fnk -> text();
    QString strS=stringa.simplified();
    if(strS.contains(".nk"))
        stringa=strS.section('.', 0, -2);
    NANK[9]=stringa;
    out << stringa << "\n";
    //media
    ci=ui->comboBox_average -> currentIndex();
    par[35][1]=ci;//N# std spectrum
    //standard spectra for mean value
    QFile file1(fStdSpect);
    if (!file1.open (QIODevice::ReadOnly | QIODevice::Text)){
        msgErrLoad("SaveSetting-fStdSpect",fStdSpect);
        occupyPF=0;
        return;
    }
    QTextStream stream1 ( &file1 );
    QString line;
    for(int i=0;i<=ci;i++){
        line = stream1.readLine();
        line = stream1.readLine();
    }
    NANK[10]=line;
    out << line <<Qt::endl;
    //sampleName
    stringa=ui->lineEdit_sample-> text();
    NANK[11]=stringa;
    out << stringa << "\n";
    //ellipsometric measurements
    lab=ui->cBmis7 -> currentText();
    state=ui->checkB_mis7_1->checkState();
    if(lab.isEmpty() || state==Qt::Unchecked)
        out << "mate/aa999.9" << "\n";
    else{
        NANK[12]=stringa+"."+lab;
        out << stringa+"."+lab << "\n";
        if(pm[86][1]<1.)
            ui->dSB_PM_86_1-> setValue(ui->cBteE1 -> currentText().toDouble());
        else
            ui->dSB_PM_86_1-> setValue(pm[86][1]);
    }
    lab=ui->cBmis9 -> currentText();
    state=ui->checkB_mis9_1->checkState();
    if(lab.isEmpty() || state==Qt::Unchecked)
        out << "mate/aa999.9" << "\n";
    else{
        NANK[13]=stringa+"."+lab;
        out << stringa+"."+lab << "\n";
        if(pm[87][1]<1.)
            ui->dSB_PM_87_1-> setValue(ui->cBteE2 -> currentText().toDouble());
        else
            ui->dSB_PM_87_1-> setValue(pm[87][1]);
    }
    lab=ui->cBmis11 -> currentText();
    state=ui->checkB_mis11_1->checkState();
    if(lab.isEmpty() || state==Qt::Unchecked)
        out << "mate/aa999.9" << "\n";
    else{
        NANK[14]=stringa+"."+lab;
        out << stringa+"."+lab << "\n";
        if(pm[88][1]<1.)
            ui->dSB_PM_88_1-> setValue(ui->cBteE3 -> currentText().toDouble());
        else
            ui->dSB_PM_88_1-> setValue(pm[88][1]);
    }
    lab=ui->cBmis13 -> currentText();
    state=ui->checkB_mis13_1->checkState();
    if(lab.isEmpty() || state==Qt::Unchecked)
        out << "mate/aa999.9" << "\n";
    else{
        NANK[15]=stringa+"."+lab;
        out << stringa+"."+lab << "\n";
        if(pm[89][1]<1.)
            ui->dSB_PM_89_1-> setValue(ui->cBteE4 -> currentText().toDouble());
        else
            ui->dSB_PM_89_1-> setValue(pm[89][1]);
    }
    //RXY
    state=ui->checkB_RXY_25_3 -> checkState();
    if(state==Qt::Unchecked){
        L1E2=1;
        rxy[25][3]=1;
    }else{
        L1E2=2;
        rxy[25][3]=2;
    }
    state=ui->checkB_RXY_25_4 -> checkState();
    if(state==Qt::Unchecked)
        rxy[25][4]=1;
    else
        rxy[25][4]=2;
    int DP0cosDtanP1=ui->comboBox_DeltaPsiScale->currentIndex();
    for(int i=1;i<=30;i++){
        for(int j=1;j<=4;j++){
            if(idToLineEdit.contains("DP_RXY_"+QString::number(i)+"_"+QString::number(j))){
                svalue=idToLineEdit["DP_RXY_"+QString::number(i)+"_"+QString::number(j)] -> text();
                if(i==7 && DP0cosDtanP1==1)
                    rxy[19][j]=svalue.toDouble();
                else if(i==8 && DP0cosDtanP1==1)
                    rxy[22][j]=svalue.toDouble();
                else
                    rxy[i][j]=svalue.toDouble();
            }
        }
    }
    //CNK
    for(int i=1;i<=15;i++){
        CNK[i][1]=idToComboBox["cB_cnk"+QString::number(i)+"a"] -> currentIndex();
        CNK[i][2]=idToLineEdit["LEcnk"+QString::number(i)+"_2"] -> text().toDouble();
        CNK[i][3]=idToLineEdit["LEcnk"+QString::number(i)+"_3"] -> text().toDouble();
        CNK[i][4]=idToComboBox["cB_cnk"+QString::number(i)+"b"] -> currentIndex();
        CNK[i][5]=idToLineEdit["LEcnk"+QString::number(i)+"_5"] -> text().toDouble();
        CNK[i][6]=idToLineEdit["LEcnk"+QString::number(i)+"_6"] -> text().toDouble();
        f2=idToDoubleSpinBox["dSB_PM_"+QString::number(70+i)+"_1"] -> value();
        pm[70+i][1]=f2;
    }
    //PAR
    //  setting measurements
    par[22][2]=0.;//reset enabled measures
    par[23][1]=-1.0; //no SF
    par[23][2]=-1.0; //no ELI
    par[24][1]=1.0;  //IUVIR=1
    int JJ=1;
    for(int i=1;i<=14;i++){
        state =idToCheckBox["checkB_mis"+QString::number(i)+"_1"] -> checkState ();
        state1=idToCheckBox["checkB_mis"+QString::number(i)+"_2"] -> checkState ();
        state2=idToCheckBox["checkB_mis"+QString::number(i)+"_3"] -> checkState ();
        if(state==Qt::Unchecked){
            idToCheckBox["checkB_mis"+QString::number(i)+"_3"] -> setCheckState(Qt::Unchecked);
            idToCheckBox["checkB_mis"+QString::number(i)+"_3"] -> setEnabled(false);
        }else
            idToCheckBox["checkB_mis"+QString::number(i)+"_3"] -> setEnabled(true);
        if(state==Qt::Unchecked && state1==Qt::Unchecked){
            par[i][4]=0;
            DATO[i]=0;
        }
        else if(state==Qt::Unchecked && state1==Qt::Checked){
            par[i][4]=-1;
            DATO[i]=-1;
        }
        else if(state==Qt::Checked && state2==Qt::Unchecked ){
            par[i][4]=1;
            DATO[i]=1;
        }
        else if(state==Qt::Checked && state2==Qt::Checked){
            par[i][4]=2;
            DATO[i]=2;
            par[22][2]++;
        }
        if(i<=6 && state==Qt::Checked && par[23][1]<0.)
            par[23][1]=1.0;
        if(i >6 && state==Qt::Checked && par[23][2]<0.)
            par[23][2]=1.0;
        lab=idToComboBox["cBmis"+QString::number(i)] ->currentText();
        if(i<=6 && state==Qt::Checked){
            if(lab.at(0).toLatin1() == 'i') par[24][1]=2.0;//SF IR
        }
        par[20+JJ][4]=(lab.at(1)).digitValue();// in Qt6 digitValue() deve essere sostituito da unicode()
        //printf("Savesetting JJ=%d (lab.at(1)).digitValue()= %f\n",JJ,par[20+JJ][4]);
        if(JJ>=7){
            i++;
            state1=idToCheckBox["checkB_mis"+QString::number(i)+"_2"] -> checkState ();
            state2=idToCheckBox["checkB_mis"+QString::number(i)+"_3"] -> checkState ();
            if(state==Qt::Unchecked){
                idToCheckBox["checkB_mis"+QString::number(i)+"_3"] -> setCheckState(Qt::Unchecked);
                idToCheckBox["checkB_mis"+QString::number(i)+"_3"] -> setEnabled(false);
            }else
                idToCheckBox["checkB_mis"+QString::number(i)+"_3"] -> setEnabled(true);
            if(state==Qt::Unchecked && state1==Qt::Unchecked){
                par[i][4]=0;
                DATO[i]=0;
            }
            else if(state==Qt::Unchecked && state1==Qt::Checked){
                par[i][4]=-1;
                DATO[i]=-1;
            }
            else if(state==Qt::Checked && state2==Qt::Unchecked ){
                par[i][4]=1;
                DATO[i]=1;
            }
            else if(state==Qt::Checked && state2==Qt::Checked){
                par[i][4]=2;
                DATO[i]=2;
                par[22][2]++;
            }
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
    state=ui->cBox_PAR_10_1 -> checkState ();
    if(state==Qt::Unchecked)
        par[10][1]=0;
    else
        par[10][1]=1;//plot eps1 eps2
    state=ui->checkBox_3points -> checkState ();
    if(state==Qt::Checked)
        par[10][2]=0;
    else
        par[10][2]=1;//no 3 points constrain in ELI resampling
    par[10][3]=ui->comboB_PAR_10_3 -> currentIndex();
    int expA=ui->spinBox_expAENS->value();
    double valA=ui->doubleSpinBox_valAENS->value();
    par[11][1]=valA*pow(10.,expA);
    par[14][1]=ui->cBteE1 -> currentText().toDouble();
    par[15][1]=ui->cBteE2 -> currentText().toDouble();
    par[16][1]=ui->cBteE3 -> currentText().toDouble();
    par[17][1]=ui->cBteE4 -> currentText().toDouble();
    state=ui->checkBox_deltaConnect->checkState();
    par[13][2]=0.;
    iDelCon=0;
    if(state==Qt::Checked){
        par[13][2]=1.;
        iDelCon=1;
    }
    pm[86][1]=ui->dSB_PM_86_1 -> value();
    pm[87][1]=ui->dSB_PM_87_1 -> value();
    pm[88][1]=ui->dSB_PM_88_1 -> value();
    pm[89][1]=ui->dSB_PM_89_1 -> value();
    state=ui->checkB_par_18_1 -> checkState ( );//verboso
    if(state==Qt::Unchecked)
        par[18][1]=0;
    else
        par[18][1]=1;
    state=ui->checkB_PAR_22_1 -> checkState ();//moltiplicazione per Rrif
    if(state==Qt::Unchecked)
        par[22][1]=0;
    else
        par[22][1]=1;
    par[27][2]=ui->cB_PAR_27_2 -> currentIndex();//S1P2
    par[27][2]++;
    par[28][1]=ui->sB_PAR_28_1 -> value();
    par[28][2]=ui->sB_PAR_28_2 -> value();
    NeV=par[28][2];
    par[29][1]=ui->sB_PAR_29_1 -> value();
    state=ui->checkBox_logScale -> checkState();
    if(state==Qt::Unchecked)
        par[31][5]=0.;
    else{
        par[31][5]=1.;
        if(rxy[17][1]<=0.)
            rxy[17][1]=rxy[17][2]/10000.;
    }
    par[33][5]=ui->comboBox_DeltaPsiScale->currentIndex();
    if(fabs(par[21][3]) < 0.003*fabs(rxy[16][2]-rxy[16][1]))
        par[21][3]=0.003*fabs(rxy[16][2]-rxy[16][1]);
    if(fabs(par[22][3]) < 0.003*fabs(rxy[17][2]-rxy[17][1]))
        par[22][3]=0.003*fabs(rxy[17][2]-rxy[17][1]);
    par[35][4]=CNK[1][1];//used in simulation
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
            if(NINT(par[50+i][1])==1) par[53][2]=i;//puntatore strato incognito
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
    par[35][5]=nPar;
    par[51][2]=nlayer;

    pm[100][1]=ui->cB_cnk1a -> currentIndex();//fit option

    //save rxy
    for(int i=1;i<=30;i++){
        for(int j=1;j<=4;j++){
            out << QString::number(rxy[i][j],'g',7) << "\t";
        }
        out << "\n";
    }
    //save cnk
    for(int i=1;i<=15;i++){
        for(int j=1;j<=6;j++){
            out << QString::number(CNK[i][j],'g',7) << "\t";
        }
        out << "\n";
    }
    //save par
    for(int i=1;i<=60;i++){
        for(int j=1;j<=5;j++){
            out << QString::number(par[i][j],'g',7) << "\t";
        }
        out << "\n";
    }
    //save pm
    for(int i=1;i<=200;i++){
        for(int j=1;j<=5;j++){
            out << QString::number(pm[i][j],'g',7) << "\t";
        }
        out << "\n";
    }
    //save pf
    for(int i=1;i<=21;i++){
        for(int j=1;j<=7;j++){
            out << QString::number(pf[j][i],'g',7) <<"\t";
        }
        out << "\n";
    }
    file.close();
    printf("exit from SaveSetting \n");
    occupyPF=0;
}


void ksemawc::msgErrLoad(QString where,QString fnERR){
    fflush(stdout);
    QMessageBox msgBox;
    msgBox.setText(where+": error loading file:\n"+fnERR);
    msgBox.setStandardButtons(QMessageBox::Ok);
    msgBox.exec();
}

void ksemawc::LoadFilenk(){
    double NKNEW[1000][6];
    QString fnam=ui->lineEdit_Fnk -> text();
    fnam=fnam.simplified();
    if(fnam.contains("mate/aa999.9")){
        fnFnk = QFileDialog::getOpenFileName(
                    this,
                    "Choose a file-nk", //title of window
                    pathroot,           //initial directory
                    "file-nk (*.nk)");  //extension file
    }
    else
        fnFnk=pathroot+fnam;
    printf("LoadFilenk-> load file nk-solutions= %s\n",fnFnk.toStdString().c_str());
    QFile file(fnFnk);
    if (!file.open (QIODevice::ReadOnly | QIODevice::Text)){
        printf("\t ERROR!!!\n");
        msgErrLoad("LoadFilenk",fnFnk);
        return;
    }
    QTextStream stream ( &file );
    QString line,pezzo;
    //display filename and save in NANK[9]
    line=fnFnk.section(pathroot, 1, 1);
    ui->lineEdit_Fnk -> setText(line);
    if(line.contains(".nk"))
        line=line.section('.', 0, -2);
    NANK[9]=line;
    //display info
    line = stream.readLine();
    ui->lineEdit_infoFnk -> setText(line.simplified());
    printf("\tInfo: %s\n",(line.simplified()).toStdString().c_str());
    //read Ndata
    line = stream.readLine();
    line=line.simplified();
    QStringList List0;
    List0 =line.split(" ");
    pezzo=List0.at(0).toLocal8Bit().constData();
    int N=pezzo.toInt();
    double a2nm=0.1;
    if(line.contains("nm"))
        a2nm=1.;
    printf("\tN. (L,n,k)= %d a2nm=%f\n",N,a2nm);
    if(N>999){
        msgErrLoad("LoadFilenk","Ndata > 999, MAXIMUM limit");
        return;
    }
    rxy[16][3]=1000.;
    rxy[16][4]=-1000.;
    rxy[17][3]=1000.;
    rxy[17][4]=-1000.;
    SOL[1][1]=N;
    //load data
    for(int L=2;L<=N+1;L++){
        line=stream.readLine();
        line=line.simplified();
        //printf("line=%s\n",line.toStdString().c_str());
        QStringList List;
        List =line.split(" ");
        int nV=List.count();
        if(nV!=5){
            List=line.split("\t");
            nV=List.count();
            if(nV!=5){
                List=line.split(" ");
                nV=List.count();
                if(nV!=5){
                    printf("nV=%d line=%s\n",nV,line.toStdString().c_str());
                    printf("LoadFilenk-> ERROR reading file-nk= %s: separator char invalid\n",fnam.toStdString().c_str());
                    QMessageBox msgBox;
                    msgBox.setText("LoadFilenk-> ERROR reading file-nk="+fnam+"\n: separator char invalid or ERRn and ERRk missing");
                    msgBox.setStandardButtons(QMessageBox::Ok);
                    msgBox.exec();
                    return;
                    continue;
                }
            }
        }
        for(int kk=1;kk<=nV;kk++)
            SOL[L][kk]=List.at(kk-1).toDouble();
        //stream>>SOL[L][1]>>SOL[L][2]>>SOL[L][3]>>SOL[L][4]>>SOL[L][5];
        SOL[L][1]=SOL[L][1]*a2nm;
        //rxy[20][3]=min(rxy[20][3],SOL[L][1]);
        //rxy[20][4]=max(rxy[20][4],SOL[L][1]);
        rxy[16][3]=min(rxy[16][3],SOL[L][2]-SOL[L][4]);
        rxy[16][4]=max(rxy[16][4],SOL[L][2]+SOL[L][4]);
        rxy[17][3]=min(rxy[17][3],SOL[L][3]-SOL[L][5]);
        rxy[17][4]=max(rxy[17][4],SOL[L][3]+SOL[L][5]);
    }
    file.close();
    ui->DP_RXY_16_3->setText(QString::number(rxy[16][3]));
    ui->DP_RXY_16_4->setText(QString::number(rxy[16][4]));
    ui->DP_RXY_17_3->setText(QString::number(rxy[17][3]));
    ui->DP_RXY_17_4->setText(QString::number(rxy[17][4]));
    int NORD=2;
    for(int I=2;I<=N+1;I++){
        double WWM=1.E30;
        for(int H=2;H<=N+1;H++){
            if(SOL[H][1]>=.0 && SOL[H][1]<WWM){
                WWM=SOL[H][1];
                NORD=H;
            }
        }
        for(int k=1;k<=5;k++){
            NKNEW[I][k]=SOL[NORD][k];
        }
        SOL[NORD][1]=-1.;
    }
    int i=2;
    int inew=2;
    while(i<=N+1){
        if(NKNEW[i][4]<=.0)
            NKNEW[i][4]=0.001*NKNEW[i][2];//set minimal n-error
        if(SOL[i][5]<=0.){//set minimal k-error
            if(NKNEW[i][3]>0.)
                NKNEW[i][5]=0.001*NKNEW[i][3];
            else
                NKNEW[i][5]=0.001*(rxy[17][4]-rxy[17][3]);
        }
        for(int k=1;k<=5;k++)
            SOL[inew][k]=NKNEW[i][k];
        SOL[inew][6]=1.;//enabled for fit
        inew++;
        i++;
    }
    SOL[1][1]=inew-2;
    N=inew-2;
    printf("\tRange: Lmin= %f Lmax= %f\n",SOL[2][1],SOL[N+1][1]);
    rxy[20][3]=min(rxy[20][3],SOL[2][1]);
    rxy[20][4]=max(rxy[20][4],SOL[N+1][1]);
}


void ksemawc::ClrFnk(){
    ui->lineEdit_Fnk -> setText("mate/aa999.9");
    ui->lineEdit_infoFnk ->setText("");
    NANK[9]="mate/aa999.9";
}


void ksemawc::SaveFnk(){
    printf("-> SaveFnk with lastAction=%s\n",lastAction.toStdString().c_str());
    QString line;
    int iCase=0;
    int N=NINT(SOL[1][1]);
    if(N==0){
        QMessageBox msgBox;
        msgBox.setText("ATTENTION: launch a best-fit procedure before to save nk-BestFit!\nOtherwise goto Simulation to save nk-Simulated!");
        msgBox.setStandardButtons(QMessageBox::Ok);
        msgBox.exec();
        return;
    }
    QString sample=ui->lineEdit_sample-> text();
    QString suggestedName=sample+"_"+lastAction+".nk";
    fnFnk = QFileDialog::getSaveFileName(
                this,
                "Filename to save",
                pathroot+"/"+suggestedName,
                "nk file (*.nk)");
    if(fnFnk.isEmpty())
        return;
    if(!fnFnk.contains(".nk"))
        fnFnk=fnFnk+".nk";
    QInputDialog pippo;
    pippo.setLabelText("Please set a comment!");
    pippo.setInputMode(QInputDialog::TextInput);
    //pippo.setTextValue(QString::number(expTime));
    pippo.setOkButtonText("Yes");
    //pippo.setCancelButtonText("No, thanks!");
    bool ok=pippo.exec();
    QString text =pippo.textValue();
    if (ok && !text.isEmpty()){
        line=text;
    }
    printf("SaveFnk-> save file.nk= %s\n",fnFnk.toStdString().c_str());
    QFile file(fnFnk);
    if(!file.open(QIODevice::WriteOnly | QIODevice::Text)){
        msgErrLoad("SaveFnk",fnFnk);
        printf("IONK-> ERROR opening file= %s\n",fnFnk.toStdString().c_str());
    }
    int iok=1;
    if( iCase==1 && file.exists() ){
        QMessageBox msgBox;
        msgBox.setText("The file already exist!");
        msgBox.setInformativeText("Do you want to save anyway?");
        msgBox.setStandardButtons(QMessageBox::Save | QMessageBox::No);
        msgBox.setDefaultButton(QMessageBox::Save);
        int ret = msgBox.exec();
        switch (ret) {
        case QMessageBox::Save:
            // Save was clicked
            iok=1;
            break;
        case QMessageBox::No:
            // Cancel was clicked
            iok=0;
            file.close();
            break;
        }
    }
    if(iok==1){
        //line=ui->lineEdit_infoFnk -> text();
        QTextStream stream (&file);
        stream<<line<<"\n";
        stream<<N<<" nm"<<"\n";
        for(int L=2;L<=N+1;L++)
            stream<<SOL[L][1]<<"\t"<<SOL[L][2]<<"\t"<<SOL[L][3]<<"\t"<<SOL[L][4]<<"\t"<<SOL[L][5]<<"\n";
        file.close();
    }
}


void ksemawc::Setnk(int ifile){
    //if(ifn==0){
    if(fnk[ifile].contains("mate/aa999") || fnk[ifile].isEmpty()){
        fnk[ifile]=QFileDialog::getOpenFileName(
                    this,
                    "Choose a file-nk", //title of window
                    pathroot,           //initial directory
                    "file-nk (*.nk)");  //extension file
    }
    QFile file(fnk[ifile]);
    if (!file.open (QIODevice::ReadOnly | QIODevice::Text)){
        msgErrLoad("Setn",fnk[ifile]);
        return;
    }
    printf("->Setnk(%d) with fnk[%d]=%s\n",ifile,ifile,fnk[ifile].toStdString().c_str());
    double val1,val2;
    QTextStream stream ( &file );
    QString line,line2,pezzo;
    line = stream.readLine();//info line
    idToLineEdit["lineEdit_infoNK_"+QString::number(ifile)] -> setText(line.simplified());
    line2=fnk[ifile].section(pathroot, 1, 1);
    idToLineEdit["lineEdit"+QString::number(ifile)] -> setText(line2.section('.', 0, -2));//setText(line2.section('.', 0, 1))
    line = stream.readLine();//Ndat & "nm" for new file.nk
    double a2nm=0.1;
    if(line.contains("nm"))
        a2nm=1.;
    line = stream.readLine();
    line=line.simplified();
    val1=line.section(' ', 0, 0).toDouble();
    line2=line;
    do {
        line=line2;
        line2 = stream.readLine();
    } while (!line2.isNull());
    line=line.simplified();
    val2=line.section(' ', 0, 0).toDouble();
    val1=val1*a2nm;
    val2=val2*a2nm;
    printf("\tval1=%f val2=%f\n",val1,val2);
    if(val1<val2){
        idToLineEdit["WLminNK"+QString::number(ifile)] -> setText(QString::number(NINT(val1)));
        idToLineEdit["WLmaxNK"+QString::number(ifile)] -> setText(QString::number(NINT(val2)));
    }
    else{
        idToLineEdit["WLminNK"+QString::number(ifile)] -> setText(QString::number(NINT(val2)));
        idToLineEdit["WLmaxNK"+QString::number(ifile)] -> setText(QString::number(NINT(val1)));
    }
    file.close();
    MCRange();
    line2=idToLineEdit["lineEdit"+QString::number(ifile)] -> text();
    updateMatenk(7+ifile,line2);
}

void ksemawc::mDwUp(int iLayer, int Dw1UpM1){
    double tmp1=par[50+iLayer][1];
    double tmp3=par[50+iLayer][3];
    par[50+iLayer][1]=par[50+iLayer+Dw1UpM1][1];
    par[50+iLayer][3]=par[50+iLayer+Dw1UpM1][3];
    par[50+iLayer+Dw1UpM1][1]=tmp1;
    par[50+iLayer+Dw1UpM1][3]=tmp3;
    tmp1=par[90+iLayer][1];
    tmp3=par[90+iLayer][3];
    par[90+iLayer][1]=par[90+iLayer+Dw1UpM1][1];
    par[90+iLayer][3]=par[90+iLayer+Dw1UpM1][3];
    par[90+iLayer+Dw1UpM1][1]=tmp1;
    par[90+iLayer+Dw1UpM1][3]=tmp3;
    for(int j=0;j<=6;j++){
        for(int k=1;k<=5;k++){
            double tmp=pm[10*j+iLayer][k];
            pm[10*j+iLayer][k]=pm[10*j+iLayer+Dw1UpM1][k];
            pm[10*j+iLayer+Dw1UpM1][k]=tmp;
        }
    }
    SetModel(nlayer);
}


void ksemawc::Clrfn(){
    if(occupyPF==10)
        return;
    occupyPF=2;
    NANK[11]="mate/aa999";
    ui->lineEdit_sample-> setText(NANK[11]);
    fnSample=pathroot+NANK[11];
    ifn=0;
    for(int imis=1;imis<=13;imis++){
        if(imis==8 || imis==10 || imis==12)
            imis++;
        idToComboBox["cBmis"+QString::number(imis)] -> clear();
        idToCheckBox["checkB_mis"+QString::number(imis)+"_1"] -> setCheckState(Qt::Unchecked);
        idToCheckBox["checkB_mis"+QString::number(imis)+"_3"] -> setCheckState(Qt::Unchecked);
    }
    ui->lineEdit_Tn->setText("");
    ui->lineEdit_Tp->setText("");
    ui->lineEdit_Rn->setText("");
    ui->lineEdit_Rp->setText("");
    ui->lineEdit_R1->setText("");
    ui->lineEdit_Apds->setText("");
    ui->lineEdit_E1->setText("");
    ui->lineEdit_E2->setText("");
    ui->lineEdit_E3->setText("");
    ui->lineEdit_E4->setText("");
    occupyPF=0;
    SaveSetting(0);
    MCRange();
}



void ksemawc::Clrnk(int ifile){
    fnk[ifile]="mate/aa999.9";
    idToLineEdit["lineEdit"+QString::number(ifile)] -> setText("mate/aa999.9");
    idToLineEdit["lineEdit_infoNK_"+QString::number(ifile)] ->setText("");
    idToLineEdit["WLminNK"+QString::number(ifile)] -> setText("");
    idToLineEdit["WLmaxNK"+QString::number(ifile)] -> setText("");
    MCRange();
    updateMatenk(7+ifile,"nk-"+QString::number(ifile));
}

void ksemawc::callSetSample(){
    iNewSample=1;
    ifn=0;
    setSample();
}

void ksemawc::setSample(){
    QString line2;
    if(ifn==0){
        fnSample=QFileDialog::getOpenFileName(
                    this,
                    "Choose a sample", //title of window
                    pathroot,          //initial directory
                    "file-exp (*.tn *.tp *.rn *.rp *.r1 *.an *.el)"); //file extensions
    }
    fnSample=fnSample.section('.', 0, -3);
    line2=fnSample.section(pathroot, 1, 1);
    ui->lineEdit_sample -> setText(line2);
    if(fnSample.isEmpty())
        Clrfn();
}

void ksemawc::listMeas(const QString &){
    if(occupyPF==2)
        return;
    if(fnSample.isEmpty())
        return;
    printf("->listMeas (occupyPF=%d) with fnSample= %s\n",occupyPF,fnSample.toStdString().c_str());
    int i,j;
    QString lab,estens[7];
    estens[1]=".tn";
    estens[2]=".tp";
    estens[3]=".rn";
    estens[4]=".rp";
    estens[5]=".r1";
    estens[6]=".an";

    occupyPF++;
    ui->lineEdit_Tn->clear();
    ui->lineEdit_Tp->clear();
    ui->lineEdit_Rn->clear();
    ui->lineEdit_Rp->clear();
    ui->lineEdit_R1->clear();
    ui->lineEdit_Apds->clear();
    for(int imis=1;imis<=6;imis++){
        idToCheckBox["checkB_mis"+QString::number(imis)+"_1"] ->setCheckState(Qt::Unchecked);
        idToComboBox["cBmis"+QString::number(imis)] -> clear();
        idToLineEdit["WLmin"+QString::number(imis)] -> clear();
        idToLineEdit["WLmax"+QString::number(imis)] -> clear();
        for(j=0;j<2;j++){
            lab="v";
            if(j==1) lab="i";
            for(i=0;i<10;i++){
                QFile file(fnSample+"."+lab+QString::number(i)+estens[imis]);
                if(file.exists()) {
                    idToComboBox["cBmis"+QString::number(imis)] -> addItem(lab+QString::number(i));
                    //printf("addItem %s%d\n",lab.toStdString().c_str(),i);
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
    if(iNewSample==1){
        //assign incremental theta if same measurement
        for(int j=2;j<=4;j++){
            int nELIpre=idToComboBox["cBmis"+QString::number(5+2*(j-1))]->currentIndex();
            int nEli=idToComboBox["cBmis"+QString::number(5+2*j)]->currentIndex();
            printf("iNewSample=%d nELIpre=%d nEli=%d\n",iNewSample,nELIpre,nEli);
            if(nELIpre>=0 && nEli==nELIpre){
                int nTheta=idToComboBox["cBteE"+QString::number(j)] -> count();
                int iThetaPre=idToComboBox["cBteE"+QString::number(j-1)] -> currentIndex();
                int iTheta=idToComboBox["cBteE"+QString::number(j)] -> currentIndex();
                printf("\tnTheta=%d iThetaPre=%d iTheta=%d\n",nTheta,iThetaPre,iTheta);
                if(iTheta<=iThetaPre && iThetaPre<nTheta-1)
                    idToComboBox["cBteE"+QString::number(j)] -> setCurrentIndex(iThetaPre+1);
            }
        }
        iNewSample=0;
    }

    occupyPF--;
}

void ksemawc::pwTn(const int &){
    if(fnSample.contains("mate/aa") || fnSample.isEmpty())
        return;
    double wmin=0.,wmax=0.;
    QString lab,info;
    lab=ui->cBmis1 -> currentText();
    if(lab.isEmpty())
        return;
    fnTn=fnSample+"."+lab+".tn";
    previewFile(fnTn,lab,info,wmin,wmax);
    ui->lineEdit_Tn -> setText(info);
    if(wmin>0. && wmax>0.){
        ui->WLmin1 -> setText(QString::number(wmin));
        ui->WLmax1 -> setText(QString::number(wmax));
    }
    else{
        ui->WLmin1 -> setText("");
        ui->WLmax1 -> setText("");
    }
    MCRange();
}

void ksemawc::pwTp(const int &){
    if(fnSample.contains("mate/aa") || fnSample.isEmpty())
        return;
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
    ui->dSB_PAR_6_1 -> setValue(par[6][1]);
    ui->dSB_PAR_6_1bis -> setValue(par[6][1]);
    ui->dSB_PAR_6_1tris -> setValue(par[6][1]);
    int s1p2=NINT(par[27][2]);
    ui->cB_PAR_27_2 -> setCurrentIndex(s1p2-1);
    MCRange();
}

void ksemawc::pwRn(const int &){
    if(fnSample.contains("mate/aa") || fnSample.isEmpty())
        return;
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

void ksemawc::pwRp(const int &){
    if(fnSample.contains("mate/aa") || fnSample.isEmpty())
        return;
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

void ksemawc::pwR1(const int &){
    if(fnSample.contains("mate/aa") || fnSample.isEmpty())
        return;
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

void ksemawc::pwApds(const int &){
    if(fnSample.contains("mate/aa") || fnSample.isEmpty())
        return;
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


void ksemawc::pwEj(int j){
    printf("pwEj(%d) occupyPF=%d ...",j,occupyPF);
    fflush(stdout);
    if(fnSample.contains("mate/aa") || fnSample.isEmpty())
        return;
    int ntel=0,ndat=0,i,index;
    double theta=0;
    double wmin=0.,wmax=0.;
    idToComboBox["cBteE"+QString::number(j)]-> clear();
    QString lab,info,line,pezzo;
    index=idToComboBox["cBmis"+QString::number(5+2*j)]->currentIndex();
    if(index>=0){
        lab=idToComboBox["cBmis"+QString::number(5+2*j)] -> currentText();
        QString fnEli=fnSample+"."+lab+".el";
        if(j==1)
            fnE1=fnEli;
        else if(j==2)
            fnE2=fnEli;
        else if(j==3)
            fnE3=fnEli;
        else if(j==4)
            fnE4=fnEli;
        QFile file(fnEli);
        if (!file.open (QIODevice::ReadOnly | QIODevice::Text)){
            msgErrLoad("pwEj("+QString::number(j)+")",fnEli);
            return;
        }
        else if(file.exists()) {
            QTextStream stream ( &file );
            info = stream.readLine();
            idToLineEdit["lineEdit_E"+QString::number(j)]-> setText(info);
            line=stream.readLine();
            line=line.simplified();
            if(!line.contains("VASE")){
                QStringList List0;
                List0 =line.split(" ");
                if(List0.count()!=3){
                    QMessageBox msgBox;
                    msgBox.setText("Format ERROR in in the second line of "+fnEli);
                    msgBox.setStandardButtons(QMessageBox::Ok);
                    msgBox.exec();
                    return;
                }
                pezzo=List0.at(0).toLocal8Bit().constData();
                ntel=pezzo.toInt();
                pezzo=List0.at(1).toLocal8Bit().constData();
                ndat=pezzo.toInt();
                info=List0.at(2).toLocal8Bit().constData();
                printf("->pwEj(%d) file=%s Ntheta=%d Ndat=%d format=%s\n",j,
                   fnEli.toStdString().c_str(),ntel,ndat,info.toStdString().c_str());
                EliTab[j][0][0]=ntel;
                for(i=0;i<ntel;i++){
                    stream >> theta >> ndat >> wmin >> wmax;
                    EliTab[j][ntel][0]=theta;
                    EliTab[j][i][1]=ndat;
                    EliTab[j][i][2]=wmin;
                    EliTab[j][i][3]=wmax;
                    printf("theta=%f Ndat=%d WLmin=%f WLmax=%f\n",theta,ndat,wmin,wmax);
                    idToComboBox["cBteE"+QString::number(j)]-> addItem(QString::number(theta));
                }
                int Npez=6;
                if(info.contains("etcab"))
                    Npez=5;

                line=stream.readLine();
                line=stream.readLine();
                line=line.simplified();
                QStringList List;
                List =line.split(" ");
                int nV=List.count();
                printf("found %d columns on %d expected\n",nV,Npez);
                if(nV!=Npez){
                    QMessageBox msgBox;
                    msgBox.setText("Format ERROR in "+fnEli+":\nthe number of column mismatch the declared format "+info);
                    msgBox.setStandardButtons(QMessageBox::Ok);
                    msgBox.exec();
                    //fnE1="";
                    if(j==1)
                        fnE1="";
                    else if(j==2)
                        fnE2="";
                    else if(j==3)
                        fnE3="";
                    else if(j==4)
                        fnE4="";
                    idToLineEdit["lineEdit_E"+QString::number(j)]-> clear();
                    idToComboBox["cBmis"+QString::number(5+2*j)]-> clear();
                }
            }
            else{
                occupyPEj=1;//for disabling pwSubEj
                line=stream.readLine();
                line=stream.readLine();
                if(line.contains("eV") || line.contains("nm")){
                    int iL0E1=0;
                    if(line.contains("eV"))
                        iL0E1=1;
                    double theta2=100.,eV;
                    ntel=-1;
                    do{
                        line=stream.readLine();
                        line=line.simplified();
                        QStringList List;
                        List=line.split(" ");
                        pezzo=List.at(0).toLocal8Bit().constData();
                        if(pezzo.contains("dpolE"))//the line does not contain ellipsometric PSI and DELTA angles
                            continue;
                        pezzo=List.at(1).toLocal8Bit().constData();
                        eV=pezzo.toDouble();
                        pezzo=List.at(2).toLocal8Bit().constData();
                        theta=pezzo.toDouble();
                        //printf("theta=%f\n",theta);
                        if(theta==theta2){
                            ndat++;
                            wmin=min(wmin,eV);
                            wmax=max(wmax,eV);
                            EliTab[j][ntel][1]=ndat;
                            if(iL0E1==0){
                                EliTab[j][ntel][2]=wmin;
                                EliTab[j][ntel][3]=wmax;
                            }
                            else{//eV->nm
                                EliTab[j][ntel][2]=1240./wmax;
                                EliTab[j][ntel][3]=1240./wmin;
                            }
                        }
                        else{
                            printf("found new theta=%f eV=%f wl%f!\n",theta,eV,1240./eV);
                            ntel++;
                            EliTab[j][ntel][0]=theta;
                            EliTab[j][0][0]=ntel+1;
                            ndat=1;
                            theta2=theta;
                            wmin=eV;
                            wmax=eV;
                            idToComboBox["cBteE"+QString::number(j)]-> addItem(QString::number(theta));
                        }
                    }while(!stream.atEnd());
                    EliTab[j][0][0]=ntel+1;
                }
                occupyPEj=0;
                pwSubEj(j);
            }
            file.close();
        }
    }
    else{
        if(j==1)
            fnE1="";
        else if(j==2)
            fnE2="";
        else if(j==3)
            fnE3="";
        else if(j==4)
            fnE4="";
        idToComboBox["cBmis"+QString::number(5+2*j)]-> clear();
        idToLineEdit["lineEdit_E"+QString::number(j)]-> clear();
        idToLineEdit["WLmin"+QString::number(5+2*j)] -> setText("");
        idToLineEdit["WLmax"+QString::number(5+2*j)] -> setText("");
    }
}


void ksemawc::pwSubEj(int j){
    if(occupyPEj==1)
        return;
    int ntel=int(EliTab[j][0][0]+0.5);
    if(ntel>0){
        int index=idToComboBox["cBteE"+QString::number(j)]->currentIndex();//ui->cBteE1 -> currentIndex();
        if(index<0)
            return;
        printf("\t->pwSubEj(%d) ntel=%d index=%d\n",j,ntel,index);
        double theta=idToComboBox["cBteE"+QString::number(j)] -> currentText().toDouble();
        idToDoubleSpinBox["dSB_PM_"+QString::number(85+j)+"_1"]->setValue(theta);
        pm[85+j][1]=theta;
        double wmin=EliTab[j][index][2];
        double wmax=EliTab[j][index][3];
        idToLineEdit["WLmin"+QString::number(5+2*j)]-> setText(QString::number(wmin));
        idToLineEdit["WLmax"+QString::number(5+2*j)]-> setText(QString::number(wmax));
    }
    else{
        idToLineEdit["WLmin"+QString::number(5+2*j)] -> setText("");
        idToLineEdit["WLmax"+QString::number(5+2*j)] -> setText("");
    }
    MCRange();
}

void ksemawc::MCRange(){
    //setting the maximum common wavelenght range
    if(occupyPF!=0) return;
    //printf("-> MCrange (occupyPF=%d)\n",occupyPF);
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
        //if(vmin != 0. && vmax != 0.)
        //    printf("filenk-%d wlmin=%f wlmax=%f\n",ink,Lmin,Lmax);
    }
    if(fnSample.isEmpty()){
        for(int imis=1;imis<=14;imis++){
            for(int k=1;k<=3;k++)
                idToCheckBox["checkB_mis"+QString::number(imis)+"_"+QString::number(k)] -> setCheckState(Qt::Unchecked);
            if(imis>=7) imis++;
        }
    }
    for(int imis=1;imis<=14;imis++){
        state1 = idToCheckBox["checkB_mis"+QString::number(imis)+"_1"] -> checkState();
        int index=idToComboBox["cBmis"+QString::number(imis)]->currentIndex();
        if(index<0){
            idToCheckBox["checkB_mis"+QString::number(imis)+"_1"] -> setCheckState(Qt::Unchecked);
            state1=Qt::Unchecked;
        }
        if( state1 == Qt::Checked ) {
            if(DATO[imis]==0) DATO[imis]=1;
            if(imis>=7 && DATO[imis+1]==0) DATO[imis+1]=1;
            idToCheckBox["checkB_mis"+QString::number(imis)+"_2"] -> setCheckState(Qt::Checked);
            if(imis>=7)
                idToCheckBox["checkB_mis"+QString::number(imis+1)+"_2"] -> setCheckState(Qt::Checked);
            str=idToLineEdit["WLmin"+QString::number(imis)] -> text();
            vmin=str.toDouble();
            str=idToLineEdit["WLmax"+QString::number(imis)] -> text();
            vmax=str.toDouble();
            if(vmin != 0.) Lmin=max(Lmin,vmin);
            if(vmax != 0.) Lmax=min(Lmax,vmax);
        }
        else{
            if(DATO[imis]>0) DATO[imis]=0;
            if(imis>=7 && DATO[imis+1]>0){
                DATO[imis+1]=0;
                idToCheckBox["checkB_mis"+QString::number(imis)+"_3"] -> setCheckState(Qt::Unchecked);
                idToCheckBox["checkB_mis"+QString::number(imis+1)+"_3"] -> setCheckState(Qt::Unchecked);
            }
        }
        if(imis>=7) imis++;
    }
    state1=ui->checkB_PAR_22_1 -> checkState();
    if(state1==Qt::Checked){
        int nrif=ui->comboB_PAR_10_3 -> currentIndex();
        if(nrif!=0){
            vmin=ui->WLminRef->text().toDouble();
            vmax=ui->WLmaxRef->text().toDouble();
            if(vmin != 0.) Lmin=max(Lmin,vmin);
            if(vmax != 0.) Lmax=min(Lmax,vmax);
        }
    }
     printf("\tLmin=%f Lmax=%f\n",Lmin,Lmax);
    if(Lmin>0. && Lmax<1.E+20){
        Qt::CheckState state;
        state = ui->checkBox_msgr -> checkState();
        ui->dSB_PAR_4_1 -> setValue(Lmin);
        ui->dSB_PAR_4_2 -> setValue(Lmax);
        if(state == Qt::Unchecked){
            ui->DP_RXY_20_1 -> setText(QString::number(Lmin));
            ui->DP_RXY_20_2 -> setText(QString::number(Lmax));
        }
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

void ksemawc::updateMatenk(int index, QString newName){
    for(int i=1;i<=15;i++){
        idToComboBox["cB_cnk"+QString::number(i)+"a"]-> setItemText(index,newName);
        idToComboBox["cB_cnk"+QString::number(i)+"b"]-> setItemText(index,newName);
    }
}

void previewFile(QString filename, QString lab,QString& info,double& wmin,double& wmax){
    if(lab.isEmpty())
        return;
    printf("->previewFile %s occupyPF=%d lab=%s\n",filename.toStdString().c_str(),occupyPF,lab.toStdString().c_str());
    int i,ndati,ilinrim;
    double ridelta,dini,dfin,dmin,dmax,val;
    QString line,line0,pezzo;
    QFile file(filename);
    if (!file.open (QIODevice::ReadOnly | QIODevice::Text)){
        QMessageBox msgBox;
        msgBox.setText("previewFile: error loading file:\n"+filename);
        msgBox.setStandardButtons(QMessageBox::Ok);
        msgBox.exec();
        return;
    }
    else if(file.exists()) {
        QTextStream stream ( &file );
        info = stream.readLine();
        line0 = stream.readLine();
          //printf("preview->line0= %s \n",(line0.toStdString()).c_str());
        if(info.contains("PE UV", Qt::CaseInsensitive)){
            // file Perkin Elmer L900 - L950
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
            stream >> ndati;
            wmin=wmax+ridelta*(ndati-1);
                 printf("file Perkin Elmer L900 - L950 \n");
                 printf(" ridelta=%f ndati=%d wmin= %f wmax= %f \n",ridelta,ndati,wmin,wmax);
        }
        else if(line0.contains("#####SCALED") || line0.contains("#####SCALEDA") ||
                line0.contains("#####SCALED%")){
            // SCALED file
            //printf("\tfile SCALED %s\n",line0.toStdString().c_str());
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
        }
        else if(line0.contains("RTmethod[")){
            // file Transmittance by VASE
            line = stream.readLine();
            line = stream.readLine();
            int iL0E1=0;
            if(line.contains("eV"))
                iL0E1=1;
            line = stream.readLine();
            line=line.simplified();
            QStringList List;
            List=line.split(" ");
            pezzo=List.at(0).toLocal8Bit().constData();
            par[27][2]=1.;
            if(pezzo.contains("pT"))
                par[27][2]=2.;
            pezzo=List.at(1).toLocal8Bit().constData();
            dini=pezzo.toDouble();
            pezzo=List.at(2).toLocal8Bit().constData();//theta
            par[6][1]=pezzo.toDouble();
            printf("\tTheta_Tp=%f\n",par[6][1]);
            do{
                line = stream.readLine();
            }while(!stream.atEnd());
            line=line.simplified();
            List=line.split(" ");
            pezzo=List.at(1).toLocal8Bit().constData();
            dfin=pezzo.toDouble();
            if(iL0E1==0){
                wmin=dini;
                wmax=dfin;
            }
            else{
                wmin=1240./dini;
                wmax=1240./dfin;
            }
        }
        else if(line0.contains("##XYUNITS= W")){
            // file Lambda9 or old FTIR Perkin Elmer
            stream >> line >> dmax >> line >> val;
            stream >> line >> dmin >> line >> val;
            //    printf("old fashion file dmin= %f dmax= %f \n",dmin,dmax);
            if(lab.contains("v")){
                wmin=dmin;
                wmax=dmax;
            }
            else{
                wmin=(1.E+7)/dmax;
                wmax=(1.E+7)/dmin;
            }
        }
        else{
            line0 = stream.readLine();
            if(line0.contains("##DATA TYPE= UL")){
                //file Lambda9 new type
                for(i=0;i<8;i++) line = stream.readLine();
                stream >> line >> dmin;
                stream >> line >> dmax;
                wmin=dmin;
                wmax=dmax;
            }
            else if(line0.contains("##DATA TYPE= UV")){
                //file Lambda19 Perkin Elmer
                for(i=0;i<11;i++) line = stream.readLine();
                stream >> line >> dmax;
                stream >> line >> dmin;
                wmin=dmin;
                wmax=dmax;
            }
            else if(line0.contains("##DATA TYPE= IN")){
                //file FTIR new type
                for(i=0;i<14;i++) line = stream.readLine();
                stream >> line >> dmin;
                stream >> line >> dmax;
                wmin=(1.E+7)/dmin;
                wmax=(1.E+7)/dmax;
            }
        }
    }
    //printf("wmin= %f wmax= %f \n",wmin,wmax);
}

void ksemawc::SetModel(const int &){
    int itab=ui->tabWidget -> currentIndex();
    if(occupyPF!=0 && itab!=3) return;
    occupyPF=1;
    int nlayerOld=nlayer;
    nlayer=ui->sB_PAR_51_2 -> value();
    int Dnl=nlayer-nlayerOld;
    printf("-> SetModel with Dnl=%d\n",Dnl);
    if(Dnl>0){
        for(int i=0;i<nlayerOld;i++){
            par[50+nlayer-i][1]=par[50+nlayerOld-i][1];
            par[50+nlayer-i][3]=par[50+nlayerOld-i][3];
            for(int j=0;j<=7;j++){
                if(j==7) j=9;
                for(int k=1;k<=5;k++)
                    pm[10*j+nlayer-i][k]=pm[10*j+nlayerOld-i][k];
            }
        }
    }
    else{
        for(int i=1;i<=nlayer;i++){
            par[50+i][1]=par[50+i-Dnl][1];
            par[50+i][3]=par[50+i-Dnl][3];
            for(int j=0;j<=7;j++){
                if(j==7) j=9;
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
            if(NINT(par[50+i][3])<1 || NINT(par[50+i][3])>5) par[50+i][3]=1.;
            if(NINT(par[50+i][1])<1 || NINT(par[50+i][1])>9) par[50+i][1]=1.;
            idToComboBox["comB_PAR_5"+QString::number(i)+"_1"] -> setEnabled(true);
            idToComboBox["comB_PAR_5"+QString::number(i)+"_3"] -> setEnabled(true);
            idToComboBox["comB_PAR_5"+QString::number(i)+"_1"] -> setCurrentIndex(NINT(par[50+i][1])-1);
            idToComboBox["comB_PAR_5"+QString::number(i)+"_3"] -> setCurrentIndex(NINT(par[50+i][3])-1);
            if(NINT(par[50+i][3])==1){
                idToDoubleSpinBox["dSB_PM_"+QString::number(i)+"_1"]    -> setEnabled(true);
                idToDoubleSpinBox["dSB_PM_"+QString::number(i)+"_1"]    -> setValue(pm[i][1]*1.E-6);
                idToDoubleSpinBox["dSB_PM_"+QString::number(10+i)+"_1"] -> setEnabled(false);
                idToDoubleSpinBox["dSB_PM_"+QString::number(20+i)+"_1"] -> setEnabled(false);
                idToDoubleSpinBox["dSB_PM_"+QString::number(30+i)+"_1"] -> setEnabled(false);
                idToDoubleSpinBox["dSB_PM_"+QString::number(40+i)+"_1"] -> setEnabled(false);
                idToDoubleSpinBox["dSB_PM_"+QString::number(60+i)+"_1"] -> setEnabled(false);
                idToDoubleSpinBox["dSB_PM_"+QString::number(50+i)+"_1"] -> setEnabled(true);
                idToDoubleSpinBox["dSB_PM_"+QString::number(50+i)+"_1"] -> setValue(pm[50+i][1]);
                idToDoubleSpinBox["dSB_PM_"+QString::number(90+i)+"_1"] -> setEnabled(false);
            } else if(NINT(par[50+i][3])==2){
                idToDoubleSpinBox["dSB_PM_"+QString::number(i)+"_1"]    -> setEnabled(true);
                idToDoubleSpinBox["dSB_PM_"+QString::number(i)+"_1"]    -> setValue(pm[i][1]);
                idToDoubleSpinBox["dSB_PM_"+QString::number(10+i)+"_1"] -> setEnabled(false);
                idToDoubleSpinBox["dSB_PM_"+QString::number(20+i)+"_1"] -> setEnabled(false);
                idToDoubleSpinBox["dSB_PM_"+QString::number(30+i)+"_1"] -> setEnabled(false);
                idToDoubleSpinBox["dSB_PM_"+QString::number(40+i)+"_1"] -> setEnabled(false);
                idToDoubleSpinBox["dSB_PM_"+QString::number(60+i)+"_1"] -> setEnabled(false);
                idToDoubleSpinBox["dSB_PM_"+QString::number(50+i)+"_1"] -> setEnabled(true);
                idToDoubleSpinBox["dSB_PM_"+QString::number(50+i)+"_1"] -> setValue(pm[50+i][1]);
                idToDoubleSpinBox["dSB_PM_"+QString::number(90+i)+"_1"] -> setEnabled(true);
                idToDoubleSpinBox["dSB_PM_"+QString::number(90+i)+"_1"] -> setValue(pm[90+i][1]);
            } else if(NINT(par[50+i][3])==3){
                idToDoubleSpinBox["dSB_PM_"+QString::number(i)+"_1"]    -> setEnabled(true);
                idToDoubleSpinBox["dSB_PM_"+QString::number(10+i)+"_1"] -> setEnabled(true);
                idToDoubleSpinBox["dSB_PM_"+QString::number(20+i)+"_1"] -> setEnabled(true);
                idToDoubleSpinBox["dSB_PM_"+QString::number(30+i)+"_1"] -> setEnabled(true);
                idToDoubleSpinBox["dSB_PM_"+QString::number(40+i)+"_1"] -> setEnabled(true);
                idToDoubleSpinBox["dSB_PM_"+QString::number(50+i)+"_1"] -> setEnabled(true);
                idToDoubleSpinBox["dSB_PM_"+QString::number(60+i)+"_1"] -> setEnabled(true);
                idToDoubleSpinBox["dSB_PM_"+QString::number(90+i)+"_1"] -> setEnabled(true);
                idToDoubleSpinBox["dSB_PM_"+QString::number(i)+"_1"]    -> setValue(pm[i][1]);
                idToDoubleSpinBox["dSB_PM_"+QString::number(10+i)+"_1"] -> setValue(pm[10+i][1]);
                idToDoubleSpinBox["dSB_PM_"+QString::number(20+i)+"_1"] -> setValue(pm[20+i][1]);
                idToDoubleSpinBox["dSB_PM_"+QString::number(30+i)+"_1"] -> setValue(pm[30+i][1]);
                idToDoubleSpinBox["dSB_PM_"+QString::number(40+i)+"_1"] -> setValue(pm[40+i][1]);
                idToDoubleSpinBox["dSB_PM_"+QString::number(50+i)+"_1"] -> setValue(pm[50+i][1]);
                idToDoubleSpinBox["dSB_PM_"+QString::number(60+i)+"_1"] -> setValue(pm[60+i][1]);
                idToDoubleSpinBox["dSB_PM_"+QString::number(90+i)+"_1"] -> setValue(pm[90+i][1]);
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
            idToDoubleSpinBox["dSB_PM_"+QString::number(90+i)+"_1"] -> setEnabled(false);
        }
    }
    //refresh f2 values
    for(int i=71;i<=85;i++)
        idToDoubleSpinBox["dSB_PM_"+QString::number(i)+"_1"] -> setValue(pm[i][1]);
    occupyPF=0;
}


void ksemawc::RefreshModel(){
    if(occupyPF!=1){
        printf("-> RefreshModel\n");
        for(int i=1;i<=9;i++){
            if(i <= nlayer){
                par[50+i][3]=idToComboBox["comB_PAR_5"+QString::number(i)+"_3"] -> currentIndex();
                par[50+i][3]++;
                par[50+i][1]=idToComboBox["comB_PAR_5"+QString::number(i)+"_1"] -> currentIndex();
                par[50+i][1]++;
                if(NINT(par[50+i][3])==1){
                    idToDoubleSpinBox["dSB_PM_"+QString::number(i)+"_1"] -> setEnabled(true);
                    //     idToDoubleSpinBox["dSB_PM_"+QString::number(i)+"_1"] -> setValue(pm[i][1]*1.E-6);
                    idToDoubleSpinBox["dSB_PM_"+QString::number(10+i)+"_1"] -> setEnabled(false);
                    idToDoubleSpinBox["dSB_PM_"+QString::number(20+i)+"_1"] -> setEnabled(false);
                    idToDoubleSpinBox["dSB_PM_"+QString::number(30+i)+"_1"] -> setEnabled(false);
                    idToDoubleSpinBox["dSB_PM_"+QString::number(40+i)+"_1"] -> setEnabled(false);
                    idToDoubleSpinBox["dSB_PM_"+QString::number(60+i)+"_1"] -> setEnabled(false);
                    idToDoubleSpinBox["dSB_PM_"+QString::number(50+i)+"_1"] -> setEnabled(true);
                    idToDoubleSpinBox["dSB_PM_"+QString::number(90+i)+"_1"] -> setEnabled(false);
                    //     idToDoubleSpinBox["dSB_PM_"+QString::number(i)+"_1"] -> setValue(pm[50+i][1]);
                } else if(NINT(par[50+i][3])==2){
                    idToDoubleSpinBox["dSB_PM_"+QString::number(i)+"_1"] -> setEnabled(true);
                    //     idToDoubleSpinBox["dSB_PM_"+QString::number(i)+"_1"] -> setValue(pm[i][1]);
                    idToDoubleSpinBox["dSB_PM_"+QString::number(10+i)+"_1"] -> setEnabled(false);
                    idToDoubleSpinBox["dSB_PM_"+QString::number(20+i)+"_1"] -> setEnabled(false);
                    idToDoubleSpinBox["dSB_PM_"+QString::number(30+i)+"_1"] -> setEnabled(false);
                    idToDoubleSpinBox["dSB_PM_"+QString::number(40+i)+"_1"] -> setEnabled(false);
                    idToDoubleSpinBox["dSB_PM_"+QString::number(60+i)+"_1"] -> setEnabled(false);
                    idToDoubleSpinBox["dSB_PM_"+QString::number(50+i)+"_1"] -> setEnabled(true);
                    idToDoubleSpinBox["dSB_PM_"+QString::number(50+i)+"_1"] -> setValue(pm[50+i][1]);
                    idToDoubleSpinBox["dSB_PM_"+QString::number(90+i)+"_1"] -> setEnabled(true);
                    idToDoubleSpinBox["dSB_PM_"+QString::number(90+i)+"_1"] -> setValue(pm[90+i][1]);
                } else if(NINT(par[50+i][3])==3){
                    idToDoubleSpinBox["dSB_PM_"+QString::number(i)+"_1"] -> setEnabled(true);
                    idToDoubleSpinBox["dSB_PM_"+QString::number(10+i)+"_1"] -> setEnabled(true);
                    idToDoubleSpinBox["dSB_PM_"+QString::number(20+i)+"_1"] -> setEnabled(true);
                    idToDoubleSpinBox["dSB_PM_"+QString::number(30+i)+"_1"] -> setEnabled(true);
                    idToDoubleSpinBox["dSB_PM_"+QString::number(40+i)+"_1"] -> setEnabled(true);
                    idToDoubleSpinBox["dSB_PM_"+QString::number(50+i)+"_1"] -> setEnabled(true);
                    idToDoubleSpinBox["dSB_PM_"+QString::number(60+i)+"_1"] -> setEnabled(true);
                    idToDoubleSpinBox["dSB_PM_"+QString::number(90+i)+"_1"] -> setEnabled(true);
                    //     idToDoubleSpinBox["dSB_PM_"+QString::number(i)+"_1"] -> setValue(pm[i][1]);
                    idToDoubleSpinBox["dSB_PM_"+QString::number(10+i)+"_1"] -> setValue(pm[10+i][1]);
                    idToDoubleSpinBox["dSB_PM_"+QString::number(20+i)+"_1"] -> setValue(pm[20+i][1]);
                    idToDoubleSpinBox["dSB_PM_"+QString::number(30+i)+"_1"] -> setValue(pm[30+i][1]);
                    idToDoubleSpinBox["dSB_PM_"+QString::number(40+i)+"_1"] -> setValue(pm[40+i][1]);
                    idToDoubleSpinBox["dSB_PM_"+QString::number(50+i)+"_1"] -> setValue(pm[50+i][1]);
                    idToDoubleSpinBox["dSB_PM_"+QString::number(60+i)+"_1"] -> setValue(pm[60+i][1]);
                    idToDoubleSpinBox["dSB_PM_"+QString::number(90+i)+"_1"] -> setValue(pm[90+i][1]);
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
                idToDoubleSpinBox["dSB_PM_"+QString::number(90+i)+"_1"] -> setEnabled(false);
            }
        }
    }
}

void ksemawc::setMat(int nM){
    printf("->setMat called by cnk %d\n",nM);
    int curI=idToComboBox["cB_cnk"+QString::number(nM)+"a"] -> currentIndex();
    idToLineEdit["LEcnk"+QString::number(nM)+"_2"] -> setEnabled(curI==0);
    idToLineEdit["LEcnk"+QString::number(nM)+"_3"] -> setEnabled(curI==0);
    if(idToComboBox["cB_cnk"+QString::number(nM)+"b"]->isEnabled()){
        int curI2=idToComboBox["cB_cnk1b"] ->currentIndex();
        idToLineEdit["LEcnk"+QString::number(nM)+"_5"] -> setEnabled(curI2==0);
        idToLineEdit["LEcnk"+QString::number(nM)+"_6"] -> setEnabled(curI2==0);
    }
}

void ksemawc::setEMA(int nM){
    Qt::CheckState state;
    state=idToCheckBox["cB_EMA_"+QString::number(nM)] -> checkState();
    printf("->setEMA called by checkBox %d\n",nM);
    idToComboBox["cB_cnk"+QString::number(nM)+"b"] -> setEnabled(state==Qt::Checked);
    idToDoubleSpinBox["dSB_PM_"+QString::number(70+nM)+"_1"] -> setEnabled(state==Qt::Checked);
    if(state!=Qt::Checked)
        idToComboBox["cB_cnk"+QString::number(nM)+"b"] ->setCurrentIndex(-1);
    int curI=idToComboBox["cB_cnk"+QString::number(nM)+"b"] ->currentIndex();
    idToLineEdit["LEcnk"+QString::number(nM)+"_5"] -> setEnabled(curI==0 && state==Qt::Checked);
    idToLineEdit["LEcnk"+QString::number(nM)+"_6"] -> setEnabled(curI==0 && state==Qt::Checked);
}

void ksemawc::setKindOsc1(){
    if(occupyPF!=0) return;
    int iKind=ui->cBpm_101_1->currentIndex();
    pm[101][1]=iKind+1;
    //printf("Osc1->%d type\n",iKind+1);
    listOsc();
}

void ksemawc::setKindOsc2(){
    if(occupyPF!=0) return;
    int iKind=ui->cBpm_106_1->currentIndex();
    pm[106][1]=iKind+1;
    //printf("Osc2->%d type\n",iKind+1);
    listOsc();
}

void ksemawc::setKindOsc3(){
    if(occupyPF!=0) return;
    int iKind=ui->cBpm_111_1->currentIndex();
    pm[111][1]=iKind+1;
    //printf("Osc3->%d type\n",iKind+1);
    listOsc();
}

void ksemawc::setKindOsc4(){
    if(occupyPF!=0) return;
    int iKind=ui->cBpm_116_1->currentIndex();
    pm[116][1]=iKind+1;
    //printf("Osc4->%d type\n",iKind+1);
    listOsc();
}

void ksemawc::setKindOsc5(){
    if(occupyPF!=0) return;
    int iKind=ui->cBpm_121_1->currentIndex();
    pm[121][1]=iKind+1;
    //printf("Osc5->%d type\n",iKind+1);
    listOsc();
}

void ksemawc::setKindOsc6(){
    if(occupyPF!=0) return;
    int iKind=ui->cBpm_126_1->currentIndex();
    pm[126][1]=iKind+1;
    //printf("Osc6->%d type\n",iKind+1);
    listOsc();
}

void ksemawc::setKindOsc7(){
    if(occupyPF!=0) return;
    int iKind=ui->cBpm_131_1->currentIndex();
    pm[131][1]=iKind+1;
    //printf("Osc7->%d type\n",iKind+1);
    listOsc();
}

void ksemawc::setKindOsc8(){
    if(occupyPF!=0) return;
    int iKind=ui->cBpm_136_1->currentIndex();
    pm[136][1]=iKind+1;
    //printf("Osc8->%d type\n",iKind+1);
    listOsc();
}

void ksemawc::setKindOsc9(){
    if(occupyPF!=0) return;
    int iKind=ui->cBpm_141_1->currentIndex();
    pm[141][1]=iKind+1;
    //printf("Osc9->%d type\n",iKind+1);
    listOsc();
}

void ksemawc::setKindOsc10(){
    if(occupyPF!=0) return;
    int iKind=ui->cBpm_146_1->currentIndex();
    pm[146][1]=iKind+1;
    //printf("Osc10->%d type\n",iKind+1);
    listOsc();
}

void ksemawc::setKindOsc11(){
    if(occupyPF!=0) return;
    int iKind=ui->cBpm_151_1->currentIndex();
    pm[151][1]=iKind+1;
    //printf("Osc11->%d type\n",iKind+1);
    listOsc();
}

void ksemawc::setKindOsc12(){
    if(occupyPF!=0) return;
    int iKind=ui->cBpm_156_1->currentIndex();
    pm[156][1]=iKind+1;
    //printf("Osc12->%d type\n",iKind+1);
    listOsc();
}

void ksemawc::setKindOsc13(){
    if(occupyPF!=0) return;
    int iKind=ui->cBpm_161_1->currentIndex();
    pm[161][1]=iKind+1;
    //printf("Osc13->%d type\n",iKind+1);
    listOsc();
}

void ksemawc::setKindOsc14(){
    if(occupyPF!=0) return;
    int iKind=ui->cBpm_166_1->currentIndex();
    pm[166][1]=iKind+1;
    //printf("Osc14->%d type\n",iKind+1);
    listOsc();
}

void ksemawc::setKindOsc15(){
    if(occupyPF!=0) return;
    int iKind=ui->cBpm_171_1->currentIndex();
    pm[171][1]=iKind+1;
    //printf("Osc15->%d type\n",iKind+1);
    listOsc();
}

void ksemawc::setKindOsc16(){
    if(occupyPF!=0) return;
    int iKind=ui->cBpm_176_1->currentIndex();
    pm[176][1]=iKind+1;
    //printf("Osc16->%d type\n",iKind+1);
    listOsc();
}

void ksemawc::setKindOsc17(){
    if(occupyPF!=0) return;
    int iKind=ui->cBpm_181_1->currentIndex();
    pm[181][1]=iKind+1;
    //printf("Osc17->%d type\n",iKind+1);
    listOsc();
}

void ksemawc::setKindOsc18(){
    if(occupyPF!=0) return;
    int iKind=ui->cBpm_186_1->currentIndex();
    pm[186][1]=iKind+1;
    //printf("Osc18->%d type\n",iKind+1);
    listOsc();
}

void ksemawc::setKindOsc19(){
    if(occupyPF!=0) return;
    int iKind=ui->cBpm_191_1->currentIndex();
    pm[191][1]=iKind+1;
    //printf("Osc19->%d type\n",iKind+1);
    listOsc();
}

void ksemawc::setKindOsc20(){
    if(occupyPF!=0) return;
    int iKind=ui->cBpm_196_1->currentIndex();
    pm[196][1]=iKind+1;
    //printf("Osc20->%d type\n",iKind+1);
    listOsc();
}

void ksemawc::setOscN(int k){
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


void ksemawc::listOsc(){
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
                if(NINT(pf[ioptFit][i])==k)iok=1;
            }
            if(iok==1){
                int ifu=NINT(pm[100+1+(k-1)*5][1]);
                idToCheckBox["cBosc_"+QString::number(k)]-> setCheckState ( Qt::Checked );
                idToComboBox["cBpm_"+QString::number(100+1+(k-1)*5)+"_1"] -> setEnabled(true);
                idToComboBox["cBpm_"+QString::number(100+1+(k-1)*5)+"_1"] -> setCurrentIndex(NINT(pm[100+1+(k-1)*5][1])-1);
                //idToComboBox["cBpm_"+QString::number(100+1+(k-1)*5)+"_1"] -> setSizeAdjustPolicy(QComboBox::AdjustToContents);
                for(int iC=2;iC<=5;iC++){
                    idToLineEdit["LEpm_"+QString::number(100+iC+(k-1)*5)+"_1"] -> setEnabled(true);
                    idToLineEdit["LEpm_"+QString::number(100+iC+(k-1)*5)+"_1"] -> setText("");
                }
                if(ifu<4){//Lorentz, Q-homo, Q-inhomo
                    idToLineEdit["LEpm_"+QString::number(100+5+(k-1)*5)+"_1"] -> setEnabled(false);
                    idToLineEdit["LEpm_"+QString::number(100+5+(k-1)*5)+"_1"] -> setText("unused");
                }
                if(ifu==4){//Flat
                    for(int iC=3;iC<=5;iC++){
                        idToLineEdit["LEpm_"+QString::number(100+iC+(k-1)*5)+"_1"] -> setEnabled(false);
                        idToLineEdit["LEpm_"+QString::number(100+iC+(k-1)*5)+"_1"] -> setText("unused");
                    }
                }
                else if(ifu==5 || ifu==17){//Drude or Drude-ionized
                    idToLineEdit["LEpm_"+QString::number(100+3+(k-1)*5)+"_1"] -> setEnabled(false);
                    idToLineEdit["LEpm_"+QString::number(100+5+(k-1)*5)+"_1"] -> setEnabled(false);
                    idToLineEdit["LEpm_"+QString::number(100+3+(k-1)*5)+"_1"] -> setText("unused");
                    idToLineEdit["LEpm_"+QString::number(100+5+(k-1)*5)+"_1"] -> setText("unused");
                }
                for(int iC=2;iC<=5;iC++){
                    QString content=idToLineEdit["LEpm_"+QString::number(100+iC+(k-1)*5)+"_1"] -> text();
                    if(content!="unused")
                       idToLineEdit["LEpm_"+QString::number(100+iC+(k-1)*5)+"_1"] -> setText(QString::number(abs(pm[100+iC+(k-1)*5][1])));
                }
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


void ksemawc::rangeWL(){
    double Lmin=ui->dSB_PAR_4_1 -> value();
    double Lmax=ui->dSB_PAR_4_2 -> value();
    ui->lineEdit_eVrange->setText(QString::number(1240./Lmin,'g',4)+"<->"+QString::number(1240./Lmax,'g',4)+" eV");
    Qt::CheckState state;
    state = ui->checkBox_msgr -> checkState();
    if(state == Qt::Unchecked){
        rxy[20][1]=Lmin;
        rxy[20][2]=Lmax;
        ui->DP_RXY_20_1 ->setText(QString::number(Lmin));
        ui->DP_RXY_20_2 ->setText(QString::number(Lmax));
    }
    ui->dSB_WLsearchNK->setMinimum(Lmin);
    ui->dSB_WLsearchNK->setMaximum(Lmax);
}


void ksemawc::saveNKsim(){
    printf("->saveNKsim with lastAction=%s\n",lastAction.toStdString().c_str());
    int ierr,iok=1;
    QString command,sample,fn2s;
    sample=ui->lineEdit_sample-> text();
    if(sample.contains("mate/aa999.9"))
        sample="";
    sample=sample+"_"+lastAction+".nk";
    printf("sample= %s\n",sample.toStdString().c_str());
    fn2s = QFileDialog::getSaveFileName(
                this,
                "Filename to save NKsim",
                pathroot+"/"+sample,
                "nk file (*.nk)");
    QFile file(fn2s);
    if( file.exists() ){
        QMessageBox msgBox;
        msgBox.setText("The file already exist!");
        msgBox.setInformativeText("Do you want to save anyway?");
        msgBox.setStandardButtons(QMessageBox::Save | QMessageBox::No);
        msgBox.setDefaultButton(QMessageBox::Save);
        int ret = msgBox.exec();
        switch (ret) {
        case QMessageBox::Save:
            // Save was clicked
            iok=1;
            break;
        case QMessageBox::No:
            // Cancel was clicked
            iok=0;
            printf("as you requested file was not saved!\n");
            break;
        }
    }
    if(iok==1){
        cout<< "fNKsim= "<<fNKsim.toStdString()<<endl;
        cout<< "fn2s= "<<fn2s.toStdString()<<endl;
        //ierr=system((command.toStdString()).c_str());
        if (IS_POSIX == 1) {
            command="cp "+fNKsim+" "+fn2s;
            ierr=system((command.toStdString()).c_str());}
        else
        {
            command="copy "+fNKsim+" "+fn2s;
            QStringList List;
            List =command.split("/");
            int nV=List.count();
            //printf("\nnV=%d\n",nV);
            QString commandw;
            QString pezzo;
            for(int iv=0;iv<nV;iv++){
                pezzo=List.at(iv).toLocal8Bit().constData();
                //printf("pezzo[%d]= %s\t",iv,pezzo.toStdString().c_str());
                if(iv==0){
                    commandw=pezzo;
                }
                else{
                    commandw=commandw+"\\"+pezzo;
                }
            }
            //printf("comando=%s\n",(commandw.toStdString()).c_str());
            ierr=system((commandw.toStdString()).c_str());
        }
        if(ierr != 0)
            printf("Error copying NKsim!!!\n");
        else
            printf("NKsim saved as %s\n",(fn2s.toStdString()).c_str());
    }
}


void ksemawc::saveSim(){
    int ierr,iok=1;
    QString command, fn2s;
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
        msgBox.setStandardButtons(QMessageBox::Save | QMessageBox::No);
        msgBox.setDefaultButton(QMessageBox::Save);
        int ret = msgBox.exec();
        switch (ret) {
        case QMessageBox::Save:
            // Save was clicked
            iok=1;
            break;
        case QMessageBox::No:
            // Cancel was clicked
            iok=0;
            printf("as you requested file was not saved!\n");
            break;
        }
    }
    if(iok==1){
        cout<< "fMisSim= "<<fMisSim.toStdString()<<endl;
        cout<< "fn2s= "<<fn2s.toStdString()<<endl;
        //ierr=system((command.toStdString()).c_str());
        if (IS_POSIX == 1) {
            command="cp "+fMisSim+" "+fn2s;
            ierr=system((command.toStdString()).c_str());}
        else
        {
            command="copy "+fMisSim+" "+fn2s;
            QStringList List;
            List =command.split("/");
            int nV=List.count();
            //printf("\nnV=%d\n",nV);
            QString commandw;
            QString pezzo;
            for(int iv=0;iv<nV;iv++){
                pezzo=List.at(iv).toLocal8Bit().constData();
                //printf("pezzo[%d]= %s\t",iv,pezzo.toStdString().c_str());
                if(iv==0){
                    commandw=pezzo;
                }
                else{
                    commandw=commandw+"\\"+pezzo;
                }
            }
            //printf("comando=%s\n",(commandw.toStdString()).c_str());
            ierr=system((commandw.toStdString()).c_str());
        }
        if(ierr != 0){
            printf("Error copying MisSim.dat!!!\n");
            msgErrLoad("saveSim","Error copying MisSim.dat!!!");
        }
        else
            printf("MisSim saved as %s\n",(fn2s.toStdString()).c_str());
    }
}

void ksemawc::tabChanged(){
    QMessageBox msgBox;
    int itab=ui->tabWidget -> currentIndex();
    printf("-> tabChanged:  itab=%d lastTab=%d\n",itab,lastTab);
    // QString filepro=ui->lineEdit_P -> text();
    // if(filepro.contains("mate/aa999") && itab!=5)
    //     SaveProject();
    readRangeEli();
    AdjRoughMax();
    SaveSetting(lastTab);
    if(lastTab==0 && itab!=0){
        if(!fnSample.isEmpty()){
            int NcBc=0;
            Qt::CheckState state;
            for(int imis=1;imis<=14;imis++){
                state=idToCheckBox["checkB_mis"+QString::number(imis)+"_1"]->checkState();
                if(state == Qt::Checked)
                    NcBc++;
                if(imis>=7)
                    imis++;
            }
            if(NcBc==0){
                QMessageBox msgBox;
                msgBox.setText("None experimental measurements is enabled to be load!");
                msgBox.setInformativeText("Do you want to enable some now?");
                msgBox.setStandardButtons(QMessageBox::Yes | QMessageBox::No);
                msgBox.setDefaultButton(QMessageBox::Yes);
                int ret = msgBox.exec();
                switch (ret) {
                case QMessageBox::Yes:
                    ui->tabWidget -> setCurrentIndex(0);
                    return;
                    break;
                case QMessageBox::No:
                    printf("None experimental measurements is loaded!\n");
                    break;
                }
            }
        }
        occupyPF=1;
        SPADA();//load measurements, file-nk and solution-nk
        occupyPF=0;
    }
    if(lastTab==3)//Analysis
        setTabSim();
    if(itab==1)//model
        SetModel(nlayer);
    if((itab==1 || itab==2) && lastTab==4)//model & simula
        PlotMENK();
    if(itab==3){ //data analysis
        refreshFitPar();
        if(ifirstcall==0 && jobtot==0)
            StoreFitSet();
        if(lastTab==4)
            RefTrackG();
        printf("ifirstcall=%d\n",ifirstcall);
        if(ifirstcall==0){
            ui->lineEdit_chi2Best->setText("1.E+99");
            ui->spinBox_jobBest->setValue(0);
        }
    }
    if(itab==4){ //range
        setRangeEli();
        ui->sB_iwl2print->setMaximum(NeV);
    }
    if(lastTab==4)
        iwl2print=ui->sB_iwl2print->value();

    lastTab=itab;
}


void ksemawc::setTabSim(){
    printf("-> setTabSim\n");
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
                for(int ii=2;ii<6;ii++){
                    QString content=idToLineEdit["LEpm_"+QString::number(100+ii+(i-1)*5)+"_1"] -> text();
                    if(content!="unused")
                        pm[100+ii+(i-1)*5][1]=content.toDouble();
                }
                if(pm[100+1+(i-1)*5][1]==15 || pm[100+1+(i-1)*5][1]==16){
                    if(pm[100+4+(i-1)*5][1]>pm[100+5+(i-1)*5][1]*0.45){
                        pm[100+4+(i-1)*5][1]=pm[100+5+(i-1)*5][1]*0.45;
                        idToLineEdit["LEpm_"+QString::number(100+4+(i-1)*5)+"_1"] -> setText(QString::number(pm[100+4+(i-1)*5][1]));
                    }
                }
                j++;
                //printf("pm[%d][1]: %f %f %f %f %f\n",100+2+(i-1)*5,pm[100+1+(i-1)*5][1],
                //        pm[100+2+(i-1)*5][1],pm[100+3+(i-1)*5][1],
                //        pm[100+4+(i-1)*5][1],pm[100+5+(i-1)*5][1]);
            }
        }
        pm[71][1]=ui->dSB_PM_1_1->value();
        pf[ioptFit][1]=j-2;
    }
}


void ksemawc::refreshFitPar(){
    int icoherent,klim,ivp,n_fresh,iDecine,ilayer,iosc,Jcombo,ip;
    occupyPF=1;
    printf("-> refreshFitPar\n");
    fflush(stdout);
    npp=ui->sB_PAR_34_5 -> value();
    n_fresh=0;
    QString Lj,Lparametro;
    for(int j=1;j<=17;j++){
        idToComboBox["cBParFit_"+QString::number(j)] -> clear();
        idToComboBox["cBParFit_"+QString::number(j)] -> addItem("none");
    }
    for(int i=1;i<=nlayer;i++){
        icoherent=idToComboBox["comB_PAR_5"+QString::number(i)+"_3"] -> currentIndex();
        Lj="L"+QString::number(i);
        if(icoherent==0)
            klim=-1;
        else if(icoherent==1)
            klim=2;
        else if(icoherent==2)
            klim=7;
        else if(icoherent == 3 || icoherent == 4)
            klim=2;
        for(int j=1;j<=npp;j++){
            idToComboBox["cBParFit_"+QString::number(j)]-> addItem(Lj+ParFitLab[6]);//roughness is always included
            for(int k=0;k<=klim;k++){
                if(k!=6)
                    idToComboBox["cBParFit_"+QString::number(j)]-> addItem(Lj+ParFitLab[k]);
            }
        }
    }
    for(int i=1;i<=15;i++){
        Qt::CheckState state=idToCheckBox["cB_EMA_"+QString::number(i)]->checkState();
        if(state==Qt::Checked){
            printf("\tstate==Qt::Checked for i=%d -> add fEMA_%d\n",i,i);
            for(int j=1;j<=npp;j++)
                idToComboBox["cBParFit_"+QString::number(j)]-> addItem(ParFitLab[12]+QString::number(i));
        }
    }
    for(int i=1;i<=4;i++){
        Qt::CheckState state=idToCheckBox["checkB_mis"+QString::number(5+2*i)+"_1"]->checkState();
        if(state==Qt::Checked){
            printf("\tstate==Qt::Checked for i=%d -> add Theta_%d\n",i,i);
            for(int j=1;j<=npp;j++)
                idToComboBox["cBParFit_"+QString::number(j)]-> addItem(ParFitLab[13]+QString::number(i));
        }
    }
    for(int i=1;i<=20;i++){
        Lj="O"+QString::number(i);
        Qt::CheckState state=idToCheckBox["cBosc_"+QString::number(i)]-> checkState();
        if( state == Qt::Checked ) {
            int ifu=idToComboBox["cBpm_"+QString::number(100+1+(i-1)*5)+"_1"]->currentIndex();
            ifu++;
            for(int j=1;j<=npp;j++){
                idToComboBox["cBParFit_"+QString::number(j)]-> addItem(Lj+ParFitLab[8]);
                if(ifu!=4 && ifu!=5 && ifu!=17)
                    idToComboBox["cBParFit_"+QString::number(j)]-> addItem(Lj+ParFitLab[9]);
                if(ifu!=4){
                    idToComboBox["cBParFit_"+QString::number(j)]-> addItem(Lj+ParFitLab[10]);
                }
                if(ifu>5 && ifu!=17)
                    idToComboBox["cBParFit_"+QString::number(j)]-> addItem(Lj+ParFitLab[11]);
            }
        }
    }
    for(int j=1;j<=17;j++){
        if(j<=npp){
            //if(lastTab != itab){
            ivp=NINT(pm[ppm[j]][2]);
            ip=ppm[j];
            if(ip<=69 || (ip>=91 && ip<=99)){
                iDecine=int(double(ip)/10.);
                ilayer=ip-iDecine*10;
                if(iDecine==9)
                    iDecine=-1;
                Lparametro="L"+QString::number(ilayer)+ParFitLab[iDecine+1];
                //      printf("j=%d ip=%d iDecine=%d ilayer=%d Lparametro=%s\n",j,ip,iDecine,ilayer,(Lparametro.toStdString()).c_str());
            }
            else if(ip>=71 && ip<=85){//EMA
                Lparametro=ParFitLab[12]+QString::number(ip-70);
            }
            else if(ip>=86 && ip<=89){//ThetaEli
                Lparametro=ParFitLab[13]+QString::number(ip-85);
            }
            else{
                iosc=int(double(ip-101)/5.)+1;
                Lparametro="O"+QString::number(iosc)+ParFitLab[7+ip-101-5*(iosc-1)];
                //      printf("j=%d ip=%d iosc=%d Lparametro=%s\n",j,ip,iosc,(Lparametro.toStdString()).c_str());
            }
            printf("j=%d ivp=%d ip=%d Lparametro=%s\n",j,ivp,ip,Lparametro.toStdString().c_str());
            Jcombo=idToComboBox["cBParFit_"+QString::number(j)] -> findText(Lparametro);
            if(Jcombo>=0) idToComboBox["cBParFit_"+QString::number(j)] -> setCurrentIndex(Jcombo);
            idToLineEdit["DPparFitV_"+QString::number(j)] -> setText(QString::number(pm[ip][1]));
            idToCheckBox["chBeParFit_"+QString::number(j)] -> setEnabled(true);
            if(ivp==0){
                idToCheckBox["chBeParFit_"+QString::number(j)] -> setCheckState ( Qt::Unchecked );
                idToLineEdit["DPparFitErr_"+QString::number(j)] -> setEnabled(false);
                idToLineEdit["DPparFitGC_"+QString::number(j)] -> setEnabled(false);
            }
            else{
                idToCheckBox["chBeParFit_"+QString::number(j)] -> setCheckState ( Qt::Checked );
                idToLineEdit["DPparFitErr_"+QString::number(j)] -> setEnabled(true);
                idToLineEdit["DPparFitGC_"+QString::number(j)] -> setEnabled(true);
                n_fresh++;
                idToLineEdit["DPparFitErr_"+QString::number(j)] -> setText(QString::number(pm[ip][4]));
                idToLineEdit["DPparFitGC_"+QString::number(j)] -> setText(QString::number(pm[ip][5]));
            }
            printf("j=% ip=%d pm[%d][1]=%f\n",j,ip,ip,pm[ip][1]);
        }
        else{
            idToCheckBox["chBeParFit_"+QString::number(j)] -> setCheckState ( Qt::Unchecked );
            idToCheckBox["chBeParFit_"+QString::number(j)] -> setEnabled(false);
        }
    }
    nPar=n_fresh;
    ui->sB_PAR_35_5 -> setValue(nPar);
    occupyPF=0;
    SaveSetting(-1);
}



void ksemawc::PanFitEnable(){
    int n_fresh,n1,n2,ip;
    int itab=ui->tabWidget -> currentIndex();

    if(occupyPF==0){
        printf("-> PanFitEnable with itab=%d\n",itab);
        fflush(stdout);
        occupyPF=1;
        npp=ui->sB_PAR_34_5 -> value();
        n_fresh=0;
        Qt::CheckState state;
        //state = Qt::Checked;
        QString Lj,Valore;
        if(itab==3 && itab==lastTab){
            printf("\t apply PanFitEnable with npp=%d\n",npp);
            fflush(stdout);
            for(int j=1;j<=17;j++){
                if(j<=npp){
                    idToLineEdit["DPparFitV_"+QString::number(j)] -> setEnabled(true);
                    idToCheckBox["chBeParFit_"+QString::number(j)] -> setEnabled(true);
                    state=idToCheckBox["chBeParFit_"+QString::number(j)]-> checkState();
                    if( state == Qt::Checked ){
                        n_fresh++;
                        idToLineEdit["DPparFitErr_"+QString::number(j)] -> setEnabled(true);
                        idToLineEdit["DPparFitGC_"+QString::number(j)] -> setEnabled(true);
                    }
                    else{
                        idToLineEdit["DPparFitErr_"+QString::number(j)] -> setEnabled(false);
                        idToLineEdit["DPparFitGC_"+QString::number(j)] -> setEnabled(false);
                    }
                    Lj=idToComboBox["cBParFit_"+QString::number(j)] -> currentText();
                    ip=0;
                    if(Lj.at(0)=='L'){
                        n1=(Lj.at(1)).digitValue();
                        n2=-1;
                        do{
                            n2++;
                        } while(!Lj.contains(ParFitLab[n2],Qt::CaseSensitive) && n2<7);
                        if(n2==0)
                            ip=90+n1;
                        else
                            ip=10*(n2-1)+n1;
                    }
                    else if(Lj.at(0)=='f'){
                        ip=70+Lj.at(5).digitValue();
                    }
                    else if(Lj.at(0)=='T')
                        ip=85+Lj.at(6).digitValue();
                    else{
                        n1=(Lj.at(1)).digitValue();
                        if(Lj.at(2).isNumber()) {
                            n1=n1*10;
                            n1=n1+((Lj.at(2)).digitValue());
                        }
                        n2=6;
                        do{
                            n2++;
                        } while(!Lj.contains(ParFitLab[n2],Qt::CaseSensitive) && n2<=12);
                        ip=100+(n1-1)*5+(n2-7)+1;
                    }
                    ppm[j]=ip;
                    Valore=idToLineEdit["DPparFitV_"+QString::number(j)] -> text();
                    pm[ip][1]=Valore.toDouble();
                    if(ip>100)
                        pm[ip][1]=abs(pm[ip][1]);
                    if( state == Qt::Unchecked )
                        pm[ip][2]=0;
                    else{
                        pm[ip][2]=j;
                        pm[n_fresh][3]=ip;
                    }
                    //     printf("j=%d ip=%d pm[%d][1]=%f pm[%d][2]=%f pm[%d][3]=%f\n",j,ip,ip,pm[ip][1],ip,pm[ip][2],n_fresh,pm[n_fresh][3]);
                }
                else{
                    idToCheckBox["chBeParFit_"+QString::number(j)] -> setCheckState ( Qt::Unchecked );
                    idToCheckBox["chBeParFit_"+QString::number(j)] -> setEnabled(false);
                    idToLineEdit["DPparFitV_"+QString::number(j)] -> setEnabled(false);
                    idToLineEdit["DPparFitErr_"+QString::number(j)] -> setEnabled(false);
                    idToLineEdit["DPparFitGC_"+QString::number(j)] -> setEnabled(false);
                }
            }
            nPar=n_fresh;
            ui->sB_PAR_35_5 -> setValue(nPar);
//            printf("END PanFitEnable\n");
        }
        occupyPF=0;
    }
}


void ksemawc::PanFitPar(){
    if(occupyPF!=0) return;
    printf("-> PanFitPar\n");
    int icoherent,klim,nppNew;
    int itab=ui->tabWidget -> currentIndex();
    if(occupyPF==0){
        occupyPF=1;
        nppNew=ui->sB_PAR_34_5 -> value();
        Qt::CheckState state;
        //state = Qt::Checked;
        QString Lj;
        if(itab==3 && itab==lastTab){
            printf("\tStart PanFitPar: nppNew=%d npp=%d\n",nppNew,npp);
            if(nppNew > npp){
                for(int j=npp+1;j<=nppNew;j++){//enable and initialize new parameters
                    idToCheckBox["chBeParFit_"+QString::number(j)] -> setEnabled(true);
                    idToComboBox["cBParFit_"+QString::number(j)] -> clear();
                    idToComboBox["cBParFit_"+QString::number(j)] -> addItem("none");
                    idToLineEdit["DPparFitV_"+QString::number(j)]-> clear();
                    idToLineEdit["DPparFitErr_"+QString::number(j)]-> clear();
                    idToLineEdit["DPparFitGC_"+QString::number(j)]-> clear();
                    idToLineEdit["DPparFitV_"+QString::number(j)] -> setEnabled(true);
                }
                for(int i=1;i<=nlayer;i++){//add avilable layer parameters
                    icoherent=idToComboBox["comB_PAR_5"+QString::number(i)+"_3"] -> currentIndex();
                    Lj="L"+QString::number(i);
                    if(icoherent==0) //bulk
                        klim=-1;
                    else if(icoherent==1) //film homo
                        klim=2;
                    else if(icoherent==2) //film inhomo
                        klim=7;
                    for(int j=npp+1;j<=nppNew;j++){
                        idToComboBox["cBParFit_"+QString::number(j)]-> addItem(Lj+ParFitLab[6]);//roughness is always included
                        for(int k=0;k<=klim;k++){
                            if(k!=6)
                                idToComboBox["cBParFit_"+QString::number(j)]-> addItem(Lj+ParFitLab[k]);
                        }
                    }
                }
                for(int i=1;i<=15;i++){//add available EMAs
                    state=idToCheckBox["cB_EMA_"+QString::number(i)]->checkState();
                    if(state==Qt::Checked){
                        for(int j=npp+1;j<=nppNew;j++)
                            idToComboBox["cBParFit_"+QString::number(j)]-> addItem(ParFitLab[12]+QString::number(i));
                    }
                }
                for(int i=1;i<=4;i++){//add available THETAs
                    Qt::CheckState state=idToCheckBox["checkB_mis"+QString::number(5+2*i)+"_1"]->checkState();
                    if(state==Qt::Checked){
                        for(int j=npp+1;j<=nppNew;j++)
                            idToComboBox["cBParFit_"+QString::number(j)]-> addItem(ParFitLab[13]+QString::number(i));
                    }
                }
                for(int i=1;i<=20;i++){// add available oscillators
                    Lj="O"+QString::number(i);
                    state=idToCheckBox["cBosc_"+QString::number(i)]-> checkState();
                    if( state == Qt::Checked ) {
                        int ifu=idToComboBox["cBpm_"+QString::number(100+1+(i-1)*5)+"_1"]->currentIndex();
                        ifu++;
                        for(int j=npp+1;j<=nppNew;j++){
                            idToComboBox["cBParFit_"+QString::number(j)]-> addItem(Lj+ParFitLab[8]);     //C
                            if(ifu!=4 && ifu!=5 && ifu!=17)
                                idToComboBox["cBParFit_"+QString::number(j)]-> addItem(Lj+ParFitLab[9]); //E
                            if(ifu!=4){
                                idToComboBox["cBParFit_"+QString::number(j)]-> addItem(Lj+ParFitLab[10]);//D
                            }
                            if(ifu>5 && ifu!=17)
                                idToComboBox["cBParFit_"+QString::number(j)]-> addItem(Lj+ParFitLab[11]);//W
                        }
                    }
                }
            }
            else{//disable redundand parameters
                for(int j=nppNew+1;j<=npp;j++){
                    idToCheckBox["chBeParFit_"+QString::number(j)] -> setEnabled(false);
                    idToComboBox["cBParFit_"+QString::number(j)] -> clear();
                    idToComboBox["cBParFit_"+QString::number(j)] -> addItem("none");
                    idToCheckBox["chBeParFit_"+QString::number(j)] -> setCheckState ( Qt::Unchecked );
                    idToLineEdit["DPparFitV_"+QString::number(j)] -> clear();//setText("");
                    idToLineEdit["DPparFitErr_"+QString::number(j)] -> clear();
                    idToLineEdit["DPparFitGC_"+QString::number(j)] -> clear();
                    idToLineEdit["DPparFitV_"+QString::number(j)] -> setEnabled(false);
                    idToLineEdit["DPparFitErr_"+QString::number(j)] -> setEnabled(false);
                    idToLineEdit["DPparFitGC_"+QString::number(j)] -> setEnabled(false);
                }
            }
            //   printf("End PanFitPar\n");
            npp=nppNew;
            PanFitEnable();
        }
        occupyPF=0;
    }
}


void ksemawc::PanFitChoice(){
    int n1,n2,ip;
    int itab=ui->tabWidget -> currentIndex();
    QString Lj;
    Qt::CheckState state;
    if(occupyPF==0){
        printf("-> PanFitChoice\n");
        occupyPF=1;
        if(itab==3 && itab==lastTab){
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
                    n2=-1;
                    ip=0;
                    do{
                        n2++;
                    } while(!Lj.contains(ParFitLab[n2],Qt::CaseSensitive) && n2<=12);
                    if(n2==0)
                        ip=90+n1;
                    else if(1<=n2 && n2<=7)
                        ip=10*(n2-1)+n1;
                    else if(8<=n2 && n2<=11)
                        ip=100+(n1-1)*5+(n2-7)+1;
                    else if(n2==12)
                        ip=70+Lj.at(5).digitValue();
                    else if(n2==13)
                        ip=85+Lj.at(6).digitValue();
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

void ksemawc::PlotMENK(){
    printf("-> PlotMENK()\n");
    iColor=0;//reset to black color
    SPADA();
    PlotME();
    PlotNK(1);
}



void ksemawc::PlotME(){
    int Ndata=NeV;
    int iRD,iPSI=0,iDELTA=0;
    double Xp[Ndata],Yp[Ndata],ErrXp[Ndata],ErrYp[Ndata];
    SaveSetting(-1);
    for(int i=1;i<=14;i++){
        int iGraph=i;
        iRD=1;
        if(i==7 || i==9 || i==11 || i==13){
            iGraph=7;
            if(iDELTA>0)
                iRD=0;
            if(DATO[i]>0)
                iDELTA++;
        }
        if(i==8 || i==10 || i==12 || i==14){
            iGraph=8;
            if(iPSI>0)
                iRD=0;
            if(DATO[i]>0)
                iPSI++;
        }
        if(DATO[i]!=0)
            printf("->PlotMe: DATO[%d]=%d ifirstWarning=%d\n",i,DATO[i],ifirstWarning);
        if(DATO[i]>0){
            double YOLD,X,Y,errY;
            if(i>=7) YOLD=ELI[i-6][1][1];
            for(int L=1;L<=Ndata;L++){
                X=MIS[7][L][1];
                if(i<=6){
                    Y=MIS[i][L][1]*100.;
                    errY=MIS[i][L][2]*100.;
                }
                else if(i>=7){
                    Y=ELI[i-6][L][1];
                    errY=ELI[i-6][L][2];
                    if(iDelCon==1){
                        if(i==7||i==9||i==11||i==13){
                            if(L==1)
                                YOLD=Y;
                            Y=Y-Dperiod*NINT((Y-YOLD)/Dperiod);
                            YOLD=Y;
                        }
                    }
                }
                Xp[L-1]=X;
                Yp[L-1]=Y;
                ErrYp[L-1]=errY;
            }
            iColor=0;
            PLOTline1bar2(2,iRD,0,iGraph,Ndata,Xp,Yp,ErrXp,ErrYp);
        }
        else if(DATO[i]!=0)
            PLOTline1bar2(2,iRD,0,iGraph,0,Xp,Yp,ErrXp,ErrYp);//erase without data plot
    }

    iRD=1;
    if(DATO[1]>0 && DATO[3]>0){
        int iPlotA=0;
        rxy[18][3]=1000.;
        rxy[18][4]=-1000.;
        for(int L=1;L<=Ndata;L++){
            Xp[L-1]=MIS[7][L][1];
            Yp[L-1]=(1.-MIS[1][L][1]-MIS[3][L][1])*100.;//A=1-T-R
            ErrYp[L-1]=(MIS[1][L][2]+MIS[3][L][2])*100.;//error
            if(Yp[L-1]<0.) iPlotA++;
            rxy[18][3]=min(rxy[18][3],Yp[L-1]);
            rxy[18][4]=max(rxy[18][4],Yp[L-1]);
        }
        ui->DP_RXY_18_3 -> setText(QString::number(rxy[18][3]));
        ui->DP_RXY_18_4 -> setText(QString::number(rxy[18][4]));
        labelQwt="A_front";
        PLOTline1bar2(2,iRD,0,9,Ndata,Xp,Yp,ErrXp,ErrYp);
        if(ifirstWarning==0 && iPlotA>0){
            QMessageBox msgBox;
            msgBox.setText("ATTENTION: A = 1 - Tn - Rn < 0 !!!\nPlease check your measurements.");
            msgBox.setStandardButtons(QMessageBox::Ok);
            msgBox.exec();
            ifirstWarning=-1;
        }
    }
    else if(DATO[1]!=0 && DATO[3]!=0)
        PLOTline1bar2(2,iRD,0,9,0,Xp,Yp,ErrXp,ErrYp);

    if(DATO[1]>0 && DATO[5]>0){
        int iPlotA=0;
        rxy[15][3]=1000.;
        rxy[15][4]=-1000.;
        for(int L=1;L<=Ndata;L++){
            Xp[L-1]=MIS[7][L][1];
            Yp[L-1]=(1.-MIS[1][L][1]-MIS[5][L][1])*100.;//A=1-T-R1
            ErrYp[L-1]=(MIS[1][L][2]+MIS[5][L][2])*100.;//error
            if(Yp[L-1]<0.) iPlotA++;
            rxy[15][3]=min(rxy[15][3],Yp[L-1]);
            rxy[15][4]=max(rxy[15][4],Yp[L-1]);
        }
        ui->DP_RXY_15_3 -> setText(QString::number(rxy[15][3]));
        ui->DP_RXY_15_4 -> setText(QString::number(rxy[15][4]));
        labelQwt="A_back";
        PLOTline1bar2(2,iRD,0,16,Ndata,Xp,Yp,ErrXp,ErrYp);
        if(ifirstWarning<=0 && iPlotA>0){
            QMessageBox msgBox;
            msgBox.setText("ATTENTION: A= 1 - Tn - R1 <0 !!!\nPlease check your measurements.");
            msgBox.setStandardButtons(QMessageBox::Ok);
            msgBox.exec();
        }
        ifirstWarning=1;
    }
    else if(DATO[1]!=0 && DATO[5]!=0)
        PLOTline1bar2(2,iRD,0,16,0,Xp,Yp,ErrXp,ErrYp);


    // save to HD splined experimental measurements
    QFile file0(pathroot+"expo/MisSFexp.dat");
    if (!file0.open(QIODevice::WriteOnly | QIODevice::Text)){
        msgErrLoad("PlotMe file0","Error opening expo/MisSFexp.dat");
        printf("Error opening expo/MisSFexp.dat\n");
    }
    QTextStream stream0 (&file0);

    QFile file1(pathroot+"expo/MisSFexpErr.dat");
    if (!file1.open(QIODevice::WriteOnly | QIODevice::Text)){
        msgErrLoad("PlotMe file1","Error opening expo/MisSFexpErr.dat");
        printf("Error opening expo/MisSFexpErr.dat\n");
    }
    QTextStream stream1 (&file1);

    QFile file2(pathroot+"expo/MisELIexp.dat");
    if (!file2.open(QIODevice::WriteOnly | QIODevice::Text)){
        msgErrLoad("PlotMe file2","Error opening expo/MisELIexp.dat");
        printf("Error opening expo/MisELIexp.dat\n");
    }
    QTextStream stream2 (&file2);

    QFile file3(pathroot+"expo/MisELIexpErr.dat");
    if (!file3.open(QIODevice::WriteOnly | QIODevice::Text)){
        msgErrLoad("PlotMe file3","Error opening expo/MisELIErr.dat");
        printf("Error opening expo/MisELIErr.dat\n");
    }
    QTextStream stream3 (&file3);

    stream0<<"Tn\tTp\tRn\tRp\tR1\tAPDS\twl(nm)\tn1\tn2\tn3\tn4\tn5\tn6\tn7\tn8\tn_ibri"<<"\n";
    stream1<<"eTn\teTp\teRn\teRp\teR1\teAPDS\teLETT\tk1\tk2\tk3\tk4\tk5\tk6\tk7\tk8\tk_ibri"<<"\n";
    stream2<<"DEL_1\tPSI_1\tDEL_2\tPSI_2\tDEL_3\tPSI_3\tDEL_4\tPSI_4"<<"\n";
    stream3<<"eDEL_1\tePSI_1\teDEL_2\tePSI_2\teDEL_3\tePSI_3\teDEL_4\tePSI_4"<<"\n";
    for(int i=1;i<=NeV;i++){
        for(int j=1;j<=16;j++) stream0<<QString::number(MIS[j][i][1],'f',7)<<"\t";
        for(int j=1;j<=16;j++) stream1<<QString::number(MIS[j][i][2],'f',7)<<"\t";
        for(int j=1;j<=8;j++)  stream2<<QString::number(ELI[j][i][1],'f',7)<<"\t";
        for(int j=1;j<=8;j++)  stream3<<QString::number(ELI[j][i][2],'f',7)<<"\t";
        stream0<<"\n";
        stream1<<"\n";
        stream2<<"\n";
        stream3<<"\n";
    }
    file0.close();
    file1.close();
    file2.close();
    file3.close();
}


void ksemawc::PlotNK(int iRD){
    int Ndata=NINT(SOL[1][1]);
    //if(Ndata==0)
    //    return;
    double Xp[Ndata],Yp[Ndata],ErrXp[Ndata],ErrYp[Ndata];
    // n-data plot
    for(int i=0;i<Ndata;i++){
        Xp[i]=SOL[i+2][1];
        Yp[i]=SOL[i+2][2];
        ErrYp[i]=SOL[i+2][4];
    }
    PLOTline1bar2(2,iRD,iColor,12,Ndata,Xp,Yp,ErrXp,ErrYp);
    // k-data plot
    for(int i=0;i<Ndata;i++){
        Xp[i]=SOL[i+2][1];
        Yp[i]=SOL[i+2][3];
        ErrYp[i]=SOL[i+2][5];
    }
    PLOTline1bar2(2,iRD,iColor,13,Ndata,Xp,Yp,ErrXp,ErrYp);
    //eps1 eps2 plot
    if(NINT(par[10][1])==1){
        //eps1-data plot
        double v1,v2,v3,v4,yd,yu;
        rxy[26][3]=1000.;
        rxy[26][4]=-1000.;
        for(int i=2;i<=Ndata+1;i++){
            v1=pow(SOL[i][1+1]-SOL[i][1+3],2.)-pow(SOL[i][2+1]-SOL[i][2+3],2.);
            v2=pow(SOL[i][1+1]-SOL[i][1+3],2.)-pow(SOL[i][2+1]+SOL[i][2+3],2.);
            v3=pow(SOL[i][1+1]+SOL[i][1+3],2.)-pow(SOL[i][2+1]-SOL[i][2+3],2.);
            v4=pow(SOL[i][1+1]+SOL[i][1+3],2.)-pow(SOL[i][2+1]+SOL[i][2+3],2.);
            yd=min(v1,v2);
            yd=min(yd,v3);
            yd=min(yd,v4);
            yu=max(v1,v2);
            yu=max(yu,v3);
            yu=max(yu,v4);
            rxy[26][3]=min(rxy[26][3],yd);
            rxy[26][4]=max(rxy[26][4],yu);
            Xp[i-2]=SOL[i][1];
            Yp[i-2]=(yu+yd)/2.;
            ErrYp[i-2]=(yu-yd)/2.;
        }
        PLOTline1bar2(2,iRD,iColor,14,Ndata,Xp,Yp,ErrXp,ErrYp);
        ui->DP_RXY_26_3->setText(QString::number(rxy[26][3]));
        ui->DP_RXY_26_4->setText(QString::number(rxy[26][4]));
        //eps2-data plot
        rxy[27][3]=1000.;
        rxy[27][4]=-1000.;
        for(int i=2;i<=Ndata+1;i++){
            v1=2.*(SOL[i][1+1]-SOL[i][1+3])*(SOL[i][2+1]-SOL[i][2+3]);
            v2=2.*(SOL[i][1+1]-SOL[i][1+3])*(SOL[i][2+1]+SOL[i][2+3]);
            v3=2.*(SOL[i][1+1]+SOL[i][1+3])*(SOL[i][2+1]-SOL[i][2+3]);
            v4=2.*(SOL[i][1+1]+SOL[i][1+3])*(SOL[i][2+1]+SOL[i][2+3]);
            yd=min(v1,v2);
            yd=min(yd,v3);
            yd=min(yd,v4);
            yu=max(v1,v2);
            yu=max(yu,v3);
            yu=max(yu,v4);
            rxy[27][3]=min(rxy[27][3],yd);
            rxy[27][4]=max(rxy[27][4],yu);
            Xp[i-2]=SOL[i][1];
            Yp[i-2]=(yu+yd)/2.;
            ErrYp[i-2]=(yu-yd)/2.;
        }
        PLOTline1bar2(2,iRD,iColor,15,Ndata,Xp,Yp,ErrXp,ErrYp);
        ui->DP_RXY_27_3->setText(QString::number(rxy[27][3]));
        ui->DP_RXY_27_4->setText(QString::number(rxy[27][4]));
    }
}


void ksemawc::Simula(){
    printf("->Simula\n");
    lastAction="Simula";
    AdjRoughMax();
    SaveSetting(-1);
    double mc[15][1000];
    Qt::CheckState state,state1;
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
    for(int i=1;i<=14;i++){
        state1=idToCheckBox["checkB_mis"+QString::number(i)+"_2"] -> checkState ();
        if(state1==Qt::Checked && DATO[i]==0){
            DATO[i]=-1;
        }
    }
    int expA=ui->spinBox_expAENS->value();//attenuation coefficient of k by Fit#N
    double valA=ui->doubleSpinBox_valAENS->value();
    par[11][1]=valA*pow(10.,expA);
    CalcMis(mc);
    if(state == Qt::Checked && iok==1){
        int i1=1;
        int i2=14;
        double INSTR[21];
        for(int j=1;j<=20;j++)
            INSTR[j]=par[j][3];
        for(int ic=i1;ic<=i2;ic++){
            if(DATO[ic]!=0){
                DATO[ic]=1;
                par[ic][4]=1.;
                idToCheckBox["checkB_mis"+QString::number(ic)+"_3"] -> setEnabled(true);
                printf("%s is enabled\n",("checkB_mis"+QString::number(ic)+"_3").toStdString().c_str());
                for(int i=1;i<=NeV;i++){
                    if(ic<=6){
                        MIS[ic][i][1]=mc[ic][i];
                        double DL=MIS[7][i][2];
                        if(ic==1 || ic==2)
                            MIS[ic][i][2]=MIS[ic][i][1]*INSTR[3]+DL;
                        if(ic==3 || ic==4 || ic==5)
                            MIS[ic][i][2]=MIS[ic][i][1]*(INSTR[3]+INSTR[4])+DL;
                        if(iw==1)
                            printf("MIS[%d][%d][1]= %f\n",ic,i,MIS[ic][i][1]);
                    }
                    else if(ic>=7 && ic<=14)
                        ELI[ic-6][i][1]=mc[ic][i];
                }
            }
        }
    }
    for(int ic=1;ic<=14;ic++){
        state1=idToCheckBox["checkB_mis"+QString::number(ic)+"_2"] -> checkState ();
        if(state1==Qt::Checked && DATO[ic]>0)
            idToLineEdit["LEpar_"+QString::number(35+ic)+"_1"] -> setText(QString::number(par[35+ic][1]));
        else
            idToLineEdit["LEpar_"+QString::number(35+ic)+"_1"] -> clear();
    }
    ui->cBox_PAR_19_1 -> setCheckState(Qt::Unchecked);

    // **** computing mean experimental and simulated SF
    //    load weight of selected mean
    QString fnmean=pathroot+NANK[10].simplified();
    QFile file(fnmean);
    if(!file.open(QIODevice::ReadOnly | QIODevice::Text)){
        msgErrLoad("Simula fnmean",fnmean);
        printf("-> Simula ERROR opening file= %s\n",fnmean.toStdString().c_str());
        return;
    }
    else{
        int ipr;
        QString line,pezzo;
        QTextStream stream (&file);
        line = stream.readLine();//info line
        cout << line.toStdString()<<"\n";
        line = stream.readLine();//Ndat & "nm" for new file
        line=line.simplified();
        QStringList List0;
        List0 =line.split(" ");
        pezzo=List0.at(0).toLocal8Bit().constData();
        ipr=pezzo.toInt();
        double a2nm=0.1;
        if(line.contains("nm"))
            a2nm=1.;
        printf("Ndat= %d a2nm=%f\n",ipr,a2nm);
        double wp[ipr+1][3];
        for(int i=1;i<=ipr;i++){
            stream>>wp[i][1]>>wp[i][2];
            wp[i][1]=wp[i][1]*a2nm;
            //printf("\t %f %f\n",wp[i][1],wp[i][2]);
        }
        file.close();
        //weight
        double speso=0.;
        int jjmin=1;
        while(wp[jjmin][1] < MIS[7][1][1] && jjmin<ipr)
            jjmin++;
        int jj=jjmin;
        int jjmax=jj;
        while(wp[jj][1] >= MIS[7][1][1] && wp[jj][1] <= MIS[7][NeV][1] && jj<=ipr){
            speso=speso+wp[jj][2];
            jjmax=jj;
            jj++;
        }
        printf("jjmin= %d wl=%f\n",jjmin,wp[jjmin][1]);
        printf("jjmax= %d wl=%f\n",jjmax,wp[jjmax][1]);
        printf("SUMpeso= %f\n",speso);
        // reset adders
        double av[7],avs[7],val[7][3];
        for( int i=1;i<=6;i++){
            av[i]=.0;
            avs[i]=.0;
        }
        // mean @theta
        for(int ii=jjmin;ii<=jjmax;ii++){
            int ij=1;
            while(wp[ii][1]>MIS[7][ij][1] && ij<=NINT(par[28][2])-1)
                ij++;
            if(ij>1) ij--;
            double x=(wp[ii][1]-MIS[7][ij][1])/(MIS[7][ij+1][1]-MIS[7][ij][1]);
            //printf("x=%f wp[%d][1]=%f MIS[7][%d][1]=%f\n",x,ii,wp[ii][1],ij,MIS[7][ij][1]);
            for(int j=1;j<=2;j++){
                int ik=ij+j-1;
                if(j==1){
                    for(int iv1=1;iv1<=6;iv1++){
                        val[iv1][1]=MIS[iv1][ik][1]*(1.-x);
                        val[iv1][2]=mc[iv1][ik]*(1.-x);
                        //if(iv1==3) printf("mc[3][%d]=%f x=%f\n",ik,mc[iv1][ik],x);
                    }
                }
                else{
                    for(int iv1=1;iv1<=6;iv1++){
                        val[iv1][1]=val[iv1][1]+MIS[iv1][ik][1]*x;
                        val[iv1][2]=val[iv1][2]+mc[iv1][ik]*x;
                    }
                }
            }
            for(int j=1;j<=6;j++){
                av[j]=av[j]+val[j][1]*wp[ii][2]/speso;
                avs[j]=avs[j]+val[j][2]*wp[ii][2]/speso;
                //if(j==3) printf("val[3][2]=%f\n",val[j][2]);
            }
        }
        for(int i=1;i<=6;i++){
            if(DATO[i]!=0){
                par[35+i][2]=av[i]*100.;
                par[35+i][3]=avs[i]*100.;
                if(DATO[i]>0)
                    idToLineEdit["LEpar_"+QString::number(35+i)+"_2"] -> setText(QString::number(par[35+i][2]));
                else
                    idToLineEdit["LEpar_"+QString::number(35+i)+"_2"] -> clear();
                idToLineEdit["LEpar_"+QString::number(35+i)+"_3"] -> setText(QString::number(par[35+i][3]));
            }
            else{
                idToLineEdit["LEpar_"+QString::number(35+i)+"_2"] -> clear();
                idToLineEdit["LEpar_"+QString::number(35+i)+"_3"] -> clear();
            }
        }
        SaveSetting(-1);
    }
}

void ksemawc::CalcMis(double mc[15][1000]){
    printf("->CalcMis\n");
    nextColor();
    int s1p2u3=ui->cB_PAR_35_2 ->currentIndex();
    s1p2u3++;
    int nwl=NeV;
    double teta[6],wl,Xp[nwl],Yp[nwl],ErrXp[nwl],ErrYp[nwl],nn[nwl],kk[nwl],e1[nwl],e2[nwl],VNK[17][3],vot[6][3],delOld[4],del=0.;
    fill_n(ErrXp, nwl, 0.);
    fill_n(ErrYp, nwl, 0.);
    teta[1]=ui->dSB_PAR_6_1->value();
    teta[2]=ui->dSB_PM_86_1->value();
    teta[3]=ui->dSB_PM_87_1->value();
    teta[4]=ui->dSB_PM_88_1->value();
    teta[5]=ui->dSB_PM_89_1->value();
    Qt::CheckState state;
    state=ui->checkBox_deltaConnect->checkState();
    int iDelCon=0;
    if(state==Qt::Checked)
        iDelCon=1;
    for(int k=1;k<6;k++){
        printf("teta[%d]=%f\n",k,teta[k]);
        teta[k]=teta[k]*deg2rad;
    }
    QFile fileNK(fNKsim);
    if (!fileNK.open(QIODevice::WriteOnly | QIODevice::Text)){
        msgErrLoad("CalcMis-fileNK",fNKsim);
        return;
    }
    QTextStream streamNK(&fileNK);
    streamNK<<"nk-values of last simulation\n";
    streamNK<<QString::number(NeV)<<" nm\n";
    for(int i=1;i<=NeV;i++){
        wl=MIS[7][i][1];
        Xp[i-1]=wl;
        COSVNK(VNK,i);
        nn[i-1]=VNK[1][1];
        kk[i-1]=VNK[1][2];
        e1[i-1]=VNK[1][1]*VNK[1][1]-VNK[1][2]*VNK[1][2];
        e2[i-1]=2.*VNK[1][1]*VNK[1][2];
        streamNK<<wl<<"\t"<<VNK[1][1]<<"\t"<<VNK[1][2]<<"\t"<<VNK[1][1]*0.001<<"\t"<<VNK[1][2]*0.001<<"\n";
        if(DATO[1]!=0 || DATO[3]!=0 || DATO[5]!=0){
            ASSEMBLER(i,1,0.,vot);
            mc[1][i]=vot[1][1];
            mc[3][i]=vot[2][1];
            mc[5][i]=vot[3][1];
        }
        else{
            mc[1][i]=.0;
            mc[3][i]=.0;
            mc[5][i]=.0;
        }
        if(DATO[2]!=0 || DATO[4]!=0){
            ASSEMBLER(i,1,teta[1],vot);
            if(s1p2u3<3){
                mc[2][i]=vot[1][s1p2u3];
                mc[4][i]=vot[2][s1p2u3];
            }
            else{
                mc[2][i]=(vot[1][1]+vot[1][2])/2.;
                mc[4][i]=(vot[2][1]+vot[2][2])/2.;
            }
        }
        else{
            mc[2][i]=.0;
            mc[4][i]=.0;
        }
        if(DATO[6]!=0){
            ASSEMBLER(i,2,0.,vot);
            mc[6][i]=vot[4][1];
        }
        else
            mc[6][i]=.0;
        for(int j=1;j<=4;j++){
            if(DATO[7+2*(j-1)]!=0 || DATO[8+2*(j-1)]!=0){
                ASSEMBLER(i,2+j,teta[j+1],vot);
                del=vot[5][2];//Delta
                if(iDelCon==1){
                    if(i==1){
                        delOld[j-1]=del;
                    }
                    del=del-Dperiod*NINT((del-delOld[j-1])/Dperiod);
                    delOld[j-1]=del;
                }
                mc[7+2*(j-1)][i]=del;//+Dperiod*N180offset;//Delta
                mc[8+2*(j-1)][i]=vot[5][1];//Psi
            }
            else{
                mc[7+2*(j-1)][i]=.0;
                mc[8+2*(j-1)][i]=.0;
            }
        }
    }
    fileNK.close();
    int Nperiod;
    double SIMchi2=0.;
    double nMIS=0.;
    //Plot calcData and computing chi2
    for(int ic=1;ic<=14;ic++){
        double FM=0.;
        state=idToCheckBox["checkB_mis"+QString::number(ic)+"_2"] -> checkState ();
        if(state==Qt::Checked && DATO[ic]!=0){
            nMIS++;
            Nperiod=0;
            if(ic==7 || ic==9 || ic==11 || ic==13){
                double Diff=.0;
                for(int i=1;i<=NeV;i++){
                    Diff=Diff+(mc[ic][i]-ELI[ic-6][i][1])/NeV;
                }
                Nperiod=NINT(Diff/Dperiod);
            }
            for(int i=1;i<=NeV;i++){
                Yp[i-1]=mc[ic][i];
                if(ic<=5)
                    Yp[i-1]=Yp[i-1]*100.;
                if(ic<=6)
                    FM=FM+pow((mc[ic][i]-MIS[ic][i][1])/MIS[ic][i][2],2.)/NeV;
                else{
                    FM=FM+pow((mc[ic][i]-ELI[ic-6][i][1]-Nperiod*Dperiod)/ELI[ic-6][i][2],2.)/NeV;
                    Yp[i-1]=Yp[i-1]-Nperiod*Dperiod;
                }
            }
            par[35+ic][1]=sqrt(FM);
            if(ic<=6)
                PLOTline1bar2(1,0,iColor,ic,NeV,Xp,Yp,ErrXp,ErrYp);
            else if(ic==7 || ic==9 || ic==11 || ic==13)
                PLOTline1bar2(1,0,iColor,7,NeV,Xp,Yp,ErrXp,ErrYp);
            else if(ic==8 || ic==10 || ic==12 || ic==14)
                PLOTline1bar2(1,0,iColor,8,NeV,Xp,Yp,ErrXp,ErrYp);
            SIMchi2=SIMchi2+FM;
            printf("ic=%d FM=%f SIMchi2=%f\n",ic,FM,SIMchi2);
        }
    }
    ui->lineEdit_SIMchi2->setText("chi2= "+QString::number(SIMchi2/nMIS));
    if(DATO[1]!=0 && DATO[3]!=0){
        for(int i=1;i<=NeV;i++){
            Yp[i-1]=(1.-mc[1][i]-mc[3][i])*100.;
        }
        labelQwt="A_front";
        PLOTline1bar2(1,0,iColor,9,NeV,Xp,Yp,ErrXp,ErrYp);
    }
    if(DATO[1]!=0 && DATO[5]!=0){
        for(int i=1;i<=NeV;i++){
            Yp[i-1]=(1.-mc[1][i]-mc[5][i])*100.;
        }
        labelQwt="A_back";
        PLOTline1bar2(1,0,iColor,16,NeV,Xp,Yp,ErrXp,ErrYp);
    }
    //if(IXW[12]>0){//plot nk
        PLOTline1bar2(1,0,iColor,12,NeV,Xp,nn,ErrXp,ErrYp);
        PLOTline1bar2(1,0,iColor,13,NeV,Xp,kk,ErrXp,ErrYp);
        if(NINT(par[10][1])==1){//plot permittivity
            PLOTline1bar2(1,0,iColor,14,NeV,Xp,e1,ErrXp,ErrYp);
            PLOTline1bar2(1,0,iColor,15,NeV,Xp,e2,ErrXp,ErrYp);
        }
    //}
    QFile file(fMisSim);
    if (!file.open(QIODevice::WriteOnly | QIODevice::Text)){
        msgErrLoad("CalcMis-fMisSim",fMisSim);
        return;
    }
    QTextStream stream(&file);
    stream<<"wl\tTn\tTp\tRn\tRp\tR1\tApds\tDELTA\tPSI\n";
    for(int i=1;i<=NeV;i++){
        stream<<MIS[7][i][1];
        for(int ii=1;ii<=14;ii++)
            stream<<"\t"<<mc[ii][i];
        stream<<"\n";
    }
    file.close();
}

void ksemawc::PlotAve(){
    //computing mean value SF VS theta
    QString line,pezzo;
    int ipr;
    double vot[6][3],vt[4][3],tet[92],tau[92],rho[92],rho1[92],ErrXp[92],ErrYp[92];
    fill_n(ErrXp, 92, 0.);
    fill_n(ErrXp, 92, 0.);
    SaveSetting(-1);
    int s1p2u3=ui->cB_PAR_35_2 ->currentIndex();
    s1p2u3++;
    // kind of average
    QString fname=pathroot+NANK[10].simplified();
    QFile file(fname);
    if (!file.open(QIODevice::ReadOnly | QIODevice::Text)){
        msgErrLoad("PlotAve fname",fname);
        printf("file %s not found\n",fname.toStdString().c_str());
        return;
    }
    QTextStream stream(&file);
    line=stream.readLine();//info line
    cout<<line.toStdString()<<"\n";
    line=stream.readLine();
    line=line.simplified();
    QStringList List;
    List=line.split(" ");
    pezzo=List.at(0).toLocal8Bit().constData();
    ipr=pezzo.toInt();
    double wl2nm=1.;
    if(line.contains("nm"))
        wl2nm=0.1;
    printf("Ndati= %d\n",ipr);
    double wp[ipr+1][3];
    double speso=0.;
    for(int i=1;i<=ipr;i++){
        stream>>wp[i][1]>>wp[i][2];
        wp[i][1]=wl2nm*wp[i][1];
        cout<< wp[i][1] <<"\t"<< wp[i][2]<<"\n";
        speso=speso+wp[i][2];
    }
    file.close();
    int NT=91;
    double teta,tetar;
    double av[7];
    double ymin,ymax;
    rxy[21][1]=ui->DP_RXY_21_1->text().toDouble();
    rxy[21][2]=ui->DP_RXY_21_2->text().toDouble();
    if(rxy[21][2]<1.){
        rxy[21][1]=0.;
        rxy[21][2]=90.;
        ui->DP_RXY_21_1 -> setText(QString::number(0));
        ui->DP_RXY_21_2 -> setText(QString::number(90));
    }
    for(int I=1;I<=NT;I++){
        teta=rxy[21][1]+(rxy[21][2]-rxy[21][1])*(I-1)/double(NT-1);
        if(teta>90.)
            teta=90.;
        tetar=teta*deg2rad;
        av[1]=.0;
        av[2]=.0;
        av[3]=.0;
        ymin=1.e+36;
        ymax=-1.e+36;
        // computing mean value@teta
        int ij=1;
        double x=0.;
        for(int ii=1;ii<=ipr;ii++){
            if(wp[ii][1]<=MIS[7][1][1]){
                ij=1;
                x=0.;
                if(wp[ii][1]<MIS[7][1][1] && I==1)
                    printf("Attention wl_mean = %f <WLmin!!!\n",wp[ii][1]);
            }
            else if(wp[ii][1]>MIS[7][1][1] && wp[ii][1]<MIS[7][NeV][1]){
                ij=1;
                while(wp[ii][1]>MIS[7][ij][1])
                    ij++;
                if(ij>1)
                    ij--;
                x=(wp[ii][1]-MIS[7][ij][1])/(MIS[7][ij+1][1]-MIS[7][ij][1]);
            }
            else if(wp[ii][1]>=MIS[7][NeV][1]){
                ij=NeV-1;
                x=1.;
                if(wp[ii][1]>MIS[7][NeV][1] && I==1)
                    printf("Attention wl_mean = %f >WLmax!!!!\n",wp[ii][1]);
            }
            for(int j=1;j<=2;j++){
                int ik=ij+j-1;
                if(teta<90.)
                    ASSEMBLER(ik,1,tetar,vot);
                else{
                    vot[1][1]=.0;
                    vot[1][2]=.0;
                    vot[2][1]=1.;
                    vot[2][2]=1.;
                    vot[3][1]=1.;
                    vot[3][2]=1.;
                }
                if(j==1){
                    for(int iv2=1;iv2<=2;iv2++){
                        for(int iv1=1;iv1<=3;iv1++)
                            vt[iv1][iv2]=vot[iv1][iv2]*(1.-x);
                    }
                }
                else{
                    for(int iv2=1;iv2<=2;iv2++){
                        for(int iv1=1;iv1<=3;iv1++)
                            vt[iv1][iv2]=vt[iv1][iv2]+vot[iv1][iv2]*x;
                    }
                }
            }
            for(int j=1;j<=3;j++){
                if(s1p2u3<=2)
                    av[j]=av[j]+vt[j][s1p2u3]*wp[ii][2]/speso;
                else
                    av[j]=av[j]+(vt[j][1]+vt[j][2])/2.*wp[ii][2]/speso;
            }
        }
        for(int j=1;j<=3;j++){
            ymin=min(ymin,av[j]);
            ymax=max(ymax,av[j]);
        }
        tet[I-1]=teta;
        tau[I-1]=av[1]*100.;
        rho[I-1]=av[2]*100.;
        rho1[I-1]=av[3]*100.;
    }
    rxy[23][3]=ymin*100.;
    rxy[23][4]=ymax*100.;
    ui->DP_RXY_23_3 -> setText(QString::number(ymin*100.));
    ui->DP_RXY_23_4 -> setText(QString::number(ymax*100.));
    if(rxy[23][2]<1.){
        rxy[23][1]=ymin*100.;
        rxy[23][2]=ymax*100.;
        ui->DP_RXY_23_1 -> setText(QString::number(int(ymin*100.)));
        ui->DP_RXY_23_2 -> setText(QString::number(NINT(ymax*100.)));
    }
    iColor=1;
    PLOTline1bar2(1,1,iColor,10,91,tet,tau,ErrXp,ErrYp);
    nextColor();
    PLOTline1bar2(1,0,iColor,10,91,tet,rho,ErrXp,ErrYp);
    nextColor();
    PLOTline1bar2(1,0,iColor,10,91,tet,rho1,ErrXp,ErrYp);
    nextColor();
    SaveSetting(-1);
    // saving <T> <R> <R1> vs teta
    QString fname2=pathroot+"expo/ThetaRhoVStheta.dat";
    QFile file2(fname2);
    if (!file2.open(QIODevice::WriteOnly | QIODevice::Text)){
        msgErrLoad("PlotAve fname2",fname2);
        printf("Error opening %s\n",fname2.toStdString().c_str());
        return;
    }
    QTextStream out(&file2);
    out<<" Teta(deg)  Tau  Rho  Rho1"<<"\n";
    out<<NT;
    for(int i=1;i<=NT;i++)
        out<< tet[i]<<"\t"<<tau[i]<<"\t"<<rho[i]<<"\t"<<rho1[i]<<"\n";
    file2.close();
    printf("Data vs theta saved in: expo/ThetaRhoVStheta.dat\n");
}

void ksemawc::PlotAbsEL(){
    int LastMediumAtWl[1002];
    double wwl[1002],x[1002],y[1002],ErrXp[1002],ErrYp[1002],Trasm1[1002],vosi[6][3],denom,ABS[1002][10];
    fill_n(ErrXp, 1002, 0.);
    fill_n(ErrYp, 1002, 0.);
    fill_n(Trasm1, 1002, 1.);
    long double Ps[1002][11],Pp[1002][11];
    complex<long double> Bs[1002][11],Cs[1002][11],Bp[1002][11],Cp[1002][11];
    complex<double> irup,irdw,irup2,irdw2,mups,mupp,mdws,mdwp,pq,ir[10], NQ, MUS, DELTA;
    int icol,ivnkdw,ivnkup,s1p2u3, LastMedium,ivnk;
    double VNK[17][3];
    SaveSetting(-1);

    s1p2u3=NINT(par[35][2]);
    double teta=par[6][1]*deg2rad;
    int iFirstCoe=1;
    int Nlayers=NINT(par[51][2]);
    ivnkdw=9;//output media for SF measurements
    ivnkup=10;//refractive index of input medium for SF measurements
    int iRD=1;//redraw the Absorptance plot

    printf("plot Abs at each layer of %d at theta=%f with polarization %d\n\t Nlayers=%d ivnkdw=%d\n",NINT(par[51][2]),par[6][1],s1p2u3,Nlayers,ivnkdw);

    //set wwl array
    for(int i=1;i<=NeV;i++){
        wwl[i]=MIS[7][i][1];
        for(int j=1;j<=Nlayers;j++)
            ABS[i][j]=-1.1;
    }

    //verify LastMedium at each wl
    for(int i=1;i<=NeV;i++){
        COSVNK(VNK,i);
        irup=complex<double>(VNK[ivnkup][1],-VNK[ivnkup][2]);
        irup2=irup*irup;
        pq=pow(irup*sin(teta),2.);// (n_input*sin(teta))**2.
        LastMediumAtWl[i]=Nlayers+1;//i.e. the medium after the last layer
        for(int j=1;j<=Nlayers;j++){
            int iWorn=0;
            if(NINT(par[50+j][3])==1 && j>1)
                iWorn=1;
            ivnk=int(par[50+j][1]);//index of the material to be assigned to the ii layer
            ir[j]=complex<double>(VNK[ivnk][1],-VNK[ivnk][2]);
            NQ=ir[j]*ir[j];
            MUS=sqrt(NQ-pq);
            DELTA=2.*PIG*pm[j][1]*MUS/wwl[i];
            //printf("iwl=%d iLayer=%d ds=%f abs(real/DELTA))=%f abs(imag(DELTA))=%f\n",i,j,pm[j][1],abs(real(DELTA)),abs(imag(DELTA)));
            if(abs(imag(DELTA))>ABSmax) {
                //if(i==iwl2print)
                    //printf("fully absorbing layer=%d wl[%d]=%f d=%f nm n=%f k=%f \n",j,i, wwl[i],pm[j][1], real(ir[j]), imag(ir[j]));
                LastMediumAtWl[i]=j;
                break;
            }
            if(iWorn==1 && LastMediumAtWl[j]!=j){
                QMessageBox msgBox;
                msgBox.setText("ATTENTION!!!");
                msgBox.setInformativeText("The procedure is not valid for the actual layer model and it will be aborted!");
                msgBox.setStandardButtons(QMessageBox::Ok);
                msgBox.setDefaultButton(QMessageBox::Ok);
                msgBox.exec();
                return;
            }
        }
        if(i==iwl2print) printf("LastMediumAtWl[%d]=%d\n",i,LastMediumAtWl[i]);
    }

    //computing at each wavelength
    for(int i=1;i<=NeV;i++){
        //refractive indices @ wl
        COSVNK(VNK,i);

        //set last medium
        ivnkdw=9;
        LastMedium=LastMediumAtWl[i];
        if(LastMedium<Nlayers+1){
            ivnkdw=NINT(par[50+LastMedium][1]);
            for(int m=LastMedium+1;m<=Nlayers+1;m++){
                Ps[i][m]=0.;
                Pp[i][m]=0.;
            }
        }
        //loop on the layers
        for(int iLayer=1;iLayer<=LastMedium;iLayer++){
            if(i==iwl2print){
                printf("\n|PlotAbsEL@wl=%f> for-loop with iLayer=%d<=%d\n",MIS[7][iwl2print][1],iLayer,LastMedium);
            }
            if(iLayer==1 && (NINT(par[51][3])==1 || LastMedium==1)){//the first layer is inchoerent or completely absorbing
                irup=complex<double>(VNK[ivnkup][1],-VNK[ivnkup][2]);
                pq=pow(irup*sin(teta),2.);// (n_input*sin(teta))**2.
                if(i==iwl2print)
                    printf("\t Incoherent ivnkup=%d nup=%f kup=%f\n",ivnkup, VNK[ivnkup][1],-VNK[ivnkup][2]);
                ivnkdw=NINT(par[51][1]);//refractive index of first layer
                irup=complex<double>(VNK[ivnkup][1],-VNK[ivnkup][2]);
                irup2=irup*irup;
                pq=pow(irup*sin(teta),2.);// (n_input*sin(teta))**2.
                irdw=complex<double>(VNK[ivnkdw][1],-VNK[ivnkdw][2]);
                irdw2=irdw*irdw;
                mups=sqrt(irup2-pq);
                mupp=irup2/mups;
                mdws=sqrt(irdw2-pq);
                mdwp=irdw2/mdws;
                double atten=exp(-4.*PIG*pm[1][1]/wwl[i]*imag(-mdws));
                double RaS=pow(abs((mups-mdws)/(mups+mdws)),2.);
                double RaP=pow(abs((mupp-mdwp)/(mupp+mdwp)),2.);
                //double TaS=1.-RaS;//T=1-R because A=0 at the interface
                //double TaP=1.-RaP;
                //double R1aS=RaS;
                //double R1aP=RaP;
                if(i==iwl2print){
                    printf("The first layer is inchoerent with first interface:\n");
                    printf("nUP=%f kUP=%f nDW=%f kDW=%f d=%f\n",
                           VNK[ivnkup][1],VNK[ivnkup][2],VNK[ivnkdw][1],VNK[ivnkdw][2],pm[1][1]);
                    printf("atten=%f RaS=%f RaP=%f\n",atten,RaS,RaP);
                    printf("2nd interface: call BUILDER: iFirstLayer=%d ncoe=%d\n",2,Nlayers-1);
                }
                BUILDER(i,1,2,Nlayers-1,pq,vosi);
                double RbS=vosi[2][1];
                double RbP=vosi[2][2];
                if(i==iwl2print) printf("RbS=%f RbP=%f\n",RbS,RbP);
                //double TbS=vosi[1][1];
                //double TbP=vosi[1][2];
                //double R1bS=vosi[3][1];
                //double R1bP=vosi[3][2];
                double Abso=0.,Trans=0.;
                if(s1p2u3==1 || s1p2u3==3){
                    denom=(1.-RaS*RbS*atten*atten);
                    Abso=(1.-RaS)*(1.-atten)*(1.+RbS*atten)/denom;
                    Trans=(1.-RaS)*atten/denom;//Transmittance at the entrance of the "b" interface
                }
                if(s1p2u3==2 || s1p2u3==3){
                    denom=(1.-RaP*RbP*atten*atten);
                    Abso=Abso+(1.-RaP)*(1.-atten)*(1.+RbP*atten)/denom;
                    Trans=Trans+(1.-RaP)*atten/denom;//Transmittance at the entrance of the "b" interface
                }
                if(s1p2u3==3){
                    Abso=Abso/2.;
                    Trans=Trans/2.;
                }
                ABS[i][1]=Abso*100.;
                Trasm1[i]=Trans;
                iFirstCoe=2;
                // if(i==iwl2print)
                //     printf("\tFirst inchoerent layer, results: Abso=%f Trans=%f\n",Abso,Trans);

            }
            else{//thin layer
                if(iLayer<LastMedium){
                    if(i==iwl2print)
                        printf("->call BUILDER: iFirstLayer=%d ncoe=%d\n",iLayer,LastMedium-iLayer);
                    BUILDER(i,1,iLayer,LastMedium-iLayer,pq,vosi);
                    Bs[i][iLayer]=freCoeff[4][0];
                    Cs[i][iLayer]=freCoeff[4][1];
                    Bp[i][iLayer]=freCoeff[5][0];
                    Cp[i][iLayer]=freCoeff[5][1];
                    Ps[i][iLayer]=real(freCoeff[0][0]);
                    Pp[i][iLayer]=real(freCoeff[0][1]);
                    if(i==iwl2print)
                        printf("|PlotAbsEL> Bs[%d][%d]=%f + %f*i\n",i,iLayer,real(freCoeff[4][0]),imag(freCoeff[4][0]));
                    if(iLayer==iFirstCoe){
                        Bs[i][iFirstCoe-1]=complex<double>(vosi[2][1],0.);//Rs full stack
                        Bp[i][iFirstCoe-1]=complex<double>(vosi[2][2],0.);
                        if(i==iwl2print){
                            printf("\tiLayer=iFirstCoe=%d -> Re{Bs[%d][%d]}=%Lf\n",iLayer,i,iFirstCoe-1,real(Bs[i][iFirstCoe-1]));
                            printf("\t\tPs[%d][%d]=%Lf \n",i,iLayer,Ps[i][iLayer]);
                        }
                    }
                }
                if(iLayer==LastMedium){
                    irdw=complex<double>(VNK[ivnkdw][1],-VNK[ivnkdw][2]);
                    NQ=irdw*irdw;
                    MUS=sqrt(NQ-pq);
                    Bs[i][iLayer]=complex<double>(1.,0.);
                    Cs[i][iLayer]=MUS;
                    Ps[i][iLayer]=real(Bs[i][iLayer]*conj(Cs[i][iLayer]));//Poynting vector s-pol
                    Bp[i][iLayer]=complex<double>(1.,0.);
                    Cp[i][iLayer]=NQ/MUS;
                    Pp[i][iLayer]=real(Bp[i][iLayer]*conj(Cp[i][iLayer]));//Poynting vector p-pol
                    if(i==iwl2print){
                        printf("|PlotAbsEL> iLayer=LastMedium=%d ->B & C computing with ivnkdw=%d irdw=%f + %f*i\n",iLayer,ivnkdw,real(irdw),imag(irdw));
                        printf("                  -> Ps[%d][%d]=%Lf \n",i,iLayer,Ps[i][iLayer]);
                    }
                }
            }
        }
    }

    //Absorptance computing for thin layers
    for(int i=1;i<=NeV;i++){
        for(int iLayer=1;iLayer<=min(LastMediumAtWl[i],Nlayers);iLayer++){
            if(ABS[i][iLayer]<-1. && iLayer<=LastMediumAtWl[i]){//absorptance must be computed
                if(LastMediumAtWl[i]<iLayer)
                    break;
                ABS[i][iLayer]=0.;
                if(s1p2u3==1 || s1p2u3==3)
                    ABS[i][iLayer]=(Ps[i][iLayer]-Ps[i][iLayer+1])/Ps[i][iFirstCoe]*(1.-real(Bs[i][iFirstCoe-1]));
                if(s1p2u3==2 || s1p2u3==3)
                    ABS[i][iLayer]=ABS[i][iLayer]+(Pp[i][iLayer]-Pp[i][iLayer+1])/Pp[i][iFirstCoe]*(1.-real(Bp[i][iFirstCoe-1]));
                if(s1p2u3==3)
                    ABS[i][iLayer]=ABS[i][iLayer]/2.;
                ABS[i][iLayer]=ABS[i][iLayer]*100.*Trasm1[i];
                if(ABS[i][iLayer]<0.) ABS[i][iLayer]=0.;
                if(i==iwl2print){
                    printf("Layer=%d ASB[%d]=%f  \nPs[%d]=%Le \nPs[%d]=%Le \nPs[%d]=%Le \nRfront=real(Bs[%d])=%Lf \n",
                           iLayer,i, ABS[i][iLayer],iLayer,Ps[i][iLayer],iLayer+1,Ps[i][iLayer+1],iFirstCoe,Ps[i][iFirstCoe],iFirstCoe-1,real(Bs[i][iFirstCoe-1]));
                }
            }
            //if(i==iwl2print)printf("calc abs with iLayer=%d -> ABS[%d]=%f \n",iLayer-1,i,ABS[i][iLayer]);
        }
    }

    //plot & save
    double Amin=1000.;
    double Amax=-1000.;
    for(int j=1;j<=Nlayers;j++){
        // saving Abs to individual files
        QString fname=pathroot+"expo/Abs#"+QString::number(j)+".dat";
        QFile file(fname);
        if (!file.open(QIODevice::WriteOnly | QIODevice::Text)){
            msgErrLoad("PlotAbsEL",fname);
            printf("Error opening %s\n",fname.toStdString().c_str());
            return;
        }
        QTextStream out(&file);
        out<<"Lambda Absorptance Layer#"<<QString::number(j)<<"\n";
        out<<QString::number(NeV)<<"\n";
        int Np=0;
        for(int i=1;i<=NeV;i++){
            if(ABS[i][j]>-1.){
                x[Np]=wwl[i];
                y[Np]=ABS[i][j];
                Amin=min(Amin,ABS[i][j]);
                Amax=max(Amax,ABS[i][j]);
                out<<x[Np]<<"\t"<<y[Np]<<"\n";
                Np++;
            }
        }
        file.close();
        if(Np>0){
            labelQwt="layer_"+QString::number(j);
            icol=j-1;
            if(icol>=7) icol=icol-6;
            PLOTline1bar2(1,iRD,icol,9,Np,x,y,ErrXp,ErrYp);
            if(iRD==1)
                iRD=0;
        }
    }
    ui->DP_RXY_18_3 -> setText(QString::number(Amin));
    ui->DP_RXY_18_4 -> setText(QString::number(Amax));
}


void ksemawc::PlotTexturized(){
    double alpha=ui->doubleSpinBox_alpha->value();
    double theta1=alpha;
    double theta2=180.-3.*alpha;
    ui->lineEdit_texture->setText("theta1="+QString::number(theta1)+" theta2="+QString::number(theta2));
    printf("->Plot texturized\n");
    AdjRoughMax();
    SaveSetting(-1);
    double wl,Xp[NeV],Atext[NeV],Rtext[NeV],Ttext[NeV],ErrXp[NeV],ErrYp[NeV],VNK[17][3],vot[6][3];
    fill_n(ErrXp, NeV, 0.);
    fill_n(ErrYp, NeV, 0.);
    for(int i=1;i<=NeV;i++){
        wl=MIS[7][i][1];
        Xp[i-1]=wl;
        COSVNK(VNK,i);
        ASSEMBLER(i,1,theta1*deg2rad,vot);
        double T1s=vot[1][1];
        double T1p=vot[1][2];
        double R1s=vot[2][1];
        double R1p=vot[2][2];
        ASSEMBLER(i,1,theta2*deg2rad,vot);
        double T2s=vot[1][1];
        double T2p=vot[1][2];
        double R2s=vot[2][1];
        double R2p=vot[2][2];
        Rtext[i-1]=(R1s*R2s+R1p*R2p+R1s*T2s*T1s+R1p*T2p*T1p)/2.*100.;
        Ttext[i-1]=(T1s+R1s*T2s*R1s+T1p+R1p*T2p*R1p)/2.*100.;
        Atext[i-1]=100.-Rtext[i-1]-Ttext[i-1];
        //Rtext[i-1]=(R1s*R2s+R1p*R2p+T1s*R1p+T1p*R1s)/2.*100.;//(R1s*R2s+R1p*R2p+R1s*T2s*T3s+R1p*T2p*T3p)/2.*100.;
        //Ttext[i-1]=(T1s+T1p)/2.*100.;//(T1s+R1s*T2s*R3s+T1p+R1p*T2p*R3p)/2.*100.;
        //printf("wl=%f Rtext=%f Ttext=%f\n",wl,Rtext[i-1],Ttext[i-1]);
    }
    nextColor();
    PLOTline1bar2(1,0,iColor,1,NeV,Xp,Ttext,ErrXp,ErrYp);
    PLOTline1bar2(1,0,iColor,3,NeV,Xp,Rtext,ErrXp,ErrYp);
    PLOTline1bar2(1,0,iColor,9,NeV,Xp,Atext,ErrXp,ErrYp);
}

void ksemawc::manageLEwl(){
    int iChoiceWL=ui->comboBox_searchNK -> currentIndex();
    ui->dSB_WLsearchNK->setEnabled(iChoiceWL==0);
}

void ksemawc::searchNK(){
    //set nk search-range
    par[1][1]=rxy[16][1];
    par[1][2]=rxy[16][2];
    par[2][1]=rxy[17][1];
    par[2][2]=rxy[17][2];

    int Ndat=0,iCol=0,Ndata=NeV*2;
    double Xp[Ndata],Yp[Ndata],ErrXp[Ndata],ErrYp[Ndata],WL;
    int iChoiceWL=ui->comboBox_searchNK -> currentIndex();
    if(iChoiceWL==0)
        WL=ui->dSB_WLsearchNK -> value();
    else if(iChoiceWL==1)
        WL=ui->dSB_PAR_4_1 -> value();
    else
        WL=ui->dSB_PAR_4_2 -> value();
    int iWL=1;
    while(MIS[7][iWL][1]<WL && iWL<NeV)
        iWL++;
    WL=MIS[7][iWL][1];
    printf("->searchNK @ WL= %f\n",WL);
    SaveSetting(-1);
    CNK[1][1]=0.;//nk unknown
    PLOTline1bar2(1,1,iCol,11,0,Xp,Yp,ErrXp,ErrYp);//erase plot
    for(int imis=1;imis<=14;imis++){
        if(DATO[imis]==2){
            if(imis<=6){
                iCol=imis;
                printf("imis=%d: %f +-%f\n ",imis,MIS[imis][iWL][1],MIS[imis][iWL][2]);
                Ndat=SOLVE(imis,iWL,Xp,Yp,ErrYp);
                PLOTline1bar2(2,0,iCol,11,Ndat,Xp,Yp,ErrXp,ErrYp);
            }
            else{
                int I=imis-6;
                if(I==1 || I==3 || I==5 || I==7) iCol=1;
                if(I==2 || I==4 || I==6 || I==8) iCol=4;
                printf("imis=%d: %f +-%f\n ",imis,ELI[I][iWL][1],ELI[I][iWL][2]);
                Ndat=SOLVE(imis,iWL,Xp,Yp,ErrYp);
                PLOTline1bar2(2,0,iCol,11,Ndat,Xp,Yp,ErrXp,ErrYp);
            }
        }
    }
}


void ksemawc::RefTrackG(){
    SaveSetting(-1);
    SPADA();
    PlotMENK();
    par[1][1]=rxy[16][1];
    par[1][2]=rxy[16][2];
    par[2][1]=rxy[17][1];
    par[2][2]=rxy[17][2];
}

void ksemawc::NumericalSearch(){
    SaveSetting(-1);
    lastAction="NumSearch";
    iStop=false;
    int NmisEnab=NINT(par[22][2]);
    printf("Exhaustive numerical n-search in wl-n and wl-k spaces with N.meas=%d sampled on %d points\n",NmisEnab,NeV);
    ui->progressBar_ENS->setValue(0);
    Qt::CheckState state;
    state=ui->checkBox_relMin->checkState();
    int irm=0;
    if(state==Qt::Checked)
        irm=1;
    int calN1K2NK3=3;
    if(NmisEnab==0){
        QMessageBox msgBox;
        msgBox.setText("No measurements is enabled!!!\nPlease select one or more!");
        msgBox.setStandardButtons(QMessageBox::Ok);
        msgBox.exec();
        return;
    }
    else if(NmisEnab==1){
        QMessageBox msgBox(QMessageBox::Question,
                            "Only 1 measure is enabled!!!",  // Window title
                            "Set the unknow to be computed:", // Message
                            QMessageBox::NoButton, // Initially no standard button
                            nullptr); // Widget
        QPushButton *nButton = msgBox.addButton("n", QMessageBox::AcceptRole);
        QPushButton *kButton = msgBox.addButton("k", QMessageBox::YesRole);
        QPushButton *abortButton = msgBox.addButton("Abort", QMessageBox::RejectRole);
        int result = msgBox.exec();
        printf("results=%d\n",result);
        QAbstractButton *clickedButton = msgBox.clickedButton();
        if(clickedButton == nButton){
            calN1K2NK3=1;
        }
        else if(clickedButton == kButton){
            calN1K2NK3=2;
        }
        else if (clickedButton == abortButton){
            return;
        }
        else{
            return;
        }
    }
    else if(NmisEnab>2){
        QMessageBox msgBox;
        msgBox.setText("Please enable only 2 measures!");
        msgBox.setStandardButtons(QMessageBox::Ok);
        msgBox.exec();
        return;
    }
    printf("calN1K2NK3= %d\n",calN1K2NK3);
    double n,k,kMin=100,kMax=-100.;
    double FM,FMold,nMin,dFMold=0.;
    CNK[1][1]=0.;//nk unknown n=CNK[1][2] k=CNK[1][3]
    par[1][1]=rxy[16][1];//nMin
    par[1][2]=rxy[16][2];//nMax
    par[2][1]=rxy[17][1];//kMin
    par[2][2]=rxy[17][2];//kMax
    par[6][1]=ui->dSB_PAR_6_1->value();
    pm[86][1]=ui->dSB_PM_86_1->value();
    pm[87][1]=ui->dSB_PM_87_1->value();
    pm[88][1]=ui->dSB_PM_88_1->value();
    pm[89][1]=ui->dSB_PM_89_1->value();
    int IL1=1;
    int IL2=NeV;
    int iSol=1;
    double Dn=(par[1][2]-par[1][1])/100.;
    k=par[2][1];
    for(int IL=IL1;IL<=IL2;IL++){
        printf("\rProgress: %d%%",int(double(IL-IL1)/double(IL2-IL1)*100.));
        fflush(stdout);
        ui->progressBar_ENS->setValue(int(double(IL-IL1)/double(IL2-IL1)*100.));
        ui->progressBar_ENS->update();
        par[24][2]=IL;
        double lam=MIS[7][IL][1];
        double Trasm=1.,dt,dlim,tol;
        if(NINT(par[50+NINT(par[53][2])][3])==1){
            if(DATO[1]>0) Trasm=MIS[1][IL][1];
            if(DATO[2]>0) Trasm=MIS[2][IL][1];
            double dinco=pm[NINT(par[53][2])][1];
            if(iw==3) printf("Trasm = %f\n",Trasm);
            dt=-0.1*lam/4./3.14/dinco*log(Trasm);
            dlim=1.e-10;
        }
        else{
            dt=0.0001;//DELTA_k for thin film
            dlim=1.e-6;
        }
        tol=dt/1000.;

        FMold=-1.;
        for(int i=0;i<101;i++){
            n=par[1][1]+i*Dn;
            if(calN1K2NK3==3){//searching nk mathematical solutions
                CNK[1][2]=n;
                CNK[1][3]=k;
                double dt=0.001*fabs(rxy[17][2]-rxy[17][1]);
                FM=FindRoot(&FMER,k,dt,dlim,tol);
                k=CNK[1][3];
                //printf("n=%f k=%f FM=%f FMold=%f dFMold=%f\n",n,k,FM,FMold,dFMold);
            }
            else if(calN1K2NK3==1){//n calculation with k=k_simulation
                double VNK[17][3];
                CNK[1][1]=par[35][4];//Fit option
                COSVNK(VNK,IL);
                k=VNK[1][2];
                CNK[1][1]=0;//nk are unknown
                CNK[1][2]=n;
                CNK[1][3]=k;
                FM=FMER(k);
            }
            else if(calN1K2NK3==2){//k calculation with n=n_simulation
                double VNK[17][3];
                CNK[1][1]=par[35][4];//Fit option
                COSVNK(VNK,IL);
                n=VNK[1][1];
                //k=VNK[1][2];
                CNK[1][1]=0;//nk are unknown
                CNK[1][2]=n;
                CNK[1][3]=k;
                double dt=0.001*fabs(rxy[17][2]-rxy[17][1]);
                FM=FindRoot(&FMER,k,dt,dlim,tol);
                k=CNK[1][3];
                iSol++;
                SOL[iSol][1]=lam;
                SOL[iSol][2]=n;
                SOL[iSol][4]=n*0.01;
                SOL[iSol][3]=k;
                double Dk=k/100.;//error evaluation
                do{
                    k=k+Dk;
                    FMold=FMER(k*1.01);
                }while(FMold<1.);
                SOL[iSol][5]=abs(k-SOL[iSol][3]);
                SOL[iSol][6]=1.;//enable data
                i=101;
                continue;
            }
            if((FMold>1.||i==0) && FM<=1.){
                iSol++;
                nMin=n;
                kMin=k;
                kMax=k;
            }
            else if(i>0 && FMold<=1. && FM<=1.){
                kMin=min(kMin,k);
                kMax=max(kMax,k);
            }
            else if(i>0 && FMold<=1. && FM>1.){
                SOL[iSol][1]=lam;
                SOL[iSol][2]=(nMin+n)/2.;
                SOL[iSol][4]=(n-nMin)/2.;
                SOL[iSol][3]=(kMin+kMax)/2.;
                SOL[iSol][5]=(kMax-kMin)/2.;
                SOL[iSol][6]=1.;//enable data
                //printf("sol-nk: %f %f %f %f %f\n",SOL[iSol][1],SOL[iSol][2],SOL[iSol][3],SOL[iSol][4],SOL[iSol][5]);
            }
            else if(irm==1&& dFMold<0. && (FM-FMold)>0. && FM>1. && FMold>1.){//found a minimum
                iSol++;
                SOL[iSol][1]=lam;
                SOL[iSol][2]=n;
                SOL[iSol][3]=k;
                SOL[iSol][4]=n*0.0001;
                SOL[iSol][5]=k*0.0001;
                SOL[iSol][6]=1.;//enable data
                //printf("Found relative minimum -> add solution-point @ lam=%f\n",lam);
            }
            dFMold=FM-FMold;
            FMold=FM;
            QCoreApplication::processEvents();
            if(iStop){
                printf("The exhaustive computing was stopped as requested by the user!\n");
                return;
            }
        }
    }
    SOL[1][1]=iSol-1;
    printf("\nFound %d nk-solutions\n",iSol-1);
    if(iSol-1==0){
        QMessageBox msgBox;
        msgBox.setText("No solution has been found!\nPlease change the n-range of the search");
        msgBox.setStandardButtons(QMessageBox::Ok);
        msgBox.exec();
    }
    nextColor();
    PlotNK(0);
}

void ksemawc::selectNsol(){
    polygonF.clear();
    QMessageBox msgBox;
    msgBox.setText("Please delimit with a polygon the n-solutions to be deleted!\nLeft click to set a point; click on the first point to terminate!");
    msgBox.setStandardButtons(QMessageBox::Ok);
    msgBox.exec();
    picker_m = new QwtPickerClickPointMachine();
    d_picker = new QwtPlotPicker( QwtAxis::XBottom, QwtAxis::YLeft,
        QwtPlotPicker::CrossRubberBand, QwtPicker::AlwaysOn,
        G12_wn->canvas() );
    d_picker->setStateMachine(picker_m);
    connect(d_picker, SIGNAL(selected(QPointF)), this, SLOT(drawPolygon(QPointF)));
}


void ksemawc::drawPolygon(QPointF pos){
    double x1,y1;
    x1=pos.x();
    y1=pos.y();
    if(L1E2==2)
        x1=1240./x1;
    polygonF<<QPointF(x1,y1);
    int Np=polygonF.count();
    //printf("pt_%d: x=%f y=%f\n",Np,x1,y1);
    if(Np<=1)
        return;
    double x[Np],y[Np],ErrXp[Np],ErrYp[Np];
    for(int i=0;i<Np;i++){
        x[i]=polygonF.at(i).x();
        y[i]=polygonF.at(i).y();
        printf("x[%d]=%f y[%d]=%f\n",i,x[i],i,y[i]);
    }
    PLOTline1bar2(1,0,1,12,Np,x,y,ErrXp,ErrYp);//n plot
    double delta=pow((x[Np-1]-x[0])/(rxy[20][2]-rxy[20][1]),2.)+
                   pow((y[Np-1]-y[0])/(rxy[16][2]-rxy[16][1]),2.);
    if(sqrt(delta)>0.01)
        return;
    QPointF ptF;
    int nSol=NINT(SOL[1][1]);
    printf("nSol=%d\n",nSol);
    for(int i=2;i<=nSol+1;i++){
        ptF.setX(SOL[i][1]);
        ptF.setY(SOL[i][2]);
        if(polygonF.containsPoint(ptF,Qt::OddEvenFill)){
            //printf("i=%d wl=%f n=%f is contained! It will be canceled\n",i,x1,y1);
            SOL[i][6]=0.;//to be canceled
        }
        else
            SOL[i][6]=1.;
    }
    fflush(stdout);
    //cancel selected nk-solutions
    int newN=0;
    double NKNEW[1000][6];
    for(int I=2;I<=nSol+1;I++){
        if(SOL[I][6]>0.5){
            newN++;
            for(int k=1;k<=5;k++)
                NKNEW[newN+1][k]=SOL[I][k];
        }
    }
    SOL[1][1]=newN;
    for(int I=2;I<=newN+1;I++){
        for(int k=1;k<=5;k++)
            SOL[I][k]=NKNEW[I][k];
        SOL[I][6]=1;
    }
    printf("now the nk-solutions are %d\n",newN);
    PlotNK(1);//plot nk
    //stop motor&picker
    delete picker_m;
    picker_m=nullptr;
    d_picker->setEnabled(false);
}


void ksemawc::IbridPlotFit(){
    lastAction="IbridCurrent";
    SaveSetting(-1);
    IbridKernel("g");
}

// void ksemawc::IbridPlotIbrid(){
//     SaveSetting(-1);
//     IbridKernel("gi");
// }

void ksemawc::FitN(){
    lastAction="FitN";
    LoadFilenk();
    par[32][5]=0.;//Fit n in IbridOne
    IbridFit();
}

void ksemawc::FitNK(){
    lastAction="FitNK";
    LoadFilenk();
    par[32][5]=1.;//Fit n&K in IbridOne
    IbridFit();
}

void ksemawc::FitE1E2(){
    lastAction="FitE1E2";
    LoadFilenk();
    par[32][5]=2.;//Fit epi&epsi2 in IbridOne
    par[10][1]=1.;//Plot epsi1 and epsi2
    IbridFit();
}

void ksemawc::FitSelExpMeas(){
    lastAction="FitExpMeas";
    SaveSetting(3);
    par[32][5]=3.;//Fit Selected Experimental Measurement in IbridOne
    IbridKernel("fsem");
}

void ksemawc::IbridFit(){
    Qt::CheckState state;
    npp=ui->sB_PAR_34_5 -> value();
    for(int j=1;j<=npp;j++){
        state=idToCheckBox["chBeParFit_"+QString::number(j)]-> checkState();
        if(ppm[j]<71 && state==Qt::Checked)
            idToCheckBox["chBeParFit_"+QString::number(j)] -> setCheckState ( Qt::Unchecked );
    }
    PanFitEnable();
    SaveSetting(3);
    IbridKernel("f");
}

void ksemawc::IbridOne(){
    lastAction="IbridOne";
    SaveSetting(3);
    IbridKernel("i");
}

void ksemawc::IbridOneStore(){
    lastAction="IbridOne";
    SaveSetting(3);
    IbridKernel("M");
}

void ksemawc::ClearTempIbri(){
    SOL[1][1]=0.;
}

void ksemawc::setStop(){
    iStop=true;
    printf("-> setStop: iStop= true\n");
    fflush(stdout);
}

void ksemawc::IbridKernel(QString rc){
    PanFitEnable();
    // data number
    printf("-> IbridKernel with rc=%s\n",rc.toStdString().c_str());
    int mwl=NeV;//SF and ELI
    int mwla=NeV;//enabled wavelengths
    int md=NINT(SOL[1][1]);//nk-data
    int mda=0;//enabled nk-data
    int ial=0;
    int nparmax=17;
    int ibridok,m=0,ifit;
    int Info=0;
    double wwl[mwl],kk[mwl],nn[mwl],e1[mwl],e2[mwl],vot[6][3],tn[mwl],tp[mwl],rn[mwl],rp[mwl],r1[mwl],ErrXp[mwl],ErrYp[mwl],
            miss[6][NeV+1],elis[9][NeV+1],pst[nparmax+1][3],ncent[mwl],ndev[mwl],kcent[mwl],kdev[mwl],ps[NeV],psi[mwl],delta[mwl];
    double PsiMat[mwl][4],DelMat[mwl][4];
    double chi2ini=1.E+09;
    par[38][1]=1.e+20;
    par[38][2]=0.;
    for(int i=2;i<=md+1;i++){
        if(NINT(SOL[i][6])>0){
            if(SOL[i][1]>=par[4][1] && SOL[i][1]<=par[4][2]){
                mda=mda+1;
                par[38][1]=min(par[38][1],SOL[i][1]);
                par[38][2]=max(par[38][2],SOL[i][1]);
            }
            else
                SOL[i][6]=0.;
        }
        if(SOL[i][4]<=0.){
            if(ial==0){
                printf("ERRn = 0 @%f  => ERRn = Dn =%f\n",SOL[i][1],par[21][3]);
                ial=1;
            }
            SOL[i][4]=abs(par[21][3]);
        }
    }
    printf("WLmin=%f WLmax=%f md=%d mda=%d\n",par[4][1],par[4][2],md,mda);
    par[38][3]=1.e+20;
    par[38][4]=0.;
    for(int i=1;i<=NeV;i++){
            par[38][3]=min(par[38][3],MIS[7][i][1]);
            par[38][4]=max(par[38][4],MIS[7][i][1]);
    }
    //initialization par for fit
    int n=ui->sB_PAR_35_5 -> value();//NINT(par[35][5]);
    double p[n];
    //double* p=nullptr;
    //p = new double[n];
    for(int i=1;i<=n;i++){
        int ipm=NINT(pm[i][3]);
        p[i-1]=pm[ipm][1];
        //printf("%d ipm=%d p[%d]=%f pm[%d][2]=%f\n",i,ipm,i,p[i-1],ipm,pm[ipm][2]);
    }
    if(ifirstcall==0){// IbridKernel has been called
        ifirstcall++;
    }
    // initialization
    int ioptf=NINT(pm[100][1]);
    double chi2fin=par[55][2];
    double fredeg=par[56][2];
    // minimal initialization
    CNK[1][1]=0.;//nk unknown
    double te=0.;//incident angle
    //int imR=3;//Rn
    //int iwr=13;// Rn
    int s1p2=NINT(par[27][2]);//polarization
    int nmisure=0;
    for(int i=1;i<=5;i++){
        if(DATO[i]==2)
            nmisure++;
    }
    int jie=0;
    int jiemax=pow(2,nmisure);
    int ieon=0;//error computing is OFF
    QString r="?";
    while(r!="X"){
        // control loop of IbridOne with error
        if(ieon==1){
            jie=jie+1;
            if(jie>jiemax){//stop
                for(int ii=1;ii<=NeV;ii++){
                    for(int i=1;i<=5;i++)
                        MIS[i][ii][1]=miss[i][ii];//restore measured values
                    for(int i=1;i<=8;i++)
                        ELI[i][ii][1]=elis[i][ii];
                }
                for(int i=1;i<=n;i++){
                    int ip=NINT(pm[i][3]);
                    p[i-1]=pst[i][1];//central value (jie=0)
                    pm[ip][1]=p[i-1];
                    if(pst[i][2]>.0)
                        pm[ip][4]=sqrt(pst[i][2]);//error = deviation rms
                    else
                        pm[ip][4]=0.;
                }
                for(int i=0;i<NeV;i++){
                    //store central value of solution
                    nn[i]=ncent[i];
                    kk[i]=kcent[i];
                    MIS[16][i+1][1]=nn[i];
                    MIS[16][i+1][2]=kk[i];
                    if(ndev[i]>0.)
                        ndev[i]=sqrt(ndev[i]);
                    else
                        ndev[i]=0.;
                    if(kdev[i]>0.)
                        kdev[i]=sqrt(kdev[i]);
                    else
                        kdev[i]=0.;
                }
//                int isol=NINT(SOL[1][1]);
//                int irima=998-isol;
                printf("<<<<<<<<<< Save nk IbridOne >>>>>>>>>>>>>>\n");
                SOL[1][1]=mwl;
                printf("Saving nk-solutions in temp:\n");
//                printf(" N. solutions already stored = %d\n",isol);
//                printf(" N. solutions to add = %d\n", mwl);
//                printf(" available space = %d\n",irima);
                for(int i=0;i<mwl;i++){
                    SOL[i+2][1]=wwl[i];
                    SOL[i+2][2]=ncent[i];
                    SOL[i+2][3]=kcent[i];
                    SOL[i+2][4]=ndev[i];
                    SOL[i+2][5]=kdev[i];
                    SOL[i+2][6]=1.;//data enabled
                }
                rc="X";
                ieon=0;
                SaveFnk();
            }
        }

        if(ieon==0)
            r=rc;
        else if(ieon==1)
            r="i";
        if(r=="M" && ieon==0){//set IbridOne for error computing
            ieon=1;
            jie=0;
            if(nmisure>=1){
                r="i";//IbridOne
                for(int ii=1;ii<=NeV;ii++){
                    for(int i=1;i<=5;i++)
                        miss[i][ii]=MIS[i][ii][1];//store measured values
                    for(int i=1;i<=8;i++)
                        elis[i][ii]=ELI[i][ii][1];
                }
            }
            else
                r=" ";//stop error computing
        }

        if(r=="i" && ieon==1 && jie>0){
            // measurement perturbation with the error
            int imisura=1;
            printf("..... IbridOne with error computing: step %d of %d\n",jie,jiemax);
            for(int i=1;i<=5;i++){
                if(DATO[i]==2){
                    int iespo=int((jie+1)/imisura);
                    printf("  measure= %d (-1)^%d= %d\n",i,iespo,int(pow(-1,iespo)));
                    for(int ii=1;ii<=NeV;ii++)
                        MIS[i][ii][1]=miss[i][ii]+MIS[i][ii][2]*pow(-1.,iespo);
                    imisura++;
                }
            }
            for(int i=8;i<=14;i=i+2){
                if(DATO[i]==2){
                    int iespo=int((jie+1)/imisura);
                    printf("  measure= %d (-1)^%d= %d\n",i,iespo,int(pow(-1,iespo)));
                    for(int ii=1;ii<=NeV;ii++)
                        ELI[i-7][ii][1]=elis[i-1][ii]+ELI[i-7][ii][2]*pow(-1.,iespo);
                    imisura++;
                }
            }
        }

        if(r=="f" || r=="fsem" || r=="i"){
            // launch fit or ibridOne
            iStop=false;
            if(r=="i"){//set IbridOne computing
                ibridok=0;
                if(DATO[1]==2){
                    printf("-> launch IbridOne: compute k from Tn\n");
                    ibridok=1;
                }
                else if(DATO[2]==2){
                    printf("-> launch IbridOne: compute k from Tp teta= %f deg  S1P2= %d\n",par[6][1],s1p2);
                    ibridok=1;
                }
                else if(DATO[1]!=2 && DATO[2]!=2){
                    printf("-> ABORT IbridOne because transmittance is not enabled\n");
                    ibridok=-100;
                    r=" ";//abort
                    QMessageBox msgBox;
                    msgBox.setText("IbridOne aborted because transmittance is not enabled\n");
                    msgBox.setStandardButtons(QMessageBox::Ok);
                    msgBox.exec();
                    return;
                }
                if(DATO[3]==2){
                    printf("\tfit of Rn\n");
                    ibridok++;
                }
                if(DATO[4]==2){
                    printf("\tfit of Rp  teta= %f deg  S1P2= %d\n",par[6][1],s1p2);
                    ibridok++;
                }
                if(DATO[5]==2){
                    printf("\tfit of R1\n");
                    ibridok++;
                }
                if(DATO[8]==2){
                    printf("\tfit of PSI_1\n");
                    ibridok++;
                }
                if(DATO[10]==2){
                    printf("\tfit of PSI_2\n");
                    ibridok++;
                }
                if(DATO[12]==2){
                    printf("\tfit of PSI_3\n");
                    ibridok++;
                }
                if(DATO[14]==2){
                    printf("\tfit of PSI_4\n");
                    ibridok++;
                }
                if(ibridok<2){
                    printf("IbridOne aborted because needs T&R or T&DELTA-PSI!\n");
                    r=" ";//abort
                    QMessageBox msgBox;
                    msgBox.setText("IbridOne aborted because needs T&R or T&DELTA-PSI!");
                    msgBox.setStandardButtons(QMessageBox::Ok);
                    msgBox.exec();
                    return;
                }
                else{
                    m=(ibridok-1)*mwla;
                    printf("     ... on %d data and %d parameters\n",m,n);
                }
            }
            else if(r=="f"){// suppress p_fit ><  SELMQ
                if(NINT(par[32][5])==0)
                    m=mda;
                else
                    m=2*mda;
                int npfit=n;
                for(int ipf=npfit;ipf>=1;ipf--){
                    int ipm=NINT(pm[ipf][3]);
                    printf("ipf = %d ipm = %d\n",ipf,ipm);
                    if(ipm<100){
                        pm[ipm][2]=.0;
                        for(int iii=ipf;iii<=n;iii++){
                            pm[iii][3]=pm[iii+1][3];// move 1 step down
                            pm[iii][4]=pm[iii+1][4];
                            pm[iii][5]=pm[iii+1][5];
                        }
                        n=n-1;
                        printf("n = %d\n",n);
                    }
                    par[35][5]=n;
                }
                QString fitwhat;
                if(NINT(par[32][5])==0)
                    fitwhat="n";
                else if(NINT(par[32][5])==1)
                    fitwhat="n & k";
                else if(NINT(par[32][5])==2)
                    fitwhat="epsi1 & epsi2";
                printf("launch fit of %s\n",fitwhat.toStdString().c_str());
            }
            else if(r=="fsem"){
                int iSEM=0;
                for(int iM=1;iM<=14;iM++){
                    if(DATO[iM]==2)
                        iSEM++;
                }
                if(iSEM==0)
                    r="X";
                m=iSEM*mwl;
                printf("launch fit of Selected Experimental Measurements: Nsem=%d Ndata=%d\n",iSEM,m);
            }

            //Fit, compute jac & cov matrices, err. e corr.
            //int* iwa=nullptr;
            //iwa = new int[n];
            //double* fvec=nullptr;
            //fvec = new double[m];
            //double* wa=nullptr;
            //wa = new double[lwa];
            fredeg=double(m-n);
            par[56][2]=fredeg;
            ifit=1;
            QFile fchi2(filechi2);
            fchi2.resize(0);
            fchi2.close();
            int iwa[n];
            int lwa=m*n+5*n+m;
            double wa[lwa];
            double fvec[m];
            double FJAC[m*n];
            double tol=sqrt(dpmpar(1));
            CNK[1][1]=0;//nk unknown
            printf("n(param)= %d m(Ndata)= %d lwa= %d fredeg=%f tol=%e s1p2=%d\n",
                   n,m,lwa,fredeg,tol,NINT(par[27][2]));
            if(n==0){
                QMessageBox msgBox;
                msgBox.setText("ATTENTION: please enable at least one parameter for best-fit!!!\n");
                msgBox.setStandardButtons(QMessageBox::Ok);
                msgBox.exec();
                return;
            }
            if(m==0){
                QMessageBox msgBox;
                msgBox.setText("ATTENTION: please select at least one measurement!!!\n");
                msgBox.setStandardButtons(QMessageBox::Ok);
                msgBox.exec();
                return;
            }
            iRecChi2=1;//save chi2 values at the end of called functions
            if(r=="f" && n>0){
                Info=lmdif1(FSQ, &pTF2,m,n,p,fvec,tol,iwa,wa,lwa);
                //printf("Info_FSQ= %d\n",Info);
            }
            else if(r=="fsem" && n>0 && m>0){
                Info=lmdif1(FSEM, &pTF2,m,n,p,fvec,tol,iwa,wa,lwa);
                //printf("Info_FSEM= %d\n",Info);
            }
            else if(r=="i"){
                Info=lmdif1(FRCK, &pTF2,m,n,p,fvec,tol,iwa,wa,lwa);
                //printf("Info_FRCK= %d\n",Info);
            }
            iRecChi2=0;

            //analysis chi2
//            if((r=="f" || r=="fsem" || (r=="i" && ieon==0)) && n>0){
                QFile fchi2b(filechi2);
                if(fchi2b.open(QIODevice::ReadOnly | QIODevice::Text)){
                    QTextStream inp(&fchi2b);
                    QString linea,st1;
                    int iLine=0;
                    double chi2;
                    do{
                        linea=inp.readLine();
                        //printf("iLine=%d-> %s\n",iLine,linea.toStdString().c_str());
                        QStringList list;
                        list = linea.split("\t");
                        st1=list.at(0).toLocal8Bit().constData();
                        chi2=st1.toDouble();
                        if(iLine==0){
                            chi2ini=chi2;
                            chi2fin=chi2ini;
                        }
                        if(iLine==0 || chi2<chi2fin){
                            chi2fin=chi2;
                            for(int i=0;i<n;i++){
                                st1=list.at(1+i).toLocal8Bit().constData();
                                p[i]=st1.toDouble();
                            }
                        }
                        //printf("chi2=%f  chi2ini=%f chi2fin=%f\n",chi2,chi2ini,chi2fin);
                        iLine++;
                    }while(!inp.atEnd());
                    fchi2b.close();
                    par[27][1]=chi2ini;
                    par[55][2]=chi2fin;
                    ui->DP_PAR_27_1->setText(QString::number(chi2ini));
                    ui->DP_PAR_55_2->setText(QString::number(chi2fin));
                }
//            }

            // set BF parameters
            //fnorm=enorm(m,fvec);
            printf("BF: chi2= %.9g\n",chi2fin);
            for(int i=1;i<=n;i++){
                int ip=NINT(pm[i][3]);
                pm[ip][1]=p[i-1];
                //check value
                if((ip>0 && ip<10) || (ip>50 && ip<60) || ip>100)//thickness or roughness or oscillator
                    if(pm[ip][1]<.0)
                        pm[ip][1]=0.;
                if((ip>=71 && ip<=85) || (ip>=91 && ip<=99)){//fEMA or d_nonuni
                    if(pm[ip][1]<.0)
                        pm[ip][1]=0.;
                    if(pm[ip][1]>1.)
                        pm[ip][1]=1.;
                }
                pm[ip][4]=.0;
                printf("\tp[%d]=pm[%d][1]=%.9g\n",i-1,ip,p[i-1]);
            }
            printf("\n");
            if(ifit==1){//write Info fit
                printf("********************INFO FIT***********************\n");
                printf("Enabled point m=%d\n",m);
                printf("N parameters n=%d => degree of freedom=%d\n",n,NINT(fredeg));
    //            if(n>0){
    //                printf("Chi^2: %f -> %f\n",chi2ini,chi2fin);
    //                par[27][1]=chi2ini;
    //            }
    //            else{
    //                printf("Chi^2 =%f\n",chi2fin);
    //                par[27][1]=chi2fin;
    //            }
                //printf("Final L2 norm of the residual = %f\n",fnorm);
                if(Info<0)
                    printf("excetution terminated by user\n");
                else if(Info==0)
                    printf("Info=0: improper input parameters\n");
                else if(Info==1)
                    printf("Info=1: algorithm estimates that the relative error in the sum of squares is at most tol\n");
                else if(Info==2)
                    printf("Info=2: algorithm estimates that the relative error between x and the solution is at most tol\n");
                else if(Info==3)
                    printf("Info=3: conditions for Info = 1 and Info = 2 both hold\n");
                else if(Info==4)
                    printf("Info=4:  fvec is orthogonal to the columns of the jacobian to machine precision\n");
                else if(Info==5)
                    printf("Info=5:  number of calls to fcn has reached or exceeded %d\n",200*(n+1));
                else if(Info==6)
                    printf("Info=6: tol is too small. no further reduction in the sum of squares is possible\n");
                else if(Info==7)
                    printf("Info=7: tol is too small. no further improvement in the approximate solution x is possible\n");
                printf("*************************************************\n");
                ifit=0;
            }

            if(m>0 && ieon==0){
//                for(int i=0;i<m;i++)
//                    printf("fvec[%d]=%f\twa[%d]=%f\n",i,fvec[i],i,wa[i]);
                //computing jabobiana fjac
                for(int i=1;i<=200;i++)//save actual status
                    ps[i]=pm[i][1];
                //int iflag=2;
                double epsfcn=100.*tol;//0.;
                printf("epsfcn=%e\n",epsfcn);
                if(r=="f" && n>0)
                    fdjac2(FSQ,&pTF2,m,n,p,fvec,FJAC,m,epsfcn,wa);
                else if(r=="fsem" && n>0)
                    fdjac2(FSEM,&pTF2,m,n,p,fvec,FJAC,m,epsfcn,wa);
                else if(r=="i" && n>0)
                    fdjac2(FRCK,&pTF2,m,n,p,fvec,FJAC,m,epsfcn,wa);
                for(int i=1;i<=200;i++)//return to best fit
                    pm[i][1]=ps[i];
                for(int i=1;i<=n;i++){
                    p[i-1]=pm[NINT(pm[i][3])][1];
                }

                //computing matrix inv covar (fjac_trasposta)*(fjac)
                double **cinv;
                cinv= new double *[n];
                for(int i=0;i<n;i++)
                    cinv[i]=new double[n];
                double **cov;
                cov= new double *[n];
                for(int i=0;i<n;i++)
                    cov[i]=new double[n];
                //double cinv[n][n],cov[n][n];
                for(int i=0;i<n;i++){
//                    for(int i2=0;i2<m;i2++)
//                        printf("fjac[%d][%d]=%f\n",i2,i,FJAC[i2+i*m]);
                    for(int j=0;j<n;j++){
                        cinv[i][j]=.0;
                        for(int i2=0;i2<m;i2++)
                            cinv[i][j]=cinv[i][j]+FJAC[i2+i*m]*FJAC[i2+j*m];//fjac[i + j * ldfjac]
                        //printf("cinv[%d][%d]=%f\n",i,j,cinv[i][j]);
                    }
                }

                //compunting matrix varianza-covarianza
                printf("computing MATINV...");
                par[11][2]=0.;
                MATINV(n,n,cinv,cov);
                if(NINT(par[11][2])==0){
                    //MATINV(n,nparmax,cinv,cov);
                    printf("done!\n");
                    //                for(int i=0;i<n;i++){
                    //                    for(int j=0;j<n;j++)
                    //                        cout<<cov[i][j]<<"\n";
                    //                }

                    //computing errors at 99% of confidency
                    for(int i=1;i<=n;i++){
                        if(cov[i-1][i-1]<1.e3){
                            pm[NINT(pm[i][3])][4]=sqrt(abs(6.635*cov[i-1][i-1]));
                            //printf("pm[%d][4]=%f\n",NINT(pm[i][3]),pm[NINT(pm[i][3])][4]);
                        }
                    }
                    fflush(stdout);

                    //computing global correlation
                    if(n>1){
                        for(int i=1;i<=n;i++){
                            if(1.-1./(cov[i-1][i-1]*cinv[i-1][i-1])>0.)
                                pm[NINT(pm[i][3])][5]=sqrt(1.-1./(cov[i-1][i-1]*cinv[i-1][i-1]));
                            else
                                pm[NINT(pm[i][3])][5]=1.;
                            if(pm[NINT(pm[i][3])][5]>1.)
                                pm[NINT(pm[i][3])][5]=1.;
                            //printf("pm[%d][5]=%f\n",NINT(pm[i][3]),pm[NINT(pm[i][3])][5]);
                        }
                    }
                }
                //refresh PanFit
                Qt::CheckState state;
                for(int i=1;i<=npp;i++){
                    state=idToCheckBox["chBeParFit_"+QString::number(i)]-> checkState();
                    int ip=NINT(ppm[i]);
                    //check value
                    if((ip>0 && ip<10) || (ip>50 && ip<60) || ip>100)//thickness or roughness or oscillator
                        if(pm[ip][1]<.0)
                            pm[ip][1]=0.;
                    if((ip>=71 && ip<=85) || (ip>=91 && ip<=99)){//fEMA or d_nonuni
                        if(pm[ip][1]<.0)
                            pm[ip][1]=0.;
                        if(pm[ip][1]>1.)
                            pm[ip][1]=1.;
                    }
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
//            if(wa){
//                delete[] wa;
//                wa = nullptr;
//            }
//            if(fvec){
//                delete [] fvec;
//                fvec=nullptr;
//            }
//            if(iwa){
//                delete []  iwa;
//                iwa=nullptr;
//            }

        }

        if(r=="f" || r=="fsem" || r=="i" || r=="g" || r=="gi"){
            //refresh plots
            nextColor();
            double ymin=1.e6;
            double ymax=-1.e6;
            double chi2=.0;
            double delOld[4],del,delDiff[4]={0};
            printf("r=%s Computing n by ioptf=%d\n",r.toStdString().c_str(),ioptf);
            if(r=="fsem" || r=="f")
                SOL[1][1]=0;
            for(int i=1;i<=mwl;i++){
                par[24][2]=i;//iWL
                double wl=MIS[7][i][1];
                double eV=1240./wl;
                if(r=="f" || r=="fsem" || r=="i" || r=="g"){
                    FDISP(ioptf,eV);
                    CNK[1][2]=sqn;
                    CNK[1][3]=sqk;
                    MIS[16][i][1]=sqn;
                    if(r=="i" || r=="g"){//IbridOne, compute k with best-fit parameters
                        par[24][2]=i;
                        MIS[16][i][1]=sqn;
                        double Trasm=1.,dt,dlim,tol;
                        if(NINT(par[50+NINT(par[53][2])][3])==1){
                            if(DATO[1]>0) Trasm=MIS[1][i][1];
                            if(DATO[2]>0) Trasm=MIS[2][i][1];
                            double dinco=pm[NINT(par[53][2])][1];
                            dt=-0.1*wl/4./3.14/dinco*log(Trasm);
                            dlim=1.e-10;
                            if(iw>0)
                                printf("Trasm = %f dt=%f\n",Trasm,dt);
                        }
                        else{
                            dt=0.001;//DELTA_k for thin film
                            dlim=1.e-5;
                        }
                        tol=dt/1000.;
                        //CNK[1][3]=MIS[16][i][2];//<<<<<<<<<<<<<prima era commentata
                        double k=CNK[1][3];
                        FindRoot(&DELTAT,k,dt,dlim,tol);
                        MIS[16][i][2]=CNK[1][3];
                        //printf("wl=%f n=%f k=%f\n",wl,CNK[1][2],CNK[1][3]);
                    }
                    else if(r=="f" || r=="fsem"){
                        MIS[16][i][1]=CNK[1][2];
                        MIS[16][i][2]=CNK[1][3];
                        SOL[1][1]++;
                        SOL[i+1][1]=wl;
                        SOL[i+1][2]=CNK[1][2];
                        SOL[i+1][3]=CNK[1][3];
                        SOL[i+1][4]=0.001;
                        SOL[i+1][5]=CNK[1][3]/100.;
                        SOL[i+1][6]=1.;//data enabled
                    }
                }
                else if(r=="gi"){
                    CNK[1][2]=MIS[16][i][1];
                    CNK[1][3]=MIS[16][i][2];
                }
                wwl[i-1]=wl;
                nn[i-1]=CNK[1][2];
                kk[i-1]=CNK[1][3];
                e1[i-1]=CNK[1][2]*CNK[1][2]-CNK[1][3]*CNK[1][3];
                e2[i-1]=2.*CNK[1][2]*CNK[1][3];
                ymax=max(kk[i-1],ymax);
                ymin=min(kk[i-1],ymin);
                ErrXp[i-1]=0.;
                ErrYp[i-1]=0.;
                //computing Tn Rn R1
                if(DATO[1]==2 || DATO[3]==2 || DATO[5]==2){
                    ASSEMBLER(i,1,0.,vot);
                    tn[i-1]=vot[1][s1p2]*100.;
                    rn[i-1]=vot[2][s1p2]*100.;
                    r1[i-1]=vot[3][s1p2]*100.;
                    if(DATO[1]==2)
                        chi2=chi2+pow((MIS[1][i][1]-tn[i-1]/100.)/MIS[1][i][2],2.)/fredeg;
                    if(DATO[3]==2)
                        chi2=chi2+pow((MIS[3][i][1]-rn[i-1]/100.)/MIS[3][i][2],2.)/fredeg;
                    if(DATO[5]==2)
                        chi2=chi2+pow((MIS[5][i][1]-r1[i-1]/100.)/MIS[5][i][2],2.)/fredeg;
                }
                if(DATO[2]==2 || DATO[4]==2){
                    te=par[6][1]*deg2rad;
                    ASSEMBLER(i,1,te,vot);
                    tp[i-1]=vot[1][s1p2]*100.;
                    rp[i-1]=vot[2][s1p2]*100.;
                    if(DATO[2]==2)
                        chi2=chi2+pow((MIS[2][i][1]-rp[i-1]/100.)/MIS[2][i][2],2.)/fredeg;
                    if(DATO[4]==2)
                        chi2=chi2+pow((MIS[4][i][1]-rp[i-1]/100.)/MIS[4][i][2],2.)/fredeg;
                }
                for(int ieli=0;ieli<4;ieli++){
                    if(DATO[7+2*ieli]==2 ||DATO[8+2*ieli]==2){
                        te=pm[86+ieli][1]*deg2rad;;//par[14+ieli][1]*deg2rad;
                        ASSEMBLER(i,1,te,vot);
                        del=vot[5][2];//Delta
                        if(iDelCon==1){
                            if(i==1){
                                delOld[ieli]=del;
                            }
                            del=del-Dperiod*NINT((del-delOld[ieli])/Dperiod);
                            delOld[ieli]=del;
                        }
                        PsiMat[i-1][ieli]=vot[5][1];
                        DelMat[i-1][ieli]=del;
                        if(DATO[7+2*ieli]>0)
                            delDiff[ieli]=delDiff[ieli]+(del-ELI[1+2*ieli][i][1])/mwl;
                    }
                }
            }
            for(int ieli=0;ieli<4;ieli++){
                if(DATO[7+2*ieli]>0)
                    printf("|IbridKernel> delDiff[%d]=%f\n",ieli,delDiff[ieli]);
            }
            if(ieon==1){// save central value and RMS deviation
                if(jie==0){
                    printf("***setting ncent & kcent and initialize ndev kdev****\n");
                    for(int i=0;i<NeV;i++){
                        ncent[i]=nn[i];
                        ndev[i]=0.;
                        kcent[i]=kk[i];
                        kdev[i]=0.;
                        printf("wl=%f ncent=%f kcent=%f\n",wwl[i],ncent[i],kcent[i]);
                    }
                    for(int i=0;i<n;i++){
                        pst[i][1]=p[i];
                        pst[i][2]=0.;
                    }
                }
                else{
                    printf("***computing ndev & kdev****\n");
                    for(int i=0;i<NeV;i++){
                        ndev[i]=ndev[i]+pow(ncent[i]-nn[i],2.)/jiemax;
                        kdev[i]=kdev[i]+pow(kcent[i]-kk[i],2.)/jiemax;
                    }
                    for(int i=0;i<n;i++)
                        pst[i][2]=pst[i][2]+pow(pst[i][1]-p[i],2.)/jiemax;
                }
            }

            //refresh plot
            PLOTline1bar2(1,0,iColor,12,mwl,wwl,nn,ErrXp,ErrYp);//n plot
            PLOTline1bar2(1,0,iColor,13,mwl,wwl,kk,ErrXp,ErrYp);//k plot
            if(NINT(par[10][1])==1){
                PLOTline1bar2(1,0,iColor,14,mwl,wwl,e1,ErrXp,ErrYp);//e1 plot
                PLOTline1bar2(1,0,iColor,15,mwl,wwl,e2,ErrXp,ErrYp);//e2 plot
            }
            if(DATO[1]==2)
                PLOTline1bar2(1,0,iColor,1,mwl,wwl,tn,ErrXp,ErrYp);//Tn plot
            if(DATO[2]==2)
                PLOTline1bar2(1,0,iColor,2,mwl,wwl,rp,ErrXp,ErrYp);//Tp plot
            if(DATO[3]==2)
                PLOTline1bar2(1,0,iColor,3,mwl,wwl,rn,ErrXp,ErrYp);//Rn plot
            if(DATO[4]==2)
                PLOTline1bar2(1,0,iColor,4,mwl,wwl,rp,ErrXp,ErrYp);//Rp plot
            if(DATO[5]==2)
                PLOTline1bar2(1,0,iColor,5,mwl,wwl,r1,ErrXp,ErrYp);//R1 plot
            for(int ieli=0;ieli<4;ieli++){
                if(DATO[7+2*ieli]==2){
                    int Nperiod=NINT(delDiff[ieli]/Dperiod);
                    for(int i=0;i<NeV;i++)
                        delta[i]=DelMat[i][ieli]-Nperiod*Dperiod;
                    PLOTline1bar2(1,0,iColor,7,mwl,wwl,delta,ErrXp,ErrYp);//DELTA plot
                }
                if(DATO[8+2*ieli]==2){
                    for(int i=0;i<NeV;i++)
                        psi[i]=PsiMat[i][ieli];
                    PLOTline1bar2(1,0,iColor,8,mwl,wwl,psi,ErrXp,ErrYp);//PSI plot
                }
            }
        }

        //save Job# in /temp and set njobBest
        if(ieon==2 || (ieon==0 && (r=="f" || r=="i" || r=="fsem"))){
            StoreFitSet();
            double chi2Best= ui->lineEdit_chi2Best->text().toDouble();
            if(chi2fin<chi2Best){
                ui->lineEdit_chi2Best->setText(QString::number(chi2fin));
                ui->spinBox_jobBest->setValue(jobtot);
            }
        }

        if(ieon==0)
            r="X";
    }
//    if(p){
//        delete [] p;
//        p=nullptr;
//    }
}


void ksemawc::StoreFitSet(){
    jobtot++;
    QString fileStorecurrent=fileStore;
    fileStore=pathroot+"temp/semaw"+QString::number(jobtot)+".Spj";
    SaveSetting(3);
    fileStore=fileStorecurrent;
    CNK[1][1]=0;//because it was modified along the saving
    ui->sB_PAR_8_1->setValue(jobtot);
    ui->sB_PAR_8_2->setValue(jobtot);
}


void ksemawc::GoBest(){
    int jobBest=ui->spinBox_jobBest->value();
    if(jobBest==0)
        return;
    QString SpjName=ui->lineEdit_P -> text();
    QString infoFnk=ui->lineEdit_infoFnk-> text();
    QString ftmp=pathroot+"temp/semaw"+QString::number(jobBest)+".Spj";
    occupyPF=10;
    caller="GoBest";
    ReadSetting(ftmp);
    ui->sB_PAR_8_2->setValue(jobBest);
    ui->sB_PAR_8_1->setValue(jobtot);
    ui->lineEdit_P -> setText(SpjName);
    ui->lineEdit_infoFnk-> setText(infoFnk);
}

void ksemawc::GoPrevious(){
    int jobview=ui->sB_PAR_8_2->value();
    int jobtot=ui->sB_PAR_8_1->value();
    if(jobview<=1)
        return;
    QString SpjName=ui->lineEdit_P -> text();
    QString infoFnk=ui->lineEdit_infoFnk-> text();
    jobview--;
    QString ftmp=pathroot+"temp/semaw"+QString::number(jobview)+".Spj";
    occupyPF=10;
    caller="GoPrevious";
    ReadSetting(ftmp);
    ui->sB_PAR_8_1->setValue(jobtot);
    ui->sB_PAR_8_2->setValue(jobview);
    ui->lineEdit_P -> setText(SpjName);
    ui->lineEdit_infoFnk-> setText(infoFnk);
}

void ksemawc::GoNext(){
    int jobview=ui->sB_PAR_8_2->value();
    int jobtot=ui->sB_PAR_8_1->value();
    if(jobview>=jobtot)
        return;
    QString SpjName=ui->lineEdit_P -> text();
    QString infoFnk=ui->lineEdit_infoFnk-> text();
    jobview++;
    QString ftmp=pathroot+"temp/semaw"+QString::number(jobview)+".Spj";
    caller="GoNext";
    ReadSetting(ftmp);
    ui->sB_PAR_8_1->setValue(jobtot);
    ui->sB_PAR_8_2->setValue(jobview);
    ui->lineEdit_P -> setText(SpjName);
    ui->lineEdit_infoFnk-> setText(infoFnk);
}


void ksemawc::SPADA(){
    printf("-> SPADA ....\n");
    int N,MR,Ndati,ilinrim,NANG;
    double instr[21],X[12000],Y[12000],Z[12000],EPR[3][12000][3],
            wmin,wmax,dmin,dmax,Div,rstep,DLAM;
    QString ST1,ST2,DIS[7],EST[8],ER,SPA,SAL,fnam,line,line2,line3,pezzo;
    // matrix initialization for measures and temporary solutions
    SOL[1][1]=0.;//none solution
    for(int I=1;I<=17;I++){
        for(int J=1;J<=NeV;J++){
            MIS[I][J][1]=0.;
            if(I<=15)
                MIS[I][J][2]=1.;
            else
                MIS[I][J][2]=0.;
        }
    }
    for(int I=1;I<=8;I++){
        for(int J=1;J<=NeV;J++){
            ELI[I][J][1]=0.;
            ELI[I][J][2]=1.;
        }
    }
    rxy[19][3]=10.;//cos(Delta)
    rxy[19][4]=-10.;
    rxy[22][3]=10.e+36;//tan(Psi)
    rxy[22][4]=-10.;
    int IUVIR=NINT(par[24][1]);
    if(IUVIR==1) ST2=".v";
    if(IUVIR==2) ST2=".i";
    DIS[1]=".tn";
    DIS[2]=".tp";
    DIS[3]=".rn";
    DIS[4]=".rp";
    DIS[5]=".r1";
    DIS[6]=".an";
    for(int i=1;i<=14;i++){
        if(DATO[i]>=1. && i<=6){
            //printf("SPADA-> DATO[%d]=%d\n",i,DATO[i]);
            ST1=QString::number(NINT(par[20+i][4]));
            EST[i]=ST2+ST1+DIS[i];
        }
    }
    if(!NANK[9].contains("mate/aa999")){
        printf("\tSPADA-> call LoadFilenk with file-solutions %s\n",NANK[9].toStdString().c_str());
        LoadFilenk();
    }
    for(int i=1;i<=20;i++)
        instr[i]=par[i][3];
    if(instr[6]<.5)
        ER="y";
    else
        ER=" ";
    SPA="L";
    if(NINT(rxy[25][4])==2)
        SPA="E";
    SAL="n";
    MR=NINT(par[22][1]);
    int nrif=NINT(par[10][3]);
    dmin=par[4][1];
    dmax=par[4][2];
    int STEP=1;
    // set WL and baseline error
    for(int L=1;L<=NeV;L++){
        if(SPA=="L")
            MIS[7][L][1]=dmin+double(L-1)/(NeV-1)*(dmax-dmin);
        else
            MIS[7][L][1]=1./(1./dmin+double(L-1)/(NeV-1)*(1./dmax-1./dmin));
        MIS[6][L][1]=1.;//reflectance correction
        MIS[7][L][2]=instr[5];
        if(ER!="y"){
            if(IUVIR==1){
                if(MIS[7][L][1]>184.9 && MIS[7][L][1]<860.81)
                    MIS[7][L][2]=.0005+.838951/(1+pow((MIS[7][L][1]-183.5)/0.7,2.))+.0026/(1.+pow((MIS[7][L][1]-870.)/25.,4.));
                if(MIS[7][L][1]>3100. && MIS[7][L][1]<3200.1)
                    MIS[7][L][2]=-4.66E-04+.00476/(1+pow((MIS[7][L][1]-3500.)/200.,2.));
            }
            if(MIS[7][L][2]<instr[5]) MIS[7][L][2]=instr[5];
        }
    }
    //file-nk load
    for(int I=1;I<=8;I++){
        fnam=pathroot+NANK[I].simplified();
        if(!fnam.contains("mate/aa999.9")){
            fnam=fnam+".nk";
            printf("\tSPADA-> load file-nk %s\n",fnam.toStdString().c_str());
            QFile file(fnam);
            if(!file.open(QIODevice::ReadOnly | QIODevice::Text)){
                msgErrLoad("SPADA-fnam",fnam);
                printf("SPADA-> ERROR opening file-nk= %s\n",fnam.toStdString().c_str());
                continue;
            }
            QTextStream stream (&file);
            line = stream.readLine();//info first line
            line = stream.readLine();//Ndat and "nm" if new file.nk
            line=line.simplified();
            QStringList List0;
            List0=line.split(" ");
            pezzo=List0.at(0).toLocal8Bit().constData();
            Ndati=pezzo.toInt();
            printf("\tNdata= %d\n",Ndati);
            double A2nm=0.1;
            if(line.contains("nm"))
                A2nm=1;
            if(Ndati>12000){
                QMessageBox msgBox;
                msgBox.setText("ATTENTION: please reduce your data number below 12000!!!\n");
                msgBox.setStandardButtons(QMessageBox::Ok);
                msgBox.exec();
                return;
            }
            for(int J=1;J<=Ndati;J++){
                line=stream.readLine();
                line=line.simplified();
                QStringList List;
                List =line.split(" ");
                int nV=List.count();
//                if(J==1){
//                    printf("nV=%d\n",nV);
//                    for(int k=0;k<nV;k++)
//                        printf("%d->%f ",k,List.at(k).toDouble());
//                    printf("\n");
//                }
                if(nV!=5){
                    List =line.split("\t");
                    nV=List.count();
                    if(nV!=5){
                        printf("nV=%d line=%s\n",nV,line.toStdString().c_str());
                        printf("SPADA-> ERROR reading file-nk= %s: separator char invalid\n",fnam.toStdString().c_str());
                        QMessageBox msgBox;
                        msgBox.setText("LoadFilenk-> ERROR reading file-nk="+fnam+"\n: separator char invalid or ERRn and ERRk missing");
                        msgBox.setStandardButtons(QMessageBox::Ok);
                        msgBox.exec();
                        continue;
                    }
                }
                X[J]=List.at(0).toDouble();//wl
                X[J]=X[J]*A2nm;
                Y[J]=List.at(1).toDouble();//n
                Z[J]=List.at(2).toDouble();//k
                //cout << X[J]<<"\t"<<Y[J]<<"\t"<<Z[J]<<"\n";
            }
            file.close();
            if(X[1]<X[Ndati]){
                N=1;
                STEP=1;
            }
            else{
                N=Ndati;
                STEP=-1;
            }
            CONVER(X,Y,Ndati,N,STEP,7+I,1);
            CONVER(X,Z,Ndati,N,STEP,7+I,2);
        }
    }
    //file-SF load
    if(NINT(par[23][1])==1){
        for(int i=0;i<=6;i++){
            if(i==0){
                if(MR!=0){
                    QFile file(fRefMir);
                    if(file.open(QIODevice::ReadOnly | QIODevice::Text)){
                        QTextStream stream (&file);
                        for(int irif=1;irif<=nrif;irif++){
                            line = stream.readLine();
                            line2 = stream.readLine();
                        }
                        file.close();
                        fnam=pathroot+line2.simplified();
                    }
                    else
                        msgErrLoad("SPADA-fRefMir",fRefMir);
                }
                else
                    continue;
            }
            else if(DATO[i]>=1){
                fnam=pathroot+NANK[11].simplified();
                fnam=fnam+EST[i];
            }
            else
                continue;
            QFile file(fnam);
            if(!file.open(QIODevice::ReadOnly | QIODevice::Text)){
                msgErrLoad("SPADA-fnam bis",fnam);
                printf("SPADA-> i=%d ERROR opening file-SF= %s\n",i,fnam.toStdString().c_str());
                continue;
            }
            else{
                printf("\tSPADA-> i=%d loading fileSF %s\n",i,fnam.toStdString().c_str());
                QTextStream stream (&file);
                line = stream.readLine();
                line2 = stream.readLine();
                printf("line= %s\n",line.toStdString().c_str());
                printf("line2= %s\n",line2.toStdString().c_str());
                if(line.contains("PE UV")){// file Perkin Elmer L900 - L950
                    printf("\tfile Perkin Elmer L900 - L950\n");
                    Div=100.;
                    for(int ir=1;ir<=11;ir++)
                        line = stream.readLine();
                    line3 = stream.readLine();
                    if(line3.contains("PerkinElmer UV WinLab 5"))
                        ilinrim=66;
                    else
                        ilinrim=70;
                    for(int ir=1;ir<=ilinrim;ir++)
                        line = stream.readLine();
                    stream >> rstep;
                    STEP=NINT(rstep);
                    stream >> Ndati;
                    N=Ndati;
                    printf("Ndati=%d\n",Ndati);
                    if(Ndati>12000){
                        QMessageBox msgBox;
                        msgBox.setText("ATTENTION: please reduce your data number below 12000!!!\n");
                        msgBox.setStandardButtons(QMessageBox::Ok);
                        msgBox.exec();
                        return;
                    }
                    do{
                        line = stream.readLine();
                    }while(!line.contains("#DATA"));
                    for(int j=1;j<=Ndati;j++){
                        stream >> X[j] >> Y[j];
                        //printf("X[%d]=%f\tY[%d]=%f\n",j,X[j],j,Y[j]);
                    }
                    fflush(stdout);
                }
                else if(line2.contains("#####SCALED")){//SCALED file
                    printf("\tSCALED file\n");
                    int ncolo;
                    double xin,xfi,zz;
                    stream >> xin;
                    stream >> xfi;
                    stream >> Ndati;
                    if(Ndati>12000){
                        QMessageBox msgBox;
                        msgBox.setText("ATTENTION: please reduce your data number below 12000!!!\n");
                        msgBox.setStandardButtons(QMessageBox::Ok);
                        msgBox.exec();
                        return;
                    }
                    if(xin<xfi){
                        STEP=1;
                        N=1;
                    }
                    else{
                        STEP=-1;
                        N=Ndati;
                    }
                    if(line2.contains("#####SCALED ")){
                        printf("\ttype: #####SCALED \n");
                        cout << "Number of column to load (2,..,6)?\n";
                        cin >> ncolo;
                        if(ncolo<2) ncolo=2;
                        if(ncolo>6) ncolo=6;
                        for(int j=1;j<=Ndati;j++){
                            if(ncolo==2) stream >> X[j] >> Y[j];
                            if(ncolo==3) stream >> X[j] >> zz >> Y[j];
                            if(ncolo==4) stream >> X[j] >> zz >> zz >> Y[j];
                            if(ncolo==5) stream >> X[j] >> zz >> zz >> zz >> Y[j];
                            if(ncolo==6) stream >> X[j] >> zz >> zz >> zz >> zz >> Y[j];
                        }
                        cout << "Absolute Unit (0) or not (1)?\n";
                        cin >> ncolo;
                        if(ncolo==0)
                            Div=1;
                        else
                            Div=100.;
                    }
                    else if(line2.contains("#####SCALEDA")){
                        printf("\ttype: #####SCALEDA\n");
                        for(int j=1;j<=Ndati;j++)
                            stream >> X[j] >> Y[j];
                        Div=1.;
                    }
                    else if(line2.contains("#####SCALED%")){
                        printf("\ttype: #####SCALED%%\n");
                        for(int j=1;j<=Ndati;j++){
                            stream >> X[j] >> Y[j];
                            //cout<< j<<"\t"<<X[j]<<"\t"<<Y[j]<<"\n";
                        }
                        Div=100.;
                    }
                }
                else if(line2.contains("RTmethod[")){
                    // file Transmittance by VASE
                    printf("\ttype VASE transmittance\n");
                    line = stream.readLine();
                    line = stream.readLine();
                    int iL0E1=0;
                    if(line.contains("eV"))
                        iL0E1=1;
                    QStringList List;
                    int j=0;
                    do{
                        j++;
                        line = stream.readLine();
                        line=line.simplified();
                        List=line.split(" ");
                        pezzo=List.at(1).toLocal8Bit().constData();
                        X[j]=pezzo.toDouble();
                        if(iL0E1==1)
                            X[j]=1240./X[j];
                        pezzo=List.at(3).toLocal8Bit().constData();
                        Y[j]=pezzo.toDouble();
                        //printf("\t X[%d]=%f Y[%d]=%f\n",j,X[j],j,Y[j]);
                    }while(!stream.atEnd());
                    Ndati=j;
                    printf("\tNdat=%d\n",Ndati);
                    Div=1.;
                    STEP=1;
                    N=1;
                }
                else if(line2.contains("##XYUNITS= W")){//Lamba9 or old FTIR
                    printf("\ttype: Lambda9 or old FTIR\n");
                    line = stream.readLine();
                    line = stream.readLine();
                    stream >> line >> DLAM >> ST1>> Div;
                    printf("DLAM= %f Div= %f\n",DLAM,Div);
                    Div=Div*100.;
                    stream >> line >> Ndati>>ST1;
                    printf("Ndati= %d\n",Ndati);
                    if(Ndati>12000){
                        QMessageBox msgBox;
                        msgBox.setText("ATTENTION: please reduce your data number below 12000!!!\n");
                        msgBox.setStandardButtons(QMessageBox::Ok);
                        msgBox.exec();
                        return;
                    }
                    do{
                        line = stream.readLine();
                        line=line.simplified();
                    }while(line.isEmpty());
                    int NRIGHE=int(Ndati/14.);
                    double DE=14.*((Ndati/14.)-NRIGHE);
                    if(DE>.5) NRIGHE++;
                    int K;
                    double X1;
                    for(int J=1;J<=NRIGHE;J++){
                        K=(J-1)*14;
                        line=stream.readLine();
                        line=line.simplified();
                        QStringList List;
                        List =line.split(" ");
                        int nV=List.count();
                        //printf("\nnV=%d\n",nV);
                        QString pezzo;
                        for(int iv=0;iv<nV;iv++){
                            pezzo=List.at(iv).toLocal8Bit().constData();
                            //printf("pezzo[%d]= %s\t",iv,pezzo.toStdString().c_str());
                            if(iv==0){
                                X1=pezzo.toDouble();
                                //printf("\n");
                            }
                            else{
                                X[K+iv]=X1+(iv-1)*DLAM;
                                Y[K+iv]=pezzo.toDouble();
                                //printf("X[%d]= %f\tY[%d]= %f\n",K+iv,X[K+iv],K+iv,Y[K+iv]);
                            }
                        }
                    }
                    if(IUVIR==1){
                        N=Ndati;
                        STEP=-1;
                    }
                    else{
                        N=1;
                        STEP=1;
                    }
                }
                else{
                    line3=stream.readLine();
                    int K;
                    double xf,yf,X1;
                    if(line3.contains("##DATA TYPE= UL")){//new file Lamba9)
                        printf("\ttype: new Lambda9\n");
                        for(int J=1;J<=6;J++)
                            line = stream.readLine();
                        stream >> line >> xf;
                        stream >> line >> yf;
                        stream >> line >> wmin;
                        stream >> line >> wmax;
                        stream >> line >> Ndati;
                        if(Ndati>12000){
                            QMessageBox msgBox;
                            msgBox.setText("ATTENTION: please reduce your data number below 12000!!!\n");
                            msgBox.setStandardButtons(QMessageBox::Ok);
                            msgBox.exec();
                            return;
                        }
                        //printf("xf= %f\tyf= %f\twmin=%f\twmax=%f Ndati=%d\n",xf,yf,wmin,wmax,Ndati);
                        line = stream.readLine();
                        line = stream.readLine();
                        line = stream.readLine();
                        DLAM=(wmax-wmin)/double(Ndati-1);
                        Div=1./yf;
                        int NRIGHE=int(Ndati/10.);
                        double DE=10.*((Ndati/10.)-int(Ndati/10.));
                        if(DE>.5) NRIGHE++;
                        int K;
                        double X1;
                        for(int J=1;J<=NRIGHE;J++){
                            K=(J-1)*10;
                            line=stream.readLine();
                            line=line.simplified();
                            //cout<<line.toStdString()<<"\n";
                            QStringList List;
                            List =line.split(" ");
                            int nV=List.count();
                            //printf("\nnV=%d\n",nV);
                            QString pezzo;
                            for(int iv=0;iv<nV;iv++){
                                pezzo=List.at(iv).toLocal8Bit().constData();
                                //printf("pezzo[%d]= %s\t",iv,pezzo.toStdString().c_str());
                                if(iv==0){
                                    X1=pezzo.toDouble();
                                    //printf("\n");
                                }
                                else{
                                    X[K+iv]=X1*xf+(iv-1)*DLAM;
                                    Y[K+iv]=pezzo.toDouble();
                                    //printf("X[%d]= %f\tY[%d]= %f\n",K+iv,X[K+iv],K+iv,Y[K+iv]);
                                }
                            }
                        }
                        N=1;
                        STEP=1;
                    }
                    else if(line3.contains("##DATA TYPE= UV")){//file Lambda19
                        printf("\ttype: Lambda19\n");
                        for(int J=1;J<=9;J++)
                            line = stream.readLine();
                        stream>>line>>xf;
                        stream>>line>>yf;
                        stream>>line>>wmax;
                        stream>>line>>wmin;
                        stream>>line>>Ndati;
                        if(Ndati>12000){
                            QMessageBox msgBox;
                            msgBox.setText("ATTENTION: please reduce your data number below 12000!!!\n");
                            msgBox.setStandardButtons(QMessageBox::Ok);
                            msgBox.exec();
                            return;
                        }
                        for(int i1=1;i1<=4;i1++)
                            line = stream.readLine();
                        DLAM=(wmax-wmin)/(Ndati-1);
                        Div=100./yf;
                        int NRIGHE=int(Ndati/5.);
                        double DE=5.*((Ndati/5.)-int(Ndati/5.));
                        if(DE>.5) NRIGHE++;
                        int K;
                        double X1;
                        for(int J=1;J<=NRIGHE;J++){
                            K=(J-1)*5;
                            line=stream.readLine();
                            line=line.simplified();
                            QStringList List;
                            List =line.split(" ");
                            int nV=List.count();
                            //printf("\nnV=%d\n",nV);
                            QString pezzo;
                            for(int iv=0;iv<nV;iv++){
                                pezzo=List.at(iv).toLocal8Bit().constData();
                                //printf("pezzo[%d]= %s\t",iv,pezzo.toStdString().c_str());
                                if(iv==0){
                                    X1=pezzo.toDouble();
                                    //printf("\n");
                                }
                                else{
                                    X[K+iv]=X1*xf-(iv-1)*DLAM;
                                    Y[K+iv]=pezzo.toDouble();
                                    //printf("X[%d]= %f\tY[%d]= %f\n",K+iv,X[K+iv],K+iv,Y[K+iv]);
                                }
                            }
                        }
                        N=Ndati;
                        STEP=-1;
                    }
                    else if(line3.contains("##DATA TYPE= IN")){//new IR file
                        printf("\ttype: new IR\n");
                        for(int J=1;J<=9;J++)
                            line = stream.readLine();
                        stream>>line>>DLAM;
                        line = stream.readLine();
                        line = stream.readLine();
                        stream>>line>>xf;
                        stream>>line>>yf;
                        stream>>line>>wmin;
                        stream>>line>>wmax;
                        stream>>line>>Ndati;
                        if(Ndati>12000){
                            QMessageBox msgBox;
                            msgBox.setText("ATTENTION: please reduce your data number below 12000!!!\n");
                            msgBox.setStandardButtons(QMessageBox::Ok);
                            msgBox.exec();
                            return;
                        }
                        for(int i1=1;i1<=4;i1++)
                            line = stream.readLine();
                        Div=1./yf;
                        K=0;
                        while(K<Ndati){
                            line = stream.readLine();
                            if(line.at(0)!=" "){// FTIR1760
                                QStringList List;
                                List =line.split("+");
                                int nV=List.count();
                                QString pezzo;
                                for(int iv=0;iv<nV;iv++){
                                    pezzo=List.at(iv).toLocal8Bit().constData();
                                    K++;
                                    if(iv==0)
                                        X[K]=pezzo.toDouble()*xf;
                                    else{
                                        X[K]=X[K-1]+DLAM;
                                        Y[K]=pezzo.toDouble()*yf;
                                    }
                                }
                            }
                            else{//FTIR SPECTRUM GS
                                if(K==0) Div=Div*100.;
                                stream>>X1>>Y[K+1]>>Y[K+2]>>Y[K+3]>>Y[K+4];
                                for(int iiii=1;iiii<=5;iiii++)
                                    X[K+iiii-1]=X1*xf+(iiii-1)*DLAM;
                                K=K+4;
                            }
                        }
                        N=1;
                        STEP=1;
                    }
                }
                file.close();
                if(Ndati>12000){
                    QMessageBox msgBox;
                    msgBox.setText("ATTENTION: please reduce your data number below 12000!!!\n");
                    msgBox.setStandardButtons(QMessageBox::Ok);
                    msgBox.exec();
                    return;
                }
                for(int L=1;L<=Ndati;L++){
                    if(IUVIR==2) X[L]=1./X[L]*1.E7;//cm^-1 -> nm
                    Y[L]=Y[L]/Div;
                }
                if(i==0){
                    CONVER(X,Y,Ndati,N,STEP,6,IUVIR);
                    rxy[6][1]=0.;
                    rxy[6][2]=1.;
                    rxy[6][3]=0.;
                    rxy[6][4]=1.;
                }
                else{
                    iw=0;
                    //printf("X[1]=%f\n",X[1]);
                    CONVER(X,Y,Ndati,N,STEP,i,IUVIR);
                    idToLineEdit["DP_RXY_"+QString::number(i)+"_"+QString::number(3)]-> setText(QString::number(rxy[i][3]));
                    idToLineEdit["DP_RXY_"+QString::number(i)+"_"+QString::number(4)]-> setText(QString::number(rxy[i][4]));
                }
            }
        }
    }
    //file ELI load
    if(NINT(par[23][2])==1){
        //int DP0cosDtanP1=ui->comboBox_DeltaPsiScale->currentIndex();
        for(int J=1;J<=4;J++){
            if(DATO[5+2*J]>0 && !NANK[11+J].contains("mate/aa999.9")){
                printf("%s\n",NANK[11+J].toStdString().c_str());
                fnam=pathroot+NANK[11+J].simplified()+".el";
                printf("\tSPADA-> J=%d loading fileELI %s\n",J,fnam.toStdString().c_str());
                QFile file(fnam);
                if(!file.open(QIODevice::ReadOnly | QIODevice::Text)){
                    msgErrLoad("SPADA-fnam tris",fnam);
                    printf("SPADA-> ERROR opening file-ELI= %s\n",fnam.toStdString().c_str());
                    continue;
                }
                else{
                    int iL0E1=0;
                    QString pezzo;
                    QTextStream stream (&file);
                    line = stream.readLine();
                    cout<<line.toStdString()<<"\n";
                    line=stream.readLine();
                    line=line.simplified();
                    if(!line.contains("VASE")){
                        QStringList List0;
                        List0 =line.split(" ");
                        pezzo=List0.at(0).toLocal8Bit().constData();
                        NANG=pezzo.toInt();
                        pezzo=List0.at(1).toLocal8Bit().constData();
                        Ndati=pezzo.toInt();
                        line2=List0.at(2);
                        for(int I=1;I<=NANG;I++){
                            do{
                                line = stream.readLine();
                            }while(line.isEmpty());
                            cout<<line.toStdString()<<"\n";
                        }
                    }
                    else{
                        line2="VASE";
                        line=stream.readLine();
                        line=stream.readLine();
                        if(line.contains("eV"))
                            iL0E1=1;
                        Ndati=0;
                        NANG=EliTab[J][0][0];
                        for(int k=0;k<NANG;k++)
                            Ndati=Ndati+EliTab[J][k][1];
                    }
                    if(Ndati>12000){
                        QMessageBox msgBox;
                        msgBox.setText("ATTENTION: please reduce your data number below 12000!!!\n");
                        msgBox.setStandardButtons(QMessageBox::Ok);
                        msgBox.exec();
                        return;
                    }
                    printf("Nang= %d  NdatiEL= %d  st4= %s\n",NANG,Ndati,line2.toStdString().c_str());
                    int NDEL=0;
                    int iformat=0;
                    double TE,LAM,DE,PS,DDE,DPS;
                    if(line2.contains("tldp"))
                        iformat=1;
                    else if(line2.contains("etpd"))
                        iformat=2;
                    else if(line2.contains("ltpd"))//elia
                        iformat=3;
                    else if(line2.contains("etcab")){
                        iformat=4;
                        TE=line.section('\t', 0, 0).toDouble();
                    }
                    else if(line2.contains("VASE")){
                        iformat=5;
                    }
                    if(iformat==0){
                        cout<<"ATTENTION: set the format of file ELI!!!"<<"\n";
                        return;
                    }
                    printf("iformat=%d par[13+J][1]=%f\n",iformat,par[13+J][1]);
                    double theta2load=idToComboBox["cBteE"+QString::number(J)]->currentText().toDouble();
                    for(int I=1;I<=Ndati;I++){
                        if(iformat==1){
                            stream>>TE>>LAM>>DE>>PS>>DDE>>DPS;
                        }
                        else if(iformat==2){
                            stream>>LAM>>TE>>PS>>DE>>DPS>>DDE;
                            LAM=1240./LAM;
                        }
                        else if(iformat==3)
                            stream>>LAM>>TE>>PS>>DE>>DPS>>DDE;
                        else if(iformat==4){
                            stream>>LAM>>PS>>DE>>DPS>>DDE;
                            LAM=1240./LAM;
                            //printf("tan(PS)= %f",PS);
                            PS=atan(PS)/deg2rad;
                            //printf(" -> PS = %f(rad) = %f(deg)\n",PS*deg2rad,PS);
                            DE=acos(DE)/deg2rad;
                            DDE=0.04;//err Delta (deg)
                            DPS=0.02;//err Psi (deg)
                        }
                        else if(iformat==5){
                            line=stream.readLine();
                            line=line.simplified();
                            QStringList List;
                            List=line.split(" ");
                            pezzo=List.at(0).toLocal8Bit().constData();
                            if(pezzo.contains("depolE"))//the line does not contain DELTA and PSI angles
                                continue;
                            pezzo=List.at(1).toLocal8Bit().constData();
                            LAM=pezzo.toDouble();
                            if(iL0E1==1)
                                LAM=1240./LAM;
                            pezzo=List.at(2).toLocal8Bit().constData();
                            TE=pezzo.toDouble();
                            pezzo=List.at(3).toLocal8Bit().constData();
                            PS=pezzo.toDouble();
                            pezzo=List.at(4).toLocal8Bit().constData();
                            DE=pezzo.toDouble();
                            pezzo=List.at(5).toLocal8Bit().constData();
                            DPS=pezzo.toDouble();
                            pezzo=List.at(6).toLocal8Bit().constData();
                            DDE=pezzo.toDouble();
                        }
                        //if(abs(TE-par[13+J][1])<=1.){//0.001
                        if(abs(TE-theta2load)<=1.){
                            //cout<<TE<<"\t"<<LAM<<"\t"<<DE<<"\t"<<PS<<"\t"<<DDE<<"\t"<<DPS<<"\n";
                            NDEL=NDEL+1;
                            X[NDEL]=LAM;
                            EPR[1][NDEL][1]=DE;
                            EPR[1][NDEL][2]=DDE;
                            EPR[2][NDEL][1]=PS;
                            EPR[2][NDEL][2]=DPS;
                        }
                    }
                    file.close();
                    printf("Theta= %f  NdatEl= %d\n",theta2load,NDEL);
                    if(NDEL==0)
                        return;
                    for(int I=1;I<=2;I++){
                        for(int II=1;II<=2;II++){
                            for(int L=1;L<=NDEL;L++)
                                Y[L]=EPR[I][L][II];
                            CONVER(X,Y,NDEL,1,1,20+2*(J-1)+I,II);
                            if(II==1){
                                idToLineEdit["DP_RXY_"+QString::number(4+2*J+I)+"_"+QString::number(3)]-> setText(QString::number(rxy[4+2*J+I][3]));
                                idToLineEdit["DP_RXY_"+QString::number(4+2*J+I)+"_"+QString::number(4)]-> setText(QString::number(rxy[4+2*J+I][4]));
                            }
                        }
                    }
                }
            }
        }
        setRangeEli();
    }
    par[38][3]=par[4][1];
    par[38][4]=par[4][2];
    par[21][1]=NeV;
    par[28][2]=NeV;
    if(NINT(par[9][1])==1)
        ui->checkBox_setPsoK -> setCheckState ( Qt::Checked );
    printf("exit SPADA\n");
}




void CONVER(double X[12000],double Y[12000],int NDATI,int N1,int STEP,int I,int IUVIR){
    //int iw=1;//verbose (0 silent, 1 verbose)
    iw=NINT(par[18][1]);
    printf("CONVER: NDATI=%d N1=%d STEP=%d I=%d IUVIR=%d\n",NDATI,N1,STEP,I,IUVIR);
    int ALARM=0;
    int ns[1000];
    double xs[500],ys[500];
    double ymin=1.e+36;
    double ymax=1.e-36;
    double cosDelMin=10.;
    double cosDelMax=-10.;
    double tanPsiMin=1.e+36;
    double tanPsiMax=-1.;
    double INSTR[21];
    double WL,WLa,WLb;
    double yOld;
    for(int j=1;j<=20;j++){
        INSTR[j]=par[j][3];
    }
    for(int H=1;H<=NeV;H++){
        WL=MIS[7][H][1];
        double FACO=1.;
        if((I==3 || I==4 || I==5) && IUVIR==1){
            if(NINT(INSTR[10])>0)
                FACO=FACO*MIS[6][H][1];
        }

        // points to average
        if(H>1)
            WLa=WL-(MIS[7][H][1]-MIS[7][H-1][1])/2.;
        else
            WLa=WL-(MIS[7][H+1][1]-MIS[7][H][1])/2.;
        if(H<NeV)
            WLb=WL+(MIS[7][H+1][1]-MIS[7][H][1])/2.;
        else
            WLb=WL+(MIS[7][H][1]-MIS[7][H-1][1])/2.;
        int N=N1;
        while(X[N]<WLa && (N+STEP)<=NDATI && (N+STEP)>=1)
            N=N+STEP;
        if(iw==1){
            printf("WLa= %f  WLb= %f\n",WLa,WLb);
            printf("N= %d X[%d]= %f\n",N,N,X[N]);
            printf("Data in [WLa , WLb]:\n");
        }
        int NPM=0;
        int J=0;
        while(X[N]>=WLa && X[N]<=WLb && N<=NDATI && N>=1){
            NPM=NPM+1;
            ns[J]=N;
            xs[J]=X[N];
            ys[J]=Y[N];
            if(iw==1)
                cout<<xs[J]<<"\t"<<ys[J]<<"\n";
            J++;
            N=N+STEP;
            if(N>NDATI ||N<1) continue;
        }

        // check
        if(NPM<3){
            if(iw==1)
                printf("... NPM= %d' < 3 -> set NPM=3\n",NPM);
            if(NPM>=1)
                N=ns[0];
            else if(NPM==0){
                NPM=1;
                ns[0]=N;
            }
            if(NPM==1 && N>1 && N<NDATI && X[N+STEP]>WLb)
                N=N-STEP;
            while((N+STEP*2)<1 || (N+STEP*2)>NDATI)
                N=N-STEP;
            NPM=3;
            for(int j=0;j<NPM;j++){
                ns[j]=N;
                xs[j]=X[N];
                ys[j]=Y[N];
                if(iw==1)
                    cout<<xs[j]<<"\t"<<ys[j]<<"\n";
                N=N+STEP;
            }
        }
        else if(NPM>3 && I>=21 && NINT(par[10][2])==0){ //ellissometric measurement
            if(iw==1)
                cout<<"... ELLI meas -> set NPM=3!"<<"\n";
            while(X[N+STEP]<WL)
                N=N+STEP;
            NPM=3;
            for(int j=0;j<NPM;j++){
                ns[j]=N;
                xs[j]=X[N];
                ys[j]=Y[N];

                if(iw==1)
                    cout<<xs[j]<<"\t"<<ys[j]<<"\n";
                N=N+STEP;
            }
        }
        // check and cure ELLI meas jump
        if(I>=21){
            double ysOld;
            for(int J=0;J<NPM;J++){
                if(J==0)
                    ysOld=ys[J];
                ys[J]=ys[J]-Dperiod*NINT((ys[J]-ysOld)/Dperiod);
                ysOld=ys[J];
            }
        }

        // fit/interpolation
        double a,b,c;
        if(NPM>3){//fit
            int mFit=NPM;
            int nFit=3;
            int lwa=mFit*nFit+5*nFit+mFit;
            int iwa[nFit];
            //int* iwa=nullptr;
            //iwa = new int[nFit];
            double p[nFit];
            //double* p=nullptr;
            //p = new double[nFit];
            double fvec[mFit];
            //double* fvec=nullptr;
            //fvec = new double[mFit];
            double wa[lwa];
            //double* wa=nullptr;
            //wa = new double[lwa];
            p[2]=0.;
            for(int j=0;j<NPM;j++){
                ARSE[j][1]=xs[j];
                ARSE[j][2]=ys[j];
                p[2]=p[2]+ys[j];
            }
            p[0]=0.;
            p[1]=0.;
            p[2]=p[2]/double(NPM);
            double tol=sqrt(dpmpar(1));
            int Info=lmdif1(FPAR, &pTF2, mFit, nFit, p,fvec, tol, iwa, wa, lwa);
            if(Info<=0)
                printf("improper input parameters\n");
            c=p[2];
            b=p[1];
            a=p[0];
            if(iw==1){
                printf("coeff of BF parabolic:\na= %f  b= %f  c= %f\n",a,b,c);
                for(int iii=0;iii<NPM;iii++)
                    cout<< xs[iii]<<"\t"<<ys[iii]<<"\n";
            }
        }
        else{//interpolation on 3 points
            if(abs(xs[0]-xs[1])<.01 || abs(xs[1]-xs[2])<.01)
                printf("ATTENTION: 2 data with same abscidssa %f %f %f\n!!!",
                       xs[0],xs[1],xs[2]);
            double XD1=xs[0];
            double XD2=xs[1];
            double XD3=xs[2];
            double YD1=ys[0];
            double YD2=ys[1];
            double YD3=ys[2];
            double DETD=XD1*XD1*(XD2-XD3)-XD1*(XD2*XD2-XD3*XD3);
            DETD=DETD+XD2*XD2*XD3-XD3*XD3*XD2;
            a=YD1*(XD2-XD3)-XD1*(YD2-YD3);
            a=(a+YD2*XD3-YD3*XD2)/DETD;
            b=XD1*XD1*(YD2-YD3)-YD1*(XD2*XD2-XD3*XD3);
            b=(b+XD2*XD2*YD3-XD3*XD3*YD2)/DETD;
            c=XD1*XD1*(XD2*YD3-XD3*YD2);
            c=c-XD1*(XD2*XD2*YD3-XD3*XD3*YD2);
            c=(c+YD1*(XD2*XD2*XD3-XD3*XD3*XD2))/DETD;
            if(iw==1 ){
                printf("Coeff parabola at 3 points:\n");
                printf("X: %f %f %f\n",XD1,XD2,XD3);
                printf("Y: %f %f %f\n",YD1,YD2,YD3);
                printf("a= %f  b= %f  c= %f\n",a,b,c);
                printf("FACO=%f \n", FACO);
            }
        }
        double YFIN=(a*WL*WL+b*WL+c)*FACO;
        if(I>=7 && I<=15){
            if(YFIN<.0 && ALARM==0 && NINT(par[9][1])!=1){
                if(IUVIR==1){
                    QMessageBox msgBox;
                    msgBox.setText("n-file "+NANK[I-7].simplified()+" is <0 !!!\n Do you want impose n>=0?");
                    msgBox.setStandardButtons(QMessageBox::Yes | QMessageBox::No);
                    msgBox.setDefaultButton(QMessageBox::Yes);
                    int ret = msgBox.exec();
                    switch (ret) {
                    case QMessageBox::Yes:
                        par[9][1]=1;
                        break;
                    }
                    printf("n-file %s is <0 !!!\n",NANK[I-7].simplified().toStdString().c_str());
                }
                else if(IUVIR==2){
                    QMessageBox msgBox;
                    msgBox.setText("k-file "+NANK[I-7].simplified()+" is <0 !!!\n Do you want impose k>=0?");
                    msgBox.setStandardButtons(QMessageBox::Yes | QMessageBox::No);
                    msgBox.setDefaultButton(QMessageBox::Yes);
                    int ret = msgBox.exec();
                    switch (ret) {
                    case QMessageBox::Yes:
                        par[9][1]=1;
                        break;
                    }
                    printf("k-file %s is <0 !!!\n",NANK[I-7].simplified().toStdString().c_str());
                }
                ALARM=1;
            }
            if(int(par[9][1])==1 && YFIN<.0)
                YFIN=.0;
            MIS[I][H][IUVIR]=YFIN;
        }
        if(I<=6){
            MIS[I][H][1]=YFIN;
            double DL=MIS[7][H][2];
            if(I==1 || I==2)
                MIS[I][H][2]=MIS[I][H][1]*INSTR[3]+DL;
            if(I==3 || I==4 || I==5)
                MIS[I][H][2]=MIS[I][H][1]*(INSTR[3]+INSTR[4])+DL;
            if(iw==1)
                printf("MIS[%d][%d][1]= %f\n",I,H,MIS[I][H][1]);
        }
        if(I>=21 && I<=28){
            if(iDelCon==1&&(I==21||I==23||I==25||I==27)){
                if(H==1)
                    yOld=YFIN;
                YFIN=YFIN-Dperiod*NINT((YFIN-yOld)/Dperiod);
                yOld=YFIN;
            }
            if(I==21||I==23||I==25||I==27){
                cosDelMin=min(cosDelMin,cos(YFIN*deg2rad));
                cosDelMax=max(cosDelMax,cos(YFIN*deg2rad));
            }
            else{
                tanPsiMin=min(tanPsiMin,tan(YFIN*deg2rad));
                tanPsiMax=max(tanPsiMax,tan(YFIN*deg2rad));
            }
            ELI[I-20][H][IUVIR]=YFIN;
        }
        ymin=min(ymin,YFIN);
        ymax=max(ymax,YFIN);
        if(iw==1){
            char ch;
            cout<<"Go on?\n";
            cin>> ch;
        }
    }
    if(I>=1 && I<=6){
        rxy[I][3]=ymin*100.;
        rxy[I][4]=ymax*100.;
    }
    else if(I>=21 && I<=28 && IUVIR==1){
        rxy[I-14][3]=ymin;
        rxy[I-14][4]=ymax;
        if(I==21||I==23||I==25||I==27){
            rxy[19][3]=min(rxy[19][3],cosDelMin);
            rxy[19][4]=max(rxy[19][4],cosDelMax);
        }
        else{
            rxy[22][3]=min(rxy[22][3],tanPsiMin);
            rxy[22][4]=max(rxy[22][4],tanPsiMax);
        }
    }
}


int FPAR(void *p, int m, int n, const double *x, double *fvec, int iflag){
    //struct pointToFit2 *pTF2 = (struct pointToFit2 *)p;
    for(int i=0;i<m;i++)
        fvec[i]=ARSE[i][2]-x[0]*ARSE[i][1]*ARSE[i][1]-x[1]*ARSE[i][1]-x[2];
    return(0);
}


void PLOTline1bar2(int iL1B2,int iRD,int iCol,int ic,int Ndata,double *Xp,double *Yp,double *ErrXp,double *ErrYp){
    //int iL1B2= 1->line 2->dot & errorBar
    //int iRD= 0->overwrite 1->redraw
    //int iCol= color index 0->black
    //                      1->blue
    //                      2->cyan
    //                      3->green
    //                      4->magenta
    //                      5->red
    //                      6->yellow
    //int ic identify the graph
    //    ic= 1: G1_Tn
    //    ic= 2: G2_Tp
    //    ic= 3: G3_Rn
    //    ic= 4: G4_Rp
    //    ic= 5: G5_R1
    //    ic= 6: G6_Apds
    //    ic= 7: G7_D(elta)
    //    ic= 8: G8_P(si)
    //    ic= 9: G9_Af
    //    ic=10: G10_tTR
    //    ic=11: G11_nk
    //    ic=12: G12_wn
    //    ic=13: G13_wk
    //    ic=14: G14_we1
    //    ic=15: G15_we2
    //    ic=16: G16_Ab
    int graphDrift=40;
    double Px[Ndata],Py[Ndata],Ymin=rxy[17][1],Ymax=rxy[17][2];
    for(int i=0;i<Ndata;i++){
        Px[i]=Xp[i];
        Py[i]=Yp[i];
    }
    if(ic==13 && NINT(par[31][5]==1)){//klog plot
        if(Ymin>0.)
            Ymin=log10(Ymin);
        else
            Ymin=log10(Ymax/1.E+04);
        Ymax=log10(Ymax);
        for(int i=0;i<Ndata;i++){
            double Yup=Yp[i]+ErrYp[i];
            double Ydw=Yp[i]-ErrYp[i];
            if(Yup<Ymin)
                Yup=Ymin;
            if(Ydw<Ymin)
                Ydw=Ymin;
            Py[i]=0.5*(log10(Yup)+log10(Ydw));
            ErrYp[i]=0.5*(log10(Yup)-log10(Ydw));
        }
    }
    if(ic==7 && NINT(par[33][5]==1) ){
        for(int i=0;i<Ndata;i++){
            double Yup=cos((Yp[i]+ErrYp[i])*deg2rad);
            double Ydw=cos((Yp[i]-ErrYp[i])*deg2rad);
            Py[i]=0.5*(Yup+Ydw);
            ErrYp[i]=0.5*abs(Yup-Ydw);
        }
    }
    if(ic==8 && NINT(par[33][5]==1) ){
        for(int i=0;i<Ndata;i++){
            double Yup=tan((Yp[i]+ErrYp[i])*deg2rad);
            double Ydw=tan((Yp[i]-ErrYp[i])*deg2rad);
            Py[i]=0.5*(Yup+Ydw);
            ErrYp[i]=0.5*abs(Yup-Ydw);
        }
    }
    if(ic<1 ||ic>16)//if(Ndata==0 || ic<1 ||ic>16)
        return;
    int wGhx=int(rxy[24][1]);
    int wGhy=int(rxy[24][2]);
    printf("PLOT: IXW[%d]=%d iL1B2=%d iRD=%d iCol=%d ic=%d Ndata=%d L1E2=%d\n",ic,IXW[ic],iL1B2,iRD,iCol,ic,Ndata,L1E2);
    if(IXW[ic]<0){
        if(ic==1){
            G1_Tn=new QwtPlot();
            G1_Tn->setGeometry(nOpenGraph*graphDrift,nOpenGraph*graphDrift,wGhx,wGhy);
            nOpenGraph++;
            QwtPlotPicker *m_picker1 = new QwtPlotPicker( QwtAxis::XBottom, QwtAxis::YLeft,
                                                         QwtPlotPicker::CrossRubberBand, QwtPicker::AlwaysOn,
                                                         G1_Tn->canvas() );
        }
        if(ic==2){
            G2_Tp=new QwtPlot();
            G2_Tp->setGeometry(nOpenGraph*graphDrift,nOpenGraph*graphDrift,wGhx,wGhy);
            nOpenGraph++;
            QwtPlotPicker *m_picker2 = new QwtPlotPicker( QwtAxis::XBottom, QwtAxis::YLeft,
                                                         QwtPlotPicker::CrossRubberBand, QwtPicker::AlwaysOn,
                                                         G2_Tp->canvas() );
        }
        if(ic==3){
            G3_Rn=new QwtPlot();
            G3_Rn->setGeometry(nOpenGraph*graphDrift,nOpenGraph*graphDrift,wGhx,wGhy);
            nOpenGraph++;
            QwtPlotPicker *m_picker3 = new QwtPlotPicker( QwtAxis::XBottom, QwtAxis::YLeft,
                                                         QwtPlotPicker::CrossRubberBand, QwtPicker::AlwaysOn,
                                                         G3_Rn->canvas() );
        }
        if(ic==4){
            G4_Rp=new QwtPlot();
            G4_Rp->setGeometry(nOpenGraph*graphDrift,nOpenGraph*graphDrift,wGhx,wGhy);
            nOpenGraph++;
            QwtPlotPicker *m_picker4 = new QwtPlotPicker( QwtAxis::XBottom, QwtAxis::YLeft,
                                                         QwtPlotPicker::CrossRubberBand, QwtPicker::AlwaysOn,
                                                         G4_Rp->canvas() );
        }
        if(ic==5){
            G5_R1=new QwtPlot();
            G5_R1->setGeometry(nOpenGraph*graphDrift,nOpenGraph*graphDrift,wGhx,wGhy);
            nOpenGraph++;
            QwtPlotPicker *m_picker5 = new QwtPlotPicker( QwtAxis::XBottom, QwtAxis::YLeft,
                                                         QwtPlotPicker::CrossRubberBand, QwtPicker::AlwaysOn,
                                                         G5_R1->canvas() );
        }
        if(ic==6){
            G6_Apds=new QwtPlot();
            G6_Apds->setGeometry(nOpenGraph*graphDrift,nOpenGraph*graphDrift,wGhx,wGhy);
            nOpenGraph++;
            QwtPlotPicker *m_picker6 = new QwtPlotPicker( QwtAxis::XBottom, QwtAxis::YLeft,
                                                         QwtPlotPicker::CrossRubberBand, QwtPicker::AlwaysOn,
                                                         G6_Apds->canvas() );
        }
        if(ic==7){
            G7_D=new QwtPlot();
            G7_D->setGeometry(nOpenGraph*graphDrift,nOpenGraph*graphDrift,wGhx,wGhy);
            nOpenGraph++;
            QwtPlotPicker *m_picker7 = new QwtPlotPicker( QwtAxis::XBottom, QwtAxis::YLeft,
                                                         QwtPlotPicker::CrossRubberBand, QwtPicker::AlwaysOn,
                                                         G7_D->canvas() );
        }
        if(ic==8){
            G8_P=new QwtPlot();
            G8_P->setGeometry(nOpenGraph*graphDrift,nOpenGraph*graphDrift,wGhx,wGhy);
            nOpenGraph++;
            QwtPlotPicker *m_picker8 = new QwtPlotPicker( QwtAxis::XBottom, QwtAxis::YLeft,
                                                         QwtPlotPicker::CrossRubberBand, QwtPicker::AlwaysOn,
                                                         G8_P->canvas() );
        }
        if(ic==9){
            G9_Af=new QwtPlot();
            G9_Af->setGeometry(nOpenGraph*graphDrift,nOpenGraph*graphDrift,wGhx,wGhy);
            nOpenGraph++;
            QwtPlotPicker *m_picker9 = new QwtPlotPicker( QwtAxis::XBottom, QwtAxis::YLeft,
                                                         QwtPlotPicker::CrossRubberBand, QwtPicker::AlwaysOn,
                                                         G9_Af->canvas() );
        }
        if(ic==10){
            G10_tTR=new QwtPlot();
            G10_tTR->setGeometry(nOpenGraph*graphDrift,nOpenGraph*graphDrift,wGhx,wGhy);
            nOpenGraph++;
            QwtPlotPicker *m_picker10 = new QwtPlotPicker( QwtAxis::XBottom, QwtAxis::YLeft,
                                                         QwtPlotPicker::CrossRubberBand, QwtPicker::AlwaysOn,
                                                         G10_tTR->canvas() );
        }
        if(ic==11){
            G11_nk=new QwtPlot();
            G11_nk->setGeometry(nOpenGraph*graphDrift,nOpenGraph*graphDrift,wGhx,wGhy);
            nOpenGraph++;
            QwtPlotPicker *m_picker11 = new QwtPlotPicker( QwtAxis::XBottom, QwtAxis::YLeft,
                                                         QwtPlotPicker::CrossRubberBand, QwtPicker::AlwaysOn,
                                                         G11_nk->canvas() );
        }
        if(ic==12){
            G12_wn=new QwtPlot();
            G12_wn->setGeometry(nOpenGraph*graphDrift,nOpenGraph*graphDrift,wGhx,wGhy);
            nOpenGraph++;
            QwtPlotPicker *m_picker12 = new QwtPlotPicker( QwtAxis::XBottom, QwtAxis::YLeft,
                                                         QwtPlotPicker::CrossRubberBand, QwtPicker::AlwaysOn,
                                                         G12_wn->canvas() );
        }
        if(ic==13){
            G13_wk=new QwtPlot();
            G13_wk->setGeometry(nOpenGraph*graphDrift,nOpenGraph*graphDrift,wGhx,wGhy);
            nOpenGraph++;
            QwtPlotPicker *m_picker13 = new QwtPlotPicker( QwtAxis::XBottom, QwtAxis::YLeft,
                                                         QwtPlotPicker::CrossRubberBand, QwtPicker::AlwaysOn,
                                                         G13_wk->canvas() );
        }
        if(ic==14){
            G14_we1=new QwtPlot();
            G14_we1->setGeometry(nOpenGraph*graphDrift,nOpenGraph*graphDrift,wGhx,wGhy);
            nOpenGraph++;
            QwtPlotPicker *m_picker14 = new QwtPlotPicker( QwtAxis::XBottom, QwtAxis::YLeft,
                                                         QwtPlotPicker::CrossRubberBand, QwtPicker::AlwaysOn,
                                                         G14_we1->canvas() );
        }
        if(ic==15){
            G15_we2=new QwtPlot();
            G15_we2->setGeometry(nOpenGraph*graphDrift,nOpenGraph*graphDrift,wGhx,wGhy);
            nOpenGraph++;
            QwtPlotPicker *m_picker15 = new QwtPlotPicker( QwtAxis::XBottom, QwtAxis::YLeft,
                                                         QwtPlotPicker::CrossRubberBand, QwtPicker::AlwaysOn,
                                                         G15_we2->canvas() );
        }
        if(ic==16){
            G16_Ab=new QwtPlot();
            G16_Ab->setGeometry(nOpenGraph*graphDrift,nOpenGraph*graphDrift,wGhx,wGhy);
            nOpenGraph++;
            QwtPlotPicker *m_picker16 = new QwtPlotPicker( QwtAxis::XBottom, QwtAxis::YLeft,
                                                         QwtPlotPicker::CrossRubberBand, QwtPicker::AlwaysOn,
                                                         G16_Ab->canvas() );
        }
        IXW[ic]=1;
    }
    else{
        if(ic==1 && iRD==1){
            G1_Tn->detachItems(QwtPlotItem::Rtti_PlotCurve, true);
            G1_Tn->detachItems(QwtPlotItem::Rtti_PlotIntervalCurve,true);
        }
        if(ic==2 && iRD==1){
            G2_Tp->detachItems(QwtPlotItem::Rtti_PlotCurve, true);
            G2_Tp->detachItems(QwtPlotItem::Rtti_PlotIntervalCurve,true);
        }
        if(ic==3 && iRD==1){
            G3_Rn->detachItems(QwtPlotItem::Rtti_PlotCurve, true);
            G3_Rn->detachItems(QwtPlotItem::Rtti_PlotIntervalCurve,true);
        }
        if(ic==4 && iRD==1){
            G4_Rp->detachItems(QwtPlotItem::Rtti_PlotCurve, true);
            G4_Rp->detachItems(QwtPlotItem::Rtti_PlotIntervalCurve,true);
        }
        if(ic==5 && iRD==1){
            G5_R1->detachItems(QwtPlotItem::Rtti_PlotCurve, true);
            G5_R1->detachItems(QwtPlotItem::Rtti_PlotIntervalCurve,true);
        }
        if(ic==6 && iRD==1){
            G6_Apds->detachItems(QwtPlotItem::Rtti_PlotCurve, true);
            G6_Apds->detachItems(QwtPlotItem::Rtti_PlotIntervalCurve,true);
        }
        if(ic==7 && iRD==1){
            G7_D->detachItems(QwtPlotItem::Rtti_PlotCurve, true);
            G7_D->detachItems(QwtPlotItem::Rtti_PlotIntervalCurve,true);
        }
        if(ic==8 && iRD==1){
            G8_P->detachItems(QwtPlotItem::Rtti_PlotCurve, true);
            G8_P->detachItems(QwtPlotItem::Rtti_PlotIntervalCurve,true);
        }
        if(ic==9 && iRD==1){
            G9_Af->detachItems(QwtPlotItem::Rtti_PlotCurve, true);
            G9_Af->detachItems(QwtPlotItem::Rtti_PlotIntervalCurve,true);
        }
        if(ic==10 && iRD==1){
            G10_tTR->detachItems(QwtPlotItem::Rtti_PlotCurve, true);
            G10_tTR->detachItems(QwtPlotItem::Rtti_PlotIntervalCurve,true);
        }
        if(ic==11 && iRD==1){
            G11_nk->detachItems(QwtPlotItem::Rtti_PlotCurve, true);
            G11_nk->detachItems(QwtPlotItem::Rtti_PlotIntervalCurve,true);
        }
        if(ic==12 && iRD==1){
            G12_wn->detachItems(QwtPlotItem::Rtti_PlotCurve,true);
            G12_wn->detachItems(QwtPlotItem::Rtti_PlotIntervalCurve,true);
        }
        if(ic==13 && iRD==1){
            G13_wk->detachItems(QwtPlotItem::Rtti_PlotCurve,true);
            G13_wk->detachItems(QwtPlotItem::Rtti_PlotIntervalCurve,true);
        }
        if(ic==14 && iRD==1){
            G14_we1->detachItems(QwtPlotItem::Rtti_PlotCurve,true);
            G14_we1->detachItems(QwtPlotItem::Rtti_PlotIntervalCurve,true);
        }
        if(ic==15 && iRD==1){
            G15_we2->detachItems(QwtPlotItem::Rtti_PlotCurve,true);
            G15_we2->detachItems(QwtPlotItem::Rtti_PlotIntervalCurve,true);
        }
        if(ic==16 && iRD==1){
            G16_Ab->detachItems(QwtPlotItem::Rtti_PlotCurve, true);
            G16_Ab->detachItems(QwtPlotItem::Rtti_PlotIntervalCurve,true);
        }
    }
    double WEmin=rxy[20][1],WEmax=rxy[20][2];
    if(L1E2==2 && ic!=10 && ic!=11){
        WEmin=1240./WEmin;
        WEmax=1240./WEmax;
        for(int i=0;i<Ndata;i++)
            Px[i]=1240./Px[i];
    }
    QwtPlotCurve *dataPlot=new QwtPlotCurve(labelQwt);
    //QwtPlotCurve *dataPlot=new QwtPlotCurve("Layer_"+QString::number(iCol));
    //QwtPlotCurve *dataPlot=new QwtPlotCurve("dataPlot");
    QwtPlotIntervalCurve *range_plot = new QwtPlotIntervalCurve("range");
    dataPlot->setSamples(Px, Py, Ndata);
    int wT=int(rxy[24][3]);
    if(iL1B2==1){
        dataPlot->setPen(QPen(myColor[iCol],wT,Qt::DotLine));
    }
    else{
        dataPlot->setSymbol(new QwtSymbol(QwtSymbol::Ellipse, Qt::NoBrush, QPen(myColor[iCol]), QSize(2, 2)));
        dataPlot->setStyle(QwtPlotCurve::NoCurve);
        dataPlot->setRenderHint(QwtPlotItem::RenderAntialiased);
        /* error bars */
        QVector<QwtIntervalSample> range(Ndata);
        for(int i = 0; i < Ndata; i++) {
            range[i] = QwtIntervalSample(Px[i],Py[i]-ErrYp[i],Py[i]+ErrYp[i]);
        }
        QwtIntervalSymbol *errorbar = new QwtIntervalSymbol(QwtIntervalSymbol::Bar);
        errorbar->setPen(QPen(myColor[iColor], wT));
        errorbar->setWidth(wT);
        range_plot->setSamples(range);
        range_plot->setSymbol(errorbar);
        range_plot->setStyle(QwtPlotIntervalCurve::NoCurve);
    }
    if(ic==1){
        G1_Tn -> setAxisTitle(0,"Tnormal (%)");
        if(L1E2==1)
            G1_Tn -> setAxisTitle(2,"wavelength (nm)");
        else
            G1_Tn -> setAxisTitle(2,"photon energy (eV)");
        dataPlot->attach(G1_Tn);
        if(iL1B2==2)
            range_plot->attach(G1_Tn);
        //printf("Tn: Range X: %f %f\nRangeY: %f %f\n",WEmin,WEmax,rxy[1][1],rxy[1][2]);
        G1_Tn->setAxisScale(0,rxy[1][1],rxy[1][2],0);
        G1_Tn->setAxisScale(2,WEmin,WEmax,0);
        G1_Tn->setAutoReplot();
        G1_Tn->show();
    }
    else if(ic==2){
        G2_Tp -> setAxisTitle(0,"Tpolarised (%)");
        if(L1E2==1)
            G2_Tp -> setAxisTitle(2,"wavelength (nm)");
        else
            G2_Tp -> setAxisTitle(2,"photon energy (eV)");
        dataPlot->attach(G2_Tp);
        if(iL1B2==2)
            range_plot->attach(G2_Tp);
        G2_Tp->setAxisScale(0,rxy[2][1],rxy[2][2],0);
        G2_Tp->setAxisScale(2,WEmin,WEmax,0);
        G2_Tp->setAutoReplot();
        G2_Tp->show();
    }
    else if(ic==3){
        G3_Rn -> setAxisTitle(0,"Rnormal (%)");
        if(L1E2==1)
            G3_Rn -> setAxisTitle(2,"wavelength (nm)");
        else
            G3_Rn -> setAxisTitle(2,"photon energy (eV)");
        dataPlot->attach(G3_Rn);
        if(iL1B2==2)
            range_plot->attach(G3_Rn);
        G3_Rn->setAxisScale(0,rxy[3][1],rxy[3][2],0);
        G3_Rn->setAxisScale(2,WEmin,WEmax,0);
        G3_Rn -> setAutoReplot();
        G3_Rn->show();
    }
    else if(ic==4){
        G4_Rp -> setAxisTitle(0,"Rpolarized (%)");
        if(L1E2==1)
            G4_Rp -> setAxisTitle(2,"wavelength (nm)");
        else
            G4_Rp -> setAxisTitle(2,"photon energy (eV)");
        dataPlot->attach(G4_Rp);
        if(iL1B2==2)
            range_plot->attach(G4_Rp);
        G4_Rp->setAxisScale(0,rxy[4][1],rxy[4][2],0);
        G4_Rp->setAxisScale(2,WEmin,WEmax,0);
        G4_Rp -> setAutoReplot();
        G4_Rp->show();
    }
    else if(ic==5){
        G5_R1 -> setAxisTitle(0,"R1normal (%)");
        if(L1E2==1)
            G5_R1 -> setAxisTitle(2,"wavelength (nm)");
        else
            G5_R1 -> setAxisTitle(2,"photon energy (eV)");
        dataPlot->attach(G5_R1);
        if(iL1B2==2)
            range_plot->attach(G5_R1);
        G5_R1->setAxisScale(0,rxy[5][1],rxy[5][2],0);
        G5_R1->setAxisScale(2,WEmin,WEmax,0);
        G5_R1 -> setAutoReplot();
        G5_R1->show();
    }
    else if(ic==6){
        G6_Apds -> setAxisTitle(0,"Apds");
        if(L1E2==1)
            G6_Apds -> setAxisTitle(2,"wavelength (nm)");
        else
            G6_Apds -> setAxisTitle(2,"photon energy (eV)");
        dataPlot->attach(G6_Apds);
        if(iL1B2==2)
            range_plot->attach(G6_Apds);
        G6_Apds->setAxisScale(0,rxy[6][1],rxy[6][2],0);
        G6_Apds->setAxisScale(2,WEmin,WEmax,0);
        G6_Apds -> setAutoReplot();
        G6_Apds->show();
    }
    else if(ic==7){
        if( NINT(par[33][5]==0) )
            G7_D -> setAxisTitle(0,"Delta (deg)");
        else
            G7_D -> setAxisTitle(0,"cos(Delta)");
        if(L1E2==1)
            G7_D -> setAxisTitle(2,"wavelength (nm)");
        else
            G7_D -> setAxisTitle(2,"photon energy (eV)");
        dataPlot->attach(G7_D);
        if(iL1B2==2)
            range_plot->attach(G7_D);
        if( NINT(par[33][5]==0) )
            G7_D->setAxisScale(0,rxy[7][1],rxy[7][2],0);
        else
            G7_D->setAxisScale(0,rxy[19][1],rxy[19][2],0);
        G7_D->setAxisScale(2,WEmin,WEmax,0);
        G7_D -> setAutoReplot();
        G7_D->show();
    }
    else if(ic==8){
        if( NINT(par[33][5]==0) )
            G8_P -> setAxisTitle(0,"Psi (deg)");
        else
            G8_P -> setAxisTitle(0,"tan(Psi)");
        if(L1E2==1)
            G8_P -> setAxisTitle(2,"wavelength (nm)");
        else
            G8_P -> setAxisTitle(2,"photon energy (eV)");
        dataPlot->attach(G8_P);
        if(iL1B2==2)
            range_plot->attach(G8_P);
        if( NINT(par[33][5]==0) )
            G8_P->setAxisScale(0,rxy[8][1],rxy[8][2],0);
        else
            G8_P->setAxisScale(0,rxy[22][1],rxy[22][2],0);
        G8_P->setAxisScale(2,WEmin,WEmax,0);
        G8_P -> setAutoReplot();
        G8_P->show();
    }
    else if(ic==9){
        G9_Af -> setAxisTitle(0,"Absorptance_front (%)");
        if(L1E2==1)
            G9_Af -> setAxisTitle(2,"wavelength (nm)");
        else
            G9_Af -> setAxisTitle(2,"photon energy (eV)");
        dataPlot->attach(G9_Af);
        if(iL1B2==2)
            range_plot->attach(G9_Af);
        G9_Af->setAxisScale(0,rxy[18][1],rxy[18][2],0);
        G9_Af->setAxisScale(2,WEmin,WEmax,0);
        if(labelQwt.contains("layer"))
            G9_Af->insertLegend( new QwtLegend(), QwtPlot::RightLegend );
        //else if(G9_Af->legend()->)
        //    G9_Af->legend()->setVisible(false);
        G9_Af -> setAutoReplot();
        G9_Af->show();
    }
    else if(ic==10){
        G10_tTR -> setAxisTitle(0,"Tau Rho Rho1 (%)");
        G10_tTR -> setAxisTitle(2,"Theta (deg)");
        dataPlot->attach(G10_tTR);
        if(iL1B2==2)
            range_plot->attach(G10_tTR);
        G10_tTR->setAxisScale(0,rxy[23][1],rxy[23][2],0);
        G10_tTR->setAxisScale(2,rxy[21][1],rxy[21][2],0);
        G10_tTR->setAutoReplot();
        G10_tTR->show();
    }
    else if(ic==11){
        G11_nk -> setAxisTitle(0,"k");
        G11_nk -> setAxisTitle(2,"n");
        dataPlot->attach(G11_nk);
        if(iL1B2==2)
            range_plot->attach(G11_nk);
        G11_nk->setAxisScale(0,rxy[17][1],rxy[17][2],0);
        G11_nk->setAxisScale(2,rxy[16][1],rxy[16][2],0);
        G11_nk -> setAutoReplot();
        G11_nk->show();
    }
    else if(ic==12){
        G12_wn->setAxisTitle(0,"n");
        if(L1E2==1)
            G12_wn->setAxisTitle(2,"wavelength (nm)");
        else
            G12_wn->setAxisTitle(2,"photon energy (eV)");
        dataPlot->attach(G12_wn);
        if(iL1B2==2)
            range_plot->attach(G12_wn);
        G12_wn->setAxisScale(0,rxy[16][1],rxy[16][2],0);
        G12_wn->setAxisScale(2,WEmin,WEmax,0);
        G12_wn->setAutoReplot();
        G12_wn->show();
    }
    else if(ic==13){
        if(NINT(par[31][5]!=1))
            G13_wk->setAxisTitle(0,"k");
        else
            G13_wk->setAxisTitle(0,"log10(k)");
        if(L1E2==1)
            G13_wk->setAxisTitle(2,"wavelength (nm)");
        else
            G13_wk->setAxisTitle(2,"photon energy (eV)");
        dataPlot->attach(G13_wk);
        if(iL1B2==2)
            range_plot->attach(G13_wk);
        G13_wk->setAxisScale(0,Ymin,Ymax,0);
        G13_wk->setAxisScale(2,WEmin,WEmax,0);
        G13_wk->setAutoReplot();
        G13_wk->show();
    }
    else if(ic==14){
        G14_we1->setAxisTitle(0,"eps1");
        if(L1E2==1)
            G14_we1->setAxisTitle(2,"wavelength (nm)");
        else
            G14_we1->setAxisTitle(2,"photon energy (eV)");
        dataPlot->attach(G14_we1);
        if(iL1B2==2)
            range_plot->attach(G14_we1);
        G14_we1->setAxisScale(0,rxy[26][1],rxy[26][2],0);
        G14_we1->setAxisScale(2,WEmin,WEmax,0);
        G14_we1->setAutoReplot();
        G14_we1->show();
    }
    else if(ic==15){
        G15_we2->setAxisTitle(0,"eps2");
        if(L1E2==1)
            G15_we2->setAxisTitle(2,"wavelength (nm)");
        else
            G15_we2->setAxisTitle(2,"photon energy (eV)");
        dataPlot->attach(G15_we2);
        if(iL1B2==2)
            range_plot->attach(G15_we2);
        G15_we2->setAxisScale(0,rxy[27][1],rxy[27][2],0);
        G15_we2->setAxisScale(2,WEmin,WEmax,0);
        G15_we2->setAutoReplot();
        G15_we2->show();
    }
    else if(ic==16){
        G16_Ab -> setAxisTitle(0,"Absorptance_back (%)");
        if(L1E2==1)
            G16_Ab -> setAxisTitle(2,"wavelength (nm)");
        else
            G16_Ab -> setAxisTitle(2,"photon energy (eV)");
        dataPlot->attach(G16_Ab);
        if(iL1B2==2)
            range_plot->attach(G16_Ab);
        G16_Ab->setAxisScale(0,rxy[18][1],rxy[18][2],0);
        G16_Ab->setAxisScale(2,WEmin,WEmax,0);
        if(labelQwt.contains("layer"))
            G16_Ab->insertLegend( new QwtLegend(), QwtPlot::RightLegend );
        G16_Ab -> setAutoReplot();
        G16_Ab->show();
    }
}





void COSVNK(double VNK[17][3],int L){
    // L = index of wavelenght
    //printf("->COSVNK iWL=%d\n",L);
    double f2;
    for(int J=1;J<=16;J++){
        int i1=NINT(CNK[J][1]);
        SETVNK(i1,J,VNK,L,2);
        int i2=NINT(CNK[J][4]);
        //if(L==iwl2print) printf("-> COSVNK @ iWL=%d J=%d i1=%d i2=%d\n",L,J,i1,i2);
        if(i2>=0){
            f2=pm[70+J][1];
            double n1=VNK[J][1];
            double k1=VNK[J][2];
            SETVNK(i2,J,VNK,L,5);
            double n2=VNK[J][1];
            double k2=VNK[J][2];
            //if(L==iwl2print) printf("\t n2=%f k2=%f\n",n1,k2);
            EMA(n1,k1,n2,k2,f2);
            VNK[J][1]=Nema;
            VNK[J][2]=Kema;
        }
        //if(L==iwl2print) printf(" VNK[%d][1]=%f VNK[%d][2]=%f\n",J,VNK[J][1],J,VNK[J][2]);
    }
}


void SETVNK(int io,int J,double VNK[17][3],int L,int Kcte){
    // L = index of wavelenght
    // set VNK(J,1) and VNK(J,2) given CNK(J,k)
    if(io==0){
        VNK[J][1]=CNK[J][Kcte];
        VNK[J][2]=CNK[J][Kcte+1];
    }
    else if(io>=1 && io<=7){ // 1-7 Fit options
        double WL=MIS[7][L][1];
        double ev=1240./WL;
        FDISP(io,ev);
        VNK[J][1]=sqn;
        VNK[J][2]=sqk;
    }
    else if(io>=8 &&io<=15){
        VNK[J][1]=MIS[io][L][1];
        VNK[J][2]=MIS[io][L][2];
    }
    else if(io==16){
        VNK[J][1]=VNK[1][1];
        VNK[J][2]=VNK[1][2];
    }
    //if(L==iwl2print) printf("L= %d n=%f k=%f io=%d\n",L,VNK[J][1],VNK[J][2],io);
}


void EMA(double N,double K,double NA,double KA,double FA){
    // Effective Medium Approximation computing
    double NE,KE;
    if(FA<0.0001){
        NE=N;
        KE=K;
    }
    else{
        if(FA>=0.0001 && FA<=0.9999){
            complex<double> EA(NA*NA-KA*KA,-2.*NA*KA);
            complex<double> EB(N*N-K*K,-2.*N*K);
            double SEGNO=(1.-FA)*imag(EB)+FA*imag(EA);
            complex<double> EE=(1.-FA)*(EA-2.*EB)+FA*(EB-2.*EA);
            complex<double> E=(-EE+sqrt(EE*EE+8.*EA*EB))/4.;
            double ER=real(E);
            double EI=imag(E);
            if(EI*SEGNO<0.){
                E=(-EE-sqrt(EE*EE+8.*EA*EB))/4.;
                ER=real(E);
                EI=imag(E);
            }
            KE=sqrt((-ER+sqrt(ER*ER+EI*EI))*.5)*copysign(1.,-EI);
            NE=sqrt(ER+KE*KE);
        }
        else{
            NE=NA;
            KE=KA;
        }
    }
    Nema=NE;
    Kema=KE;
}


void FDISP(int iopt,double eV){
    // nk computing by oscillators
    //printf("->FDISP @ E=%f\n",eV);
    double cln2=0.693147181;
    double sqrln2=0.832554611;
    long double E3, Ex, Ey, W, KTi, r, Eb, Eh, El, derexp, derparab, TUx, TUy;

    // epsilon and nk by sum of oscillators
    int noscFT=NINT(pf[iopt][1]);
    long double E=eV;
    long double epr=0.;
    long double epi=0.;
    long double eprj,epij;
    for(int io=2;io<=noscFT+1;io++){
        int i=NINT(pf[iopt][io]);
        int ifu=NINT(pm[101+(i-1)*5][1]);
        long double C=abs(pm[102+(i-1)*5][1]);
        long double E0=abs(pm[103+(i-1)*5][1]);
        long double D=abs(pm[104+(i-1)*5][1]);
        long double K=abs(pm[105+(i-1)*5][1]);
        eprj=0.;
        epij=0.;
        if(ifu==1){// classic oscillator
            long double den=pow(E0*E0-E*E,2.)+pow(2.*E*D,2.);
            eprj=2.*C/PIG*(E0*E0-E*E)/den;
            epij=2.*C/PIG*2.*E*D/den;//K was deleted
        }
        else if(ifu==2){// homogeneous quantum oscillator
            long double x=(E0-E)/D;
            eprj=C/(PIG*D*E0)*x/(1+x*x);
            epij=C/(PIG*D*E0)/(1+x*x);//K was deleted
        }
        else if(ifu==3){//inomogeneous quantum oscillator
            long double x=(E0-E)/D;
            long double Kinho=pow((D/sqrln2),2.)/2*exp(-pow(E0*sqrln2/D,2.))+
                                E0*D/sqrln2/2.*sqrt(PIG)*erfcl(-E0*sqrln2/D);
            eprj=C/Kinho*DAWS(sqrln2*x)*2./sqrt(PIG);
            epij=C/Kinho*exp(-cln2*x*x);//K was deleted
        }
        else if(ifu==4){
            eprj=C*C;// only in this case C=n_cte !
            epij=0.;
        }
        else if(ifu==5){// drude model
            eprj=-2.*C/PIG/(D*D+E*E);
            epij= 2.*C/PIG*D/(pow(E,3.)+E*D*D);
        }
        else if(ifu==6){// indirect Gap Cody
            double De=K;
            double KCi=E0*De+De*De/2.;
            eprj=C/KCi*(2.*ICCODYRE(E,E0+3.*De/4.,D)-2.*ICCODYRE(E,E0+De/4.,D)+
                              16./De/De*(I2CODYRE(E,E0,E0+De/4.,D)-I2CODYRE(E,E0,E0,D)-
                                               I2CODYRE(E,E0+De/2.,E0+3.*De/4.,D)+
                                               I2CODYRE(E,E0+De/2.,E0+De/4.,D)+I2CODYRE(E,E0+De,E0+De,D)-
                                               I2CODYRE(E,E0+De,E0+3.*De/4.,D)));

            epij=C/KCi*(2.*ICCODYIM(E,E0+3.*De/4.,D)-2.*ICCODYIM(E,E0+De/4.,D)+
                              16./De/De*(I2CODYIM(E,E0,E0+De/4.,D)-I2CODYIM(E,E0,E0,D)-
                                               I2CODYIM(E,E0+De/2.,E0+3.*De/4.,D)+
                                               I2CODYIM(E,E0+De/2.,E0+De/4.,D)+I2CODYIM(E,E0+De,E0+De,D)-
                                               I2CODYIM(E,E0+De,E0+3.*De/4.,D)));
        }
        else if(ifu==7){// indirect Gap Tauc
            double De=K;
            double r=E0/De;
            double KTi=2.*(8.*r*r*log((4.*r+1.)/(4.*r)) + 8.*pow((r+1.),2.)*log((4.*r+4.)/(4.*r+3.))
                               - (1.+8.*r*(r+1.))*log((4.*r+3.)/(4.*r+1.)) );
            eprj=C/KTi*(2.*ICTAUCRE(E,E0+3.*De/4.,D)-2.*ICTAUCRE(E,E0+De/4.,D)+
                              16./De/De*(I2TAUCRE(E,E0,E0+De/4.,D)-I2TAUCRE(E,E0,E0,D)-
                                               I2TAUCRE(E,E0+De/2.,E0+3.*De/4.,D)+
                                               I2TAUCRE(E,E0+De/2.,E0+De/4.,D)+I2TAUCRE(E,E0+De,E0+De,D)-
                                               I2TAUCRE(E,E0+De,E0+3.*De/4.,D)));

            epij=C/KTi*(2.*ICTAUCIM(E,E0+3.*De/4.,D)-2.*ICTAUCIM(E,E0+De/4.,D)+
                              16./De/De*(I2TAUCIM(E,E0,E0+De/4.,D)-I2TAUCIM(E,E0,E0,D)-
                                               I2TAUCIM(E,E0+De/2.,E0+3.*De/4.,D)+
                                               I2TAUCIM(E,E0+De/2.,E0+De/4.,D)+I2TAUCIM(E,E0+De,E0+De,D)-
                                               I2TAUCIM(E,E0+De,E0+3.*De/4.,D)));
        }
        else if(ifu==8){// direct Gap Cody
            double E3=K+E0;
            double KCd=PIG/16.*K*K*(2.*E0+K);
            eprj=C/KCd*ReM0M3(E,E0,E3,D);
            epij=C/KCd*aImM0M3(E,E0,E3,D);
        }
        else if(ifu==9){//direct Gap Tauc
            double E3=K+E0;
            double KTd=PIG/2.*(2.*E0+K-2.*sqrt(E0*(E0+K)));
            eprj=C/KTd*ReTaucL(E,E0,D,E3);
            epij=C/KTd*aImTaucL(E,E0,D,E3);
        }
        else if(ifu==10){// Indirect gap Cody with Urbach tail and K-K integral
            long double E3=K+E0;
            eprj=ReIndirCodyUrb(E,C,E0,E3,D);
            epij=ImIndirCodyUrb(E,C,E0,E3,D);
            // if(eV>3.8 && eV<4.)
            //     printf("->ifu=10 E=%Lf C=%Lf E0=%Lf E3=%Lf D=%Lf,eprj=%Lf epij=%Lf\n",E,C,E0,E3,D,eprj,epij);
        }
        else if(ifu==11){// Indirect gap Tauc with Urbach tail and K-K integral
            E3=K+E0;
            W=E3-E0;
            r=E0/W;
            KTi=2.*(8.*r*r*log((4.*r+1.)/(4.*r)) + 8.*pow((r+1.),2.)*log((4.*r+4.)/(4.*r+3.))
                               - (1.+8.*r*(r+1.))*log((4.*r+3.)/(4.*r+1.)) );
            Ex=(E0+sqrt(E0*E0+E0*D*8.))/2.;
            Ey=(E3+sqrt(E3*E3-E3*D*8.))/2.;
            TUx=C*16./KTi/(W*W)/(Ex*Ex)*(Ex-E0)*(Ex-E0); //  if(E<=Ex) then chi2=TUx*exp((E-Ex)/D)
            TUy=C*16./KTi/(W*W)/(Ey*Ey)*(Ey-E3)*(Ey-E3); //  if(E>=Ey) then chi2=TUy*exp(-(E-Ey)/D)
            if(Ex>E0+W/4.){
                Eh=Ex;
                El=E0+W/4.;
                for(int j=1;j<=13;j=j+1){
                    Eb=(Eh+El)/2.; //bisection to find the energy for joining exponential and Tauc on th low energy side
                    derexp=C/KTi/(Eb*Eb)*(-(Eb-E0-W/2.)*(Eb-E0-W/2.)*16./(W*W)+2.)/D;
                    derparab=C/KTi/(Eb*Eb)*(32./W/W*(Eb-E0-W/2.)*((Eb-E0-W/2.)/Eb-1.)-4./Eb);
                    //if(abs(E-3.)<0.01) printf("E=%Lf j=%d Eb=%Lf derexp=%Lf derparab=%Lf \n",E, j, Eb,derexp, derparab);
                    if (derexp>derparab) {
                        Eh=Eb;
                    }
                    else{
                        El=Eb;
                    }
                }
                Ex=Eb;
                TUx=C/KTi/(Ex*Ex)*(-(Ex-E0-W/2.)*(Ex-E0-W/2.)*16./(W*W)+2.);
            }
            if(Ey<E0+W*3./4.){
                Eh=E0+W*3./4.;
                El=Ey;
                for(int j=1;j<=13;j=j+1){
                    Eb=(Eh+El)/2.; //bisection to find the energy for joining exponential and Tauc on th high energy side
                    derexp=-C/KTi/(Eb*Eb)*(-(Eb-E0-W/2.)*(Eb-E0-W/2.)*16./(W*W)+2.)/D;
                    derparab=C/KTi/(Eb*Eb)*(32./W/W*(Eb-E0-W/2.)*((Eb-E0-W/2.)/Eb-1.)-4./Eb);
                    //if(abs(E-3.)<0.01) printf("E=%Lf j=%d Eb=%Lf derexp=%Lf derparab=%Lf \n",E, j, Eb,derexp, derparab);
                    if (derexp<derparab) {
                        El=Eb;
                    }
                    else{
                        Eh=Eb;
                    }
                }
                Ey=Eb;
                TUy=C/KTi/(Ey*Ey)*(-(Ey-E0-W/2.)*(Ey-E0-W/2.)*16./(W*W)+2.);
            }
            eprj=ReIndirTaucUrb(E,C,E0,E3,D,Ex,Ey,TUx,TUy,KTi);
            epij=ImIndirTaucUrb(E,C,E0,E3,D,Ex,Ey,TUx,TUy,KTi);
        }
        else if(ifu==12){// Direct gap Cody with Urbach tail and K-K integral
            long double E3=K+E0;
            eprj=ReDirCodyUrb(E,C,E0,E3,D);
            epij=ImDirCodyUrb(E,C,E0,E3,D);
        }
        else if(ifu==13){// Direct gap Tauc with Urbach tail and K-K integral
            long double E3=K+E0;
            eprj=ReDirTaucUrb(E,C,E0,E3,D);
            epij=ImDirTaucUrb(E,C,E0,E3,D);
        }
        else if(ifu==14){// Direct gap Tauc with Exciton
            complex<long double> uno(1.,0.);
            complex<long double> due(2.,0.);
            complex<long double> pigc(PIG,0.);
            complex<long double> xi,cmpxOut;
            complex<long double> aim(0.,1.);
            complex<double> psi;
            gsl_sf_result re, im;
            int err;
            xi=sqrt(K/(E0-E-aim*D));
            err=gsl_sf_complex_psi_e(double(real(xi)), double(imag(xi)), &re, &im );
            psi=complex<double>(re.val, im.val);
//            if(E>1.4 && E<1.6)
//                printf("xi_re=%f xi_in=%f psi_re=%f psi_im=%f E=%Lf \n" ,double(real(xi)), double(imag(xi)),re.val,im.val,E);
            cmpxOut=due*(log(xi)-pigc/tan(pigc*xi)-(complex<long double>)psi)-uno/xi;
            xi=sqrt(K/(E0+E+aim*D));
            err=gsl_sf_complex_psi_e(double(real(xi)), double(imag(xi)), &re, &im );
            psi=complex<double>(re.val, im.val);
            cmpxOut=cmpxOut+due*(log(xi)-pigc/tan(pigc*xi)-(complex<long double>)psi)-uno/xi;
            xi=sqrt(K/E0);
            err=gsl_sf_complex_psi_e(double(real(xi)), double(imag(xi)), &re, &im );
            psi=complex<double>(re.val, im.val);
            cmpxOut=cmpxOut-due*(due*(log(xi)-pigc/tan(pigc*xi)-(complex<long double>)psi)-uno/xi);
            cmpxOut=cmpxOut*C*sqrt(K)/((E+aim*D)*(E+aim*D));
            eprj=real(cmpxOut);
            epij=imag(cmpxOut);
        }
        else if(ifu==15){// Direct gap Cody M1-M2 CP with K-K integral
            if(D>K*0.45)
                D=K*0.45;
            long double E3=K+E0;
            eprj=ReCodyM1M2(E,C,E0,E3,D);
            epij=ImCodyM1M2(E,C,E0,E3,D);
        }
        else if(ifu==16){// Direct gap Tauc M1-M2 CP with K-K integral
            if(D>K*0.45)
                D=K*0.45;
            long double E3=K+E0;
            eprj=ReTaucM1M2(E,C,E0,E3,D);
            epij=ImTaucM1M2(E,C,E0,E3,D);
        }
        if(ifu==17){// Drude ionized
            long double DOM=D;
            long double Etr=sqrt(2.*C/PIG)*0.75;
            if(E>Etr) DOM=D*pow(Etr/E,1.5);
            eprj=-2.*C/PIG/(DOM*DOM+E*E);
            epij= 2.*C/PIG*DOM/(pow(E,3.)+E*DOM*DOM);
        }
        if(ifu==18){// Lorentz-Dirac oscillator
            long double den=pow(E0*E0-E*E,2.)+pow(2.*E*D+E*E*E/K,2.);
            eprj=2.*C/PIG*(E0*E0-E*E)/den;
            epij=2.*C/PIG*(2.*E*D+E*E*E/K)/den;
// Expression as given in Prokopidis_2014
//            long double den=pow(E0*E0-E*E,2.)+pow(E*(2.*D+E0*E0/K),2.);
//            eprj=2.*C/PIG*((E0*E0-E*E)+E*E/K*(2.*D+E0*E0/K))/den;
//            epij=2.*C/PIG*(2.*E*D+E*E*E/K)/den;
        }
        epr=epr+eprj;
        epi=epi+epij;
        if(isnan(eprj) || isnan(epij))
            printf("->FDISP-ERROR: OscNumber=%d E(eV)=%f eprj=%Lf epij=%Lf\n",io-1,eV,eprj,epij);
    }
    EpsiR=epr;
    EpsiI=epi;
    sqn=sqrt(sqrt(epr*epr+epi*epi)+epr)/sqrt(2.);
    sqk=par[11][1]*sqrt(sqrt(epr*epr+epi*epi)-epr)/sqrt(2.);//par[11][1] is the attenuation coeff of k by Fit#N
    //printf("ioptf=%d eV=%f sqn=%f sqk=%f\n",iopt,eV,sqn,sqk);
}

long double I2CODYRE(long double E,long double E0,long double Ep,long double D){
    long double value=0.5/PIG*((Ep-E)*(3.*E-4.*E0+Ep)+4.*D*(E-E0)*
                      atan2(E-Ep,D)+(pow(E-E0,2.)-D*D)*log(D*D+pow(E-Ep,2.)));
    return(value);
}

long double ICCODYRE(long double E,long double Ep,long double D){
    long double value=0.5/PIG*log(D*D+pow(E-Ep,2.));
    return(value);
}

long double I2CODYIM(long double E,long double E0,long double Ep,long double D){
    long double value=1./PIG*((D*D-pow(E-E0,2.))*atan2(E-Ep,D)+
                       D*(Ep-E+(E-E0)*log(D*D+pow(E-Ep,2.))));
    return(value);
}

long double ICCODYIM(long double E,long double Ep,long double D){
    long double value=-1./PIG*atan2(E-Ep,D);
    return(value);
}

long double ICTAUCRE(long double E,long double Ep,long double D){
    long double value=(-4.*D*E*atan2(E-Ep,D)+2.*E*(D*D+E*E)/Ep+(E*E-D*D)*
                      log((D*D+pow(E-Ep,2.))/(Ep*Ep)))/(2.*PIG*pow(D*D+E*E,2.));
    return(value);
}

long double I2TAUCRE(long double E,long double E0,long double Ep,long double D){
    long double value=(4.*D*E0*(D*D+E*E-E*E0)*atan2(E-Ep,D)+(pow(D,4.)+
                       E*E*pow(E-E0,2.)+D*D*(2.*E*E-2.*E*E0-E0*E0))*log(D*D+pow(E-Ep,2.))+
                        2.*E0*(E*(D*D+E*E)*E0/Ep+(E*E*(2.*E-E0)+D*D*(2.*E+E0))*log(Ep)))/
                        (2.*PIG*pow(D*D+E*E,2.));
    return(value);
}

long double ICTAUCIM(long double E,long double Ep,long double D){
    long double value=-((E*E-D*D)*atan2(E-Ep,D)+D*((D*D+E*E)/Ep+
                       E*log((D*D+pow(E-Ep,2.))/(Ep*Ep))))/(PIG*pow(D*D+E*E,2.));
    return(value);
}

long double I2TAUCIM(long double E,long double E0,long double Ep,long double D){
    long double value=-((pow(D,4.)+E*E*pow(E-E0,2.)+D*D*(2.*E*E-2.*E*E0-E0*E0))*
                       atan2(E-Ep,D)+D*E0*(E0*(D*D+E*E)/Ep+(D*D+E*E-E*E0)*
                       log(Ep*Ep/(D*D+pow(E-Ep,2.)))))/(PIG*pow(D*D+E*E,2.));
    return(value);
}

long double DAWS(long double x){
    static int init=0;
    int NMAX=6;
    long double H=0.4;
    long double A1=2./3.;
    long double A2=0.4;
    long double A3=2./7.;
    int n0;
    long double d1,d2,e1,e2,sum,x2,xp,xx,value;
    if(init==0){//c initialization
        init=1;
        for(int i=1;i<=NMAX;i++)
            cDAW[i]=exp(-pow((2.*double(i)-1.)*H,2.));
    }
    if(abs(x)<0.2){// Use series expansion
        x2=x*x;
        value=x*(1.-A1*x2*(1.-A2*x2*(1.-A3*x2)));
    }
    else{// Use sampling theorem representation
        xx=abs(x);
        n0=2*NINT(0.5*xx/H);
        xp=xx-double(n0)*H;
        e1=exp(2.*xp*H);
        e2=e1*e1;
        d1=double(n0+1);
        d2=d1-2.;
        sum=0.;
        for(int i=1;i<=NMAX;i++){
            sum=sum+cDAW[i]*(e1/d1+1./(d2*e1));
            d1=d1+2.;
            d2=d2-2.;
            e1=e2*e1;
        }
        value=0.5641895835*copysign(exp(-xp*xp),x)*sum;// Constant is 1/sqrt(pig)
    }
    return(value);
}

long double ReFdirGap(long double E,long double E0,long double D,long double x){
    complex<long double> b(E-E0,0.);
    complex<long double> dc(D,0.);
    complex<long double> xc(x,0.);
    complex<long double> aim(0.,1.);
    long double ldv2=2.0;
    long double value=real(ldv2*sqrt(dc*xc+b)-sqrt(aim*dc-b)*atan(sqrt((dc*xc+b)/
                    (aim*dc-b)))-sqrt(b+aim*dc)*atanh(sqrt((dc*xc+b)/(aim*dc+b))));
    return(value);
}

long double aImFdirGap(long double E,long double E0,long double D,long double x){
    complex<long double> b(E-E0,0.);
    complex<long double> dc(D,0.);
    complex<long double> xc(x,0.);
    complex<long double> aim(0.,1.);
    long double value=real(aim*(-sqrt(aim*dc-b)*atan(sqrt((dc*xc+b)/
                       (aim*dc-b)))+sqrt(b+aim*dc)*atanh(sqrt((dc*xc+b)/(aim*dc+b)))));
    return(value);
}


long double ReM0M3(long double E,long double E0,long double E3,long double reD){
    complex<long double> a((E0-E)/reD,0.);
    complex<long double> b((E3-E)/reD,0.);
    complex<long double> d(reD,0.);
    complex<long double> aim(0.,1.);
    long double ldv=0.5;
    complex<long double> ctmp=ldv*d*(a+b-2.*real(sqrt(a+aim)*sqrt(b+aim)));
    return(real(ctmp));
}

long double aImM0M3(long double E,long double E0,long double E3,long double reD){
    complex<long double> a((E0-E)/reD,0.);
    complex<long double> b((E3-E)/reD,0.);
    complex<long double> d(reD,0.);
    complex<long double> aim(0.,1.);
    complex<long double> ctmp=-d*(1.-imag(sqrt(a+aim)*sqrt(b+aim)));
    return(real(ctmp));
}

long double ReTaucL(long double E,long double E0,long double D,long double E3){
    complex<long double> aim(0.,1.);
    complex<long double> cE(E ,0.);
    complex<long double> cE0(E0,0.);
    complex<long double> cD(D ,0.);
    complex<long double> cE3(E3,0.);
    long double ldv=0.5;
    long double ldv1=1.0;
    long double ldv2=2.0;
    long double ldv4=4.0;
    complex<long double> cTmp=ldv/(cD-aim*cE)*(ldv2*aim*cD+ldv2*cE-ldv2*sqrt(cE0*cE3)-(aim-ldv1)*
                               sqrt(ldv2)*sqrt(cD-aim*(cE-cE0))*sqrt(-aim*cD-cE+cE3) );
    cTmp=cTmp+ldv/(cD+aim*cE)*(-ldv2*aim*cD+ldv2*cE-ldv2*sqrt(cE0*cE3)+(aim+ldv1)*
            sqrt(ldv2)*sqrt(cD+aim*(cE-cE0))*sqrt(aim*cD-cE+cE3) );
    cTmp=cTmp+cE*pow(cD-aim*cE,2.)*sqrt(cD+aim*(cE-cE0))*sqrt(cD+aim*(cE-cE3))/
            pow(cD*cD+cE*cE,2.);
    cTmp=cTmp+cE*pow(cD+aim*cE,2.)*sqrt(cD-aim*(cE-cE0))*sqrt(cD-aim*(cE-cE3))/
            pow(cD*cD+cE*cE,2.);
    cTmp=cTmp-cE*cD*(cD*cD*(cE0+cE3)+cE*(-ldv4*cE0*cE3+cE*(cE0+cE3)))/
            (pow(cD*cD+cE*cE,2.)*sqrt(cE0*cE3));
    cTmp=ldv/cD*cTmp;
    return(real(cTmp));
}

long double aImTaucL(long double E,long double E0,long double D,long double E3){
    //COMPLEX aim,cE,cE0,cD,cE3,cTmp
    complex<long double> aim(0.,1.);
    complex<long double> cE(E ,0.);
    complex<long double> cE0(E0,0.);
    complex<long double> cD(D ,0.);
    complex<long double> cE3(E3,0.);
    long double ldv=0.5;
    long double ldv4=4.0;
    complex<long double> cTmp=-pow(cD-aim*cE,2.)*sqrt(cD+aim*(cE-cE0))*sqrt(cD+aim*(cE-cE3))/
                                pow(cD*cD+cE*cE,2.);
    cTmp=cTmp-pow(cD+aim*cE,2.)*sqrt(cD-aim*(cE-cE3))*sqrt(cD-aim*(cE-cE0))/
          pow(cD*cD+cE*cE,2.);
    cTmp=cTmp+cD*(cD*cD*(cE0+cE3)+cE*(-ldv4*cE0*cE3+cE*(cE0+cE3)))/
           pow(cD*cD+cE*cE,2.)/sqrt(cE0*cE3);
    return(real(ldv*cTmp));
}

long double ImDirCodyUrb(long double E,long double C,long double E0,long double E3,long double D){
    //printf("->DirCodyUrb E=%Lf\n",E);
    long double Ex, Ey, chi2;						  
    long double KCd=PIG/16.*(E3-E0)*(E3-E0)*(2.*E0+E3-E0);
    Ex=((D+E3+E0)-sqrt(D*D+(E3-E0)*(E3-E0)))/2.;
    Ey=((E3+E0-D)+sqrt(D*D+(E3-E0)*(E3-E0)))/2.;
    if(E<Ex){
        chi2=C/KCd*sqrt((Ex-E0)*(E3-Ex))*exp((E-Ex)/D);
    }
    else if(E<Ey){
        chi2=C/KCd*sqrt((E-E0)*(E3-E));
    }
    else{
        chi2=C/KCd*sqrt((Ey-E0)*(E3-Ey))*exp((Ey-E)/D);
    }
    return(chi2);
}

long double ReDirCodyUrb(long double E,long double C,long double E0,long double E3,long double D){
// chi1 calculation by Kramers-Kronig integration
    int nstep=50;
    long double Ex, Ey, chi1, chi1old, Emin, Emax,Eminin, Emaxin, KCd, de, x;
    int n1=0,n2=0;
    KCd=PIG/16.*(E3-E0)*(E3-E0)*(2.*E0+E3-E0);
    Ex=((D+E3+E0)-sqrt(D*D+(E3-E0)*(E3-E0)))/2.;
    Ey=((E3+E0-D)+sqrt(D*D+(E3-E0)*(E3-E0)))/2.;
    //Emin, Emax calculation so that alpha(Em)<1e-2 cm^-1
    Eminin=log(1.973e-7/(E0*C/KCd*sqrt((Ex-E0)*(E3-Ex))*exp(-Ex/D)))*D;
    Emaxin=-log(1.973e-7/(E3*C/KCd*sqrt((Ey-E0)*(E3-Ey))*exp(Ey/D)))*D;
    chi1=0.;
    do{
        chi1old=chi1;
        chi1=0.;
        Emax=Emaxin;
        Emin=Eminin;
        if (E>=Emax || E<=Emin){
            de=(Emax-Emin)/nstep;
            n1=2*int((Emax-Emin)/de/2.); //n1 must be even
            n2=0.;
        }
        else{
            de=(Emax-Emin)/nstep;
            if((Emax-E)<de*2.){
                Emax=E+de*2.;
                n1=2;
            }
            else{
                n1=2*int((Emax-E)/de/2.); //n1 must be even
            }
            de=(Emax-E)/n1;
            n2=2*int((E-Emin)/de/2.); //n2 must be even
            Emin=E-de*n2;
        }
        //if(E>4.067 && E<4.069)
        //    printf("Emin=%Lf Emax=%Lf n1=%d n2=%d nstep=%d\n",Emin,Emax,n1,n2,nstep);
        for(int j=1;j<=(n1+n2+1);j=j+2){
            if (n2!=0){
                x=E+de*(j-n2);
            }
            else{
                x=Emin+de*j;
            }
            chi1=chi1+ImDirCodyUrb(x,C,E0,E3,D)*x/(x*x-E*E);
        }
        chi1=chi1*2.*de*2./PIG;
        nstep=nstep*2;
    }while(abs((chi1-chi1old)/chi1)>0.001 && nstep<801);
    return(chi1);
}

long double ImDirTaucUrb(long double E,long double C,long double E0,long double E3,long double D){
    long double Ex, Ey, chi2, q, a1,a2,a0,r;
    complex<long double> s1,s2;
    complex<long double> aim(0.,1.);
    long double KTd=PIG/2.*(E0+E3-2.*sqrt(E0*E3));

    a2=D-E0-E3;
    a1=E0*E3-D*(E3+E0)*1.5;
    a0=D*E0*E3*2.;
    q=a1/3.-a2*a2/9.;
    r=(a1*a2-a0*3.)/6.-a2*a2*a2/27.;

//  long double dis=q*q*q+r*r;
    s1=pow((r+aim*sqrt(-(q*q*q+r*r))),1./3.);
    s2=pow((r-aim*sqrt(-(q*q*q+r*r))),1./3.);

//    Ex=real(s1+s2-a2/3.); unphysical root
    long double r32=sqrt(3.)/2.;
    long double due=2.;
//   Ex=real(aim*(s1-s2)*r32-(s1+s2)/due-a2/3.); unphysical root
    Ex=real(-aim*(s1-s2)*r32-(s1+s2)/due-a2/3.);

    a2=-D-E0-E3;
    a1=E0*E3+D*(E3+E0)*1.5;
    a0=-D*E0*E3*2.;
    q=a1/3.-a2*a2/9.;
    r=(a1*a2-a0*3.)/6.-a2*a2*a2/27.;

//  long double dis=q*q*q+r*r;
    s1=pow((r+aim*sqrt(-(q*q*q+r*r))),1./3.);
    s2=pow((r-aim*sqrt(-(q*q*q+r*r))),1./3.);

    Ey=real(s1+s2-a2/3.);
//   Ey=real(aim*(s1-s2)*r32-(s1+s2)/due-a2/3.); unphysical root
//   Ey=real(-aim*(s1-s2)*r32-(s1+s2)/due-a2/3.); unphysical root

    if(E<Ex){
        chi2=C/KTd*sqrt((Ex-E0)*(E3-Ex))*exp((E-Ex)/D)/(Ex*Ex);
    }
    else if(E<Ey){
        chi2=C/KTd*sqrt((E-E0)*(E3-E))/(E*E);
    }
    else{
        chi2=C/KTd*sqrt((Ey-E0)*(E3-Ey))*exp((Ey-E)/D)/(Ey*Ey);
    }
    return(chi2);
}


long double ReDirTaucUrb(long double E,long double C,long double E0,long double E3,long double D){
    long double Ex, Ey, chi1, chi1old, Emin, Emax, Eminin, Emaxin, de, x, q, a1,a2,a0,r;
    int n1=0,n2=0;
    complex<long double> s1,s2;
    complex<long double> aim(0.,1.);
    long double KTd=PIG/2.*(E0+E3-2.*sqrt(E0*E3));

    a2=D-E0-E3;
    a1=E0*E3-D*(E3+E0)*1.5;
    a0=D*E0*E3*2.;
    q=a1/3.-a2*a2/9.;
    r=(a1*a2-a0*3.)/6.-a2*a2*a2/27.;
    s1=pow((r+aim*sqrt(-(q*q*q+r*r))),1./3.);
    s2=pow((r-aim*sqrt(-(q*q*q+r*r))),1./3.);

//    Ex=real(s1+s2-a2/3.); unphysical root
    long double r32=sqrt(3.)/2.;
    long double due=2.;
//   Ex=real(aim*(s1-s2)*r32-(s1+s2)/due-a2/3.); unphysical root
    Ex=real(-aim*(s1-s2)*r32-(s1+s2)/due-a2/3.);

    a2=-D-E0-E3;
    a1=E0*E3+D*(E3+E0)*1.5;
    a0=-D*E0*E3*2.;
    q=a1/3.-a2*a2/9.;
    r=(a1*a2-a0*3.)/6.-a2*a2*a2/27.;
    s1=pow((r+aim*sqrt(-(q*q*q+r*r))),1./3.);
    s2=pow((r-aim*sqrt(-(q*q*q+r*r))),1./3.);

    Ey=real(s1+s2-a2/3.);
//   Ey=real(aim*(s1-s2)*r32-(s1+s2)/due-a2/3.); unphysical root
//   Ey=real(-aim*(s1-s2)*r32-(s1+s2)/due-a2/3.); unphysical root


    // chi1 calculation by Kramers-Kronig integration

    int nstep=50;
    //Emin, Emax calculation so that alpha(Em)<1e-2 cm^-1
    Eminin=log(1.973e-7/(E0*C/KTd/(Ex*Ex)*sqrt((Ex-E0)*(E3-Ex))*exp(-Ex/D)))*D;
    Emaxin=-log(1.973e-7/(E3*C/KTd/(Ey*Ey)*sqrt((Ey-E0)*(E3-Ey))*exp(Ey/D)))*D;
    //if(E>4.067 && E<4.069)
    //    printf("Ex=%Lf Ey=%Lf Eminin=%Lf Emaxin=%Lf \n",Ex,Ey, Eminin, Emaxin);

    chi1=0.;
    do{
        chi1old=chi1;
        chi1=0.;
        Emax=Emaxin;
        Emin=Eminin;
        if (E>=Emax || E<=Emin){
            de=(Emax-Emin)/nstep;
            n1=2*int((Emax-Emin)/de/2.); //n1 must be even
            n2=0.;
        }
        else{
            de=(Emax-Emin)/nstep;
            if((Emax-E)<de*2.){
                Emax=E+de*2.;
                n1=2;
            }
            else{
                n1=2*int((Emax-E)/de/2.); //n1 must be even
            }
            de=(Emax-E)/n1;
            n2=2*int((E-Emin)/de/2.); //n2 must be even
            Emin=E-de*n2;
        }
        for(int j=1;j<=(n1+n2+1);j=j+2){
            if (n2!=0){
                x=E+de*(j-n2);
            }
            else{
                x=Emin+de*j;
            }
            chi1=chi1+ImDirTaucUrb(x,C,E0,E3,D)*x/(x*x-E*E);
        }
        chi1=chi1*2.*de*2./PIG;
        nstep=nstep*2;
    }while(abs((chi1-chi1old)/chi1)>0.001 && nstep<801);
    //if(E>3.8 && E<4.) printf("E=%Lf chi1=%Lf chi1old=%Lf nstep=%d \n",E,chi1,chi1old, nstep);
    return(chi1);
}

long double ImIndirCodyUrb(long double E,long double C,long double E0,long double E3,long double D){
    long double Ex, Ey, chi2, W, KCi;
    W=E3-E0;
    KCi=W*(E0+W*0.5);
    if(2.*D<=W/4.){
        Ex=E0+D*2.;
        Ey=E3-D*2.;
        if(E<=Ex){
            chi2=C*16./KCi/(W*W)*(Ex-E0)*(Ex-E0)*exp((E-Ex)/D);
        }
        else if(E>Ex && E<=(E0+W/4.)){
            chi2=C*16./KCi/(W*W)*(E-E0)*(E-E0);
        }
        else if(E>(E0+W/4.) && E<=(E0+W*3./4.)){
            chi2=C/KCi*(-(E-E0-W/2.)*(E-E0-W/2.)*16./(W*W)+2.);
        }
        else if(E>(E0+W*3./4.) && E<Ey){
            chi2=C*16./KCi/(W*W)*(E-E3)*(E-E3);
        }
        else{
            chi2=C*16./KCi/(W*W)*(Ey-E3)*(Ey-E3)*exp(-(E-Ey)/D);
        }
    }
    else{
        Ex=E0+W/2.+D-sqrt(D*D+W*W/8.);
        Ey=E0+W/2.-D+sqrt(D*D+W*W/8.);
        if(E<=Ex){
            chi2=C/KCi*(-(Ex-E0-W/2.)*(Ex-E0-W/2.)*16./(W*W)+2.)*exp((E-Ex)/D);
        }
        else if(E<=Ey){
            chi2=C/KCi*(-(E-E0-W/2.)*(E-E0-W/2.)*16./(W*W)+2.);
        }
        else{
            chi2=C/KCi*(-(Ey-E0-W/2.)*(Ey-E0-W/2.)*16./(W*W)+2.)*exp(-(E-Ey)/D);
        }
    }
    return(chi2);
}

long double ReIndirCodyUrb(long double E,long double C,long double E0,long double E3,long double D){
    long double Ex, Ey, chi1, chi1old, Emin, Emax,Eminin, Emaxin,W, KCi, de, x;
    int n1=0,n2=0;
    W=E3-E0;
    KCi=W*(E0+W*0.5);
    if(2.*D<=W/4.){
        Ex=E0+D*2.;
        Ey=E3-D*2.;
    }
    else{
        Ex=E0+W/2.+D-sqrt(D*D+W*W/8.);
        Ey=E0+W/2.-D+sqrt(D*D+W*W/8.);
    }
    //Emin, Emax calculation so that alpha(Em)<1e-2 cm^-1
    Eminin=log(1.973e-7/(E0*C*16./KCi/(W*W)*(Ex-E0)*(Ex-E0)*exp(-Ex/D)))*D;
    Emaxin=-log(1.973e-7/(E3*C*16./KCi/(W*W)*(Ey-E3)*(Ey-E3)*exp(Ey/D)))*D;

    // chi1 calculation by Kramers-Kronig integration

    int nstep=50;
    chi1=0.;
    do{
        chi1old=chi1;
        chi1=0.;
        Emax=Emaxin;
        Emin=Eminin;
        if (E>=Emax || E<=Emin){
            de=(Emax-Emin)/nstep;
            n1=2*int((Emax-Emin)/de/2.); //n1 must be even
            n2=0.;
        }
        else{
            de=(Emax-Emin)/nstep;
            if((Emax-E)<de*2.){
                Emax=E+de*2.;
                n1=2;
            }
            else{
                n1=2*int((Emax-E)/de/2.); //n1 must be even
            }
            de=(Emax-E)/n1;
            n2=2*int((E-Emin)/de/2.); //n2 must be even
            Emin=E-de*n2;
        }
        for(int j=1;j<=(n1+n2+1);j=j+2){
            if (n2!=0){
                x=E+de*(j-n2);
            }
            else{
                x=Emin+de*j;
            }
            chi1=chi1+ImIndirCodyUrb(x,C,E0,E3,D)*x/(x*x-E*E);
        }
        chi1=chi1*2.*de*2./PIG;
        nstep=nstep*2;
    }while(abs((chi1-chi1old)/chi1)>0.001 && nstep<801);
    return(chi1);
}

long double ImIndirTaucUrb(long double E,long double C,long double E0,long double E3,long double D,long double Ex,long double Ey,long double TUx,long double TUy, long double KTi){
    long double chi2, W;
    W=E3-E0;
    if(E<=Ex){
        chi2=TUx*exp((E-Ex)/D);
    }
    else if(E>Ex && E<=(E0+W/4.)){
        chi2=C*16./KTi/(W*W)/(E*E)*(E-E0)*(E-E0);
        //if(abs(E-3.)<0.005) printf("E=%Lf Ex=%Lf E0+W/4.=%Lf \n",E, Ex,E0+W/4.);
    }
    else if(E>max(Ex,(E0+W/4.)) && E<=min(Ey,(E0+W*3./4.))){
        chi2=C/KTi/(E*E)*(-(E-E0-W/2.)*(E-E0-W/2.)*16./(W*W)+2.);
    }
    else if(E>(E0+W*3./4.) && E<Ey){
        chi2=C*16./KTi/(W*W)/(E*E)*(E-E3)*(E-E3);
    }
    else{
        chi2=TUy*exp(-(E-Ey)/D);
    }
    return(chi2);
}

long double ReIndirTaucUrb(long double E,long double C,long double E0,long double E3,long double D,long double Ex,long double Ey,long double TUx,long double TUy, long double KTi){
    long double chi1, chi1old, Eminin, Emaxin, Emin, Emax, de, x;
    int n1=0,n2=0;
    //Emin, Emax calculation so that alpha(Em)<1e-2 cm^-1
    Eminin=log(1.973e-7/(E0*TUx*exp(-Ex/D)))*D;
    Emaxin=-log(1.973e-7/(E3*TUy*exp(Ey/D)))*D;

    // chi1 calculation by Kramers-Kronig integration

    int nstep=50;
    chi1=0.;
    do{
        chi1old=chi1;
        chi1=0.;
        Emax=Emaxin;
        Emin=Eminin;
        if (E>=Emax || E<=Emin){
            de=(Emax-Emin)/nstep;
            n1=2*int((Emax-Emin)/de/2.); //n1 must be even
            n2=0.;
        }
        else{
            de=(Emax-Emin)/nstep;
            if((Emax-E)<de*2.){
                Emax=E+de*2.;
                n1=2;
            }
            else{
                n1=2*int((Emax-E)/de/2.); //n1 must be even
            }
            de=(Emax-E)/n1;
            n2=2*int((E-Emin)/de/2.); //n2 must be even
            Emin=E-de*n2;
        }
        for(int j=1;j<=(n1+n2+1);j=j+2){
            if (n2!=0){
                x=E+de*(j-n2);
            }
            else{
                x=Emin+de*j;
            }
            chi1=chi1+ImIndirTaucUrb(x,C,E0,E3,D,Ex,Ey,TUx,TUy,KTi)*x/(x*x-E*E);
        }
        chi1=chi1*2.*de*2./PIG;
        nstep=nstep*2;
    }while(abs((chi1-chi1old)/chi1)>0.001 && nstep<801);
    return(chi1);
}

long double ImCodyM1M2(long double E,long double C,long double E0,long double E3,long double D){
    //printf("->DirCodyM1M2 E=%Lf C=%Lf E0=%Lf E3=%Lf D=%Lf E3-E0-D*2.=%Lf\n",E,C,E0,E3,D,E3-E0-D*2.);
    long double chi2,K1, EM, KCM1M2;
    K1=(E3-E0)*sqrt(E3-E0-D*2.)*sqrt(2.);
    EM=(E0+E3)/2.;
    KCM1M2=(E3*E3-E0*E0)/6.*(E0*E3-(E3*E3+E0*E0)/2.)
    +(E3*E3-E0*E0)*K1/2.*sqrt((E3-E0)/2.-D)
    -K1*2./15.*sqrt(EM-D-E0)*((EM-D)*(EM-D)*2.-E0*E0*3.+E0*(EM-D))
    -K1*2./15.*sqrt(-EM-D+E3)*(-(EM+D)*(EM+D)*2.+E3*E3*3.-E3*(EM+D));
    chi2=(E-E0)*(E-E3)+K1*sqrt((E3-E0)/2. - D);
    if(E>E3 || E<E0) {
        chi2=0;
    }
    else if (abs(E-EM)>=D) {
        chi2=chi2-K1*sqrt(abs(E-EM)-D);
    }
    //printf("chi2=%Lf KCM1M2=%Lf\n",chi2,KCM1M2);
    chi2=chi2*C/KCM1M2;
    return(chi2);
}


long double ReCodyM1M2(long double E,long double C,long double E0,long double E3,long double D){
// chi1 calculation by Kramers-Kronig integration
    int nstep=50;
    long double chi1, chi1old, Emin, Emax,Eminin, Emaxin, de, x;
    int n1=0,n2=0;
    Eminin=E0;
    Emaxin=E3;
    chi1=0.;
    do{
        chi1old=chi1;
        chi1=0.;
        Emax=Emaxin;
        Emin=Eminin;
        if (E>=Emax || E<=Emin){
            de=(Emax-Emin)/nstep;
            n1=2*int((Emax-Emin)/de/2.); //n1 must be even
            n2=0.;
        }
        else{
            de=(Emax-Emin)/nstep;
            if((Emax-E)<de*2.){
                Emax=E+de*2.;
                n1=2;
            }
            else{
                n1=2*int((Emax-E)/de/2.); //n1 must be even
            }
            de=(Emax-E)/n1;
            n2=2*int((E-Emin)/de/2.); //n2 must be even
            Emin=E-de*n2;
        }
        for(int j=1;j<=(n1+n2+1);j=j+2){
            if (n2!=0){
                x=E+de*(j-n2);
            }
            else{
                x=Emin+de*j;
            }
            chi1=chi1+ImCodyM1M2(x,C,E0,E3,D)*x/(x*x-E*E);
        }
//        if(E>1.998)
//            printf("chi1=%Lf\n",chi1);
        chi1=chi1*2.*de*2./PIG;
        nstep=nstep*2;
    }while(abs((chi1-chi1old)/chi1)>0.001 && nstep<801);
    return(chi1);
}

long double ImTaucM1M2(long double E,long double C,long double E0,long double E3,long double D){
    long double chi2,K1, EM, KTM1M2;
    K1=(E3-E0)*sqrt(E3-E0-D*2.)*sqrt(2.);
    EM=(E0+E3)/2.;
    KTM1M2=-(E3*E3-E0*E0)/2.+log(E3/E0)*(E3*E0+K1*sqrt((E3-E0)/2.-D))
    +K1*2.*(sqrt(EM-D-E0)-sqrt(EM-D)*log((sqrt(EM-D)+sqrt(EM-D-E0))/sqrt(E0)))
    -K1*2.*(sqrt(-EM-D+E3)+sqrt(EM+D)*(asin(sqrt((EM+D)/E3))-PIG/2.));
    chi2=(E-E0)*(E-E3)+K1*sqrt((E3-E0)/2. - D);
    if(E>E3 || E<E0) {
        chi2=0;
    }
    else if (abs(E-EM)>=D) {
        chi2=chi2-K1*sqrt(abs(E-EM)-D);
    }
    chi2=chi2*C/KTM1M2/E/E;
    return(chi2);
}


long double ReTaucM1M2(long double E,long double C,long double E0,long double E3,long double D){
// chi1 calculation by Kramers-Kronig integration
    int nstep=50;
    long double chi1, chi1old, Emin, Emax,Eminin, Emaxin, de, x;
    int n1=0,n2=0;
    Eminin=E0;
    Emaxin=E3;
    chi1=0.;
    do{
        chi1old=chi1;
        chi1=0.;
        Emax=Emaxin;
        Emin=Eminin;
        if (E>=Emax || E<=Emin){
            de=(Emax-Emin)/nstep;
            n1=2*int((Emax-Emin)/de/2.); //n1 must be even
            n2=0.;
        }
        else{
            de=(Emax-Emin)/nstep;
            if((Emax-E)<de*2.){
                Emax=E+de*2.;
                n1=2;
            }
            else{
                n1=2*int((Emax-E)/de/2.); //n1 must be even
            }
            de=(Emax-E)/n1;
            n2=2*int((E-Emin)/de/2.); //n2 must be even
            Emin=E-de*n2;
        }
        for(int j=1;j<=(n1+n2+1);j=j+2){
            if (n2!=0){
                x=E+de*(j-n2);
            }
            else{
                x=Emin+de*j;
            }
            chi1=chi1+ImTaucM1M2(x,C,E0,E3,D)*x/(x*x-E*E);
        }
        chi1=chi1*2.*de*2./PIG;
        nstep=nstep*2;
    }while(abs((chi1-chi1old)/chi1)>0.001 && nstep<801);
    return(chi1);
}

void ASSEMBLER(int iwl, int ikind, double teta, double vot[6][3]){
    int ivnkdw,ivnkup;
    double VNK[17][3],vosi[6][3];
    double alfa=1.;
    complex<double> irup,irdw,irup2,irdw2,mups,mupp,mdws,mdwp,pq,rr;

/* Subroutine for assembling the not coherent contributions arising from the interfaces of the multilayer
   Single or multiple attigue thin films are reduced to a simple interface by the BUILDER subroutine
      iwl: wavelength index

      IKIND kind of measurement
          1 SF
          2 PDS
          3 ELI_1
          .......
          6 ELI_4

      vot(5,2): vector with computed quantities
               1     2
          1    Ts    Tp
          2    Rs    Rp
          3    R1s   R1p
          4    Apds  /
          5    PSI   DELTA
*/    
    double wl=MIS[7][iwl][1];//wavelength

    //vot inizialization
    vot[1][1]=1.;
    vot[1][2]=1.;
    for(int j=2;j<=5;j++){
        vot[j][1]=0.;
        vot[j][2]=0.;
    }

    //complex refractive indices and initialization entrance medium
    COSVNK(VNK,iwl);
    ivnkup=9+ikind;// entrance medium; it depends of the kind of measurement
    irup=complex<double>(VNK[ivnkup][1],-VNK[ivnkup][2]);
    pq=pow(irup*sin(teta),2.);// (N_entrance*sin(teta))**2.
    // if(wl==iwl2print){
    //     printf("pq=%f+i%f lambda=%f iwl=%d \n",real(pq), imag(pq), wl, iwl);
    //     printf("->ASSEMBLER theta=%f ivnkup=%d n=%f k=%f\n",teta/deg2rad,ivnkup,VNK[ivnkup][1],VNK[ivnkup][2]);
    // }
    //*** multilayer parameters
    int Nlayer=NINT(par[51][2]);//Nlayer

    //*** loop on the interfaces = Nlayer+1
    int i=1;
    while(i<=Nlayer+1){
        int ncoe=0;
        // if(iwl==iwl2print)
        //     printf("|ASS> while loop with i=%d of %d interfaces\n",i,Nlayer+1);
        //*** computing R T at i interface
        if(NINT(par[50+i][3])>1 || pm[50+i][1]>.0){//film or roughness
            int ifst=i;
            //thin film
            if(NINT(par[50+i][3])>1){//film!
                ncoe=1;
                while(ifst+ncoe<=Nlayer && NINT(par[50+ifst+ncoe][3])>1){//multilayer
                    ncoe=ncoe+1;
                }
                BUILDER(iwl,ikind,ifst,ncoe,pq,vosi);
                //save the physical parameters affering to the single interface
                vot[4][1]=vosi[4][1];// save Apds
                vot[5][1]=vosi[5][1];// save PSI
                vot[5][2]=vosi[5][2];// salva DELTA
                // if(iwl==iwl2print){
                //     printf("|ASS> found N=%d film \n",ncoe);
                //     printf("\t vosi[1][1]=%f vosi[2][1]=%f vosi[3][1]=%f\n",vosi[1][1],vosi[2][1],vosi[3][1]);
                // }
            }

            //bulk with roughness
            else if(pm[50+i][1]>.0){ //rough bulk!
                // if(iwl==iwl2print)
                //     printf("|ASS> Found rough bulk @i=%d\n",i);
                ncoe=1;//need to use builder
                BUILDER(iwl,ikind,ifst,ncoe,pq,vosi);//pq was theta
                //save the physical parameters affering to the single interface
                vot[4][1]=vosi[4][1];// Apds
                vot[5][1]=vosi[5][1];// PSI
                vot[5][2]=vosi[5][2];// DELTA
                ncoe=0;//reset
            }
        }

        //ideal interface between two media
        else{
            if(i==1)
                ivnkup=9+ikind;// input medium; it depends of the kind of measurement
            else
                ivnkup=NINT(par[50+i-1][1]);
            if(i<=Nlayer)
                ivnkdw=NINT(par[50+i][1]);
            else //the next medium is the exit one for SF
                ivnkdw=9;//16
            irup=complex<double>(VNK[ivnkup][1],-VNK[ivnkup][2]);
            irup2=irup*irup;
            irdw=complex<double>(VNK[ivnkdw][1],-VNK[ivnkdw][2]);
            irdw2=irdw*irdw;
            mups=sqrt(irup2-pq);
            mupp=irup2/mups;
            mdws=sqrt(irdw2-pq);
            mdwp=irdw2/mdws;
            // if(iwl==iwl2print){
            //     printf("|ASS> Ideal interface\n\t ivnkup=%d n=%f k=%f\n",ivnkup,VNK[ivnkup][1],VNK[ivnkup][2]);
            //     printf("\t ivnkdw=%d n=%f k=%f\n",ivnkdw,VNK[ivnkdw][1],VNK[ivnkdw][2]);
            // }
            //symmetric coating on back surface
            if(i==Nlayer+1 && NINT(par[50+i-1][3])==1 && NINT(par[52][2])==1){ //symmetric multistrate on the rear face last layer
                // if(iwl==iwl2print)
                //     printf("|ASS> found symmetric coating on back surface\n");
                for(int ii=1;ii<=2;ii++){
                    vosi[1][ii]=vot[1][ii];
                    vosi[2][ii]=vot[3][ii];
                    vosi[3][ii]=vot[2][ii];
                }
            }
            else{
                // if(iwl==iwl2print)
                //     printf("|ASS> ideal interface between two media\n");
                //vosi[1][1]=4.*real(mups)*real(mdws)/pow(abs(mups+mdws),2.);//requires transparent in out media
                //vosi[1][2]=4.*real(mupp)*real(mdwp)/pow(abs(mupp+mdwp),2.);
                vosi[2][1]=pow(abs((mups-mdws)/(mups+mdws)),2.);
                vosi[2][2]=pow(abs((mupp-mdwp)/(mupp+mdwp)),2.);
                vosi[1][1]=1.-vosi[2][1];//T=1-R because A=0 at the interface
                vosi[1][2]=1.-vosi[2][2];
                vosi[3][1]=vosi[2][1];
                vosi[3][2]=vosi[2][2];
                if(i==1 && Nlayer==1){ //BARE SUBSTRATE
                    //*** Apds
                    vot[4][1]=4.*real(irup*(conj(mdws)-irdw))/pow(abs(irup+mdws),2.);
                    //*** PSI e DELTA
                    rr=((mupp-mdwp)/(mupp+mdwp))/((mups-mdws)/(mups+mdws));
                    vot[5][1]=atan(abs(rr))/deg2rad;//PSI
                    vot[5][2]=arg(rr)/deg2rad;//DELTA
                    // if(iwl==iwl2print)
                    //     printf("|ASS> This is the first interface -> DELTA=%f PSI=%f\n",vot[5][2],vot[5][1]);
                }
            }
        }


        // combination with R T of the previous interface
        if(i==1)
            alfa=1.;
        else{
            ivnkup=NINT(par[50+i-1][1]);
            irup=complex<double>(VNK[ivnkup][1],-VNK[ivnkup][2]);
            irup2=irup*irup;
            mups=sqrt(irup2-pq);
            alfa=exp(-4.*PIG*pm[i-1][1]/wl*imag(-mups));
        }
        // if(iwl==iwl2print){
        //     printf("|ASS> alfa=%f\n",alfa);
        //     printf("|ASS> combining vosi[1][1]=%f vosi[2][1]=%f vosi[3][1]=%f with previous interface\n \tvot[1][1]=%f vot[2][1]=%f vot[3][1]=%f ...\n",
        //            vosi[1][1],vosi[2][1],vosi[3][1],vot[1][1],vot[2][1],vot[3][1]);
        // }
        double uguales=alfa/(1.-vot[3][1]*vosi[2][1]*alfa*alfa);
        double ugualep=alfa/(1.-vot[3][2]*vosi[2][2]*alfa*alfa);
        double ts=vot[1][1]*vosi[1][1]*uguales;
        double tp=vot[1][2]*vosi[1][2]*ugualep;
        double rs=vot[2][1]+pow(vot[1][1],2.)*vosi[2][1]*alfa*uguales;
        double rp=vot[2][2]+pow(vot[1][2],2.)*vosi[2][2]*alfa*ugualep;
        double r1s=vosi[3][1]+pow(vosi[1][1],2.)*vot[3][1]*alfa*uguales;
        double r1p=vosi[3][2]+pow(vosi[1][2],2.)*vot[3][2]*alfa*ugualep;
        vot[1][1]=ts;
        vot[1][2]=tp;
        vot[2][1]=rs;
        vot[2][2]=rp;
        vot[3][1]=r1s;
        vot[3][2]=r1p;
        // if(iwl==iwl2print){
        //     printf("|ASS> ..... ts= %f rs=%f r1s=%f\n",ts,rs,r1s);
        // }
        i=i+ncoe+1;
        fflush(stdout);
    }
}



void BUILDER(int iwl,int ikind,int ifst,int ncoe, complex<double> pq,double vosi[6][3]){
    int ivnk,NFA,iv[10]={},irougfa[11],irougms[11],iv2[10]={};
    int nu[10][3];
    //           [ i][0]=Nlayer with thickness nonuniform
    //           [ i][1]=NFAmin
    //           [ i][1]=NFAmax
    double ds[10],gn[10],cn[10],gk[10],ck[10],ru[10],nmedio,kmedio,al[10],d[999],VNK[17][3],
            dz,BpN,zz,ApN,ApK,BpK,CpK,CpN,Rs,R1s,Rp,R1p,Ts,Tp,As,Ap,ra,Ps,Pp,dsnu[10],DELTA,PSI;
    complex<double> ir[999],NQ,mus,mup,rhoS,rhopS,tauS,rhoP,rhopP,tauP,nq1,mus1,mup1,out[10][3],Bs,Cs,Bp,Cp,denoS,denoP;

    // Subroutine for build-up the indicated coherent part of the multilayer and/or roughness
    double wl=MIS[7][iwl][1];//wavelength

    //if(iwl==iwl2print)
    //    printf("->Start BUILDER with ifst=%d ncoe=%d \n",ifst,ncoe);

    //layer parameters
    int imax=NINT(par[51][2]);// N. layers
    int nino=NINT(par[29][1]);//discretization inhomogeneity
    int N=NINT(par[28][1]);//   discretization roughness
    for(int i=1;i<=9;i++){  //9 is N_layer max
        ds[i]=pm[i][1];     // layer thickness
        al[i]=par[50+i][3]; // kind of layer
        gn[i]=pm[10+i][1]+1240./wl*pm[60+i][1]; //Dn/<n>
        cn[i]=pm[20+i][1];   //curv_n
        gk[i]=pm[30+i][1];   //Dk/<k>
        ck[i]=pm[40+i][1];   //curv_k
        ru[i]=pm[50+i][1];   //roughness
        dsnu[i]=pm[90+i][1]; //nonuniform thickness = (d_max-d_min)/(d_max+d_min)
        if(ru[i]>0.){//reduce the thickness of the layers which interface is rough for 3*sigma
            ds[i]=ds[i]-3.*ru[i];
            if(i>1)
                ds[i-1]=ds[i-1]-3.*ru[i];
        }
    }

    //refractive indices
    COSVNK(VNK,iwl);
    int ivnk1=9+ikind;
    ir[1]=complex<double>(VNK[ivnk1][1],-VNK[ivnk1][2]);
    d[1]=0.;
    if(ifst>1){
        ivnk1=NINT(par[50+ifst-1][1]);
        ir[1]=complex<double>(VNK[ivnk1][1],-VNK[ivnk1][2]);
    }
    NFA=1;
    //if(iwl==iwl2print)printf("\tiNFA=%d n=%f k=%f d=%f\n",1,real(ir[1]),imag(ir[1]),d[1]);

    // set index vector and refractive indices
    int ims=ifst;
    int nroug=0;
    int Nnouni=0;
    int iDWroug=0;
    while(ims<=(ifst+ncoe-1) || (ims==ifst+ncoe && ru[ims]>0. && iDWroug==0)){
        if(ds[ims]<0.00001){
            if(iwl==iwl2print) printf("Skip layer %d because has zero thickness\n",ims);
            ims++;
            continue;
        }
        // roughness
        if(ru[ims]>0.01){// add double layer for roughness
            nroug=nroug+1;//one more rough-interface
            NFA++;
            ir[NFA]=ir[NFA-1];// same refractive index of the previous fase
            d[NFA]=3.*ru[ims];//first half rough layer
            NFA++;
            irougfa[nroug]=NFA;//NFA assigned to the new layer
            irougms[nroug]=ims;//index of the concerning model-layer
            ivnk=int(par[50+ims][1]);//to assign the complex refractive index
            ir[NFA]=complex<double>(VNK[ivnk][1],-VNK[ivnk][2]);
            d[NFA]=3.*ru[ims];//second half rough thickness
            //if(iwl==iwl2print) printf("|BUILD> found ru[%d]=%e => additional double layer:\n NFA=%d d[%d]=%f\n NFA=%d d[%d]=%f\n",
            //                    ims,ru[ims],NFA-1,NFA-1,d[NFA-1],NFA,NFA,d[NFA]);
            if(ims==ifst+ncoe){
                iDWroug++;
                continue;
            }
        }
        if(ds[ims]>0.00001){
            if(NINT(al[ims])==3){ //inomogeneous layer
                //if(iwl==iwl2print) printf("|BUILD> inhomogeneous layer\n");
                ivnk=NINT(par[50+ims][1]);
                nmedio=VNK[ivnk][1];
                kmedio=VNK[ivnk][2];
                dz=1./real(nino);
                ApN=6.*cn[ims]*nmedio;
                ApK=6.*ck[ims]*kmedio;
                BpN=gn[ims]*nmedio-ApN;
                BpK=gk[ims]*kmedio-ApK;
                CpN=nmedio-ApN/3.-BpN/2.;
                CpK=kmedio-ApK/3.-BpK/2.;
                for(int i=1;i<=nino;i++){
                    NFA=NFA+1;
                    zz=1.+dz/2.-real(i)*dz;
                    ir[NFA]=complex<double>(ApN*zz*zz+BpN*zz+CpN,-(ApK*zz*zz+BpK*zz+CpK));
                    d[NFA]=dz*ds[ims];
                }
                if(ru[ims]>0.)//roughness
                    ir[NFA-nino]=ir[NFA-nino+1];//nk for roughness
                if(dsnu[ims]>0.){//nonuniform thickness
                    Nnouni++;
                    nu[Nnouni][0]=ims;//N layer
                    nu[Nnouni][2]=NFA;//NFA_max
                    nu[Nnouni][1]=NFA-nino+1;//NFA_min
                    //if(iwl==iwl2print) printf("|BUILD> found d-nonuniform @ Nlayer=%d dsnu[%d]=%f NFA_min=%d NFA_max=%d\n",ims,ims,dsnu[ims],nu[Nnouni][1],nu[Nnouni][2]);
                }
            }
            else if(NINT(al[ims])==2){// homogeneous layer
                //if(iwl==iwl2print) printf("|BUILD> homogeneous layer\n");
                NFA=NFA+1;
                ivnk=NINT(par[50+ims][1]);
                ir[NFA]=complex<double>(VNK[ivnk][1],-VNK[ivnk][2]);
                d[NFA]=ds[ims];
                //if(iwl==iwl2print)printf("CNK[ivnk][1]=%f\n",CNK[ivnk][1]);
                //if(iwl==iwl2print)printf("\tNFA=%d ivnk=%d n=%f k=%f d=%f par[50+%d][1]=%f\n ",NFA,ivnk,real(ir[NFA]),imag(ir[NFA]),d[NFA],ims,par[50+ims][1]);
                if(dsnu[ims]>0.){//nonuniform thickness
                    Nnouni++;
                    nu[Nnouni][0]=ims;//N layer
                    nu[Nnouni][2]=NFA;//NFA_max
                    nu[Nnouni][1]=NFA;//NFA_min
                    //if(iwl==iwl2print) printf("|BUILD> found d-nonuniform @ Nlayer=%d dsnu[%d]=%f NFA_min=%d NFA_max=%d\n",ims,ims,dsnu[ims],nu[Nnouni][1],nu[Nnouni][2]);
                }
            }
        }
        ims++;
    }
    // substrate
    NFA=NFA+1;
    //if(iwl==iwl2print) printf("|BUILD> ifst=%d ncoe=%d ims=%d imax=%d NFA=%d \n", ifst, ncoe, ims, imax, NFA);
    if(ims<=imax)
        ivnk=NINT(par[50+ims][1]);
    else if(imax==1 && NINT(al[1])==1)
        ivnk=NINT(par[50+imax][1]);
    else{ //without substrate
        ivnk=9;//16
    }
    ir[NFA]=complex<double>(VNK[ivnk][1],-VNK[ivnk][2]);
    d[NFA]=0.;
    // if(iwl==iwl2print){
    //     printf("|BUILDER> exit medium: ims=%d NFA=%d imax=%d\n",ims,NFA,imax);
    //     for(int i=1;i<=NFA;i++)
    //         printf("\tiNFA=%d n=%f k=%f d=%f\n",i,real(ir[i]),imag(ir[i]),d[i]);
    //     printf("|BUILD> nroug=%d\n",nroug);
    // }

    //VOSI & freCoeff initialization for d-nonuniform averaging
    for(int i=0;i<6;i++){
        for(int k=0;k<3;k++){
            vosi[i][k]=0.;
            if(k<2)
                freCoeff[i][k]=complex<double>(0.,0.);
        }
    }
    int M=2*N+1;//effective discretization of rough-layer thickness & Niteration for d-nonuniform
    int Nit2=pow(M,Nnouni);
    double Dit2=double(Nit2);
    for(int it2=0;it2<Nit2;it2++){//iteration over combinations of d-nonuniform
        for(int Jnu=0;Jnu<Nnouni;Jnu++){
            if(Jnu+1<Nnouni)
                iv2[Jnu]=int(it2/pow(M,Nnouni-Jnu-1));
            else{
                if(it2>0)
                    iv2[Jnu]++;
            }
            if(iv2[Jnu]+1>M)
                iv2[Jnu]=iv[Jnu]-M*int(iv2[Jnu]/M);
            int NFAmin=nu[Jnu+1][1];
            int NFAmax=nu[Jnu+1][2];
            double Dd=1./double(NFAmax-NFAmin+1);
            for(int NFAi=NFAmin;NFAi<=NFAmax;NFAi++){
                d[NFAi]=ds[nu[Jnu+1][0]]*Dd*(1.+dsnu[nu[Jnu+1][0]]*(iv2[Jnu]-N)/N);
                //if(iwl==iwl2print) printf("|BUILD> indexLoop d-nonuni it2=%d iv2[%d]=%d NFAmin=%d NFAmax=%d NFAi=%d d[%d]=%f\n",it2,Jnu,iv2[Jnu],NFAmin,NFAmax,NFAi,NFAi,d[NFAi]);
            }
        }
        // managment of roughness and Fresnel's coefficient computing
        rhoS=complex<double>(.0,.0);
        rhopS=complex<double>(.0,.0);
        tauS=complex<double>(.0,.0);
        rhoP=complex<double>(.0,.0);
        rhopP=complex<double>(.0,.0);
        tauP=complex<double>(.0,.0);
        Rs=0.;
        R1s=0.;
        Rp=0.;
        R1p=0.;
        Ts=0.;
        Tp=0.;
        As=0.;
        Ap=0.;
        Ps=0.;
        Pp=0.;
        Bs=complex<double>(.0,.0);
        Cs=complex<double>(.0,.0);
        Bp=complex<double>(.0,.0);
        Cp=complex<double>(.0,.0);
        denoS=complex<double>(.0,.0);
        denoP=complex<double>(.0,.0);
        DELTA=0.;
        PSI=0.;

        for(int it=0;it<pow(M,nroug);it++){//iteration over roughness combinations
            double wtot=1.;
            for(int j=0;j<nroug;j++){//index assignation to each rough layer and each nonuniform thickness
                if(j+1<nroug)
                    iv[j]=int(it/pow(M,nroug-j-1));
                else{
                    if(it>0)
                        iv[j]++;
                }
                if(iv[j]+1>M)
                    iv[j]=iv[j]-M*int(iv[j]/M);
                //printf("iv[%d]=%d ",j,iv[j]);
                ra=pow(-1.,iv[j])*dx[N][int((iv[j]+1)/2.+1)];//normalized thickness variation
                wtot=wtot*w[N][int((iv[j]+1)/2.+1)];         //cumulative weight
                d[irougfa[j+1]-1]=(3.+ra)*ru[irougms[j+1]];   //thickness first half rough layer
                d[irougfa[j+1]]=(3.-ra)*ru[irougms[j+1]];//thickness seconf half rough layer
                // if(iwl==iwl2print) {
                //     printf("|BUILD> indexLoop rough it=%d j=%d iv[%d]=%d d[%d]=%f d[%d]=%f wtot=%f\n",
                //            it,j,j,iv[j],irougfa[j+1]-1,d[irougfa[j+1]-1],irougfa[j+1],d[irougfa[j+1]],wtot);
                //     //iw=1;
                // }
            }
            // if(iwl==iwl2print){
            //     printf("callCALFRE with wl=%f pq=%f +%f*i\n",wl,real(pq),imag(pq));
            //     for(int k=1;k<=NFA;k++)
            //         printf("callCALFRE iNFA=%d n=%f k=%f d=%f\n",k,real(ir[k]),imag(ir[k]),d[k]);
            // }
            CALFRE(NFA,wl,pq,ir,d,out);
            // if(iwl==iwl2print){
            //     for(int k=0;k<=7;k++)
            //         printf("   out[%d][1]=%f + %f*i\n",k,real(out[k][1]),imag(out[k][1]));
            // }
            iw=0;
            tauS=tauS+out[1][1]*wtot;
            rhoS=rhoS+out[2][1]*wtot;
            rhopS=rhopS+out[3][1]*wtot;
            Ts=Ts+real(out[4][1])*wtot;
            Rs=Rs+real(out[5][1])*wtot;
            R1s=R1s+real(out[6][1])*wtot;
            As=As+real(out[7][1])*wtot;
            tauP=tauP+out[1][2]*wtot;
            rhoP=rhoP+out[2][2]*wtot;
            rhopP=rhopP+out[3][2]*wtot;
            Tp=Tp+real(out[4][2])*wtot;
            Rp=Rp+real(out[5][2])*wtot;
            R1p=R1p+real(out[6][2])*wtot;
            Ap=Ap+real(out[7][2])*wtot;
            Bs=Bs+out[0][0]*wtot;
            Cs=Cs+out[1][0]*wtot;
            Bp=Bp+out[2][0]*wtot;
            Cp=Cp+out[3][0]*wtot;
            Ps=Ps+real(out[0][1])*wtot;
            Pp=Pp+real(out[0][2])*wtot;
            PSI=PSI+atan(abs(rhoP/rhoS))/deg2rad*wtot;
            DELTA=DELTA+arg(rhoP/rhoS)/deg2rad*wtot;
        }
        //save Fresnell coefficients
        freCoeff[0][0]=freCoeff[0][0]+Ps/Dit2;//Poynting s-pol
        freCoeff[0][1]=freCoeff[0][1]+Pp/Dit2;//Poynting p-pol
        freCoeff[1][0]=freCoeff[1][0]+rhoS/Dit2;
        freCoeff[1][1]=freCoeff[1][1]+rhoP/Dit2;
        freCoeff[2][0]=freCoeff[2][0]+rhopS/Dit2;
        freCoeff[2][1]=freCoeff[2][1]+rhopP/Dit2;
        freCoeff[3][0]=freCoeff[3][0]+ir[NFA]/Dit2;
        freCoeff[4][0]=freCoeff[4][0]+Bs/Dit2;
        freCoeff[4][1]=freCoeff[4][1]+Cs/Dit2;
        freCoeff[5][0]=freCoeff[5][0]+Bp/Dit2;
        freCoeff[5][1]=freCoeff[5][1]+Cp/Dit2;

        //save to VOSI
        if(NINT(par[54][2])==0){ //specular SF measurements
            NQ=ir[NFA]*ir[NFA];
            mus=sqrt(NQ-pq);
            mup=NQ/mus;
            nq1=ir[1]*ir[1];
            mus1=sqrt(nq1-pq);
            mup1=nq1/mus1;

            vosi[1][1]=vosi[1][1]+real(mus)/real(mus1)*pow(abs(tauS),2.)/Dit2;
            vosi[1][2]=vosi[1][2]+real(mup)/real(mup1)*pow(abs(tauP),2.)/Dit2;
            vosi[2][1]=vosi[2][1]+pow(abs(rhoS),2.)/Dit2;
            vosi[2][2]=vosi[2][2]+pow(abs(rhoP),2.)/Dit2;
            vosi[3][1]=vosi[3][1]+pow(abs(rhopS),2.)/Dit2;
            vosi[3][2]=vosi[3][2]+pow(abs(rhopP),2.)/Dit2;
        }
        else{ // hemispherical SF measurements
            vosi[1][1]=vosi[1][1]+Ts/Dit2;
            vosi[1][2]=vosi[1][2]+Tp/Dit2;
            vosi[2][1]=vosi[2][1]+Rs/Dit2;
            vosi[2][2]=vosi[2][2]+Rp/Dit2;
            vosi[3][1]=vosi[3][1]+R1s/Dit2;
            vosi[3][2]=vosi[3][2]+R1p/Dit2;
        }
        vosi[4][1]=vosi[4][1]+As/Dit2;
        vosi[4][2]=vosi[4][2]+Ap/Dit2;
        // vosi[5][1]=vosi[5][1]+atan(abs(rhoP/rhoS))/deg2rad/Dit2;//PSI
        // vosi[5][2]=vosi[5][2]+arg(rhoP/rhoS)/deg2rad/Dit2;//DELTA
        vosi[5][1]=vosi[5][1]+PSI/Dit2;//PSI
        vosi[5][2]=vosi[5][2]+DELTA/Dit2;//DELTA
        // if(iwl==1wl2print){
        //     printf("|BUILD>final output:\n DELTA=%f  PSI=%f\n",vosi[5][2],vosi[5][1]);
        //     printf(" tauP= %f +i%f Tp=%f\n tauS= %f +i%f Ts=%f\n",
        //        real(tauP),imag(tauP),real(mup)/real(mup1)*pow(abs(tauP),2.),real(tauS),imag(tauS),real(mus)/real(mus1)*pow(abs(tauS),2.));
        //     printf(" rhoP= %f +i%f Rp=%f\n rhoS= %f +i%f Rs=%f\n",
        //        real(rhoP),imag(rhoP),pow(abs(rhoP),2.),real(rhoP),imag(rhoP),pow(abs(rhoP),2.));
        //     printf(" rho1P= %f +i%f R1p=%f\n rho1S= %f +i%f R1s=%f\n",
        //        real(rhopP),imag(rhopP),pow(abs(rhopP),2.),real(rhopS),imag(rhopS),pow(abs(rhopS),2.));
        //     fflush(stdout);
        // }
    }
}


void CALFRE(int NFAin,double wl,complex<double> pq,complex<double> IR[999],double d[999],complex<double> out[10][3]){
    complex<double> SC[3][3],PC[3][3],MUS,MUP,NQ,CI,DELTA,C,S,B,BP,CP,deno,nq1,mus1,mup1,
            rhoS,rhopS,tauS,rhoP,rhopP,tauP;
    complex<double> a11,a12,a21,a22,b11,b12,b21,b22;
    CI=complex<double>(0.,1.);

    //if(iw==1) printf("-> |CALFRE> NFAin=%d wl=%f\n",NFAin,wl);

//check if a layer absorbs completely the light. In that case this layer becomes the exit medium
    int NFA=NFAin;
    for(int I=2;I<=NFAin-1;I++){
        NQ=IR[I]*IR[I];
        MUS=sqrt(NQ-pq);
        MUP=NQ/MUS;
        DELTA=2.*PIG*d[I]*MUS/wl;
        if(abs(imag(DELTA))>ABSmax) {
            NFA=I;
            break;
        }
//        if(wl>11820. && wl<11850.){
//            printf("pq=%f+i%f theta[%d]=%f k[%d]=%f k_enhanced[%d]=%f \n",real(pq), imag(pq), I,atan(sqrt(real(pq))/real(MUS))/deg2rad, I, imag(IR[I]),I, imag(MUS));
//            printf("n*sin(theta)[%d]=%f+i%f \n",I,real(IR[I]*sin(acos(MUS/IR[I]))), imag(IR[I]*sin(acos(MUS/IR[I]))));
//        }
    }

    //Computing the characteristic matrix and Fresnel's coefficients
    //initialization of SC and PC
    SC[1][1]=complex<double>(1.,0.);
    SC[2][2]=complex<double>(1.,0.);
    PC[1][1]=complex<double>(1.,0.);
    PC[2][2]=complex<double>(1.,0.);
    SC[1][2]=complex<double>(0.,0.);
    SC[2][1]=complex<double>(0.,0.);
    PC[1][2]=complex<double>(0.,0.);
    PC[2][1]=complex<double>(0.,0.);
    // characteristic matrix
    for(int I=2;I<=NFA-1;I++){
        NQ=IR[I]*IR[I];
        MUS=sqrt(NQ-pq);//IR[I]*cos(theta_I)
        MUP=NQ/MUS;
        DELTA=2.*PIG*d[I]*MUS/wl;
        C=cos(DELTA);
        S=sin(DELTA)*CI;

        //s-pol
        a11=SC[1][1];
        a12=SC[1][2];
        a21=SC[2][1];
        a22=SC[2][2];
        b11=C;
        b12=S/MUS;
        b21=S*MUS;
        b22=C;
        SC[1][1]=a11*b11+a12*b21;
        SC[1][2]=a11*b12+a12*b22;
        SC[2][1]=a21*b11+a22*b21;
        SC[2][2]=a21*b12+a22*b22;

        //p-pol
        a11=PC[1][1];
        a12=PC[1][2];
        a21=PC[2][1];
        a22=PC[2][2];
        //b11=C;
        b12=S/MUP;
        b21=S*MUP;
        //b22=C;
        PC[1][1]=a11*b11+a12*b21;
        PC[1][2]=a11*b12+a12*b22;
        PC[2][1]=a21*b11+a22*b21;
        PC[2][2]=a21*b12+a22*b22;
        // if(abs(SC[1][1])>1.E+06 || abs(SC[1][2])>1.E+06){
        //     printf("CALFRE @wl=%f & I=%d SC[1][1]=%f +%f*i SC[1][2]=%f +%f*i\n",wl,I,real(SC[1][1]),imag(SC[1][1]),real(SC[1][2]),imag(SC[1][2]));
        //     printf("|NQ|=%f |MUS|=%f |DELTA|=%f\n",abs(NQ),abs(MUS),abs(imag(DELTA)));
        // }
    }
    NQ=IR[NFA]*IR[NFA];
    MUS=sqrt(NQ-pq);
    MUP=NQ/MUS;
    nq1=IR[1]*IR[1];
    mus1=sqrt(nq1-pq);
    mup1=nq1/mus1;

    //Fresnel's coefficient and R R1 T A
//  polarization-s
    B=SC[1][1]+SC[1][2]*MUS;
    C=SC[2][1]+SC[2][2]*MUS;
    //BP=-SC[1][1]+SC[1][2]*MUS;
    //CP=-SC[2][1]+SC[2][2]*MUS;
    BP=SC[2][2]+SC[1][2]*mus1;//Macleod's formula for inverse path
    CP=SC[2][1]+SC[1][1]*mus1;
    deno=mus1*B+C;
    rhoS=(mus1*B-C)/deno;
    tauS=2.*mus1/deno;
    rhopS=(MUS*BP-CP)/(MUS*BP+CP);
    double Rs=pow(abs(rhoS),2.);
    double R1s=pow(abs(rhopS),2.);
    double Ts=real(MUS)/real(mus1)*pow(abs(tauS),2.);
    double As=4.*real(mus1)*real(B*conj(C)-MUS)/pow(abs(deno),2.);
    out[0][0]=B;
    out[1][0]=C;
    out[0][1]=real(B*conj(C));//Poynting vector s-pol

//  polarization-p
    B=PC[1][1]+PC[1][2]*MUP;
    C=PC[2][1]+PC[2][2]*MUP;
    //BP=-PC[1][1]+PC[1][2]*MUP;
    //CP=-PC[2][1]+PC[2][2]*MUP;
    BP=PC[2][2]+PC[1][2]*mup1;//Macleod's formula for inverse path
    CP=PC[2][1]+PC[1][1]*mup1;
    deno=mup1*B+C;
    rhoP=(mup1*B-C)/deno;
    tauP=2.*mup1/deno;
    rhopP=(MUP*BP-CP)/(MUP*BP+CP);
    double Rp=pow(abs(rhoP),2.);
    double R1p=pow(abs(rhopP),2.);
    double Tp=real(MUP)/real(mup1)*pow(abs(tauP),2.);
    double Ap=4.*real(mup1)*real(B*conj(C)-MUP)/pow(abs(deno),2.);
    out[2][0]=B;
    out[3][0]=C;
    out[0][2]=real(B*conj(C));//Poynting vector p-pol

    // PSI & DELTA
    double PSI=atan(abs(rhoP/rhoS))/deg2rad;
    double DEL=atan2(-imag(rhoP/rhoS),-real(rhoP/rhoS))/deg2rad;//periodicity 2*pig

    // save to OUT
    out[1][1]=tauS;
    out[2][1]=rhoS;
    out[3][1]=rhopS;
    out[4][1]=Ts;
    out[5][1]=Rs;
    out[6][1]=R1s;
    out[7][1]=As;
    out[8][1]=PSI;
    out[1][2]=tauP;
    out[2][2]=rhoP;
    out[3][2]=rhopP;
    out[4][2]=Tp;
    out[5][2]=Rp;
    out[6][2]=R1p;
    out[7][2]=Ap;
    out[8][2]=DEL;

    if (NFA<NFAin) {
        out[1][1]=0.;
        out[4][1]=0.;
        out[7][1]=1.-Rs;
        out[1][2]=0.;
        out[4][2]=0.;
        out[7][2]=1.-Rp;
    }
    // if(iw==1){
    //     printf(" tauP= %f +i%f Tp=%f\n tauS= %f +i%f Ts=%f\n",
    //        real(tauP),imag(tauP),real(MUP)/real(mup1)*pow(abs(tauP),2.),real(tauS),imag(tauS),real(MUS)/real(mus1)*pow(abs(tauS),2.));
    //     printf(" rhoP= %f +i%f Rp=%f\n rhoS= %f +i%f Rs=%f\n",
    //        real(rhoP),imag(rhoP),pow(abs(rhoP),2.),real(rhoP),imag(rhoP),pow(abs(rhoP),2.));
    //     printf(" rho1P= %f +i%f R1p=%f\n rho1S= %f +i%f R1s=%f\n",
    //        real(rhopP),imag(rhopP),pow(abs(rhopP),2.),real(rhopS),imag(rhopS),pow(abs(rhopS),2.));
    //     fflush(stdout);
    // }
}



int SOLVE(int imis,int iWL,double *Xp,double *Yp,double *ErrYp){
    int Ndat=0;
    int ikind=0,ivot=0,IPD=0;
    double tetar=0.,KCL,vot[6][3];
    //   iWL= wl index

    // set type and input medium
    if(imis==1){
        ikind=1;
        tetar=.0;
        ivot=1;
    }
    else if(imis==2){
        ikind=1;
        tetar=par[6][1];
        ivot=1;
    }
    else if(imis==3){
        ikind=1;
        tetar=.0;
        ivot=2;
    }
    else if(imis==4){
        ikind=1;
        tetar=par[6][1];
        ivot=2;
    }
    else if(imis==5){
        ikind=1;
        tetar=.0;
        ivot=3;
    }
    else if(imis==6){
        ikind=2;
        tetar=.0;
        ivot=4;
    }
    else if(imis==7 || imis==8){
        ikind=3;
        tetar=pm[86][1];
        ivot=5;
        IPD=imis-int(imis/2.)*2;
    }
    else if(imis==9 || imis==10){
        ikind=4;
        tetar=pm[87][1];
        ivot=5;
        IPD=imis-int(imis/2.)*2;
    }
    else if(imis==11 ||imis==12){
        ikind=5;
        tetar=pm[88][1];
        ivot=5;
        IPD=imis-int(imis/2.)*2;
    }
    else if(imis==13||imis==14){
        ikind=6;
        tetar=pm[89][1];
        ivot=5;
        IPD=imis-int(imis/2.)*2;
    }
    int s1p2=NINT(par[27][2]);
    double LAM=MIS[7][iWL][1];
    printf("-> SOLVE: imis=%d ikind=%d iWL=%d s1p2=%d teta=%f\n",imis,ikind,iWL,s1p2,tetar);
    tetar=tetar*deg2rad;
    double Dn=(par[1][2]-par[1][1])/100.;//n-step
    double Dk=(par[2][2]-par[2][1])/500.;//k-step
    double k1,FM,FMold=1.E+10;
    CNK[1][1]=0;//nk are unknown
    for(int iN=0;iN<101;iN++){
        CNK[1][2]=par[1][1]+iN*Dn;
        FMold=1.E+10;
        for(int iK=0;iK<501;iK++){
            KCL=par[2][1]+iK*Dk;
            CNK[1][3]=KCL;
            ASSEMBLER(iWL,ikind,tetar,vot);
            if(imis<=6)
                FM=abs(MIS[imis][iWL][1]-vot[ivot][s1p2])/MIS[imis][iWL][2];
            else
                FM=abs(sin(deg2rad*(ELI[imis-6][iWL][1]-vot[ivot][IPD+1])/2.))/(deg2rad*ELI[imis-6][iWL][2]);
            if(FM<=1. && FMold>1.){
                k1=KCL;
                //printf("k1=%f ",KCL);
            }
            else if(FM>1. && FMold<=1.){
                Xp[Ndat]=CNK[1][2];
                Yp[Ndat]=0.5*(k1+KCL);
                ErrYp[Ndat]=0.5*(KCL-k1);
                if(ErrYp[Ndat]<Dk/4.)
                    ErrYp[Ndat]=Dk/4.;
                //printf("Ndat=%d Xp=%f Yp=%f ErrYp=%f ",Ndat,Xp[Ndat],Yp[Ndat],ErrYp[Ndat]);
                Ndat++;
            }
            //printf("n=%f k=%f FM=%f FMold=%f\n",CNK[1][2],CNK[1][3],FM,FMold);
            FMold=FM;
        }
    }
    return(Ndat);
}


double FMER(double k){
    int iWL=NINT(par[24][2]);
    int s1p2=NINT(par[27][2]);
    double teta[6];
    CNK[1][3]=k;
    teta[1]=par[6][1];
    teta[2]=pm[86][1];
    teta[3]=pm[87][1];
    teta[4]=pm[88][1];
    teta[5]=pm[89][1];
    for(int i=1;i<=5;i++)
        teta[i]=teta[i]*deg2rad;
    double vot[6][3];
    double wl=MIS[7][iWL][1];
    double FM=.0;
    int NMIS=0;
    if(DATO[1]==2 || DATO[3]==2 ||DATO[5]==2){
        ASSEMBLER(iWL,1,0.,vot);
        for(int i=1;i<=3;i++){
            int j=2*i-1;
            if(DATO[j]==2){
                FM=FM+pow((vot[i][s1p2]-MIS[j][iWL][1])/MIS[j][iWL][2],2.);
                NMIS++;
            }
        }
    }
    if(DATO[2]==2 || DATO[4]==2){
        ASSEMBLER(iWL,1,teta[1],vot);
        for(int i=1;i<=2;i++){
            int j=2*i;
            if(DATO[j]==2){
                FM=FM+pow((vot[i][s1p2]-MIS[j][iWL][1])/MIS[j][iWL][2],2.);
                NMIS++;
            }
        }
    }
    if(DATO[6]==2){
        ASSEMBLER(iWL,2,0.,vot);
        FM=FM+pow(vot[4][1]-MIS[6][iWL][1],2.);
        NMIS++;
    }
    for(int I=1;I<=4;I++){
        if(DATO[5+2*I]==2 || DATO[6+2*I]==2){
            ASSEMBLER(iWL,I+2,teta[1+I],vot);
            double DE=ELI[5+2*I-6][iWL][1];
            double DC=vot[5][2];//Delta;
            double DDE=ELI[5+2*I-6][iWL][2]*deg2rad;
            double PE=ELI[6+2*I-6][iWL][1];
            double PC=vot[5][1];//PSI
            double DP=ELI[6+2*I-6][iWL][2];
            double DFMD=pow((cos(DC*deg2rad)-cos(DE*deg2rad))/DDE,2.);
            double DFMP=pow((PC-PE)/DP,2.);
            if(DATO[5+2*I]==2){
                FM=FM+DFMD;
                NMIS++;
            }
            if(DATO[6+2*I]==2){
                FM=FM+DFMP;
                NMIS++;
            }
        }
    }
    if(NMIS>0)
        FM=FM/NMIS;
    return(FM);
}


int FSQ(void *p, int m, int n, const double *x, double *fvec, int iflag){
    /* calculate the functions at x and return the values in fvec[0] through fvec[m-1] */
    //struct pointToFit2 *pTF2 = (struct pointToFit2 *)p;
    int ioptf=NINT(pm[100][1]);// FIT#
    // addressing fit-parameters to PM
    for(int i=0;i<n;i++){
        int ip=NINT(pm[i+1][3]);
        pm[ip][1]=x[i];
    }

    // computing fvec for n fit
    int md=NINT(SOL[1][1]);
    int j=-1;
    double chi2=0.;
    double fredeg=real(m-n);
    double erynMIN=(rxy[16][2]-rxy[16][1])/1.e4;
    double erykMIN=(rxy[17][2]-rxy[17][1])/1.e4;
    double ere1MIN=(rxy[26][2]-rxy[26][1])/1.e4;
    double ere2MIN=(rxy[27][2]-rxy[27][1])/1.e4;
    double wl,ev,y,ery;
    for(int i=1;i<=md;i++){
        if(NINT(SOL[i][6])==1){//the data is enabled
            wl=SOL[i][1];
            ev=1240./wl;
            FDISP(ioptf,ev);
            if(NINT(par[32][5])<=1){
                j++;
                y=SOL[i][2];
                ery=SOL[i][4];
                if(ery<erynMIN) ery=erynMIN;
                fvec[j]=(y-sqn)/ery;
                chi2=chi2+pow(fvec[j],2.)/fredeg;
            }
            if(NINT(par[32][5])==1){
                j++;
                y=SOL[i][3];
                ery=SOL[i][5];
                if(ery<erykMIN) ery=erykMIN;
                fvec[j]=(y-sqk)/ery;
                chi2=chi2+pow(fvec[j],2.)/fredeg;
            }
            if(NINT(par[32][5])==2){
                j++;
                ery=2.*(SOL[i][2]*SOL[i][4]+SOL[i][3]*SOL[i][5]);
                if(ery<ere1MIN) ery=ere1MIN;
                fvec[j]=((SOL[i][2]*SOL[i][2]-SOL[i][3]*SOL[i][3])-EpsiR)/ery;
                chi2=chi2+pow(fvec[j],2.)/fredeg;
                j++;
                ery=2.*(SOL[i][4]*SOL[i][3]+SOL[i][2]*SOL[i][5]);
                if(ery<ere2MIN) ery=ere2MIN;
                fvec[j]=(2.*SOL[i][2]*SOL[i][3]-EpsiI)/ery;;
                chi2=chi2+pow(fvec[j],2.)/fredeg;
            }
        }
    }
    // monitor values
    printf("|FSQ> chi2=%.9g\t",chi2);
//    for(int i=0;i<n;i++)
//        printf("p[%d]=%.9g\t",i,x[i]);
    printf("\n");
    if(iRecChi2==1){
        QFile fchi2(filechi2);
        if(fchi2.open(QIODevice::WriteOnly | QIODevice::Text | QIODevice::Append)){
            QTextStream out(&fchi2);
            out<<QString::number(chi2,'g',9);
            for(int i=0;i<n;i++)
                out<<"\t"<<QString::number(x[i],'g',9);
            out<<"\n";
            fchi2.close();
        }
    }
    QCoreApplication::processEvents();
    if(iStop){
        iflag=-1;
        return(iflag);
    }
    return(0);
}


int FSEM(void *p, int m, int n, const double *x, double *fvec, int iflag){
    //***** Function to fit selected experimental measurements
    //****
    int ioptf=NINT(pm[100][1]);// set FIT#
    int mwl=NeV;
    int s1p2=NINT(par[27][2]);
    double vot[6][3];
    double te=0.;
    // addressing fit parameter onto PM
    int ip;
    for(int i=1;i<=n;i++){
        ip=NINT(pm[i][3]);
        pm[ip][1]=x[i-1];
        if(iw>0) printf("|FSEM> p[%d]=%f\n",i,x[i-1]);
    }
    double chi2=.0;
    double fredeg=double(m-n);
    //computing fvec
    int ivec=0;
    double wl,ev;
    for(int i=1;i<=mwl;i++){
        wl=MIS[7][i][1];
        ev=1240./wl;
        par[24][2]=i;
        if(ioptf<1 || ioptf>7){
            FDISP(ioptf,ev);
            MIS[16][i][1]=sqn;
            CNK[1][2]=sqn;
            MIS[16][i][2]=sqk;
            CNK[1][3]=sqk;
        }
        CNK[1][1]=ioptf;
        //computing Tn Rn R1
        if(DATO[1]==2 || DATO[3]==2 || DATO[5]==2){
            ASSEMBLER(i,1,0.,vot);
            if(DATO[1]==2){
                fvec[ivec]=(MIS[1][i][1]-vot[1][s1p2])/MIS[1][i][2];
                chi2=chi2+pow(fvec[ivec],2.)/fredeg;
                ivec++;
            }
            if(DATO[3]==2){
                fvec[ivec]=(MIS[3][i][1]-vot[2][s1p2])/MIS[3][i][2];
                chi2=chi2+pow(fvec[ivec],2.)/fredeg;
                ivec++;
            }
            if(DATO[5]==2){
                fvec[ivec]=(MIS[5][i][1]-vot[3][s1p2])/MIS[5][i][2];
                chi2=chi2+pow(fvec[ivec],2.)/fredeg;
                ivec++;
            }
        }
        //computing Tp Rp
        if(DATO[2]==2 || DATO[4]==2){
            te=par[6][1]*deg2rad;
            ASSEMBLER(i,1,te,vot);
            if(DATO[2]==2){
                fvec[ivec]=(MIS[2][i][1]-vot[1][s1p2])/MIS[2][i][2];
                chi2=chi2+pow(fvec[ivec],2.)/fredeg;
                ivec++;
            }
            if(DATO[4]==2){
                fvec[ivec]=(MIS[4][i][1]-vot[2][s1p2])/MIS[4][i][2];
                chi2=chi2+pow(fvec[ivec],2.)/fredeg;
                ivec++;
            }
        }
        //computing PSI DELTA
        for(int j=1;j<=4;j++){
            if(DATO[7+2*(j-1)]==2 || DATO[8+2*(j-1)]==2){
                te=pm[85+j][1]*deg2rad;
                ASSEMBLER(i,2+j,te,vot);
                if(DATO[7+2*(j-1)]==2){//DELTA
                    //double Der=abs(cos((ELI[1+2*(j-1)][i][1]+ELI[1+2*(j-1)][i][2])*deg2rad)-cos((ELI[1+2*(j-1)][i][1]-ELI[1+2*(j-1)][i][2])*deg2rad));
                    //fvec[ivec]=(cos(ELI[1+2*(j-1)][i][1]*deg2rad)-cos(vot[5][2]*deg2rad))/Der;
                    double offset=NINT((ELI[1+2*(j-1)][i][1]-vot[5][2])/Dperiod)*Dperiod;
                    fvec[ivec]=(ELI[1+2*(j-1)][i][1]-vot[5][2]-offset)/ELI[1+2*(j-1)][i][2];
                    chi2=chi2+pow(fvec[ivec],2.)/fredeg;
                    ivec++;
                }
                if(DATO[8+2*(j-1)]==2){//PSI
                    fvec[ivec]=(ELI[2+2*(j-1)][i][1]-vot[5][1])/ELI[2+2*(j-1)][i][2];
                    chi2=chi2+pow(fvec[ivec],2.)/fredeg;
                    ivec++;
                }
            }
        }
    }
    // monitor values
    printf("|FSEM> ioptf=%d chi2=%.9g\t",ioptf,chi2);
//    for(int i=0;i<n;i++)
//        printf("p[%d]=%.9g\t",i,x[i]);
    printf("\n");
    if(iRecChi2==1){
        QFile fchi2(filechi2);
        if (fchi2.open(QIODevice::WriteOnly | QIODevice::Text | QIODevice::Append)){
            QTextStream out(&fchi2);
            out<<QString::number(chi2,'g',9);
            for(int i=0;i<n;i++)
                out<<"\t"<<QString::number(x[i],'g',9);
            out<<"\n";
            fchi2.close();
        }
    }
    QCoreApplication::processEvents();
    if(iStop){
        iflag=-1;
        return(iflag);
    }
    return(0);
}


int FRCK(void *p, int m, int n, const double *x, double *fvec, int iflag){
    //***** Subroutine of IbridOne: computing k from T given n
    //****
    int ioptf=NINT(pm[100][1]);// set FIT#
    int mwl=NeV;
    int s1p2=NINT(par[27][2]);
    double vot[6][3];
    double te=0.;
    double FM;
    if(DATO[4]==2)
        te=par[6][1]*deg2rad;
    else if(DATO[8]==2)
        te=pm[86][1]*deg2rad;
    else if(DATO[10]==2)
        te=pm[87][1]*deg2rad;
    else if(DATO[12]==2)
        te=pm[88][1]*deg2rad;
    else if(DATO[14]==2)
        te=pm[89][1]*deg2rad;
    // addressing fit parameter onto PM
    int ip;
    for(int i=1;i<=n;i++){
        ip=NINT(pm[i][3]);
        pm[ip][1]=x[i-1];
        if(iw>0) printf("|FRCK> p[%d]=%f\n",i,x[i-1]);
    }
    double chi2=.0;
    double fredeg=double(m-n);
    //computing fvec
    int ivec=0;
    double wl,ev;
    for(int i=1;i<=mwl;i++){
        wl=MIS[7][i][1];
        ev=1240./wl;
        par[24][2]=i;
        FDISP(ioptf,ev);
        MIS[16][i][1]=sqn;
        CNK[1][2]=sqn;
        CNK[1][3]=MIS[16][i][2];
        // computing k from T
        //DELTA_k
        double Trasm=1.,dt,dlim,tol;;
        if(NINT(par[50+NINT(par[53][2])][3])==1){//incoherent layer
            if(DATO[1]>0) Trasm=MIS[1][i][1];
            if(DATO[2]>0) Trasm=MIS[2][i][1];
            double dinco=pm[NINT(par[53][2])][1];
            if(iw==1) printf("Trasm = %f\n",Trasm);
            dt=-0.1*wl/4./3.14/dinco*log(Trasm);
            dlim=1.e-10;//dlim
        }
        else{
            dt=0.001;
            dlim=1.e-5;//dlim
        }
        tol=dt/1000.;//tol
        double k=CNK[1][3];
        FM=FindRoot(&DELTAT,k,dt,dlim,tol);
        //printf("k_ini= %f dt=%f tol=%f k_fin=%f\n",k,dt,tol,CNK[1][3]);
        //printf("FM= %f\n",FM);
        MIS[16][i][2]=CNK[1][3];
        //computing R
        if(DATO[3]==2 || DATO[5]==2){
            ASSEMBLER(i,1,0.,vot);
            if(DATO[3]==2){
                fvec[ivec]=(MIS[3][i][1]-vot[2][s1p2])/MIS[3][i][2];
                chi2=chi2+pow(fvec[ivec],2.)/fredeg;
                ivec++;
            }
            if(DATO[5]==2){
                fvec[ivec]=(MIS[5][i][1]-vot[3][s1p2])/MIS[5][i][2];
                chi2=chi2+pow(fvec[ivec],2.)/fredeg;
                ivec++;
            }
        }
        if(DATO[4]==2){
            ASSEMBLER(i,1,te,vot);
            fvec[ivec]=(MIS[4][i][1]-vot[2][s1p2])/MIS[4][i][2];
            chi2=chi2+pow(fvec[ivec],2.)/fredeg;
            ivec++;
        }
        if(DATO[8]==2){
            ASSEMBLER(i,1,te,vot);
            fvec[ivec]=(ELI[2][i][1]-vot[5][1])/ELI[2][i][2];
            chi2=chi2+pow(fvec[ivec],2.)/fredeg;
            ivec++;
        }
        if(DATO[10]==2){
            ASSEMBLER(i,1,te,vot);
            fvec[ivec]=(ELI[4][i][1]-vot[5][1])/ELI[4][i][2];
            chi2=chi2+pow(fvec[ivec],2.)/fredeg;
            ivec++;
        }
        if(DATO[12]==2){
            ASSEMBLER(i,1,te,vot);
            fvec[ivec]=(ELI[6][i][1]-vot[5][1])/ELI[6][i][2];
            chi2=chi2+pow(fvec[ivec],2.)/fredeg;
            ivec++;
        }
        if(DATO[14]==2){
            ASSEMBLER(i,1,te,vot);
            fvec[ivec]=(ELI[8][i][1]-vot[5][1])/ELI[8][i][2];
            chi2=chi2+pow(fvec[ivec],2.)/fredeg;
            ivec++;
        }
        //printf("wl=%f n=%f k=%f chi2=%f\n",wl,CNK[1][2],CNK[1][3],chi2);
        //printf("continue?\n");
        //char rispo;
        //cin>>rispo;
    }
    // monitor values
    printf("|FRCK> chi2=%g\t",chi2);
    for(int i=0;i<n;i++)
        printf("p[%d]=%g\t",i,x[i]);
    printf("\n");
    if(iRecChi2==1){
        QFile fchi2(filechi2);
        if (fchi2.open(QIODevice::WriteOnly | QIODevice::Text | QIODevice::Append)){
            QTextStream out(&fchi2);
            out<<QString::number(chi2,'g',9);
            for(int i=0;i<n;i++)
                out<<"\t"<<QString::number(x[i],'g',9);
            out<<"\n";
            fchi2.close();
        }
    }
    QCoreApplication::processEvents();
    if(iStop){
        iflag=-1;
        return(iflag);
    }
    return(0);
}

double DELTAT(double k){
    //subroutine for computing the difference Texp and Tcalc
    int s1p2=1;
    double texp=0.,tcalc,te=0.,vot[6][3];
    int iwl=NINT(par[24][2]);
    if(DATO[1]==2){
        texp=MIS[1][iwl][1];
        te=0.;
        s1p2=1;
    }
    else if(DATO[2]==2){
        texp=MIS[2][iwl][1];
        te=par[6][1]*deg2rad;
        s1p2=NINT(par[27][2]);
    }
    CNK[1][2]=sqn;
    CNK[1][3]=k;
    ASSEMBLER(iwl,1,te,vot);
    tcalc=vot[1][s1p2];
    double dT=tcalc-texp;
    //        printf("|DELTAT> wl = %f  iwl = %d  Texp = %f\n",wl,iwl,texp);
    //        printf("\t n = %f  k = %f\n",CNK[1][2],k);
    //        printf("\t Tcalc = %f  dT=%f\n",tcalc,dT);
    return(dT);
}


double FindRoot(double (*FT)(double),double t,double dt,double dlim,double tol){
    int iWarning=0;
    double t0,d1,d0;
    //char ch;

    d1=FT(t);
    t0=t;
    d0=d1;
    int istop=0;
    int irif=0;
    int inul=0;
    int imi=0;
    int ilim=1000;
    if(iWarning==1)
        printf("\n FindRoot: t=%e dt=%e dlim=%e\n",t,dt,dlim);
    if(abs(d1) <= dlim) istop=1;
    while(istop == 0){
        t=t+dt;
        d1=FT(t);
        if(iWarning==1){
            printf("t: %e -> %e dt=%e\td: %e -> %e DELTAd=%e\n",t0,t,dt,d0,d1,d1-d0);
            //           cout << "To continue press a key & enter: ";
            //           cin  >> ch;
        }
        if(abs(d1) <= dlim){
            istop=1;
            //    printf("  STOP per raggiunto limite precisione!\n");
        }
        else if(d0*d1 <.0 && abs(d1) >= (10.*abs(d0))){
            t=t0;
            d1=d0;
            dt=dt/10.;
            //    printf("dt decimato = %e\n",dt);
            if(abs(dt) <= dlim){
                istop=1;
                //      printf(" STOP per raggiunto limite precisione!\n");
            }
        }
        else if(d0*d1 <.0 && abs(d1) < (10.*abs(d0))) {// Brent
            double a=t;
            double b=t0;
            double c=b;
            double fa=d1;
            double fb=d0;
            double fc=d0;
            int itmax=1000;
            double eps=3.e-9;
            double xm,ss,pp,q,rr,e=0.,d=0.,tol1;
            int iter=0;
            while(istop == 0.){
                iter=iter+1;
                if((fb>0. && fc>0.) || (fb<0. && fc<0.)){
                    c=a; //Rename a, b, c and adjust bounding interval d.
                    fc=fa;
                    d=b-a;
                    e=d;
                }
                if(abs(fc) < abs(fb)){
                    a=b;
                    b=c;
                    c=a;
                    fa=fb;
                    fb=fc;
                    fc=fa;
                }
                //      printf("\nBrent: N_iteration=%d \n",iter);
                //   "a,b,c=",a,b,c
                //   "fa,fb,fc=",fa,fb,fc
                tol1=2.*eps*abs(b)+0.5*tol; // Convergence check.
                xm=.5*(c-b);
                if(abs(xm)<=tol1 || fb==0. || iter>=itmax){
                    //	printf("eps=%e tol=%e tol1=%e \n xm=%e\n",eps,tol,tol1,xm);
                    t=b;
                    d1=fb;
                    istop=1;
                    //        if(iter<itmax)
                    //          printf("STOP per raggiunto limite di precisione\n");
                    //        else
                    //          printf("STOP per raggiunto N.max iterazioni\n");
                }
                if(istop==0){
                    if(abs(e)>=tol1 && abs(fa)>abs(fb)){
                        ss=fb/fa; //Attempt inverse quadratic interpolation.
                        if(a==c){
                            pp=2.*xm*ss;
                            q=1.-ss;
                        }
                        else{
                            q=fa/fc;
                            rr=fb/fc;
                            pp=ss*(2.*xm*q*(q-rr)-(b-a)*(rr-1.));
                            q=(q-1.)*(rr-1.)*(ss-1.);
                        }
                        if(pp>0.) q=-q; //Check whether in bounds.
                        pp=abs(pp);
                        if(2.*pp<min(3.*xm*q-abs(tol1*q),abs(e*q))){
                            e=d; //Accept interpolation.
                            d=pp/q;
                            //            printf("...quadratic interpolation...\n");
                        }
                        else{
                            //            printf("...bisection..\n");
                            d=xm; //Interpolation failed, use bisection.
                            e=d;
                        }
                    }
                    else{ //Bounds decreasing too slowly, use bisection.
                        //          printf("...bisection to speed up..\n");
                        d=xm;
                        e=d;
                    }
                    a=b; //Move last best guess to a.
                    fa=fb;
                    if(abs(d)>tol1) //Evaluate new trial root.
                        b=b+d;
                    else{
                        if(xm>0)
                            b=b+tol1;
                        else
                            b=b-tol1;
                    }
                    t=b;
                    fb=FT(t);
                    //        printf("'b=%e,fb=%e \n",b,fb);
                }
            }
        }
        else if(abs(d1)<abs(d0)) {
            //      printf("  -> better\n");
            t0=t;
            d0=d1;
            irif=0;
            imi=imi+1;
        }
        else if(abs(d1)>abs(d0) && d0*d1>0.){
            //      printf("  -> worst\n");
            t=t0;
            dt=-dt;
            irif=irif+1;
            imi=0;
            if(irif>2){
                if(abs(d0-d1)>0.3*dlim){
                    dt=dt/2.;
                    //          printf("step->step/2 !\n");
                }
                else
                    irif=ilim+10;
            }
            d1=d0;
        }
        else
            inul=inul+1;
        if(irif>ilim || imi>ilim || inul>ilim){
            //      printf("C0: irif o imi o inul > 1000! -> stop\n");
            istop=1;
        }
    }
    d1=FT(t);
    if(iWarning==1)
        printf("Root -> t=%e d1=%e\n",t,d1);
    //   std::cout << "To continue press a key & enter: ";
    //   std::cin  >> ch;
    return(d1);
}



int NINT(double x){
    int n;
    if(x>=0.)
        n=int(x+0.5);
    else
        n=int(x-0.5);
    return(n);
}


void nextColor(){
    iColor++;
    if(iColor>=7)
        iColor=1;
}



void MATINV(int n, int np, double **as, double **b){
    /* subroutine for computing the inverse matrix by subroutines given by Numerical Recipe
      n = dimension of the matrix to be inverted
      np = max dimension of matrices used in called subroutine
      as(n,n) = matrix to be inverted
      b(n,n) = inverted matrix
    */
    //printf("-> MATINV\n");
    int indx[np];
    double y[np][np];
    double **a;
    a=new double *[np];
    for(int i=0;i<np;i++)
        a[i]=new double[np];

    // copy as-matrix into a-matrix that will be lost
    for(int i=0;i<n;i++){
        for(int j=0;j<n;j++)
            a[i][j]=as[i][j];
    }

    //Set up identity matrix.
    for(int i=0;i<n;i++){
        for(int j=0;j<n;j++)
            y[i][j]=0.;
        y[i][i]=1.;
    }
    double d=0;
    ludcmp(np,a,n,indx,d);// Decompose the matrix just once.
    if(NINT(par[11][2])==1)
        return;
    double *yCol;
    yCol=new double [n];
    for(int j=0;j<n;j++){//  Find inverse by columns.
        for(int j2=0;j2<n;j2++)
            yCol[j2]=y[j2][j];
        lubksb(np,a,n,indx,yCol);
        for(int j2=0;j2<n;j2++)
            y[j2][j]=yCol[j2];
    }
    // Note that FORTRAN stores two-dimensional matrices by column,
    // so y(1,j) is the address of the jth column of y.

    // save results into b-matrix
    for(int ii=0;ii<n;ii++){
        for(int jj=0;jj<n;jj++){
            b[ii][jj]=y[ii][jj];
        }
    }
    //printf("\tMATINV completed\n");
    //fflush(stdout);
    return;
}




void lubksb(int np,double **a,int n, int *indx,double *b){
    /*  Solves the set of n linear equations A * X = B. Here a is input,
*** not as the matrix A but rather as its LU decomposition, determined
*** by the routine ludcmp. indx is input as the permutation vector
*** returned by ludcmp. b(1:n) is input as the right-hand side vector B,
*** and returns with the solution vector X. a, n, np, and indx are not
*** modified by this routine and can be left in place for successive
*** calls with different right-hand sides b. This routine takes into
*** account the possibility that b will begin with many zero elements,
*** so it is e cient for use in matrix inversion.
*/
    //printf("-> lubksb\n");
    int i,ii,j,ll;
    double sum;
    ii=0;
    /* When ii is set to a positive value, it will become the index
    of the first nonvanishing element of b. We now do
    the forward substitution, equation (2.3.6). The only new
    wrinkle is to unscramble the permutation as we go.
    */
    for(i=0;i<n;i++){
        ll=indx[i];
        sum=b[ll];
        b[ll]=b[i];
        if(ii!=0){
            for(j=ii;j<i-1;j++)
                sum=sum-a[i][j]*b[j];
        }
        else if(sum!=0.)
            ii=i;// A nonzero element was encountered, so from now on we
        // will have to do the sums in the loop above.
        b[i]=sum;
    }
    for(i=n-1;i>=0;i--){//Now we do the backsubstitution, equation (2.3.7).
        sum=b[i];
        for(j=i+1;j<n;j++)
            sum=sum-a[i][j]*b[j];
        b[i]=sum/a[i][i];//Store a component of the solution vector X.
    }
    return;//All done!
}



void ludcmp(int np,double **a,int n,int *indx,double d){
    int NMAX=500;
    double TINY=1.0e-20;
    /*  From "Numerical recipe"
*** Largest expected n, and a small number.
*** Given a matrix a(1:n,1:n), with physical dimension np by np,
*** this routine replaces it by the LU decomposition of a rowwise
*** permutation of itself. a and n are input.a is output, arranged
*** as in equation (2.3.14) above; indx(1:n) is an output vector
*** that records the row permutation effected by the partial pivoting;
*** d is output as +-1 depending on whether the number of row interchanges
*** was even or odd,respectively. This routine is used in combination with
*** lubksb to solve linear equations or invert a matrix.
*/
    //printf("-> ludcmp\n");
    int i,imax,j,k;
    double aamax,dum,sum,vv[NMAX];// vv stores the implicit scaling of each row.
    d=1.;// No row interchanges yet.
    for(i=0;i<n;i++){// Loop over rows to get the implicit scaling information.
        aamax=0.;
        for(j=0;j<n;j++)
            if(abs(a[i][j])>aamax) aamax=abs(a[i][j]);
        if(aamax==0.){
            printf("singular matrix in ludcmp!!!!\n");
            par[11][2]=1.;
            return;
        }
        vv[i]=1./aamax;//Save the scaling.
    }
    for(j=0;j<n;j++){// This is the loop over columns of Crout's method.
        for(i=0;i<j-1;i++){// This is equation (2.3.12) except for i=j.
            sum=a[i][j];
            for(k=0;k<i-1;k++)
                sum=sum-a[i][k]*a[k][j];
            a[i][j]=sum;
        }
        aamax=0.;// Initialize for the search for largest pivot element.
        for(i=j;i<n;i++){// This is i=j of equation (2.3.12) and i=j+1...N of equation (2.3.13).
            sum=a[i][j];
            for(k=0;k<j-1;k++)
                sum=sum-a[i][k]*a[k][j];
            a[i][j]=sum;
            dum=vv[i]*abs(sum);// Figure of merit for the pivot.
            if(dum>=aamax){// Is it better than the best so far?
                imax=i;
                aamax=dum;
            }
        }
        if(j!=imax){// Do we need to interchange rows?
            for(k=0;k<n;k++){// Yes, do so...
                dum=a[imax][k];
                a[imax][k]=a[j][k];
                a[j][k]=dum;
            }
            d=-d;         //         ...and change the parity of d
            vv[imax]=vv[j];//         Also interchange the scale factor.
        }
        indx[j]=imax;
        if(a[j][j]==0.) a[j][j]=TINY;
        /* If the pivot element is zero the matrix is singular (at least to the precision of the algorithm).
   For some applications on singular matrices, it is desirable to substitute TINY for zero.
   */
        if(j!=n){//         Now,finally, divide by the pivot element.
            dum=1./a[j][j];
            for(i=j+1;i<n;i++)
                a[i][j]=a[i][j]*dum;
        }
    }//! Go back for the next column in the reduction.
    return;
}



/*void AbsLayInco(double atten, double Ra, double R1a, double Rb){
    //absorptance of a not coherent layer
    double uguale=(1.-R1a*Rb*atten*atten);
    Vservice[0]=(1.-Ra)*(1.-atten)*(1.+Rb*atten)/uguale;//Absorptance
    Vservice[1]=(1.-Ra)*atten/uguale;//Transmittance at the entrance of the "b" interface
    return;
}
*/


//save oscillator parameters at any value change
void ksemawc::on_LEpm_102_1_textChanged(const QString &arg1) {SaveSetting(2);}
void ksemawc::on_LEpm_103_1_textChanged(const QString &arg1) {SaveSetting(2);}
void ksemawc::on_LEpm_104_1_textChanged(const QString &arg1) {SaveSetting(2);}
void ksemawc::on_LEpm_105_1_textChanged(const QString &arg1) {SaveSetting(2);}
void ksemawc::on_LEpm_107_1_textChanged(const QString &arg1) {SaveSetting(2);}
void ksemawc::on_LEpm_108_1_textChanged(const QString &arg1) {SaveSetting(2);}
void ksemawc::on_LEpm_109_1_textChanged(const QString &arg1) {SaveSetting(2);}
void ksemawc::on_LEpm_110_1_textChanged(const QString &arg1) {SaveSetting(2);}

void ksemawc::on_LEpm_112_1_textChanged(const QString &arg1) {SaveSetting(2);}
void ksemawc::on_LEpm_113_1_textChanged(const QString &arg1) {SaveSetting(2);}
void ksemawc::on_LEpm_114_1_textChanged(const QString &arg1) {SaveSetting(2);}
void ksemawc::on_LEpm_115_1_textChanged(const QString &arg1) {SaveSetting(2);}
void ksemawc::on_LEpm_117_1_textChanged(const QString &arg1) {SaveSetting(2);}
void ksemawc::on_LEpm_118_1_textChanged(const QString &arg1) {SaveSetting(2);}
void ksemawc::on_LEpm_119_1_textChanged(const QString &arg1) {SaveSetting(2);}
void ksemawc::on_LEpm_120_1_textChanged(const QString &arg1) {SaveSetting(2);}

void ksemawc::on_LEpm_122_1_textChanged(const QString &arg1) {SaveSetting(2);}
void ksemawc::on_LEpm_123_1_textChanged(const QString &arg1) {SaveSetting(2);}
void ksemawc::on_LEpm_124_1_textChanged(const QString &arg1) {SaveSetting(2);}
void ksemawc::on_LEpm_125_1_textChanged(const QString &arg1) {SaveSetting(2);}
void ksemawc::on_LEpm_127_1_textChanged(const QString &arg1) {SaveSetting(2);}
void ksemawc::on_LEpm_128_1_textChanged(const QString &arg1) {SaveSetting(2);}
void ksemawc::on_LEpm_129_1_textChanged(const QString &arg1) {SaveSetting(2);}
void ksemawc::on_LEpm_130_1_textChanged(const QString &arg1) {SaveSetting(2);}

void ksemawc::on_LEpm_132_1_textChanged(const QString &arg1) {SaveSetting(2);}
void ksemawc::on_LEpm_133_1_textChanged(const QString &arg1) {SaveSetting(2);}
void ksemawc::on_LEpm_134_1_textChanged(const QString &arg1) {SaveSetting(2);}
void ksemawc::on_LEpm_135_1_textChanged(const QString &arg1) {SaveSetting(2);}
void ksemawc::on_LEpm_137_1_textChanged(const QString &arg1) {SaveSetting(2);}
void ksemawc::on_LEpm_138_1_textChanged(const QString &arg1) {SaveSetting(2);}
void ksemawc::on_LEpm_139_1_textChanged(const QString &arg1) {SaveSetting(2);}
void ksemawc::on_LEpm_140_1_textChanged(const QString &arg1) {SaveSetting(2);}

void ksemawc::on_LEpm_142_1_textChanged(const QString &arg1) {SaveSetting(2);}
void ksemawc::on_LEpm_143_1_textChanged(const QString &arg1) {SaveSetting(2);}
void ksemawc::on_LEpm_144_1_textChanged(const QString &arg1) {SaveSetting(2);}
void ksemawc::on_LEpm_145_1_textChanged(const QString &arg1) {SaveSetting(2);}
void ksemawc::on_LEpm_147_1_textChanged(const QString &arg1) {SaveSetting(2);}
void ksemawc::on_LEpm_148_1_textChanged(const QString &arg1) {SaveSetting(2);}
void ksemawc::on_LEpm_149_1_textChanged(const QString &arg1) {SaveSetting(2);}
void ksemawc::on_LEpm_150_1_textChanged(const QString &arg1) {SaveSetting(2);}

void ksemawc::on_LEpm_152_1_textChanged(const QString &arg1) {SaveSetting(2);}
void ksemawc::on_LEpm_153_1_textChanged(const QString &arg1) {SaveSetting(2);}
void ksemawc::on_LEpm_154_1_textChanged(const QString &arg1) {SaveSetting(2);}
void ksemawc::on_LEpm_155_1_textChanged(const QString &arg1) {SaveSetting(2);}
void ksemawc::on_LEpm_157_1_textChanged(const QString &arg1) {SaveSetting(2);}
void ksemawc::on_LEpm_158_1_textChanged(const QString &arg1) {SaveSetting(2);}
void ksemawc::on_LEpm_159_1_textChanged(const QString &arg1) {SaveSetting(2);}
void ksemawc::on_LEpm_160_1_textChanged(const QString &arg1) {SaveSetting(2);}

void ksemawc::on_LEpm_162_1_textChanged(const QString &arg1) {SaveSetting(2);}
void ksemawc::on_LEpm_163_1_textChanged(const QString &arg1) {SaveSetting(2);}
void ksemawc::on_LEpm_164_1_textChanged(const QString &arg1) {SaveSetting(2);}
void ksemawc::on_LEpm_165_1_textChanged(const QString &arg1) {SaveSetting(2);}
void ksemawc::on_LEpm_167_1_textChanged(const QString &arg1) {SaveSetting(2);}
void ksemawc::on_LEpm_168_1_textChanged(const QString &arg1) {SaveSetting(2);}
void ksemawc::on_LEpm_169_1_textChanged(const QString &arg1) {SaveSetting(2);}
void ksemawc::on_LEpm_170_1_textChanged(const QString &arg1) {SaveSetting(2);}

void ksemawc::on_LEpm_172_1_textChanged(const QString &arg1) {SaveSetting(2);}
void ksemawc::on_LEpm_173_1_textChanged(const QString &arg1) {SaveSetting(2);}
void ksemawc::on_LEpm_174_1_textChanged(const QString &arg1) {SaveSetting(2);}
void ksemawc::on_LEpm_175_1_textChanged(const QString &arg1) {SaveSetting(2);}
void ksemawc::on_LEpm_177_1_textChanged(const QString &arg1) {SaveSetting(2);}
void ksemawc::on_LEpm_178_1_textChanged(const QString &arg1) {SaveSetting(2);}
void ksemawc::on_LEpm_179_1_textChanged(const QString &arg1) {SaveSetting(2);}
void ksemawc::on_LEpm_180_1_textChanged(const QString &arg1) {SaveSetting(2);}

void ksemawc::on_LEpm_182_1_textChanged(const QString &arg1) {SaveSetting(2);}
void ksemawc::on_LEpm_183_1_textChanged(const QString &arg1) {SaveSetting(2);}
void ksemawc::on_LEpm_184_1_textChanged(const QString &arg1) {SaveSetting(2);}
void ksemawc::on_LEpm_185_1_textChanged(const QString &arg1) {SaveSetting(2);}
void ksemawc::on_LEpm_187_1_textChanged(const QString &arg1) {SaveSetting(2);}
void ksemawc::on_LEpm_188_1_textChanged(const QString &arg1) {SaveSetting(2);}
void ksemawc::on_LEpm_189_1_textChanged(const QString &arg1) {SaveSetting(2);}
void ksemawc::on_LEpm_190_1_textChanged(const QString &arg1) {SaveSetting(2);}

void ksemawc::on_LEpm_192_1_textChanged(const QString &arg1) {SaveSetting(2);}
void ksemawc::on_LEpm_193_1_textChanged(const QString &arg1) {SaveSetting(2);}
void ksemawc::on_LEpm_194_1_textChanged(const QString &arg1) {SaveSetting(2);}
void ksemawc::on_LEpm_195_1_textChanged(const QString &arg1) {SaveSetting(2);}
void ksemawc::on_LEpm_197_1_textChanged(const QString &arg1) {SaveSetting(2);}
void ksemawc::on_LEpm_198_1_textChanged(const QString &arg1) {SaveSetting(2);}
void ksemawc::on_LEpm_199_1_textChanged(const QString &arg1) {SaveSetting(2);}
void ksemawc::on_LEpm_200_1_textChanged(const QString &arg1) {SaveSetting(2);}

//save model parameters at any value change
void ksemawc::on_dSB_PM_1_1_valueChanged(const QString &arg1) {SaveSetting(1);}
void ksemawc::on_dSB_PM_2_1_valueChanged(const QString &arg1) {SaveSetting(1);}
void ksemawc::on_dSB_PM_3_1_valueChanged(const QString &arg1) {SaveSetting(1);}
void ksemawc::on_dSB_PM_4_1_valueChanged(const QString &arg1) {SaveSetting(1);}
void ksemawc::on_dSB_PM_5_1_valueChanged(const QString &arg1) {SaveSetting(1);}
void ksemawc::on_dSB_PM_6_1_valueChanged(const QString &arg1) {SaveSetting(1);}
void ksemawc::on_dSB_PM_7_1_valueChanged(const QString &arg1) {SaveSetting(1);}
void ksemawc::on_dSB_PM_8_1_valueChanged(const QString &arg1) {SaveSetting(1);}
void ksemawc::on_dSB_PM_9_1_valueChanged(const QString &arg1) {SaveSetting(1);}

void ksemawc::on_dSB_PM_51_1_valueChanged(const QString &arg1) {SaveSetting(1);}
void ksemawc::on_dSB_PM_52_1_valueChanged(const QString &arg1) {SaveSetting(1);}
void ksemawc::on_dSB_PM_53_1_valueChanged(const QString &arg1) {SaveSetting(1);}
void ksemawc::on_dSB_PM_54_1_valueChanged(const QString &arg1) {SaveSetting(1);}
void ksemawc::on_dSB_PM_55_1_valueChanged(const QString &arg1) {SaveSetting(1);}
void ksemawc::on_dSB_PM_56_1_valueChanged(const QString &arg1) {SaveSetting(1);}
void ksemawc::on_dSB_PM_57_1_valueChanged(const QString &arg1) {SaveSetting(1);}
void ksemawc::on_dSB_PM_58_1_valueChanged(const QString &arg1) {SaveSetting(1);}
void ksemawc::on_dSB_PM_59_1_valueChanged(const QString &arg1) {SaveSetting(1);}

void ksemawc::on_dSB_PM_91_1_valueChanged(const QString &arg1) {SaveSetting(1);}
void ksemawc::on_dSB_PM_92_1_valueChanged(const QString &arg1) {SaveSetting(1);}
void ksemawc::on_dSB_PM_93_1_valueChanged(const QString &arg1) {SaveSetting(1);}
void ksemawc::on_dSB_PM_94_1_valueChanged(const QString &arg1) {SaveSetting(1);}
void ksemawc::on_dSB_PM_95_1_valueChanged(const QString &arg1) {SaveSetting(1);}
void ksemawc::on_dSB_PM_96_1_valueChanged(const QString &arg1) {SaveSetting(1);}
void ksemawc::on_dSB_PM_97_1_valueChanged(const QString &arg1) {SaveSetting(1);}
void ksemawc::on_dSB_PM_98_1_valueChanged(const QString &arg1) {SaveSetting(1);}
void ksemawc::on_dSB_PM_99_1_valueChanged(const QString &arg1) {SaveSetting(1);}

void ksemawc::on_dSB_PM_11_1_valueChanged(const QString &arg1) {SaveSetting(1);}
void ksemawc::on_dSB_PM_12_1_valueChanged(const QString &arg1) {SaveSetting(1);}
void ksemawc::on_dSB_PM_13_1_valueChanged(const QString &arg1) {SaveSetting(1);}
void ksemawc::on_dSB_PM_14_1_valueChanged(const QString &arg1) {SaveSetting(1);}
void ksemawc::on_dSB_PM_15_1_valueChanged(const QString &arg1) {SaveSetting(1);}
void ksemawc::on_dSB_PM_16_1_valueChanged(const QString &arg1) {SaveSetting(1);}
void ksemawc::on_dSB_PM_17_1_valueChanged(const QString &arg1) {SaveSetting(1);}
void ksemawc::on_dSB_PM_18_1_valueChanged(const QString &arg1) {SaveSetting(1);}
void ksemawc::on_dSB_PM_19_1_valueChanged(const QString &arg1) {SaveSetting(1);}

void ksemawc::on_dSB_PM_61_1_valueChanged(const QString &arg1) {SaveSetting(1);}
void ksemawc::on_dSB_PM_62_1_valueChanged(const QString &arg1) {SaveSetting(1);}
void ksemawc::on_dSB_PM_63_1_valueChanged(const QString &arg1) {SaveSetting(1);}
void ksemawc::on_dSB_PM_64_1_valueChanged(const QString &arg1) {SaveSetting(1);}
void ksemawc::on_dSB_PM_65_1_valueChanged(const QString &arg1) {SaveSetting(1);}
void ksemawc::on_dSB_PM_66_1_valueChanged(const QString &arg1) {SaveSetting(1);}
void ksemawc::on_dSB_PM_67_1_valueChanged(const QString &arg1) {SaveSetting(1);}
void ksemawc::on_dSB_PM_68_1_valueChanged(const QString &arg1) {SaveSetting(1);}
void ksemawc::on_dSB_PM_69_1_valueChanged(const QString &arg1) {SaveSetting(1);}

void ksemawc::on_dSB_PM_21_1_valueChanged(const QString &arg1) {SaveSetting(1);}
void ksemawc::on_dSB_PM_22_1_valueChanged(const QString &arg1) {SaveSetting(1);}
void ksemawc::on_dSB_PM_23_1_valueChanged(const QString &arg1) {SaveSetting(1);}
void ksemawc::on_dSB_PM_24_1_valueChanged(const QString &arg1) {SaveSetting(1);}
void ksemawc::on_dSB_PM_25_1_valueChanged(const QString &arg1) {SaveSetting(1);}
void ksemawc::on_dSB_PM_26_1_valueChanged(const QString &arg1) {SaveSetting(1);}
void ksemawc::on_dSB_PM_27_1_valueChanged(const QString &arg1) {SaveSetting(1);}
void ksemawc::on_dSB_PM_28_1_valueChanged(const QString &arg1) {SaveSetting(1);}
void ksemawc::on_dSB_PM_29_1_valueChanged(const QString &arg1) {SaveSetting(1);}

void ksemawc::on_dSB_PM_31_1_valueChanged(const QString &arg1) {SaveSetting(1);}
void ksemawc::on_dSB_PM_32_1_valueChanged(const QString &arg1) {SaveSetting(1);}
void ksemawc::on_dSB_PM_33_1_valueChanged(const QString &arg1) {SaveSetting(1);}
void ksemawc::on_dSB_PM_34_1_valueChanged(const QString &arg1) {SaveSetting(1);}
void ksemawc::on_dSB_PM_35_1_valueChanged(const QString &arg1) {SaveSetting(1);}
void ksemawc::on_dSB_PM_36_1_valueChanged(const QString &arg1) {SaveSetting(1);}
void ksemawc::on_dSB_PM_37_1_valueChanged(const QString &arg1) {SaveSetting(1);}
void ksemawc::on_dSB_PM_38_1_valueChanged(const QString &arg1) {SaveSetting(1);}
void ksemawc::on_dSB_PM_39_1_valueChanged(const QString &arg1) {SaveSetting(1);}

void ksemawc::on_dSB_PM_41_1_valueChanged(const QString &arg1) {SaveSetting(1);}
void ksemawc::on_dSB_PM_42_1_valueChanged(const QString &arg1) {SaveSetting(1);}
void ksemawc::on_dSB_PM_43_1_valueChanged(const QString &arg1) {SaveSetting(1);}
void ksemawc::on_dSB_PM_44_1_valueChanged(const QString &arg1) {SaveSetting(1);}
void ksemawc::on_dSB_PM_45_1_valueChanged(const QString &arg1) {SaveSetting(1);}
void ksemawc::on_dSB_PM_46_1_valueChanged(const QString &arg1) {SaveSetting(1);}
void ksemawc::on_dSB_PM_47_1_valueChanged(const QString &arg1) {SaveSetting(1);}
void ksemawc::on_dSB_PM_48_1_valueChanged(const QString &arg1) {SaveSetting(1);}
void ksemawc::on_dSB_PM_49_1_valueChanged(const QString &arg1) {SaveSetting(1);}

