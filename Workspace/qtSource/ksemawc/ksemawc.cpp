/*Main author: Marco Montecchi
             ENEA (Italy)
             email: marco.montecchi@enea.it
Porting to Windows and advanced oscillators by
             Alberto Mittiga
             ENEA (Italy)
             email: alberto.mittiga@enea.it


kSEMAWc — C++/Qt5 port of kSEMAW
Workspace for the analysis of
Spectrophotometric (SP), Ellipsometric (ELI) and
Photothermal Deflection Spectroscopy (PDS) measurements

Milestones of version history:
  v0.9.6 First public release: kernel wrote in FORTRAN and GUI based on C++/Qt4 (Montecchi)
  v1.0.0 Complete rewrite in C++ (Montecchi) and MS Windows executable (Mittiga)
  v2.x   Additional features (Montecchi) and new physical oscillators (Mittiga)
  v3.x   Improvements of loading and processing of ellipsometric spectra
  v4.0   Modern C++ refactoring and new pseudo-oscillator Volumetric-scattering
  v5.0   With the help of ClaudeCode, materials extended to 99 (user 1-99, special media 100-103),
         NK files extended to 99, [EXTRA_NK_V2] + [EXTRA_CNK_V2] file
         blocks, sB_nmat/sB_nnk computed from loaded data


   Copyright (C) 2026  Marco Montecchi

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
//C     [21][1],[21][2]....... / , /
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
//C     [35][1-    4]... standard-spectrum_for_weighting , S1P2U3 , / , /
//C
//C     [51-60][1]      pointer to VNK (material index for each layer)
//C         [51][2]     Nlayer [max N_LAYER_MAX=20, hard max N_LAYER_HARD_MAX=99]
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
//c **** per-layer parameter arrays (split from the original flat pm[201][6]) ****
//c
//c  Each array has shape [N_LAYER_HARD_MAX+1][6] (1-based, index 0 unused):
//c    pmD [i]  thickness            layer i
//c    pmDn[i]  Dn/<n>  n-gradient   layer i
//c    pmNc[i]  n-curvature          layer i
//c    pmDk[i]  Dk/<k>  k-gradient   layer i
//c    pmKc[i]  k-curvature          layer i
//c    pmRg[i]  roughness            layer i
//c    pmSd[i]  slope_Dn/<n> in 1/eV layer i
//c    pmNu[i]  non-uniformity (%)   layer i
//c    pmFe[j]  f-EMA fraction       material j  (index 1..N_CNK_HARD_MAX)
//c    pmTe[j]  ThetaEli             angle j     (index 1..4)
//c    pmOs[k]  oscillator params    (index 0..N_OSC_MAX*5)
//c
//c  All arrays store [col][1]=value, [col][2]=fit-link, [col][3]=ppm-ptr,
//c                   [col][4]=error,  [col][5]=global-corr
//c
//c  pmAt(ip): uniform flat ip index → pointer into named arrays.
//c  Used by the fit-parameter system (ppm[]). ip encoding (v5.2+):
//c
//c  Layer parameters — all 99 layers, 8 properties:
//c    ip = 1 + prop*99 + (layer-1),   layer=1..99,  prop=0..7
//c    prop  0=D      ip   1.. 99     pmD [1..99]
//c    prop  1=Dn     ip 100..198     pmDn[1..99]
//c    prop  2=Nc     ip 199..297     pmNc[1..99]
//c    prop  3=Dk     ip 298..396     pmDk[1..99]
//c    prop  4=Kc     ip 397..495     pmKc[1..99]
//c    prop  5=Rg     ip 496..594     pmRg[1..99]
//c    prop  6=Sd     ip 595..693     pmSd[1..99]
//c    prop  7=Nu     ip 694..792     pmNu[1..99]
//c
//c  Material / oscillator params:
//c    ip 793..891  f-EMA fraction  pmFe[1..99]     (ip = 792 + material)
//c    ip 892..895  ThetaEli 1..4   pmTe[1..4]      (ip = 891 + angle)
//c    ip 896       Fit#            pmOs[0]
//c    ip 897..996  oscillator params pmOs[1..100]   (ip = 896 + pmOs-index)
//c
//c  migrateIp(old_ip): converts any old encoding (stride-10, v4 layers 10-20,
//c  v5.1 layers 21-99) to the current uniform encoding.
//C
//C
//C
//C   ppm[nparmax=17] pointer pm stored in par:
//C                       PPM[i]=par[35+i][5] i=1,17
//C                       n_parameters_fit=par[35][5]
//C
//C
//C
//C  Experimental measurements are stored in
//C      ms.measures[i].value[iwl]
//C      ms.measures[i].error[iwl]
//C
//C       i=0: Tnormal
//C       i=1: Tpolarised
//C       i=2: Rnormal
//C       i=3: Rpolarised
//C       i=4: R1_normal
//C       i=5: Apds
//C       i=6: Delta1
//C       i=7: Psi1
//C       i=8: Delta2
//C       i=9: Psi2
//C       i=10: Delta3
//C       i=11: Psi3
//C       i=12: Delta4
//C       i=13: Psi4
//C
//C  NK-known data are stored in
//C      nkMaterials[k].n[iwl] with k=0,...,N_CNK_USER_MAX-1 (0..98): one slot per nk file
//C      nkMaterials[k].k[iwl]                  (sized N_CNK_USER_MAX=99 in SPADA)
//C      nk file #f (f=1..99) lives in nkMaterials[f-1]; its nSrc code is 7+f, so
//C      the slot is nkMaterials[nSrc-8]. This range (nSrc 8..106) is disjoint from the
//C      special input/output media cnk[100..103] (output/input SF, input PDS, input ELI)
//C      and from the ELI CONVER dispatch codes 107..114, so there is no overlap.
//C
//C  with the common wavelenght set stored in ms.lambda[iwl]
//C
//C
//C  The filenames are stored as follow:
//C   NANK[1,..,8] : nome "mate/aa999.9" nk-known
//C   NANK[9]      : nome "mate/aa999.9" nk-solution
//C   NANK[10]     : nome "mate/aa999  " custom weights for mean computing
//C   NANK[11]     : nome "mate/aa999" SF spectra
//C   NANK[12,.,15]: nome "mate/aa999" ELI spectra
//C   NANK[16]     : nome "mate/aa999.9" ksemaw project
//C
//C
//C  The instructions on how to assign n,k values to each layer are in:
//C   cnk[J]  (struct MaterialCNK, see ksemawc.h)
//C   J = 1 .. N_CNK_USER_MAX  (1..99)  : user materials
//C   J = N_CNK_SPEC_START .. N_CNK_MAX (100..103) : special input/output media:
//C       J=100 <-> output medium SF
//C       J=101 <-> input  medium SF
//C       J=102 <-> input  medium PDS
//C       J=103 <-> input  medium ELI
//C
//C   VNK[J][0..2] is built from cnk[J] by SETVNK for J=1..N_CNK_MAX:
//C       VNK[J][0] unused
//C       VNK[J][1] = n (real part)
//C       VNK[J][2] = k (imaginary part)
//C
//C   cnk[J].nSrc selects the n,k source (= combo index in the GUI):
//C       0        => constant: n=cnk[J].nConst, k=cnk[J].kConst
//C       1..7     => oscillator model (Fit# 1-7)
//C       8..7+N   => NK file: nkMaterials[nSrc-8]  (N = N_CNK_USER_MAX = 99,
//C                   so valid nSrc range for NK files is 8..106)
//C
//C   cnk[J].emaSrc < 0  => no EMA
//C   cnk[J].emaSrc >= 0 => EMA mixture: second material = emaSrc,
//C                          fill fraction = pmFe[J][1]
//C                          For J=1..20:  pmFe[J][1] stored in pmOs[70+J][1]
//C                          For J=21..N_CNK_MAX: pmFe[J][1] stored directly
//C
//C   cnk[J].forceMode:
//C       0 => no override
//C       1 => force n = cnk[J].nForced
//C       2 => force k = cnk[J].kForced
//C       3 => force both n and k
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
*/


#include "ksemawc.h"
#include <QListView>
#include <QDebug>
#include <cstdio>
#include <cstdlib>
#define _USE_MATH_DEFINES
#include <cmath>
#include <vector>
#include <qfile.h>
#include <qtextstream.h>
#include <complex>
#include <QFileDialog>
#include <QInputDialog>
#include <QMessageBox>
#include <QTimer>
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
#include <qwt_picker_machine.h>
#include <qwt_picker.h>
#include <qwt_plot_picker.h>
#include <qwt_plot_zoomer.h>
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

#ifdef __unix__
    static constexpr bool IS_POSIX = true;
#else
    static constexpr bool IS_POSIX = false;
#endif

// shortcut
inline int nint(double val) {
    return static_cast<int>(std::lround(val));
}

#ifndef M_PI
    #define M_PI 3.14159265358979323846
#endif
const double PIG=M_PI;
const double deg2rad=PIG/180.;
const double Dperiod=180.;//period of DELTA

//global variables
QString fnproject,pathroot,fileStore,fileStore0,fStdSpect,fRefMir,fNKsim,fMisSim,filechi2,lastAction,fileStoreSpjName,caller;
QString info,fnk[N_CNK_USER_MAX+1],fnSample,fnTn,fnTp,fnRn,fnRp,fnR1,fnApds,fnE1,fnE2,fnE3,fnE4,fnFnk,ParFitLab[14],NANK[17];

QPolygonF polygonF;
QColor myColor[7]={Qt::black,Qt::blue,Qt::cyan,Qt::green,Qt::magenta,Qt::red,Qt::yellow};
QString labelQwt;//label to use in absorptance plot
int iColor=0;
int iSelected=0;//used to stop recursive n-selection by polygon

int ink,lastIndex,ifn,npp,ppm[N_FIT_MAX+1],nPar,nlayer,lastTab,DATO[15],iw=0,L1E2=2,iRecChi2=0;
int ifirstcall=0;//used to initialize fit
int ifirstWarning=0;//used to warn about A= 1-T-R<0
int jobtot=0;
int nOpenGraph=3;
int NeV; // number of interpolated points in Ev or lambda
int iDelCon=0;
int iNewSample=0;
int iwl2print=0;//index of WL at which print when it is due

//attention, the index [0] is not used (FORTRAN-style 1-based indexing)
double pf[8][22],par[50+N_LAYER_HARD_MAX+1][6],rxy[31][5];
MaterialCNK cnk[N_CNK_HARD_MAX+1]; // 1-based: cnk[1..N_CNK_MAX] user+special media

// pm[] split into named per-property arrays (index 1..N_LAYER_HARD_MAX)
double pmD [N_LAYER_HARD_MAX+1][6]; // thickness
double pmDn[N_LAYER_HARD_MAX+1][6]; // Dn/<n>  gradient
double pmNc[N_LAYER_HARD_MAX+1][6]; // n-curvature
double pmDk[N_LAYER_HARD_MAX+1][6]; // Dk/<k>
double pmKc[N_LAYER_HARD_MAX+1][6]; // k-curvature
double pmRg[N_LAYER_HARD_MAX+1][6]; // roughness
double pmSd[N_LAYER_HARD_MAX+1][6]; // slope_Dn
double pmNu[N_LAYER_HARD_MAX+1][6]; // non-uniformity
double pmFe[N_CNK_HARD_MAX+1][6]; // f-EMA fraction  (index 1..N_CNK_HARD_MAX)
double pmTe[5][6];              // ThetaEli         (index 1..4)
double pmOs[N_OSC_MAX*5+1][6]; // oscillators: pmOs[iC+(k-1)*5] for osc k, param iC

// Uniform ip encoding — all layers 1..N_LAYER_HARD_MAX treated identically:
//   Layer params (8 props × 99 layers): ip = 1 + prop*99 + (layer-1),  ip 1..792
//     prop 0=D, 1=Dn, 2=Nc, 3=Dk, 4=Kc, 5=Rg, 6=Sd, 7=Nu
//   EMA fraction material 1..99:        ip = 792 + material,           ip 793..891
//   ThetaEli angle 1..4:                ip = 891 + angle,              ip 892..895
//   Fit# (pmOs[0]):                     ip = 896
//   Oscillator params (20×5=100):       ip = 896 + pmOs-index,         ip 897..996
static inline double* pmAt(int ip) {
    if (ip >= 1 && ip <= 792) {
        int prop  = (ip-1) / 99;
        int layer = (ip-1) % 99 + 1;
        switch(prop) {
            case 0: return pmD[layer];
            case 1: return pmDn[layer];
            case 2: return pmNc[layer];
            case 3: return pmDk[layer];
            case 4: return pmKc[layer];
            case 5: return pmRg[layer];
            case 6: return pmSd[layer];
            case 7: return pmNu[layer];
        }
    }
    if (ip >= 793 && ip <= 891) return pmFe[ip-792];   // EMA mat 1..99
    if (ip >= 892 && ip <= 895) return pmTe[ip-891];   // ThetaEli 1..4
    if (ip == 896)               return pmOs[0];         // Fit#
    if (ip >= 897 && ip <= 996)  return pmOs[ip-896];   // oscillator params
    static double dummy[6]{};
    return dummy;
}

// Migrate ip from any old encoding to the current uniform encoding.
// Called once at project load for each ppm[] value stored in par[35+i][5].
static int migrateIp(int old_ip) {
    if (old_ip <= 0) return 0;
    // stride-10 layers 1-9
    if (old_ip >= 1  && old_ip <= 9)   return 1 + 0*99 + (old_ip-1);    // D
    if (old_ip >= 11 && old_ip <= 19)  return 1 + 1*99 + (old_ip-11);   // Dn
    if (old_ip >= 21 && old_ip <= 29)  return 1 + 2*99 + (old_ip-21);   // Nc
    if (old_ip >= 31 && old_ip <= 39)  return 1 + 3*99 + (old_ip-31);   // Dk
    if (old_ip >= 41 && old_ip <= 49)  return 1 + 4*99 + (old_ip-41);   // Kc
    if (old_ip >= 51 && old_ip <= 59)  return 1 + 5*99 + (old_ip-51);   // Rg
    if (old_ip >= 61 && old_ip <= 69)  return 1 + 6*99 + (old_ip-61);   // Sd
    if (old_ip >= 71 && old_ip <= 85)  return 792 + (old_ip-70);         // EMA mat 1-15
    if (old_ip >= 86 && old_ip <= 89)  return 891 + (old_ip-85);         // ThetaEli 1-4
    if (old_ip == 90)                   return 792 + 20;                   // EMA mat 20
    if (old_ip >= 91 && old_ip <= 99)  return 1 + 7*99 + (old_ip-91);   // Nu layers 1-9
    if (old_ip == 100)                  return 896;                        // Fit#
    if (old_ip >= 101 && old_ip <= 200) return 896 + (old_ip-100);        // oscillators
    if (old_ip >= 201 && old_ip <= 288) {                                  // layers 10-20
        int off=old_ip-201, prop=off/11, lyr=off%11+10;
        return 1 + prop*99 + (lyr-1);
    }
    if (old_ip >= 289 && old_ip <= 920) {                                  // layers 21-99
        int off=old_ip-289, prop=off/79, lyr=off%79+21;
        return 1 + prop*99 + (lyr-1);
    }
    return old_ip; // already in current encoding or unknown
}

// The pm[1..200] block in .spj files always uses this fixed legacy encoding,
// regardless of file age. Layers 1-9 + EMA 1-15/20 + ThetaEli + oscillators.
// Layers 10-99 are NOT covered here — [STACK_V2] is authoritative for those.
static inline double* pmAtLegacy(int ip) {
    if (ip >= 1  && ip <=  9) return pmD[ip];
    if (ip >= 11 && ip <= 19) return pmDn[ip-10];
    if (ip >= 21 && ip <= 29) return pmNc[ip-20];
    if (ip >= 31 && ip <= 39) return pmDk[ip-30];
    if (ip >= 41 && ip <= 49) return pmKc[ip-40];
    if (ip >= 51 && ip <= 59) return pmRg[ip-50];
    if (ip >= 61 && ip <= 69) return pmSd[ip-60];
    if (ip >= 71 && ip <= 85) return pmFe[ip-70];   // EMA mat 1..15
    if (ip >= 86 && ip <= 89) return pmTe[ip-85];   // ThetaEli 1..4
    if (ip == 90)              return pmFe[20];       // EMA mat 20
    if (ip >= 91 && ip <= 99) return pmNu[ip-90];   // Nu layers 1-9
    if (ip == 100)             return pmOs[0];        // Fit#
    if (ip >= 101 && ip <= 200) return pmOs[ip-100]; // oscillators (20×5)
    static double dummy[6]{};
    return dummy;
}

// Build a Stack from the current global pm-arrays and par[50+i] for i=1..n.
// Called from ReadSetting after all pm data and the [EXTRA] block are loaded.
static Stack stackFromPm(int n){
    n = std::min(n, N_LAYER_HARD_MAX);
    Stack s;
    s.layers.resize(n);
    auto loadLP = [](LayerParam& lp, double (*arr)[6], int row){
        lp.value      = arr[row][1];
        lp.error      = arr[row][2];
        lp.globalCorr = arr[row][3];
        lp.fitIndex   = static_cast<int>(arr[row][4]);
    };
    for(int i = 1; i <= n; i++){
        Layer& l        = s.layers[i-1];
        l.materialIndex = nint(par[50+i][1]);
        l.type          = nint(par[50+i][3]) - 1;
        int mi          = l.materialIndex;
        l.materialEMA   = (mi >= 1 && mi <= N_CNK_MAX) ? nint(cnk[mi].emaSrc) : -1;
        loadLP(l.thickness,      pmD,  i);
        loadLP(l.nGrad,          pmDn, i);
        loadLP(l.nCurv,          pmNc, i);
        loadLP(l.kGrad,          pmDk, i);
        loadLP(l.kCurv,          pmKc, i);
        loadLP(l.roughness,      pmRg, i);
        loadLP(l.slopeNGrad,     pmSd, i);
        loadLP(l.nonUniformity,  pmNu, i);
        if(mi >= 1 && mi <= N_CNK_MAX)
            loadLP(l.fEMA, pmFe, mi);
    }
    return s;
}

// Sync pm-arrays and par[50+i] from the Stack (inverse of stackFromPm).
static void stackToPm(const Stack& s){
    auto saveLP = [](double (*arr)[6], int row, const LayerParam& lp){
        arr[row][1] = lp.value;
        arr[row][2] = lp.error;
        arr[row][3] = lp.globalCorr;
        arr[row][4] = static_cast<double>(lp.fitIndex);
    };
    for(int i = 1; i <= std::min(s.nLayers(), N_LAYER_HARD_MAX); i++){
        const Layer& l = s.layers[i-1];
        par[50+i][1] = l.materialIndex;
        par[50+i][3] = l.type + 1;
        saveLP(pmD,  i, l.thickness);
        saveLP(pmDn, i, l.nGrad);
        saveLP(pmNc, i, l.nCurv);
        saveLP(pmDk, i, l.kGrad);
        saveLP(pmKc, i, l.kCurv);
        saveLP(pmRg, i, l.roughness);
        saveLP(pmSd, i, l.slopeNGrad);
        saveLP(pmNu, i, l.nonUniformity);
    }
}

// Value-only inverse of stackToPm(): copy the per-layer scalar VALUES from the
// legacy pm* arrays back into the modern Stack. The fit-parameter editor in the
// Data Analysis tab writes new layer values into the pm* arrays via pmAt(), but
// the simulation engine (CalcMis) reads layer thickness/roughness/gradients from
// stack.layers[], so without this sync the edited layer fit values are ignored
// on "Simulate!" (oscillator params live in pmOs[], read directly by the engine,
// hence they already worked). Only column [1]=value is copied: columns [2..4]
// carry incompatible meanings between the two subsystems (fit-link in the pm/fit
// convention vs error/globalCorr/fitIndex in stackToPm's saveLP).
static void pmValuesToStack(Stack& s){
    for(int i = 1; i <= std::min(s.nLayers(), N_LAYER_HARD_MAX); i++){
        Layer& l = s.layers[i-1];
        l.thickness.value     = pmD [i][1];
        l.nGrad.value         = pmDn[i][1];
        l.nCurv.value         = pmNc[i][1];
        l.kGrad.value         = pmDk[i][1];
        l.kCurv.value         = pmKc[i][1];
        l.roughness.value     = pmRg[i][1];
        l.slopeNGrad.value    = pmSd[i][1];
        l.nonUniformity.value = pmNu[i][1];
    }
}

std::vector<std::array<double,6>> Sol;// solutions from the search algorithm
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
bool iShemispherical=false;

//invoked functions
void previewFile(QString filename, QString lab,QString& info,double& wmax,double& wmin);
int FPAR(void *p, int m, int n, const double *x, double *fvec, int iflag);
EmaResult EMA(double N,double K,double NA,double KA,double FA);
FdispResult FDISP(int iopt,double eV);
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
void CALFRE(int NFA,double wl,complex<double> pq,const std::vector<std::complex<double>>& ir,std::vector<double> d,complex<double> out[10][3]);
bool MATINV(int n, int np,double **as,double **b);
void lubksb(int np,double **a,int n,int *indx,double *b);
bool ludcmp(int np,double **a,int n,int *indx,double& d);
void nextColor();

// ── Layer widget name helpers ────────────────────────────────────────────────
// Layers 1-9  : old stride-10 scheme  (dSB_PM_${base+i}_1, comB_PAR_5${i}_N)
// Layers 10-20: semantic scheme       (dSB_L${i}_${prop}_1, comB_L${i}_N)
namespace {
    QString wD  (int i){ return i<=9 ? "dSB_PM_"+QString::number(i)+"_1"     : "dSB_L"+QString::number(i)+"_D_1";  }
    QString wNu (int i){ return i<=9 ? "dSB_PM_"+QString::number(90+i)+"_1"  : "dSB_L"+QString::number(i)+"_Nu_1"; }
    QString wRg (int i){ return i<=9 ? "dSB_PM_"+QString::number(50+i)+"_1"  : "dSB_L"+QString::number(i)+"_Rg_1"; }
    QString wDn (int i){ return i<=9 ? "dSB_PM_"+QString::number(10+i)+"_1"  : "dSB_L"+QString::number(i)+"_Dn_1"; }
    QString wNc (int i){ return i<=9 ? "dSB_PM_"+QString::number(20+i)+"_1"  : "dSB_L"+QString::number(i)+"_Nc_1"; }
    QString wDk (int i){ return i<=9 ? "dSB_PM_"+QString::number(30+i)+"_1"  : "dSB_L"+QString::number(i)+"_Dk_1"; }
    QString wKc (int i){ return i<=9 ? "dSB_PM_"+QString::number(40+i)+"_1"  : "dSB_L"+QString::number(i)+"_Kc_1"; }
    QString wSd (int i){ return i<=9 ? "dSB_PM_"+QString::number(60+i)+"_1"  : "dSB_L"+QString::number(i)+"_Sd_1"; }
    QString wCB1(int i){ return i<=9 ? "comB_PAR_5"+QString::number(i)+"_1"  : "comB_L"+QString::number(i)+"_1";   }
    QString wCB3(int i){ return i<=9 ? "comB_PAR_5"+QString::number(i)+"_3"  : "comB_L"+QString::number(i)+"_3";   }
    QString wBtn(const char* dir, int i){ return QString(dir)+QString::number(i); }
    // Copy src over dst with Qt only (no shell). The previous shell "copy"/"cp"
    // route via system() broke on paths containing '+' (cmd.exe COPY reads '+' as the
    // file-concatenation operator) or spaces. QFile::copy won't overwrite, so remove
    // the destination first. Returns true on success.
    bool copyFileOverwrite(const QString& src, const QString& dst){
        if(QFile::exists(dst) && !QFile::remove(dst))
            return false;
        return QFile::copy(src, dst);
    }
    // EMA spinbox key: materials 1-15 have static .ui widgets "dSB_PM_{70+i}_1"
    // (dSB_PM_71..85); materials 16+ are created dynamically and MUST use
    // "dSB_EMA_{i}_1". Reusing "dSB_PM_{70+i}_1" for i>=16 would clash with existing
    // static widgets: ThetaEli angles dSB_PM_86..89 (i=16..19) and layer-Nu
    // dSB_PM_91..99 (i=21..29). The clash silently overwrote the idToDoubleSpinBox
    // entry for the ThetaEli widgets, so saving the model zeroed the fitted angles.
    QString emaKey(int i){ return i<=15 ? "dSB_PM_"+QString::number(70+i)+"_1" : "dSB_EMA_"+QString::number(i)+"_1"; }

    // ip-encoding for PanFit: n2=0(Nu),1(D),2(Dn),3(Nc),4(Dk),5(Kc),6(Rg),7(Sd)
    // Uniform encoding — all layers 1..N_LAYER_HARD_MAX treated identically
    inline int ipLayer(int n1, int n2){
        int prop = (n2==0) ? 7 : (n2-1);   // D=0,Dn=1,Nc=2,Dk=3,Kc=4,Rg=5,Sd=6,Nu=7
        return 1 + prop*99 + (n1-1);
    }
}

ksemawc::ksemawc(QWidget *parent) :
    QWidget(parent),
    ui(new Ui::ksemawc)
{
    //setupUi(this); // this sets up GUI
    ui->setupUi(this);
    // default cnk[i].emaSrc=-1 means "no EMA" for all rows
    for(int i=1;i<=N_CNK_MAX;i++) cnk[i].emaSrc=-1.;

//    const QByteArray value = qgetenv("USER");
//    QString uName=QString::fromLocal8Bit(value);
//    cout << "current user = " << uName.toStdString() <<endl;

    // signals/slots mechanism in action
    connect( ui->pushButton_about,SIGNAL( clicked() ), this, SLOT(about()));
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
    connect( ui->sB_nnk, QOverload<int>::of(&QSpinBox::valueChanged), this, &ksemawc::updateNkRows);
    connect( ui->pBsetSample,SIGNAL( clicked() ), this, SLOT(callSetSample()));
    connect( ui->pBclearFn,SIGNAL( clicked() ), this, SLOT(Clrfn()));
    connect( ui->lineEdit_sample,SIGNAL(textChanged(QString)),this, SLOT(listMeas()));
    connect( ui->cBmis1,SIGNAL(currentIndexChanged(int)),this,SLOT(pwTn()));
    connect( ui->cBmis2,SIGNAL(currentIndexChanged(int)),this,SLOT(pwTp()));
    connect( ui->cBmis3,SIGNAL(currentIndexChanged(int)),this,SLOT(pwRn()));
    connect( ui->cBmis4,SIGNAL(currentIndexChanged(int)),this,SLOT(pwRp()));
    connect( ui->cBmis5,SIGNAL(currentIndexChanged(int)),this,SLOT(pwR1()));
    connect( ui->cBmis6,SIGNAL(currentIndexChanged(int)),this,SLOT(pwApds()));
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
    connect(ui->pushButton_saveNKmodel,SIGNAL(clicked()), this, SLOT(saveNKmodel()));
    connect(ui->pushButton_loadNKmodel,SIGNAL(clicked()), this, SLOT(loadNKmodel()));
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

    // ── Register layer 1-9 buttons in idToPushButton ─────────────────────────
    idToPushButton["pBmDw1"]=ui->pBmDw1; idToPushButton["pBmUp1"]=ui->pBmUp1;
    idToPushButton["pBmDw2"]=ui->pBmDw2; idToPushButton["pBmUp2"]=ui->pBmUp2;
    idToPushButton["pBmDw3"]=ui->pBmDw3; idToPushButton["pBmUp3"]=ui->pBmUp3;
    idToPushButton["pBmDw4"]=ui->pBmDw4; idToPushButton["pBmUp4"]=ui->pBmUp4;
    idToPushButton["pBmDw5"]=ui->pBmDw5; idToPushButton["pBmUp5"]=ui->pBmUp5;
    idToPushButton["pBmDw6"]=ui->pBmDw6; idToPushButton["pBmUp6"]=ui->pBmUp6;
    idToPushButton["pBmDw7"]=ui->pBmDw7; idToPushButton["pBmUp7"]=ui->pBmUp7;
    idToPushButton["pBmDw8"]=ui->pBmDw8; idToPushButton["pBmUp8"]=ui->pBmUp8;
    idToPushButton["pBmDw9"]=ui->pBmDw9; idToPushButton["pBmUp9"]=ui->pBmUp9;

    connect(ui->pushButton_saveSim, SIGNAL( clicked() ), this, SLOT(saveSim()));
    connect(ui->pushButton_saveNKsim,SIGNAL( clicked() ), this, SLOT(saveNKsim()));
    connect(ui->pushButton_PlotAbsEL,SIGNAL( clicked() ),this, SLOT(PlotAbsEL()));
    connect(ui->comboBox_searchNK,SIGNAL(currentIndexChanged(int)),this,SLOT(manageLEwl()));
    connect(ui->pushButton_selectN,SIGNAL( clicked() ), this, SLOT(selectNsol()));
    connect(ui->pushButton_reset,SIGNAL( clicked() ), this, SLOT(reset()));

    for (int N = 1; N <= 20; N++) {
        int pmIdx = 101 + 5 * (N - 1);
        QString name = QString("cBpm_%1_1").arg(pmIdx);
        QComboBox* w = findChild<QComboBox*>(name);
        if (w)
            connect(w, qOverload<int>(&QComboBox::currentIndexChanged),
                    this, [this, N]() { setKindOsc(N); });
    }

    connect(ui->comboBox_DeltaPsiScale,SIGNAL(currentIndexChanged(int)),this,SLOT(setRangeEli()));

    connect(ui->pushButton_setFont, SIGNAL( clicked() ), this, SLOT(setFontDia()));

    const QList<int> dSB_indices = {
        1,2,3,4,5,6,7,8,9,
        11,12,13,14,15,16,17,18,19,
        21,22,23,24,25,26,27,28,29,
        31,32,33,34,35,36,37,38,39,
        41,42,43,44,45,46,47,48,49,
        51,52,53,54,55,56,57,58,59,
        61,62,63,64,65,66,67,68,69,
        91,92,93,94,95,96,97,98,99
    };
    for (int idx : dSB_indices) {
        QString name = QString("dSB_PM_%1_1").arg(idx);
        QDoubleSpinBox* w = findChild<QDoubleSpinBox*>(name);
        if (w)
            connect(w, &QDoubleSpinBox::textChanged,
                    this, [this]() { SaveSetting(1); });
        else
            qWarning() << "widget non trovato:" << name;
    }

    // LEpm_*_1 → SaveSetting(2)
    const QList<int> LE_indices = {
        160,162,163,164,165,167,168,169,170,
        172,173,174,175,177,178,179,180,
        182,183,184,185,187,188,189,190,
        192,193,194,195,197,198,199,200
    };
    for (int idx : LE_indices) {
        QString name = QString("LEpm_%1_1").arg(idx);
        QLineEdit* w = findChild<QLineEdit*>(name);
        if (w)
            connect(w, &QLineEdit::textChanged,
                    this, [this]() { SaveSetting(2); });
        else
            qWarning() << "widget non trovato:" << name;
    }

    setGuiMute(true);
    qApp->installEventFilter(this);

    // parameter initialization
    QDir dir;  //current directory
    dir.cdUp();//cd ..
    dir.cdUp();//cd ..
    if(!IS_POSIX)
        dir.cdUp();
    pathroot=dir.absolutePath()+"/";
    fRefMir=pathroot+"qtSource/ksemawc/referenceMirrors.txt";
    fStdSpect=pathroot+"qtSource/ksemawc/standardSpectra.txt";
    fileStore=pathroot+"temp/defau.1.Spj";
    fileStore0=pathroot+"temp/defau.0.Spj";
    fileStoreSpjName=pathroot+"temp/SpjName.txt";
    fNKsim=pathroot+"expo/NKsim.dat";
    fMisSim=pathroot+"expo/MisSim.dat";
    filechi2=pathroot+"temp/ksemawc.log";

    qDebug()<<"pathroot= "<<pathroot;

    // stylesheet
    qDebug() << "styleSheet= " << QApplication::style()->metaObject()->className();
    qDebug() << "Current StyleSheet:" << qApp->styleSheet().toStdString().c_str();
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
    for(int i=1;i<=N_CNK_MAX;i++){
        cnk[i].forceMode=0.;//no forced
        cnk[i].nSrc=.0;//option "constant nk"
        cnk[i].nConst=1.;//n1 value
        cnk[i].kConst=.0;//k1 value
        cnk[i].emaSrc=-1.;//No EMA
        cnk[i].n2Const=1.;//n2 value
        cnk[i].k2Const=.0;//k2 value
        cnk[i].nForced=1.;//n forced
        cnk[i].kForced=.0;//k forced
    }
    // initialization of model parameters pm (uniform encoding, ip 1..996)
    for(int i=1;i<=996;i++){
        for(int j=1;j<=5;j++){
            pmAt(i)[j]=0.;
        }
    }
    pmTe[1][1]=-1.;
    pmTe[2][1]=-1.;
    pmTe[3][1]=-1.;
    pmTe[4][1]=-1.;
    pmOs[0][1]=1;    //option FT#1
    pmOs[1][1]=2.;   //quantistic homogeneous oscillator in UV
    pmOs[2][1]=0.;   //C1    "          "
    pmOs[3][1]=5.0;  //E1    "          "
    pmOs[4][1]=0.04; //D1    "          "
    pmOs[5][1]=0.;   //K1    "          "
    pmOs[6][1]=2.;   //quantistic homogeneous oscillator in IR
    pmOs[7][1]=0.;   //C2
    pmOs[8][1]=0.3;  //E2    "          "
    pmOs[9][1]=0.005;//D2    "          "
    pmOs[10][1]=0.;   //K2    "          "
    pmOs[11][1]=4.;   //FLAT term
    pmOs[12][1]=1.5;  //C3    "          "
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
    // store the fit-parameter grid for ensureFitRow (dynamic rows 18+)
    m_fitGrid = findChild<QGridLayout*>("gridLayout_10");
    ui->sB_PAR_34_5->setMaximum(N_FIT_MAX);
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

    // ── Create widget rows for layers 10..N_LAYER_MAX ─────────────────────────
    {
        QGridLayout* layerGrid = findChild<QGridLayout*>("gridLayout_8");
        QComboBox*   kindRef   = idToComboBox["comB_PAR_51_3"];
        QComboBox*   matRef    = idToComboBox["comB_PAR_51_1"];
        QWidget*     gw        = layerGrid ? layerGrid->parentWidget() : this;
        for(int i=10; i<=N_LAYER_MAX; i++){
            // Layer number label (col 0 = first)
            QLabel* numLbl = new QLabel(QString::number(i), gw);
            numLbl->setAlignment(Qt::AlignCenter);
            numLbl->setSizePolicy(QSizePolicy::Minimum, QSizePolicy::MinimumExpanding);
            if(layerGrid) layerGrid->addWidget(numLbl, i, 0);

            // Kind combo (col 1)
            QComboBox* cbK = new QComboBox(gw);
            cbK->setObjectName(wCB3(i)); cbK->setEnabled(false);
            cbK->setSizePolicy(QSizePolicy::Minimum, QSizePolicy::MinimumExpanding);
            for(int j=0; j<kindRef->count(); j++) cbK->addItem(kindRef->itemText(j));
            if(layerGrid) layerGrid->addWidget(cbK, i, 1);

            // Up/Down buttons (col 2, 3)
            QPushButton* btnUp = new QPushButton("up", gw);
            btnUp->setObjectName("pBmUp"+QString::number(i));
            btnUp->setEnabled(false);
            btnUp->setSizePolicy(QSizePolicy::Minimum, QSizePolicy::MinimumExpanding);
            btnUp->setMaximumWidth(60);
            if(layerGrid) layerGrid->addWidget(btnUp, i, 2);

            QPushButton* btnDw = new QPushButton("dw", gw);
            btnDw->setObjectName("pBmDw"+QString::number(i));
            btnDw->setEnabled(false);
            btnDw->setSizePolicy(QSizePolicy::Minimum, QSizePolicy::MinimumExpanding);
            btnDw->setMaximumWidth(60);
            if(layerGrid) layerGrid->addWidget(btnDw, i, 3);

            // Spinboxes
            auto mkDSB = [&](const QString& nm, int col, int dec, double mn, double mx) -> QDoubleSpinBox* {
                QDoubleSpinBox* sb = new QDoubleSpinBox(gw);
                sb->setObjectName(nm); sb->setEnabled(false);
                sb->setSizePolicy(QSizePolicy::Minimum, QSizePolicy::MinimumExpanding);
                sb->setDecimals(dec); sb->setMinimum(mn); sb->setMaximum(mx);
                if(layerGrid) layerGrid->addWidget(sb, i, col);
                return sb;
            };
            QDoubleSpinBox* sbD  = mkDSB(wD(i),  5,  2,     0.0,  1e9);
            QDoubleSpinBox* sbNu = mkDSB(wNu(i), 6,  4,     0.0,  1.0);
            QDoubleSpinBox* sbRg = mkDSB(wRg(i), 7,  1,     0.0,  1000.0);
            QDoubleSpinBox* sbDn = mkDSB(wDn(i), 8,  4,   -10.0,  10.0);
            QDoubleSpinBox* sbSd = mkDSB(wSd(i), 9,  3, -1000.0,  1000.0);
            QDoubleSpinBox* sbNc = mkDSB(wNc(i), 10, 4,   -10.0,  10.0);
            QDoubleSpinBox* sbDk = mkDSB(wDk(i), 11, 4,   -10.0,  10.0);
            QDoubleSpinBox* sbKc = mkDSB(wKc(i), 12, 4,   -10.0,  10.0);

            // Material combo (col 13)
            QComboBox* cbM = new QComboBox(gw);
            cbM->setObjectName(wCB1(i)); cbM->setEnabled(false);
            cbM->setSizePolicy(QSizePolicy::Minimum, QSizePolicy::MinimumExpanding);
            for(int j=0; j<matRef->count(); j++) cbM->addItem(matRef->itemText(j));
            if(layerGrid) layerGrid->addWidget(cbM, i, 13);

            // Register in maps
            idToComboBox[wCB3(i)]             = cbK;
            idToComboBox[wCB1(i)]             = cbM;
            idToPushButton["pBmUp"+QString::number(i)] = btnUp;
            idToPushButton["pBmDw"+QString::number(i)] = btnDw;
            idToDoubleSpinBox[wD(i)]  = sbD;
            idToDoubleSpinBox[wNu(i)] = sbNu;
            idToDoubleSpinBox[wRg(i)] = sbRg;
            idToDoubleSpinBox[wDn(i)] = sbDn;
            idToDoubleSpinBox[wSd(i)] = sbSd;
            idToDoubleSpinBox[wNc(i)] = sbNc;
            idToDoubleSpinBox[wDk(i)] = sbDk;
            idToDoubleSpinBox[wKc(i)] = sbKc;

            // Connect signals
            connect(cbK,   QOverload<int>::of(&QComboBox::currentIndexChanged),    this, &ksemawc::RefreshModel);
            connect(sbRg,  QOverload<double>::of(&QDoubleSpinBox::valueChanged),  this, &ksemawc::AdjRoughMax);
            connect(btnDw, &QPushButton::clicked, this, [this,i](){ mDwUp(i,  1); });
            connect(btnUp, &QPushButton::clicked, this, [this,i](){ mDwUp(i, -1); });
        }
        // Layer-number labels: header in row 0, then one per row 1-9 (col 0 = first)
        if(layerGrid){
            QLabel* hdr = new QLabel("#", gw);
            hdr->setAlignment(Qt::AlignCenter);
            hdr->setSizePolicy(QSizePolicy::Minimum, QSizePolicy::MinimumExpanding);
            layerGrid->addWidget(hdr, 0, 0);
            for(int i=1; i<=9; i++){
                QLabel* lbl = new QLabel(QString::number(i), gw);
                lbl->setAlignment(Qt::AlignCenter);
                lbl->setSizePolicy(QSizePolicy::Minimum, QSizePolicy::MinimumExpanding);
                layerGrid->addWidget(lbl, i, 0);
            }
        }
        // Spacer in row N_LAYER_MAX+1: assorbe lo spazio residuo mantenendo le righe alla loro altezza naturale
        if(layerGrid)
            layerGrid->setRowStretch(N_LAYER_MAX+1, 1);
        m_layerGrid   = layerGrid;
        m_maxLayerRow = N_LAYER_MAX;
    }
    ui->sB_PAR_51_2->setMaximum(N_LAYER_HARD_MAX);

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
    // Register nk rows 1-8 buttons in idToPushButton for uniform access in updateNkRows
    for(int i=1;i<=8;i++){
        idToPushButton["pBnk"+QString::number(i)]      = findChild<QPushButton*>("pBnk"+QString::number(i));
        idToPushButton["pBclearnk"+QString::number(i)] = findChild<QPushButton*>("pBclearnk"+QString::number(i));
    }
    // Create nk rows 9-20 programmatically
    {
        QGridLayout* nkGrid = findChild<QGridLayout*>("gridLayout_4");
        QWidget*     gw     = nkGrid ? nkGrid->parentWidget() : this;
        for(int i=9; i<=N_CNK_USER_MAX; i++){
            QString si = QString::number(i);
            QPushButton* btnLoad = new QPushButton("Load nk_"+si, gw);
            btnLoad->setObjectName("pBnk"+si);
            btnLoad->setSizePolicy(QSizePolicy::Minimum, QSizePolicy::MinimumExpanding);
            if(nkGrid) nkGrid->addWidget(btnLoad, i-1, 0);

            QLineEdit* lePath = new QLineEdit("mate/aa999.9", gw);
            lePath->setObjectName("lineEdit"+si);
            lePath->setSizePolicy(QSizePolicy::Minimum, QSizePolicy::MinimumExpanding);
            lePath->setMinimumWidth(200);
            if(nkGrid) nkGrid->addWidget(lePath, i-1, 1);

            QPushButton* btnClr = new QPushButton("Clear", gw);
            btnClr->setObjectName("pBclearnk"+si);
            btnClr->setSizePolicy(QSizePolicy::Minimum, QSizePolicy::MinimumExpanding);
            if(nkGrid) nkGrid->addWidget(btnClr, i-1, 2);

            QLineEdit* leInfo = new QLineEdit(gw);
            leInfo->setObjectName("lineEdit_infoNK_"+si);
            leInfo->setSizePolicy(QSizePolicy::Expanding, QSizePolicy::MinimumExpanding);
            leInfo->setReadOnly(true);
            if(nkGrid) nkGrid->addWidget(leInfo, i-1, 3);

            QLineEdit* leMin = new QLineEdit(gw);
            leMin->setObjectName("WLminNK"+si);
            leMin->setSizePolicy(QSizePolicy::Minimum, QSizePolicy::MinimumExpanding);
            leMin->setReadOnly(true);
            if(nkGrid) nkGrid->addWidget(leMin, i-1, 4);

            QLineEdit* leMax = new QLineEdit(gw);
            leMax->setObjectName("WLmaxNK"+si);
            leMax->setSizePolicy(QSizePolicy::Minimum, QSizePolicy::MinimumExpanding);
            leMax->setReadOnly(true);
            if(nkGrid) nkGrid->addWidget(leMax, i-1, 5);

            idToPushButton["pBnk"+si]          = btnLoad;
            idToPushButton["pBclearnk"+si]      = btnClr;
            idToLineEdit["lineEdit"+si]          = lePath;
            idToLineEdit["lineEdit_infoNK_"+si]  = leInfo;
            idToLineEdit["WLminNK"+si]           = leMin;
            idToLineEdit["WLmaxNK"+si]           = leMax;

            connect(btnLoad, &QPushButton::clicked, this, [this,i](){ Setnk(i); });
            connect(btnClr,  &QPushButton::clicked, this, [this,i](){ Clrnk(i); });
        }
        // Spacer absorbs leftover space keeping row heights uniform
        if(nkGrid) nkGrid->setRowStretch(N_CNK_USER_MAX, 1);
    }
    // The .ui statically appends a trailing "Material#1" item to cB_cnk2a..15a (it lets a
    // material reuse material #1's n,k). Strip it here so the nk-9..99 slots land at
    // contiguous indices 8..106 in every combo; it is re-appended uniformly at the trailing
    // index NSRC_MATREF1 for all material #2+ source combos below.
    for(int i=1;i<=15;i++){
        QComboBox* ca=idToComboBox["cB_cnk"+QString::number(i)+"a"];
        int mref=ca->findText("Material#1");
        if(mref>=0) ca->removeItem(mref);
    }
    // Add nk-9..nk-N_CNK_USER_MAX to all static cB_cnk combos 1-15 (currently have items 0..15)
    for(int i=1;i<=15;i++){
        for(int j=9;j<=N_CNK_USER_MAX;j++){
            idToComboBox["cB_cnk"+QString::number(i)+"a"]->addItem("nk-"+QString::number(j));
            idToComboBox["cB_cnk"+QString::number(i)+"b"]->addItem("nk-"+QString::number(j));
        }
    }
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
        idToComboBox["cBpm_"+QString::number(100+1+(k-1)*5)+"_1"]->addItem("Volume-Scattering");     //#19
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
            qDebug()<<"-> RefMir # -> "<<iB103<<" "<<line;
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
            qDebug()<< "-> StdSpectrum #"<<iStdS<<line;
            line = stream1.readLine();
            qDebug()<<" -> "<<line;
            iStdS++;
        }while(!stream1.atEnd());
        file1.close();
    }

    //label
    ParFitLab[0]="_dNonUni";
    ParFitLab[1]="_d";
    ParFitLab[2]="_grad-n";
    ParFitLab[3]="_curv-n";
    ParFitLab[4]="_grad-k";
    ParFitLab[5]="_curv-k";
    ParFitLab[6]="_sigma_rough";
    ParFitLab[7]="_slopeGrad-n";
    ParFitLab[8]="_C";
    ParFitLab[9]="_E";
    ParFitLab[10]="_D";
    ParFitLab[11]="_W";
    ParFitLab[12]="fEMA_";
    ParFitLab[13]="Theta_";

    //Graph initialization
    for(int i = 0; i <= N_GRAPHS; ++i){
        m_graphs[i]  = nullptr;
        m_pickers[i] = nullptr;
    }

    //set the latest Size & font if settings exist
    QSettings settings("ksemawc","ksemawc");
    QString usnam=settings.fileName();
    qDebug()<<"setting-file exists! "<<usnam;
    if (settings.status() == QSettings::NoError && settings.status() != QSettings::FormatError && settings.status() !=QSettings::AccessError){
        restoreGeometry(settings.value("geometry").toByteArray());
        QString fontFamily=settings.value("fontFamily",QString()).toString();
        int fontSize = settings.value("pointSize", 0).toInt();
        //apply the stored font only if the user actually chose one via setFontDia,
        //otherwise the widgets keep the fonts coming from the .ui
        if(!fontFamily.isEmpty() && fontSize > 0){
            QFont font;
            font.setFamily(fontFamily);
            font.setPointSize(fontSize);
            font.setBold(settings.value("bold", false).toBool());
            font.setItalic(settings.value("italic", false).toBool());
            const QWidgetList allWidgets = QApplication::allWidgets();
            for (QWidget *widget : allWidgets){
                widget->setFont(font);
                widget->update();
            }
            QCoreApplication::processEvents();
        }
    }

    // Hide nk rows 2-N_CNK_USER_MAX — updateNkRows will reveal them
    for(int i=2; i<=N_CNK_USER_MAX; i++){
        idToPushButton["pBnk"+QString::number(i)]->setVisible(false);
        idToPushButton["pBclearnk"+QString::number(i)]->setVisible(false);
        idToLineEdit["lineEdit"+QString::number(i)]->setVisible(false);
        idToLineEdit["lineEdit_infoNK_"+QString::number(i)]->setVisible(false);
        idToLineEdit["WLminNK"+QString::number(i)]->setVisible(false);
        idToLineEdit["WLmaxNK"+QString::number(i)]->setVisible(false);
    }
    // Hide all layer rows except row 1 — SetModel will reveal them as needed
    for(int i=2; i<=N_LAYER_MAX; i++){
        idToComboBox[wCB3(i)]->setVisible(false);
        idToComboBox[wCB1(i)]->setVisible(false);
        idToDoubleSpinBox[wD(i)]->setVisible(false);
        idToDoubleSpinBox[wDn(i)]->setVisible(false);
        idToDoubleSpinBox[wNc(i)]->setVisible(false);
        idToDoubleSpinBox[wDk(i)]->setVisible(false);
        idToDoubleSpinBox[wKc(i)]->setVisible(false);
        idToDoubleSpinBox[wRg(i)]->setVisible(false);
        idToDoubleSpinBox[wSd(i)]->setVisible(false);
        idToDoubleSpinBox[wNu(i)]->setVisible(false);
        idToPushButton["pBmUp"+QString::number(i)]->setVisible(false);
        idToPushButton["pBmDw"+QString::number(i)]->setVisible(false);
        if(m_layerGrid)
            if(auto* it = m_layerGrid->itemAtPosition(i, 0))
                if(auto* lbl = qobject_cast<QLabel*>(it->widget()))
                    lbl->setVisible(false);
    }

    // Add Material #10-N_CNK_USER_MAX to all wCB1 layer material combos (items #1-9 in .ui)
    for(int i=1; i<=N_LAYER_MAX; i++){
        QComboBox* cb=idToComboBox[wCB1(i)];
        for(int m=10; m<=N_CNK_USER_MAX; m++)
            cb->addItem("Material #"+QString::number(m));
    }

    // Create CNK rows 16-20 (user materials) inside scrollArea_mater,
    // rows 21-24 (special media) in a separate panel always visible below it
    {
        QGridLayout* cnkGrid = findChild<QGridLayout*>("gridLayout_11");
        QWidget* gwScroll = cnkGrid ? cnkGrid->parentWidget() : this;
        QComboBox* aRef = idToComboBox["cB_cnk1a"];

        // Build the special-media panel and insert it after scrollArea_mater
        QScrollArea* sa_mater = findChild<QScrollArea*>("scrollArea_mater");
        QGridLayout* specGrid = nullptr;
        QWidget* specWidget = nullptr;
        if(sa_mater && sa_mater->parentWidget()){
            specWidget = new QWidget(sa_mater->parentWidget());
            specGrid = new QGridLayout(specWidget);
            specGrid->setContentsMargins(0,0,0,0);
            if(cnkGrid){
                specGrid->setHorizontalSpacing(cnkGrid->horizontalSpacing());
                specGrid->setVerticalSpacing(cnkGrid->verticalSpacing());
            }
            auto* parentBox = qobject_cast<QBoxLayout*>(sa_mater->parentWidget()->layout());
            if(parentBox){
                // btnSpec is defined statically in the .ui as the last item of the
                // Sample Model tab layout (right after scrollArea_mater); reuse it and
                // insert the special-media panel right after it.
                QPushButton* btnSpec = findChild<QPushButton*>("pushButton_setIOmedia");
                int btnIdx = -1;
                if(btnSpec){
                    for(int k=0; k<parentBox->count(); k++){
                        if(parentBox->itemAt(k)->widget() == btnSpec){ btnIdx=k; break; }
                    }
                }
                if(btnIdx >= 0)
                    parentBox->insertWidget(btnIdx+1, specWidget);
                else
                    parentBox->addWidget(specWidget);
                specWidget->setVisible(false);
                if(btnSpec)
                    connect(btnSpec, &QPushButton::clicked, this, [specWidget](){
                        specWidget->setVisible(!specWidget->isVisible());
                    });
            }
        }

        // Shared lambda: create all widgets for one CNK material row
        // (aRef already has nk-9..N_CNK_USER_MAX items at this point)
        auto makeCnkRow = [&](int i, QGridLayout* grid, QWidget* gw, int gridRow, const QString& label){
            QString si = QString::number(i);
            QLabel* lbl = new QLabel(label, gw);
            lbl->setMinimumWidth(72); lbl->setMaximumWidth(100);
            if(grid) grid->addWidget(lbl, gridRow, 0);
            QComboBox* ca = new QComboBox(gw);
            ca->setObjectName("cB_cnk"+si+"a");
            ca->setSizePolicy(QSizePolicy::MinimumExpanding, QSizePolicy::MinimumExpanding);
            ca->setMinimumWidth(120);
            if(aRef) for(int j=0;j<aRef->count();j++) ca->addItem(aRef->itemText(j));
            if(grid) grid->addWidget(ca, gridRow, 1);
            idToComboBox["cB_cnk"+si+"a"] = ca;
            connect(ca, qOverload<int>(&QComboBox::currentIndexChanged), this, [this,i](){ setMat(i); });
            auto makeLe = [&](const QString& name, const QString& defVal, int col, bool enabled=true) -> QLineEdit* {
                QLineEdit* le = new QLineEdit(defVal, gw);
                le->setObjectName(name);
                le->setSizePolicy(QSizePolicy::Minimum, QSizePolicy::MinimumExpanding);
                le->setMinimumWidth(60); le->setMaximumWidth(100);
                le->setEnabled(enabled);
                if(grid) grid->addWidget(le, gridRow, col);
                return le;
            };
            idToLineEdit["LEcnk"+si+"_2"] = makeLe("LEcnk"+si+"_2","1",2);
            idToLineEdit["LEcnk"+si+"_3"] = makeLe("LEcnk"+si+"_3","0",3);
            QCheckBox* ema = new QCheckBox("EMA", gw);
            ema->setObjectName("cB_EMA_"+si);
            ema->setChecked(false);
            ema->setSizePolicy(QSizePolicy::Minimum, QSizePolicy::MinimumExpanding);
            ema->setMaximumWidth(100);
            if(grid) grid->addWidget(ema, gridRow, 4);
            idToCheckBox["cB_EMA_"+si] = ema;
            connect(ema, &QCheckBox::stateChanged, this, [this,i](){ setEMA(i); });
            QDoubleSpinBox* dsb = new QDoubleSpinBox(gw);
            dsb->setObjectName(emaKey(i));
            dsb->setEnabled(false);
            dsb->setMaximum(1.0); dsb->setSingleStep(0.01); dsb->setDecimals(3);
            dsb->setSizePolicy(QSizePolicy::Minimum, QSizePolicy::MinimumExpanding);
            if(grid) grid->addWidget(dsb, gridRow, 5);
            idToDoubleSpinBox[emaKey(i)] = dsb;
            QComboBox* cb = new QComboBox(gw);
            cb->setObjectName("cB_cnk"+si+"b");
            cb->setEnabled(false);
            cb->setSizePolicy(QSizePolicy::Minimum, QSizePolicy::MinimumExpanding);
            if(aRef) for(int j=0;j<aRef->count();j++) cb->addItem(aRef->itemText(j));
            cb->setCurrentIndex(-1);
            if(grid) grid->addWidget(cb, gridRow, 6);
            idToComboBox["cB_cnk"+si+"b"] = cb;
            connect(cb, qOverload<int>(&QComboBox::currentIndexChanged), this, [this,i](){ setEMA(i); });
            idToLineEdit["LEcnk"+si+"_5"] = makeLe("LEcnk"+si+"_5","1",7,false);
            idToLineEdit["LEcnk"+si+"_6"] = makeLe("LEcnk"+si+"_6","0",8,false);
        };
        // User materials 16..N_CNK_USER_MAX in cnkGrid (row = i-1)
        QWidget* gwSc = gwScroll ? gwScroll : this;
        for(int i=16; i<=N_CNK_USER_MAX; i++)
            makeCnkRow(i, cnkGrid, gwSc, i-1, "Material #"+QString::number(i));
        if(cnkGrid) cnkGrid->setRowStretch(N_CNK_USER_MAX, 1);
        m_cnkGrid = cnkGrid;
        // Special media N_CNK_SPEC_START..N_CNK_MAX in specGrid (rows 0-3)
        if(specGrid && specWidget){
            static const QStringList specLabels = {"output SF","input SF","input PDS","input ELI"};
            for(int idx=0; idx<4; idx++)
                makeCnkRow(N_CNK_SPEC_START+idx, specGrid, specWidget, idx, specLabels[idx]);
        }
    }
    // Append the "Material#1" source option to every material #2+ source combo so a
    // material can reuse material #1's (the unknown's) n,k. It lands at the trailing index
    // NSRC_MATREF1 (=107), well clear of the nk-file slots, so it never collides with them.
    for(int i=2;i<=N_CNK_MAX;i++){
        QComboBox* ca=idToComboBox.value("cB_cnk"+QString::number(i)+"a",nullptr);
        if(ca) ca->addItem("Material#1");
    }
    connect(ui->sB_nmat, QOverload<int>::of(&QSpinBox::valueChanged), this, &ksemawc::updateMatRows);
    // NK items for combos 16..N_CNK_MAX already included via aRef copy in makeCnkRow

    // scrollArea_layers and scrollArea_mater both carry Expanding size policy with
    // verstretch=1 in the .ui, so they share the tab's vertical space equally and grow
    // with the window. (Previously scrollArea_layers was capped via setMaximumHeight to
    // show only 8 layer rows, which prevented it from adapting to the available space.)

    // CNK rows 1-15 come from the .ui, where the EMA-source combo holds items and therefore
    // starts on index 0, while cnk[i].emaSrc is initialised to -1 (= no EMA). SaveSetting
    // reads emaSrc back from that combo, so an untouched row would switch EMA on by itself.
    // Put the static rows in the same "no EMA" state the dynamic ones get in makeCnkRow.
    for(int i=1;i<=15;i++){
        idToComboBox["cB_cnk"+QString::number(i)+"b"] -> setCurrentIndex(-1);
        idToComboBox["cB_cnk"+QString::number(i)+"b"] -> setEnabled(false);
        idToCheckBox["cB_EMA_"+QString::number(i)]    -> setChecked(false);
        idToDoubleSpinBox[emaKey(i)]                  -> setEnabled(false);
    }

    setGuiMute(false);

    if(qApp->platformName()=="offscreen"){
        QTimer::singleShot(200, this, [this](){
            caller="initialization";
            ReadSetting(pathroot+"ZnSe/ZnSe.2.Spj");
            SPADA();
            Simula();
            qApp->quit();
        });
        return;
    }

    //make your choice
    QMessageBox msgBox;
    msgBox.setWindowTitle("Starting ksemawc project");
    msgBox.setText("Please make a choice:");
    QPushButton* continueButton = msgBox.addButton("Continue the previous session", QMessageBox::AcceptRole);
    QPushButton* newButton = msgBox.addButton("Start a new project", QMessageBox::RejectRole);
    QPushButton* openButton = msgBox.addButton("Open an existing project", QMessageBox::DestructiveRole);
    msgBox.setDefaultButton(newButton);
    int ret=msgBox.exec();
    QString SpjName;
    QFile file(fileStoreSpjName);
    switch (ret) {

    case QMessageBox::AcceptRole:
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
        SetModel(ui->sB_PAR_51_2->value());
        updateMatRows(ui->sB_nmat->value());
        updateNkRows(ui->sB_nnk->value());
        break;

    case QMessageBox::DestructiveRole:
        LoadProject();
        break;
    }

    SaveSetting(-1);//save
}


void ksemawc::about(){
    QMessageBox::about(this, "kSEMAWc",
                       "<b>kSEMAWc</b> <br>"
                       "Spectro-Ellipsometric Measurement Analysis Workbench<br>"
                       "Version 5.0 — 15 July 2026<br><br>"
                       "Main author: Marco Montecchi, ENEA (Italy)<br>"
                       "email: marco.montecchi@enea.it<br><br>"
                       "Porting to Windows and advanced oscillators by<br>"
                       "Alberto Mittiga, ENEA (Italy)<br>"
                       "email: alberto.mittiga@enea.it"
                       );
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
    // Use QFile::copy (cross-platform, no shell): the shell "copy"/"cp" route broke on
    // paths containing '+' (cmd.exe COPY treats '+' as the concat operator) or spaces.
    int ierr = copyFileOverwrite(fileStore, fileStore0) ? 0 : 1;
    if(ierr != 0){
        msgErrLoad("SaveProject",fnproject);
        qWarning()<<"Error copying project!!!";
    }
    else
        qInfo()<<"Project saved as "+fnproject;
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
        //save the choice
        QSettings settings("ksemawc","ksemawc");
        settings.setValue("fontFamily", font.family());
        settings.setValue("pointSize", font.pointSize());
        settings.setValue("bold", font.bold());
        settings.setValue("italic", font.italic());
    }
    // if the user canceled the dialog nothing is applied nor stored
}


void ksemawc::ensureLayerRow(int i){
    if(!m_layerGrid || i <= m_maxLayerRow) return;
    QComboBox* kindRef = idToComboBox[wCB3(1)];
    QComboBox* matRef  = idToComboBox[wCB1(1)];
    QWidget*   gw      = m_layerGrid->parentWidget();
    for(int r = m_maxLayerRow+1; r <= i; r++){
        QLabel* numLbl = new QLabel(QString::number(r), gw);
        numLbl->setAlignment(Qt::AlignCenter);
        numLbl->setSizePolicy(QSizePolicy::Minimum, QSizePolicy::MinimumExpanding);
        m_layerGrid->addWidget(numLbl, r, 0);

        QComboBox* cbK = new QComboBox(gw);
        cbK->setObjectName(wCB3(r)); cbK->setEnabled(false);
        cbK->setSizePolicy(QSizePolicy::Minimum, QSizePolicy::MinimumExpanding);
        for(int j=0; j<kindRef->count(); j++) cbK->addItem(kindRef->itemText(j));
        m_layerGrid->addWidget(cbK, r, 1);

        QPushButton* btnUp = new QPushButton("up", gw);
        btnUp->setObjectName("pBmUp"+QString::number(r));
        btnUp->setEnabled(false);
        btnUp->setSizePolicy(QSizePolicy::Minimum, QSizePolicy::MinimumExpanding);
        btnUp->setMaximumWidth(60);
        m_layerGrid->addWidget(btnUp, r, 2);

        QPushButton* btnDw = new QPushButton("dw", gw);
        btnDw->setObjectName("pBmDw"+QString::number(r));
        btnDw->setEnabled(false);
        btnDw->setSizePolicy(QSizePolicy::Minimum, QSizePolicy::MinimumExpanding);
        btnDw->setMaximumWidth(60);
        m_layerGrid->addWidget(btnDw, r, 3);

        auto mkDSB = [&](const QString& nm, int col, int dec, double mn, double mx) -> QDoubleSpinBox* {
            QDoubleSpinBox* sb = new QDoubleSpinBox(gw);
            sb->setObjectName(nm); sb->setEnabled(false);
            sb->setSizePolicy(QSizePolicy::Minimum, QSizePolicy::MinimumExpanding);
            sb->setDecimals(dec); sb->setMinimum(mn); sb->setMaximum(mx);
            m_layerGrid->addWidget(sb, r, col);
            return sb;
        };
        QDoubleSpinBox* sbD  = mkDSB(wD(r),  5,  2,     0.0,  1e9);
        QDoubleSpinBox* sbNu = mkDSB(wNu(r), 6,  4,     0.0,  1.0);
        QDoubleSpinBox* sbRg = mkDSB(wRg(r), 7,  1,     0.0,  1000.0);
        QDoubleSpinBox* sbDn = mkDSB(wDn(r), 8,  4,   -10.0,  10.0);
        QDoubleSpinBox* sbSd = mkDSB(wSd(r), 9,  3, -1000.0,  1000.0);
        QDoubleSpinBox* sbNc = mkDSB(wNc(r), 10, 4,   -10.0,  10.0);
        QDoubleSpinBox* sbDk = mkDSB(wDk(r), 11, 4,   -10.0,  10.0);
        QDoubleSpinBox* sbKc = mkDSB(wKc(r), 12, 4,   -10.0,  10.0);

        QComboBox* cbM = new QComboBox(gw);
        cbM->setObjectName(wCB1(r)); cbM->setEnabled(false);
        cbM->setSizePolicy(QSizePolicy::Minimum, QSizePolicy::MinimumExpanding);
        for(int j=0; j<matRef->count(); j++) cbM->addItem(matRef->itemText(j));
        { int nm = ui->sB_nmat->value();
          for(int j = 0; j < cbM->count(); j++) if(auto* lv=qobject_cast<QListView*>(cbM->view())) lv->setRowHidden(j, j >= nm); }
        m_layerGrid->addWidget(cbM, r, 13);

        idToComboBox[wCB3(r)]                      = cbK;
        idToComboBox[wCB1(r)]                      = cbM;
        idToPushButton["pBmUp"+QString::number(r)] = btnUp;
        idToPushButton["pBmDw"+QString::number(r)] = btnDw;
        idToDoubleSpinBox[wD(r)]  = sbD;
        idToDoubleSpinBox[wNu(r)] = sbNu;
        idToDoubleSpinBox[wRg(r)] = sbRg;
        idToDoubleSpinBox[wDn(r)] = sbDn;
        idToDoubleSpinBox[wSd(r)] = sbSd;
        idToDoubleSpinBox[wNc(r)] = sbNc;
        idToDoubleSpinBox[wDk(r)] = sbDk;
        idToDoubleSpinBox[wKc(r)] = sbKc;

        connect(cbK,   QOverload<int>::of(&QComboBox::currentIndexChanged),   this, &ksemawc::RefreshModel);
        connect(sbRg,  QOverload<double>::of(&QDoubleSpinBox::valueChanged),  this, &ksemawc::AdjRoughMax);
        connect(btnDw, &QPushButton::clicked, this, [this,r](){ mDwUp(r,  1); });
        connect(btnUp, &QPushButton::clicked, this, [this,r](){ mDwUp(r, -1); });

        m_layerGrid->setRowStretch(r,   0);
        m_layerGrid->setRowStretch(r+1, 1);
    }
    m_maxLayerRow = i;
}


void ksemawc::ensureFitRow(int j){
    if(!m_fitGrid || j <= m_maxFitRow) return;
    QWidget* gw = m_fitGrid->parentWidget();
    for(int r = m_maxFitRow+1; r <= j; r++){
        QString sr = QString::number(r);

        // col 0: enable checkbox
        QCheckBox* chk = new QCheckBox(gw);
        chk->setEnabled(false);
        chk->setSizePolicy(QSizePolicy::Minimum, QSizePolicy::MinimumExpanding);
        idToCheckBox["chBeParFit_"+sr] = chk;
        m_fitGrid->addWidget(chk, r, 0);
        connect(chk, &QCheckBox::stateChanged, this, &ksemawc::PanFitEnable);

        // col 1: parameter selector combobox
        QComboBox* cb = new QComboBox(gw);
        cb->addItem("none");
        cb->setEnabled(true);
        cb->setSizePolicy(QSizePolicy::MinimumExpanding, QSizePolicy::MinimumExpanding);
        cb->setMinimumWidth(140);
        idToComboBox["cBParFit_"+sr] = cb;
        m_fitGrid->addWidget(cb, r, 1);
        connect(cb, qOverload<int>(&QComboBox::currentIndexChanged), this, &ksemawc::PanFitChoice);

        // col 2: value
        QLineEdit* leV = new QLineEdit(gw);
        leV->setEnabled(false);
        leV->setSizePolicy(QSizePolicy::Minimum, QSizePolicy::MinimumExpanding);
        idToLineEdit["DPparFitV_"+sr] = leV;
        m_fitGrid->addWidget(leV, r, 2);

        // col 3: error
        QLineEdit* leE = new QLineEdit(gw);
        leE->setEnabled(false);
        leE->setReadOnly(true);
        leE->setSizePolicy(QSizePolicy::Minimum, QSizePolicy::MinimumExpanding);
        idToLineEdit["DPparFitErr_"+sr] = leE;
        m_fitGrid->addWidget(leE, r, 3);

        // col 4: global correlation
        QLineEdit* leG = new QLineEdit(gw);
        leG->setEnabled(false);
        leG->setReadOnly(true);
        leG->setSizePolicy(QSizePolicy::Minimum, QSizePolicy::MinimumExpanding);
        idToLineEdit["DPparFitGC_"+sr] = leG;
        m_fitGrid->addWidget(leG, r, 4);
    }
    m_maxFitRow = j;
}


bool ksemawc::eventFilter(QObject* obj, QEvent* event){
    if(event->type() == QEvent::Show)
        if(auto* mb = qobject_cast<QMessageBox*>(obj)){
            mb->setStyleSheet(
                "QMessageBox        { background-color: #ffffc8; }"
                "QMessageBox QLabel { background-color: #ffffc8; }"
            );
            mb->raise();
            //Wayland forbids an application to activate its own windows: requestActivate()
            //would just print a warning. Modal dialogs get the focus from the compositor anyway.
            if(!QApplication::platformName().startsWith("wayland"))
                mb->activateWindow();
        }
    return QWidget::eventFilter(obj, event);
}


void ksemawc::setGuiMute(bool mute) {
    if(mute)
        qDebug()<<"=> setGuiMute(true)";
    else
        qDebug()<<"=> setGuiMute(false)";
    m_guiIsMuted = mute;
    // Trova tutti i widget figli (pulsanti, spinbox, checkbox, ecc.)
    // e silenzia i loro segnali in un colpo solo
    const auto widgets = this->findChildren<QWidget*>();
    for (auto w : widgets) {
        w->blockSignals(mute);
    }
}


void ksemawc::LoadProject(){
    qDebug()<<"-> LoadProject";
    fnproject = QFileDialog::getOpenFileName(
                this,
                "Choose a SEMAW project", //window title
                pathroot,                 //initial directory
                "Semaw Project (*.Spj)"); //file extension
    if(fnproject.isEmpty())
        return;

    ui->tabWidget -> setCurrentIndex(0);
    setGuiMute(true);
    //thetaEli reset
    pmTe[1][1]=0.;
    pmTe[2][1]=0.;
    pmTe[3][1]=0.;
    pmTe[4][1]=0.;
    // reset of pm - model parameters (uniform encoding, ip 1..996)
    for(int i=1;i<=996;i++){
        for(int j=1;j<=5;j++){
            pmAt(i)[j]=0.;
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
    for(int j=1;j<=m_maxFitRow;j++){
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
    setGuiMute(false);
    ReadSetting(fnproject);
    //SPADA();//load file-nk & measurements
    SaveSetting(-1);
    ifirstWarning=0;
}


void ksemawc::ReadSetting(QString filename){
    qDebug()<<"-> ReadSetting from "<<filename;
    setGuiMute(true);
    for(int i=1;i<=N_CNK_USER_MAX;i++)
        Clrnk(i);
    ClrFnk();
    Clrfn();
    Qt::CheckState state;
    QFileInfo info(filename);
    QDateTime dtTm=info.lastModified();
    QString as = dtTm.toString("yyyy-MM-dd HH:mm:ss");
    qInfo()<<"date: "<<as;
    QDateTime tLim(QDate(2022, 1, 1), QTime(0, 00, 0));
    QDateTime tLim2(QDate(2023, 4, 11), QTime(0, 00, 0));
    //QDateTime tLim3(QDate(2025, 10, 25), QTime(0, 00, 0));
    int old0new1=0;
    if(dtTm<tLim)
        qInfo()<<"The project is older than January 2022";
    else if(dtTm>=tLim && dtTm<tLim2){
        qInfo()<<"The project is newer than January 2022 and older than 11 April 2023";
        old0new1=1;
    }
    else{
        qInfo()<<"The project version is up-to-date";
        old0new1=2;
    }
    double a2nm=0.1;
    QFile file(filename);
    if (!file.open (QIODevice::ReadOnly | QIODevice::Text)){
        setGuiMute(false);
        return;
    }
    QTextStream stream ( &file );
    QString line,line2,lab;
    int irx=30;
    // "ufip" in the header marks a file whose fit-parameter pointers (ppm / par[35+i][5])
    // are already stored in the uniform ip encoding (v5.2+); such values must NOT be run
    // through migrateIp(), whose legacy ranges (e.g. 289..920 = old layers 21-99) would
    // mangle them. Files without the marker keep the legacy best-effort migration.
    bool uniformIpFile=false;
    line = stream.readLine();
    if(line.contains("iVspj=1")){
        if(line.contains("nm"))
            a2nm=1;
        if(line.contains("ufip"))
            uniformIpFile=true;
        line = stream.readLine();
    }
    else{
        irx=25;
        qInfo()<<"-> Spj old version";
    }
    ui->lineEdit_infoP -> setText(line.simplified());
    line2=fnproject.section(pathroot, 1, 1);
    ui->lineEdit_P -> setText(line2);

    int nnkLoaded=1;
    for(int i=1;i<=8;i++){
        line = stream.readLine();
        NANK[i]=line.simplified();
        idToLineEdit["lineEdit"+QString::number(i)] -> setText(NANK[i]);
        fnk[i]=pathroot+line.simplified()+".nk";
        ifn=0;
        if(!NANK[i].contains("mate/aa999.9")){
            ifn=1;
            nnkLoaded=i;
            Setnk(i);
        }
    }
    ui->sB_nnk->setValue(nnkLoaded);
    updateNkRows(nnkLoaded);

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
    line = stream.readLine();
    NANK[12]=line.simplified();
    fnE1=pathroot+line.simplified()+".el";
    line2=fnE1.section('.', 1, 1);//N. mis el
    par[27][4]=line2.toDouble();
    //qDebug("from fnE1 -> par[27][4]=%f",par[27][4]);
    line = stream.readLine();
    NANK[13]=line.simplified();
    fnE2=pathroot+line.simplified()+".el";
    line2=fnE2.section('.', 1, 1);//N. mis el
    par[28][4]=line2.toDouble();
    line = stream.readLine();
    NANK[14]=line.simplified();
    fnE3=pathroot+line.simplified()+".el";
    line2=fnE3.section('.', 1, 1);//N. mis el
    par[29][4]=line2.toDouble();
    line = stream.readLine();
    NANK[15]=line.simplified();
    fnE4=pathroot+line.simplified()+".el";
    line2=fnE4.section('.', 1, 1);//N. mis el
    par[30][4]=line2.toDouble();
    // The ellipsometric measurement number derived from the file name (above) is the
    // reliable source for BOTH old and new projects: the save routine always writes
    // NANK[12..15] = sample+"."+misNumber. Preserve it here because the PAR-block read
    // below overwrites par[27..30][4]: in old projects (pre-oscillator-fit format) that
    // slot held -1, which clobbered the right value and pushed the ELI measurements into
    // the wrong rows (ELI-3/ELI-4) leaving ELI-1/ELI-2 unselectable.
    int nMisEfn[5];
    for(int k=1;k<=4;k++)
        nMisEfn[k]=nint(par[26+k][4]);

    // setGuiMute(false);
    if(fnSample.isEmpty() || fnSample.contains("mate/aa999")){
        Clrfn();
        ifn=0;
    }
    else{
        ifn=1;
        listMeas();
        //qDebug("After listMeas par[27][4]=%f",par[27][4]);
    }
    // setGuiMute(true);

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
    if(nint(rxy[25][3])==1)
        ui->checkB_RXY_25_3 -> setCheckState ( Qt::Unchecked );
    else
        ui->checkB_RXY_25_3 -> setCheckState ( Qt::Checked );
    if(nint(rxy[25][4])==1)
        ui->checkB_RXY_25_4 -> setCheckState ( Qt::Unchecked );
    else
        ui->checkB_RXY_25_4 -> setCheckState ( Qt::Checked );

    //CNK
    for(int i=1;i<=15;i++){
        line = stream.readLine();
        line=line.simplified();
        QStringList List0;
        List0 =line.split(" ");
        int nV=List0.count();
        fflush(stdout);
        for(int j=1;j<=nV;j++){
            pezzo=List0.at(j-1).toLocal8Bit().constData();
            cnk[i].col(j)=pezzo.toDouble();
        }
        if(nV==3){//old project
            if(cnk[i].nSrc<=17){
                cnk[i].emaSrc=-1.;
            }
            else{
                int i1=nint(cnk[i].nSrc/1000.);
                int i2=nint((cnk[i].nSrc-i1*1000.)/10.);
                cnk[i].nSrc=i1;
                cnk[i].emaSrc=i2;
            }
        }
    }

    //PAR
    for(int i=1;i<=60;i++){
        for(int j=1;j<=5;j++){
            stream>> par[i][j];
        }
    }
    //qDebug("After read par: par[27][4]=%f",par[27][4]);
    // Restore the file-name-derived ELI measurement numbers clobbered by the PAR read.
    for(int k=1;k<=4;k++)
        par[26+k][4]=nMisEfn[k];
    par[4][1]=par[4][1]*a2nm;//WLmin in nm
    par[4][2]=par[4][2]*a2nm;//WLmax in nm
    double WLmin=par[4][1];
    double WLmax=par[4][2];
    // Select the combo item matching the measurement number (findText, not index n-1:
    // the combo only lists existing .el files, so number != position in general).
    // For unused rows (fnE = mate/aa999) deselect explicitly: listMeas() leaves every
    // ELI combo defaulted to index 0, so without this the unused rows would still feed
    // a measurement to pwEj() and load it into ELI-3/ELI-4.
    const QString cBmisEli[5]={"","cBmis7","cBmis9","cBmis11","cBmis13"};
    const QString fnEli[5]={"",fnE1,fnE2,fnE3,fnE4};
    for(int k=1;k<=4;k++){
        if(fnEli[k].contains("mate/aa999")){
            idToComboBox[cBmisEli[k]]->setCurrentIndex(-1);
            continue;
        }
        int idx=idToComboBox[cBmisEli[k]]->findText(QString::number(nMisEfn[k]),Qt::MatchExactly);
        idToComboBox[cBmisEli[k]]->setCurrentIndex(idx);
    }

    int JJ=1,Jitem=-1,nMis;
    for(int i=1;i<=14;i++){
        if(nint(par[i][4])==-1){
            idToCheckBox["checkB_mis"+QString::number(i)+"_1"] -> setCheckState ( Qt::Unchecked );
            idToCheckBox["checkB_mis"+QString::number(i)+"_2"] -> setCheckState ( Qt::Checked );
            idToCheckBox["checkB_mis"+QString::number(i)+"_3"] -> setCheckState ( Qt::Unchecked );
            idToCheckBox["checkB_mis"+QString::number(i)+"_3"] -> setEnabled(false);
            DATO[i]=-1;
        }
        else if(nint(par[i][4])==0){
            idToCheckBox["checkB_mis"+QString::number(i)+"_1"] -> setCheckState ( Qt::Unchecked );
            idToCheckBox["checkB_mis"+QString::number(i)+"_2"] -> setCheckState ( Qt::Unchecked );
            idToCheckBox["checkB_mis"+QString::number(i)+"_3"] -> setCheckState ( Qt::Unchecked );
            idToCheckBox["checkB_mis"+QString::number(i)+"_3"] -> setEnabled(false);
            DATO[i]=0;
        }
        else if(nint(par[i][4])==1){
            idToCheckBox["checkB_mis"+QString::number(i)+"_1"] -> setCheckState ( Qt::Checked );
            idToCheckBox["checkB_mis"+QString::number(i)+"_2"] -> setCheckState ( Qt::Checked );
            idToCheckBox["checkB_mis"+QString::number(i)+"_3"] -> setCheckState ( Qt::Unchecked );
            idToCheckBox["checkB_mis"+QString::number(i)+"_3"] -> setEnabled(true);
            DATO[i]=1;
        }
        else if(nint(par[i][4])==2){
            idToCheckBox["checkB_mis"+QString::number(i)+"_1"] -> setCheckState ( Qt::Checked );
            idToCheckBox["checkB_mis"+QString::number(i)+"_2"] -> setCheckState ( Qt::Checked );
            idToCheckBox["checkB_mis"+QString::number(i)+"_3"] -> setCheckState ( Qt::Checked );
            idToCheckBox["checkB_mis"+QString::number(i)+"_3"] -> setEnabled(true);
            DATO[i]=2;
        }
        nMis=nint(par[20+JJ][4]);
        qDebug("\ti=%d nMis=par[20+%d][4]=%d",i,JJ,nMis);
        if(nMis>=0){
            int iCount=idToComboBox["cBmis"+QString::number(i)] -> count();
            for(int ii=0;ii<iCount;ii++){
                QString Lcb=idToComboBox["cBmis"+QString::number(i)] -> itemText(ii);
                // item is "<v|i><digit>" (SP) or "<number>" (ELI): match the trailing number
                // EXACTLY. contains() matched substrings, so nMis=1 wrongly hit "10","11",...
                int p=Lcb.size();
                while(p>0 && Lcb.at(p-1).isDigit()) p--;
                if(p<Lcb.size() && Lcb.mid(p).toInt()==nMis)
                    idToComboBox["cBmis"+QString::number(i)] -> setCurrentIndex(ii);
            }
            //idToComboBox["cBmis"+QString::number(i)] -> setCurrentIndex(nMis);
        }
        if(i==1)
            pwTn();
        else if(i==2)
            pwTp();
        else if(i==3)
            pwRn();
        else if(i==4)
            pwRp();
        else if(i==5)
            pwR1();
        else if(i==6)
            pwApds();
        if(i==7 || i==9 || i==11 || i==13){//if(JJ>=7)//set angle and Psi checkBox
            if(i==7)
                pwEj(1);
            else if(i==9)
                pwEj(2);
            else if(i==11)
                pwEj(3);
            else if(i==13)
                pwEj(4);
            Jitem=idToComboBox["cBteE"+QString::number(JJ-6)] -> findText(QString::number(par[13+JJ-6][1]),Qt::MatchExactly);
            if(Jitem>=0) idToComboBox["cBteE"+QString::number(JJ-6)] ->setCurrentIndex(Jitem);
            // Signals are muted during load, so changing cBteE above does not fire
            // pwSubEj via its connection; call it explicitly so the simulation-tab
            // angle (dSB_PM_8x_1) and WL range follow the selected theta instead of
            // staying on the first theta set by the earlier pwEj() call.
            pwSubEj(JJ-6);
            i++;//PSI
            if(nint(par[i][4])==-1){
                idToCheckBox["checkB_mis"+QString::number(i)+"_2"] -> setCheckState ( Qt::Checked );
                idToCheckBox["checkB_mis"+QString::number(i)+"_3"] -> setCheckState ( Qt::Unchecked );
                idToCheckBox["checkB_mis"+QString::number(i)+"_3"] -> setEnabled(false);
                DATO[i]=-1;
            }
            else if(nint(par[i][4])==0){
                idToCheckBox["checkB_mis"+QString::number(i)+"_2"] -> setCheckState ( Qt::Unchecked );
                idToCheckBox["checkB_mis"+QString::number(i)+"_3"] -> setCheckState ( Qt::Unchecked );
                idToCheckBox["checkB_mis"+QString::number(i)+"_3"] -> setEnabled(false);
                DATO[i]=0;
            }
            else if(nint(par[i][4])==1){
                idToCheckBox["checkB_mis"+QString::number(i)+"_2"] -> setCheckState ( Qt::Checked );
                idToCheckBox["checkB_mis"+QString::number(i)+"_3"] -> setCheckState ( Qt::Unchecked );
                idToCheckBox["checkB_mis"+QString::number(i)+"_3"] -> setEnabled(true);
                DATO[i]=1;
            }
            else if(nint(par[i][4])==2){
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
        reset();
    }
    if(nint(par[9][1])==1)
        ui->checkBox_setPsoK -> setCheckState ( Qt::Checked );
    else
        ui->checkBox_setPsoK -> setCheckState ( Qt::Unchecked );
    if(nint(par[10][1])==0)
        ui->cBox_PAR_10_1 -> setCheckState (Qt::Unchecked);
    else
        ui->cBox_PAR_10_1 -> setCheckState (Qt::Checked);
    if(nint(par[10][2])==0)
        ui->checkBox_3points -> setCheckState (Qt::Checked);
    else
        ui->checkBox_3points -> setCheckState (Qt::Unchecked);
    if(old0new1==0){
        QMessageBox msgBox;
        msgBox.setText("The project kind is old => k=extinction_coefficient=0");
        msgBox.setInformativeText("Would you like to update to k evaluated from epsilon?");
        msgBox.setStandardButtons(QMessageBox::Yes | QMessageBox::No);
        msgBox.setDefaultButton(QMessageBox::No);
        int ret = msgBox.exec();
        switch (ret) {
        case QMessageBox::Yes:
            par[11][1]=1.;
            break;
        case QMessageBox::No:
            break;
        }
    }
    int expA=0;
    if(par[11][1]>0.)
        expA=nint(log10(par[11][1]));//attenuation coefficient of k by Fit#N
    else
        expA=0;
    double valA=par[11][1]/pow(10.,expA);
    ui->spinBox_expAENS->setValue(expA);
    ui->doubleSpinBox_valAENS->setValue(valA);
    int irif=nint(par[10][3]);
    ui->comboB_PAR_10_3 -> setCurrentIndex(irif);
    if(nint(par[13][2])==0){
        iDelCon=0;
        ui->checkBox_deltaConnect->setCheckState(Qt::Unchecked);
    }
    else{
        iDelCon=1;
        ui->checkBox_deltaConnect->setCheckState(Qt::Checked);
    }
    iw=nint(par[18][1]);
    if(iw==1)
        ui->checkB_par_18_1 -> setCheckState ( Qt::Checked );
    else
        ui->checkB_par_18_1 -> setCheckState ( Qt::Unchecked );
    int imultiply=nint(par[22][1]);
    if(imultiply==0)
        ui->checkB_PAR_22_1 -> setCheckState ( Qt::Unchecked );
    else{
        ui->checkB_PAR_22_1 -> setCheckState ( Qt::Checked );
        ui->comboB_PAR_10_3 -> setEnabled(true);
    }
    ui->DP_PAR_27_1  -> setText(QString::number(par[27][1]));
    int s1p2=nint(par[27][2]);
    ui->cB_PAR_27_2 -> setCurrentIndex(s1p2-1);
    ui->sB_PAR_28_1 -> setValue(nint(par[28][1]));
    ui->sB_PAR_28_2 -> setValue(nint(par[28][2]));
    ui->sB_PAR_29_1 -> setValue(nint(par[29][1]));
    if(nint(par[31][5])==0)
        ui->checkBox_logScale->setCheckState(Qt::Unchecked);
    else
        ui->checkBox_logScale->setCheckState(Qt::Checked);
    ui->comboBox_DeltaPsiScale ->setCurrentIndex(nint(par[33][5]));
    nlayer=nint(par[51][2]);//N, layer
    int isimmetry=nint(par[52][2]);
    if(isimmetry==0)
        ui->checkB_PAR_52_2 -> setCheckState ( Qt::Unchecked );
    else
        ui->checkB_PAR_52_2 -> setCheckState ( Qt::Checked );
    int ihemi=nint(par[54][2]);
    if(ihemi==0)
        ui->checkB_PAR_54_2 -> setCheckState ( Qt::Unchecked );
    else
        ui->checkB_PAR_54_2 -> setCheckState ( Qt::Checked );
    ui->DP_PAR_55_2 -> setText(QString::number(par[55][2]));
    npp=nint(par[34][5]);//N par fit
    if(npp>N_FIT_MAX) npp=N_FIT_MAX;//clamp legacy/hand-edited files: ppm[] is sized N_FIT_MAX+1 and the fit kernel caps at nparmax=N_FIT_MAX
    ui->cB_PAR_35_2 -> setCurrentIndex(nint(par[35][2]-1.));
    nPar=nint(par[35][5]);//N par enabled for fit
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

    ui->comboBox_average->setCurrentIndex(nint(par[35][1]));//set std spectrum

    for(int i=1;i<=npp;i++)
        ppm[i]= uniformIpFile ? nint(par[35+i][5]) : migrateIp(nint(par[35+i][5]));
    ensureFitRow(npp);   // create widget rows if saved project had more than current m_maxFitRow

    //PM — legacy encoding: layers 1-9, EMA 1-15/20, ThetaEli, oscillators
    for(int i=1;i<=200;i++){
        for(int j=1;j<=5;j++){
            stream>> pmAtLegacy(i)[j];
        }
    }
    if(a2nm<1.){//Angstrom -> nm (layers 1-9; layers 10-99 handled in [EXTRA]/[STACK_V2])
        for(int i=1;i<=9;i++){
            pmD[i][1]*=a2nm;
            pmRg[i][1]*=a2nm;
        }
    }
    if(old0new1<2){//renormalization of oscillator coeff for old project
        for(int i=1;i<=20;i++){
            if(pmOs[1+(i-1)*5][1]==1){//Lorentz
                pmOs[2+(i-1)*5][1]=pmOs[2+(i-1)*5][1]*PIG/2.;
                pmOs[4+(i-1)*5][1]=pmOs[4+(i-1)*5][1]/2.;
            }
            else if(pmOs[1+(i-1)*5][1]==2){//Quantum
                pmOs[2+(i-1)*5][1]=pmOs[2+(i-1)*5][1]*PIG*pmOs[4+(i-1)*5][1]*pmOs[3+(i-1)*5][1];
            }
            else if(pmOs[1+(i-1)*5][1]==3){//Inhomo
                double sqrln2=0.832554611;
                long double E0=std::abs(pmOs[3+(i-1)*5][1]);
                long double D=std::abs(pmOs[4+(i-1)*5][1]);
                long double Kinho=pow((D/sqrln2),2.)/2*exp(-pow(E0*sqrln2/D,2.))+
                                    E0*D/sqrln2/2.*sqrt(PIG)*erfcl(-E0*sqrln2/D);
                pmOs[2+(i-1)*5][1]=pmOs[2+(i-1)*5][1]*Kinho;
            }
            else if(pmOs[1+(i-1)*5][1]==5){//Drude
                double E0=std::abs(pmOs[3+(i-1)*5][1]);
                pmOs[2+(i-1)*5][1]=E0*E0*PIG/2.;
                pmOs[3+(i-1)*5][1]=0.;
            }
            else if(pmOs[1+(i-1)*5][1]==6){//Indirect-Cody
                double E0=std::abs(pmOs[3+(i-1)*5][1]);
                double De=std::abs(pmOs[5+(i-1)*5][1]);
                double KCi=E0*De+De*De/2.;
                pmOs[2+(i-1)*5][1]=pmOs[2+(i-1)*5][1]*KCi;
            }
            else if(pmOs[1+(i-1)*5][1]==7){//Indirect-Tauc
                double E0=std::abs(pmOs[3+(i-1)*5][1]);
                double De=std::abs(pmOs[5+(i-1)*5][1]);
                double r=E0/De;
                double KTi=2.*(8.*r*r*log((4.*r+1.)/(4.*r)) + 8.*pow((r+1.),2.)*log((4.*r+4.)/(4.*r+3.))
                           - (1.+8.*r*(r+1.))*log((4.*r+3.)/(4.*r+1.)) );
                pmOs[2+(i-1)*5][1]=pmOs[2+(i-1)*5][1]*KTi;
            }
            else if(pmOs[1+(i-1)*5][1]==8){//Direct-Cody
                double E0=std::abs(pmOs[3+(i-1)*5][1]);
                double K=std::abs(pmOs[5+(i-1)*5][1]);
                double KCd=PIG/16.*K*K*(2.*E0+K);
                pmOs[2+(i-1)*5][1]=pmOs[2+(i-1)*5][1]*KCd*PIG;
            }
            else if(pmOs[1+(i-1)*5][1]==9){//Direct-Tauc
                double E0=std::abs(pmOs[3+(i-1)*5][1]);
                double K=std::abs(pmOs[5+(i-1)*5][1]);
                double KTd=PIG/2.*(2.*E0+K-2.*sqrt(E0*(E0+K)));
                pmOs[2+(i-1)*5][1]=pmOs[2+(i-1)*5][1]*KTd;
            }
        }
    }

    //PF
    for(int i=1;i<=21;i++){
        for(int j=1;j<=7;j++){
            stream>> pf[j][i];
        }
    }
    //EXTRA — optional block for layers 10..N_LAYER_MAX (absent in v1 files)
    {
        QString tag;
        while(!stream.atEnd()){
            tag=stream.readLine().simplified();
            if(!tag.isEmpty()) break;
        }
        // The extended blocks below ([EXTRA_NK]/[EXTRA_NK_V2]/[EXTRA_CNK*]/[STACK_V2]) are
        // written by BOTH modern files (which have NO leading [EXTRA]) and legacy files
        // (which store layers 10-20 in an [EXTRA] block first). Enter for either tag: read
        // the legacy layers only when [EXTRA] is actually present, then fall through to the
        // extended-block parsing. (Gating all of this on [EXTRA] dropped file-nk 9+, CNK
        // 16+ and [STACK_V2] from every modern project on reload.)
        if(tag=="[EXTRA]" || tag=="[EXTRA_NK]"){
            QString tag2;
            if(tag=="[EXTRA]"){
            for(int i=10;i<=N_LAYER_MAX;i++){
                for(int j=1;j<=5;j++) stream>> par[50+i][j];
            }
            for(int i=10;i<=N_LAYER_MAX;i++){
                for(int j=1;j<=5;j++) stream>> pmD[i][j];
            }
            for(int i=10;i<=N_LAYER_MAX;i++){
                for(int j=1;j<=5;j++) stream>> pmDn[i][j];
            }
            for(int i=10;i<=N_LAYER_MAX;i++){
                for(int j=1;j<=5;j++) stream>> pmNc[i][j];
            }
            for(int i=10;i<=N_LAYER_MAX;i++){
                for(int j=1;j<=5;j++) stream>> pmDk[i][j];
            }
            for(int i=10;i<=N_LAYER_MAX;i++){
                for(int j=1;j<=5;j++) stream>> pmKc[i][j];
            }
            for(int i=10;i<=N_LAYER_MAX;i++){
                for(int j=1;j<=5;j++) stream>> pmRg[i][j];
            }
            for(int i=10;i<=N_LAYER_MAX;i++){
                for(int j=1;j<=5;j++) stream>> pmSd[i][j];
            }
            for(int i=10;i<=N_LAYER_MAX;i++){
                for(int j=1;j<=5;j++) stream>> pmNu[i][j];
            }
            if(a2nm<1.){
                for(int i=10;i<=N_LAYER_MAX;i++){
                    pmD[i][1]*=a2nm;
                    pmRg[i][1]*=a2nm;
                }
            }
                // legacy [EXTRA] held layers 10-20; advance to the next tag (the first
                // extended-block header, the same one a modern file starts at directly).
                while(!stream.atEnd()){
                    tag2=stream.readLine().simplified();
                    if(!tag2.isEmpty()) break;
                }
            }
            else
                tag2=tag;   // modern file: tag is already the first extended block
            {
                bool hasExtraCnk=false;
                if(tag2=="[EXTRA_NK]"){
                    int nnk=1;
                    { QString ln=stream.readLine().simplified(); QTextStream ls(&ln); ls >> nnk; }
                    ui->sB_nnk->setValue(nnk);
                    updateNkRows(nnk);
                    // NK files 9-20 (always present in EXTRA_NK for backward compat)
                    for(int i=9;i<=20;i++){
                        QString line=stream.readLine().simplified();
                        fnk[i]=pathroot+line+".nk";
                        idToLineEdit["lineEdit"+QString::number(i)]->setText(line);
                        if(!line.contains("mate/aa999"))
                            Setnk(i);
                    }
                    // Read next tag (may be [EXTRA_NK_V2], [EXTRA_CNK], [EXTRA_CNK_V2])
                    QString tag3;
                    while(!stream.atEnd()){
                        tag3=stream.readLine().simplified();
                        if(!tag3.isEmpty()) break;
                    }
                    // Optional [EXTRA_NK_V2]: NK files 21..N_CNK_USER_MAX
                    if(tag3=="[EXTRA_NK_V2]"){
                        int nkV2=0; stream>>nkV2; stream.readLine();
                        for(int i=21; i<21+nkV2 && i<=N_CNK_USER_MAX; i++){
                            QString line=stream.readLine().simplified();
                            fnk[i]=pathroot+line+".nk";
                            idToLineEdit["lineEdit"+QString::number(i)]->setText(line);
                            if(!line.contains("mate/aa999"))
                                Setnk(i);
                        }
                        while(!stream.atEnd()){
                            tag3=stream.readLine().simplified();
                            if(!tag3.isEmpty()) break;
                        }
                    }
                    if(tag3=="[EXTRA_CNK]"){
                        hasExtraCnk=true;
                        // v5.0 format: 9 entries for materials 16-24 (special at 21-24)
                        for(int i=16;i<=24;i++){
                            for(int j=1;j<=6;j++) stream >> cnk[i].col(j);
                            stream.readLine();
                        }
                        // v5.0→v5.1 migration: special media 21-24 → 100-103
                        for(int j=1;j<=6;j++){
                            cnk[100].col(j)=cnk[21].col(j);
                            cnk[101].col(j)=cnk[22].col(j);
                            cnk[102].col(j)=cnk[23].col(j);
                            cnk[103].col(j)=cnk[24].col(j);
                        }
                        for(int i=21;i<=24;i++) cnk[i]=MaterialCNK(); // slots 21-24 now user materials: reset
                        // remap layer material indices 21-24 → 100-103
                        for(int i=1;i<=N_LAYER_MAX;i++){
                            int m=nint(par[50+i][1]);
                            if(m>=21&&m<=24) par[50+i][1]=m+79; // 21→100,22→101,23→102,24→103
                        }
                    } else if(tag3=="[EXTRA_CNK_V2]"){
                        hasExtraCnk=true;
                        // v5.1 format: count + entries starting at 16
                        int cnt=0; stream>>cnt; stream.readLine();
                        for(int i=16; i<16+cnt && i<=N_CNK_MAX; i++){
                            for(int j=1;j<=6;j++) stream >> cnk[i].col(j);
                            stream.readLine();
                        }
                    }
                    // [STACK_V2] or legacy [EXTRA2] detection
                    if(hasExtraCnk){
                        QString tagE2;
                        while(!stream.atEnd()){
                            tagE2=stream.readLine().simplified();
                            if(!tagE2.isEmpty()) break;
                        }
                        if(tagE2=="[EXTRA2]"){
                            // legacy: stores only layers N_LAYER_MAX+1..n
                            int nExtra2=0; stream >> nExtra2;
                            int iEnd=std::min(N_LAYER_MAX+nExtra2, N_LAYER_HARD_MAX);
                            for(int i=N_LAYER_MAX+1;i<=iEnd;i++){
                                for(int j=1;j<=5;j++) stream>> par[50+i][j];
                            }
                            double (*pmPtrs[8])[6]={pmD,pmDn,pmNc,pmDk,pmKc,pmRg,pmSd,pmNu};
                            for(int a=0;a<8;a++){
                                for(int i=N_LAYER_MAX+1;i<=iEnd;i++){
                                    for(int j=1;j<=5;j++) stream>> pmPtrs[a][i][j];
                                }
                            }
                        } else if(tagE2=="[STACK_V2]"){
                            // v5.0+: unified block for ALL layers 1..N
                            int N=0; stream >> N;
                            N=std::min(N, N_LAYER_HARD_MAX);
                            double (*pmPtrs[8])[6]={pmD,pmDn,pmNc,pmDk,pmKc,pmRg,pmSd,pmNu};
                            for(int i=1;i<=N;i++){
                                int matIdx=0, ltype=0;
                                stream >> matIdx >> ltype;
                                par[50+i][1]=matIdx;
                                par[50+i][3]=ltype;
                                for(int a=0;a<8;a++)
                                    for(int j=1;j<=4;j++)
                                        stream >> pmPtrs[a][i][j];
                            }
                        }
                    }
                }
                if(!hasExtraCnk){
                    // very old format: CNK[9..12] held special media → migrate to CNK[100..103]
                    for(int j=1;j<=6;j++){
                        cnk[100].col(j)=cnk[9].col(j);
                        cnk[101].col(j)=cnk[10].col(j);
                        cnk[102].col(j)=cnk[11].col(j);
                        cnk[103].col(j)=cnk[12].col(j);
                    }
                    // clear old slots (now available as user materials 9-12)
                    for(int i=9;i<=15;i++) cnk[i]=MaterialCNK();
                    // remap layer material indices 9-12 → 100-103
                    for(int i=1;i<=N_LAYER_MAX;i++){
                        int m=nint(par[50+i][1]);
                        if(m>=9 && m<=12) par[50+i][1]=m+91; // 9→100,10→101,11→102,12→103
                    }
                }
                if(!hasExtraCnk){
                    // old format: CNK[9..12] held special media → migrate to CNK[21..24]
                    for(int j=1;j<=6;j++){
                        cnk[21].col(j)=cnk[9].col(j);
                        cnk[22].col(j)=cnk[10].col(j);
                        cnk[23].col(j)=cnk[11].col(j);
                        cnk[24].col(j)=cnk[12].col(j);
                    }
                    // clear old slots (now user materials 9-12)
                    for(int i=9;i<=15;i++)
                        for(int j=1;j<=6;j++) cnk[i].col(j)=0.;
                    // remap layer material indices that used old special media slots
                    for(int i=1;i<=N_LAYER_MAX;i++){
                        int m=nint(par[50+i][1]);
                        if(m>=9 && m<=12) par[50+i][1]=m+12; // 9→21,10→22,11→23,12→24
                    }
                }
            }
        }
    }
    file.close();

    // Phase 1 (v5.0): populate Stack in parallel with pm-arrays; assert consistency.
    stack = stackFromPm(nlayer);
    Q_ASSERT(stack.nLayers() == nlayer);
    for(int i = 0; i < nlayer; i++){
        Q_ASSERT(stack.layers[i].thickness.value     == pmD [i+1][1]);
        Q_ASSERT(stack.layers[i].roughness.value     == pmRg[i+1][1]);
        Q_ASSERT(stack.layers[i].nonUniformity.value == pmNu[i+1][1]);
        Q_ASSERT(stack.layers[i].materialIndex       == nint(par[51+i][1]));
    }

    //inserimento valori cnk
    double f2=0.;
    for(int i=1;i<=N_CNK_MAX;i++){
        idToComboBox["cB_cnk"+QString::number(i)+"a"] -> setCurrentIndex(nint(cnk[i].nSrc));
        if(nint(cnk[i].nSrc)==0){
            idToLineEdit["LEcnk"+QString::number(i)+"_2"] -> setEnabled(true);
            idToLineEdit["LEcnk"+QString::number(i)+"_3"] -> setEnabled(true);
            idToLineEdit["LEcnk"+QString::number(i)+"_2"] -> setText(QString::number(cnk[i].nConst));
            idToLineEdit["LEcnk"+QString::number(i)+"_3"] -> setText(QString::number(cnk[i].kConst));
        }
        else{
            idToLineEdit["LEcnk"+QString::number(i)+"_2"] -> setEnabled(false);
            idToLineEdit["LEcnk"+QString::number(i)+"_3"] -> setEnabled(false);
        }
        if(nint(cnk[i].emaSrc)<0.){
            idToCheckBox["cB_EMA_"+QString::number(i)] -> setCheckState ( Qt::Unchecked );
            idToComboBox["cB_cnk"+QString::number(i)+"b"] -> setCurrentIndex(-1);
            idToComboBox["cB_cnk"+QString::number(i)+"b"] -> setEnabled(false);
            idToLineEdit["LEcnk"+QString::number(i)+"_5"] -> setEnabled(false);
            idToLineEdit["LEcnk"+QString::number(i)+"_6"] -> setEnabled(false);
            idToDoubleSpinBox[emaKey(i)] -> setEnabled(false);
        }
        else{
            idToCheckBox["cB_EMA_"+QString::number(i)] -> setCheckState ( Qt::Checked );
            idToComboBox["cB_cnk"+QString::number(i)+"b"] -> setEnabled(true);
            idToComboBox["cB_cnk"+QString::number(i)+"b"] -> setCurrentIndex(nint(cnk[i].emaSrc));
            if(nint(cnk[i].emaSrc)==0){
                idToLineEdit["LEcnk"+QString::number(i)+"_5"] -> setEnabled(true);
                idToLineEdit["LEcnk"+QString::number(i)+"_6"] -> setEnabled(true);
                idToLineEdit["LEcnk"+QString::number(i)+"_5"] -> setText(QString::number(cnk[i].n2Const));
                idToLineEdit["LEcnk"+QString::number(i)+"_6"] -> setText(QString::number(cnk[i].k2Const));
            }
            else{
                idToLineEdit["LEcnk"+QString::number(i)+"_5"] -> setEnabled(false);
                idToLineEdit["LEcnk"+QString::number(i)+"_6"] -> setEnabled(false);
            }
            f2=pmFe[i][1];
            idToDoubleSpinBox[emaKey(i)] -> setEnabled(true);
            idToDoubleSpinBox[emaKey(i)] -> setValue(f2);
        }
    }
    // Compute nmat from cnk data: highest material index (2-N_CNK_USER_MAX) with non-default settings.
    // Default = constant nk (nSrc==0), n=1, k=0, no EMA.  Works for both old and new files.
    {
        int nmat = 1;
        for(int i = 2; i <= N_CNK_USER_MAX; i++){
            const bool isDefault = (nint(cnk[i].nSrc) == 0 &&
                                    std::abs(cnk[i].nConst - 1.0) < 1e-9 &&
                                    std::abs(cnk[i].kConst)       < 1e-9 &&
                                    cnk[i].emaSrc < 0.0);
            if(!isDefault) nmat = i;
        }
        ui->sB_nmat->setValue(nmat);
        updateMatRows(nmat);
    }
    //aggiornamento valori PanFit
    for(int i=1;i<=npp;i++){
        state=idToCheckBox["chBeParFit_"+QString::number(i)]-> checkState();
        int ip=nint(ppm[i]);
        if(ip>=897 && ip<=996)//oscillator params (uniform ip encoding v5.2)
            pmAt(ip)[1]=std::abs(pmAt(ip)[1]);
        idToLineEdit["DPparFitV_"+QString::number(i)] -> setEnabled(true);
        idToLineEdit["DPparFitErr_"+QString::number(i)] -> setEnabled(true);
        idToLineEdit["DPparFitGC_"+QString::number(i)] -> setEnabled(true);
        idToLineEdit["DPparFitV_"+QString::number(i)] -> setText(QString::number(pmAt(ip)[1]));
        if( state == Qt::Checked ){
            idToLineEdit["DPparFitErr_"+QString::number(i)] -> setText(QString::number(pmAt(ip)[4]));
            idToLineEdit["DPparFitGC_"+QString::number(i)] -> setText(QString::number(pmAt(ip)[5]));
        }
        else{
            idToLineEdit["DPparFitErr_"+QString::number(i)] -> setText(QString::number(0));
            idToLineEdit["DPparFitGC_"+QString::number(i)] -> setText(QString::number(0));
        }
    }
    for(int i=npp+1;i<=m_maxFitRow;i++){
        idToLineEdit["DPparFitV_"+QString::number(i)] -> setEnabled(false);
        idToLineEdit["DPparFitErr_"+QString::number(i)] -> setEnabled(false);
        idToLineEdit["DPparFitGC_"+QString::number(i)] -> setEnabled(false);
    }
    ui->cB_cnk1a -> setCurrentIndex(nint(pmOs[0][1]));//fit option
    ui->sB_PAR_51_2 -> setValue(nlayer);
    setGuiMute(false);
    qDebug()<<"ReadSetting completed. Now call setRifMir MCRange, SetModel, ListOsc, PanFitPar";

    setRifMir();
    MCRange();
    SetModel(nlayer);
    listOsc();
    PanFitPar();

    //set the WL range as stored in the project
    ui->dSB_PAR_4_1 -> setValue(WLmin);//Lmin
    ui->dSB_PAR_4_2 -> setValue(WLmax);//Lmax

    qDebug()<<"exit from ReadSetting";
}


void ksemawc::setRifMir(){
    Qt::CheckState state;
    state=ui->checkB_PAR_22_1 -> checkState();
    ui->comboB_PAR_10_3-> setEnabled(state==Qt::Checked);
    if(state==Qt::Checked){
        int nrif=ui->comboB_PAR_10_3 -> currentIndex();
        qDebug()<<"->setRifMir nrif="<<nrif;
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
    for(int i=1;i<=N_LAYER_MAX;i++){
        int iKind=idToComboBox[wCB3(i)] ->currentIndex();
        if(iKind>0){
            d=idToDoubleSpinBox[wD(i)] ->value();
            dUP=d;
            if(i>1)
               dUP=idToDoubleSpinBox[wD(i-1)] ->value();
            rghMax=min(d/3.,dUP/3.);
            idToDoubleSpinBox[wRg(i)] ->setMaximum(rghMax);
        }
        else
            idToDoubleSpinBox[wRg(i)] ->setMaximum(100.);
    }
}


void ksemawc::setRangeEli(){
    qDebug()<<"-> setRangeEli";
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
    qDebug()<<"-> readRangeEli";
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
    QString subfnproject,sample;
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
        // QFile::copy (no shell) — safe with '+'/spaces/long names in fnproject.
        ierr = copyFileOverwrite(fileStore, fnproject) ? 0 : 1;
        if(ierr != 0){
            msgErrLoad("SaveProject",fnproject);
            qWarning()<<"Error copying project!!!";
        }
        else
            qInfo()<<"Project saved as "<<fnproject;
    }
}


void ksemawc::SaveSetting(int iCall){
    if(isGuiMuted())
        return;
    qDebug()<<"-> SaveSetting option_iCall="<<iCall<<" to "<<fileStore;
    setGuiMute(true);
    Qt::CheckState state,state1,state2;
    QString stringa,lab,svalue;
    int itab;
    if(iCall<0)
        itab=ui->tabWidget -> currentIndex();
    else
        itab=iCall;
    qDebug()<<"itab="<<itab;
    if(itab==1){//Model TAB
        state=ui->checkBox_setPsoK -> checkState ();
        if( state == Qt::Checked )
            par[9][1]=1.;
        else
            par[9][1]=0.;
        for(int i=1;i<=stack.nLayers();i++){
            Layer& l        = stack.layers[i-1];
            l.materialIndex = idToComboBox[wCB1(i)]->currentIndex() + 1;
            l.type          = idToComboBox[wCB3(i)]->currentIndex();
            l.roughness.value     = idToDoubleSpinBox[wRg(i)]->value();
            l.nonUniformity.value = idToDoubleSpinBox[wNu(i)]->value();
            if(l.type == 0){
                l.thickness.value = idToDoubleSpinBox[wD(i)]->value() / 1.E-6; // mm→nm
            } else {
                l.thickness.value = idToDoubleSpinBox[wD(i)]->value();
            }
            if(l.type == 2){
                l.nGrad.value      = idToDoubleSpinBox[wDn(i)]->value();
                l.nCurv.value      = idToDoubleSpinBox[wNc(i)]->value();
                l.kGrad.value      = idToDoubleSpinBox[wDk(i)]->value();
                l.kCurv.value      = idToDoubleSpinBox[wKc(i)]->value();
                l.slopeNGrad.value = idToDoubleSpinBox[wSd(i)]->value();
            }
        }
        stackToPm(stack);
        // The dSB_PM_71..89 widgets are named in the legacy pm convention
        // (fEMA materials 1-15 → pmFe, ThetaEli 1-4 → pmTe); use pmAtLegacy so
        // they land in the right arrays, not in pmD[71..89] as plain pmAt would.
        for(int i=71;i<=89;i++)
            pmAtLegacy(i)[1]=idToDoubleSpinBox["dSB_PM_"+QString::number(i)+"_1"] -> value();//save f2 and ThetaEli values
        //refresh PanFit values
        for(int i=1;i<=npp;i++){
            state=idToCheckBox["chBeParFit_"+QString::number(i)]-> checkState();
            int ip=nint(ppm[i]);
            if(ip>=897 && ip<=996)//oscillator params (uniform ip encoding v5.2)
                pmAt(ip)[1]=std::abs(pmAt(ip)[1]);
            idToLineEdit["DPparFitV_"+QString::number(i)] -> setText(QString::number(pmAt(ip)[1]));
            if( state == Qt::Checked ){
                idToLineEdit["DPparFitErr_"+QString::number(i)] -> setText(QString::number(pmAt(ip)[4]));
                idToLineEdit["DPparFitGC_"+QString::number(i)] -> setText(QString::number(pmAt(ip)[5]));
            }
            else{
                idToLineEdit["DPparFitErr_"+QString::number(i)] -> setText(QString::number(0));
                idToLineEdit["DPparFitGC_"+QString::number(i)] -> setText(QString::number(0));
            }
        }
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
        //qDebug()<<"\t apply savePanFit";
        int n,n1,n2,ip,jpf;
        QString Lj,Valore;
        Qt::CheckState state;
        { int nMax=ui->sB_PAR_34_5->value(); n=nMax; npp=nMax; }
        jpf=0;
        for(int j=1;j<=npp;j++){
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
                if(n2==0 || (1<=n2 && n2<=7))//film parameters
                    ip=ipLayer(n1,n2);
                else if(8<=n2 && n2<=11)//oscillator parameters
                    ip = 896 + (n1-1)*5 + (n2-7) + 1;
                else if(n2==12)//fEMA
                    ip = 792 + Lj.mid(5).toInt();
                else if(n2==13)//Theta
                    ip = 891 + Lj.at(6).digitValue();

                ppm[j]=ip;
                Valore=idToLineEdit["DPparFitV_"+QString::number(j)] -> text();
                pmAt(ip)[1]=Valore.toDouble();

                //check value
                if((ip>=1 && ip<=99) || (ip>=496 && ip<=594) || (ip>=897 && ip<=996))//D, Rg, oscillators
                    if(pmAt(ip)[1]<.0)
                        pmAt(ip)[1]=0.;
                if(ip>=694 && ip<=891){//Nu (694-792) + EMA fractions (793-891)
                    if(pmAt(ip)[1]<.0)
                        pmAt(ip)[1]=0.;
                    if(pmAt(ip)[1]>1.)
                        pmAt(ip)[1]=1.;
                }
                idToLineEdit["DPparFitV_"+QString::number(j)] -> setText(QString::number(pmAt(ip)[1]));//eventually correct the value

                //qDebug()<<"\tSavePanFit:  j="<<j<<" n1="<<n1<<" n2="<<n2<<" ip="<<ip<<" pm["<<ip<<"][1]="<<pmAt(ip)[1];
                if( state == Qt::Unchecked ){
                    pmAt(ip)[2]=0;
                    n--;
                }
                else{
                    jpf++;
                    pmAt(ip)[2]=jpf;
                    pmAt(jpf)[3]=ip;
                }
            }
        }
        int nCeck=ui->sB_PAR_34_5 -> value();
        if(nCeck!=npp){
            QMessageBox msgBox;
            msgBox.setText("Please select each individual fit-parameter before fitting!!!");
            msgBox.setStandardButtons(QMessageBox::Ok);
            msgBox.exec();
            setGuiMute(false);
            return;
        }
        par[34][5]=npp;//N par fit
        par[35][5]=n;//N par enabled
        for(int i=1;i<=npp;i++){
            par[35+i][5]=ppm[i];
        }
        ui->sB_PAR_35_5 -> setValue(n);
        ui->sB_PAR_34_5 -> setValue(npp);
        // Propagate the just-edited layer fit values (pmD/pmRg/pmDn/...) into
        // stack.layers[] so the simulation engine actually uses them; SetModel
        // then refreshes the Model-tab widgets from the updated stack.
        pmValuesToStack(stack);
        SetModel(nlayer);
        listOsc();
    }
    int ci;
    double f2;
    QFile file(fileStore);
    if (!file.open(QIODevice::WriteOnly | QIODevice::Text)){
        msgErrLoad("SaveSetting-fileStore",fileStore);
        setGuiMute(false);
        return;
    }
    QTextStream out(&file);
    out << "iVspj=1nm ufip" << "\n";// "ufip" = ppm stored in uniform ip encoding (see ReadSetting)
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
        setGuiMute(false);
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
        par[27][4]=lab.toDouble();
        // Use the stored ThetaEli only if it is a refined value of the selected file
        // angle; if unset (<1) or far from it (stale/corrupt pmTe, e.g. projects saved
        // with the angle widget wrongly at 1), fall back to the reliable combo angle.
        double th1=ui->cBteE1 -> currentText().toDouble();
        if(pmTe[1][1]<1. || std::abs(pmTe[1][1]-th1)>5.)
            ui->dSB_PM_86_1-> setValue(th1);
        else
            ui->dSB_PM_86_1-> setValue(pmTe[1][1]);
    }
    lab=ui->cBmis9 -> currentText();
    state=ui->checkB_mis9_1->checkState();
    if(lab.isEmpty() || state==Qt::Unchecked)
        out << "mate/aa999.9" << "\n";
    else{
        NANK[13]=stringa+"."+lab;
        out << stringa+"."+lab << "\n";
        par[28][4]=lab.toDouble();
        double th2=ui->cBteE2 -> currentText().toDouble();
        if(pmTe[2][1]<1. || std::abs(pmTe[2][1]-th2)>5.)
            ui->dSB_PM_87_1-> setValue(th2);
        else
            ui->dSB_PM_87_1-> setValue(pmTe[2][1]);
    }
    lab=ui->cBmis11 -> currentText();
    state=ui->checkB_mis11_1->checkState();
    if(lab.isEmpty() || state==Qt::Unchecked)
        out << "mate/aa999.9" << "\n";
    else{
        NANK[14]=stringa+"."+lab;
        out << stringa+"."+lab << "\n";
        par[29][4]=lab.toDouble();
        double th3=ui->cBteE3 -> currentText().toDouble();
        if(pmTe[3][1]<1. || std::abs(pmTe[3][1]-th3)>5.)
            ui->dSB_PM_88_1-> setValue(th3);
        else
            ui->dSB_PM_88_1-> setValue(pmTe[3][1]);
    }
    lab=ui->cBmis13 -> currentText();
    state=ui->checkB_mis13_1->checkState();
    if(lab.isEmpty() || state==Qt::Unchecked)
        out << "mate/aa999.9" << "\n";
    else{
        NANK[15]=stringa+"."+lab;
        out << stringa+"."+lab << "\n";
        par[30][4]=lab.toDouble();
        double th4=ui->cBteE4 -> currentText().toDouble();
        if(pmTe[4][1]<1. || std::abs(pmTe[4][1]-th4)>5.)
            ui->dSB_PM_89_1-> setValue(th4);
        else
            ui->dSB_PM_89_1-> setValue(pmTe[4][1]);
    }
    //qDebug("saving par[27][4]=%f",par[27][4]);
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
    for(int i=1;i<=N_CNK_MAX;i++){
        cnk[i].nSrc=idToComboBox["cB_cnk"+QString::number(i)+"a"] -> currentIndex();
        cnk[i].nConst=idToLineEdit["LEcnk"+QString::number(i)+"_2"] -> text().toDouble();
        cnk[i].kConst=idToLineEdit["LEcnk"+QString::number(i)+"_3"] -> text().toDouble();
        // The EMA checkbox is the switch: only when it is ticked does the source combo
        // decide the 2nd material. Reading the combo unconditionally let a leftover index
        // (never cleared, e.g. a row the user never touched) turn EMA on behind his back.
        const bool emaOn=(idToCheckBox["cB_EMA_"+QString::number(i)] -> checkState()==Qt::Checked);
        cnk[i].emaSrc=emaOn ? idToComboBox["cB_cnk"+QString::number(i)+"b"] -> currentIndex() : -1.;
        cnk[i].n2Const=idToLineEdit["LEcnk"+QString::number(i)+"_5"] -> text().toDouble();
        cnk[i].k2Const=idToLineEdit["LEcnk"+QString::number(i)+"_6"] -> text().toDouble();
        f2=idToDoubleSpinBox[emaKey(i)] -> value();
        pmFe[i][1]=f2;
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
        if(i<=6){
            if(lab.size() >= 2)
                par[20+JJ][4] = (lab.at(1)).digitValue();// in Qt6 digitValue() deve essere sostituito da unicode()
            else
                par[20+JJ][4] = 1.;  // valore di default
        }
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
    pmTe[1][1]=ui->dSB_PM_86_1 -> value();
    pmTe[2][1]=ui->dSB_PM_87_1 -> value();
    pmTe[3][1]=ui->dSB_PM_88_1 -> value();
    pmTe[4][1]=ui->dSB_PM_89_1 -> value();
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
    if(std::abs(par[21][3]) < 0.003*std::abs(rxy[16][2]-rxy[16][1]))
        par[21][3]=0.003*std::abs(rxy[16][2]-rxy[16][1]);
    if(std::abs(par[22][3]) < 0.003*std::abs(rxy[17][2]-rxy[17][1]))
        par[22][3]=0.003*std::abs(rxy[17][2]-rxy[17][1]);
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
    for(int i=1;i<=N_LAYER_MAX;i++){
        if(i<=nlayer){
            par[50+i][1]=idToComboBox[wCB1(i)] -> currentIndex();
            par[50+i][1]++;
            par[50+i][3]=idToComboBox[wCB3(i)] -> currentIndex();
            par[50+i][3]++;
            if(nint(par[50+i][1])==1) par[53][2]=i;//puntatore strato incognito
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

    pmOs[0][1]=ui->cB_cnk1a -> currentIndex();//fit option

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
            out << QString::number(cnk[i].col(j),'g',7) << "\t";
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
    //save pm — legacy encoding: layers 1-9, EMA 1-15/20, ThetaEli, oscillators
    for(int i=1;i<=200;i++){
        for(int j=1;j<=5;j++){
            out << QString::number(pmAtLegacy(i)[j],'g',7) << "\t";
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
    //save extra nk files 9-20 (v5.1: 9-20 in EXTRA_NK, 21-99 in EXTRA_NK_V2)
    out << "[EXTRA_NK]\n";
    out << ui->sB_nnk->value() << "\n";
    for(int i=9;i<=20;i++){
        out << idToLineEdit["lineEdit"+QString::number(i)]->text() << "\n";
    }
    // NK files 21-N_CNK_USER_MAX in separate block for backward compat
    out << "[EXTRA_NK_V2]\n" << (N_CNK_USER_MAX-20) << "\n";
    for(int i=21;i<=N_CNK_USER_MAX;i++){
        out << idToLineEdit["lineEdit"+QString::number(i)]->text() << "\n";
    }
    // save all CNK rows 16-N_CNK_MAX (user 16-99 + special 100-103)
    int cnkSaveN = N_CNK_MAX - 16 + 1;  // 88
    out << "[EXTRA_CNK_V2]\n" << cnkSaveN << "\n";
    for(int i=16;i<=N_CNK_MAX;i++){
        for(int j=1;j<=6;j++)
            out << QString::number(cnk[i].col(j),'g',7) << "\t";
        out << "\n";
    }
    // [STACK_V2]: unified clean storage for all layers (v5.0+, replaces legacy [EXTRA2])
    {
        int N = std::min(stack.nLayers(), N_LAYER_HARD_MAX);
        out << "[STACK_V2]\n" << N << "\n";
        double (*pmPtrs[8])[6] = {pmD, pmDn, pmNc, pmDk, pmKc, pmRg, pmSd, pmNu};
        for(int i = 1; i <= N; i++){
            out << stack.layers[i-1].materialIndex << "\t" << (stack.layers[i-1].type + 1);
            for(int a = 0; a < 8; a++)
                for(int j = 1; j <= 4; j++)
                    out << "\t" << QString::number(pmPtrs[a][i][j], 'g', 7);
            out << "\n";
        }
    }
    file.close();
    qDebug()<<"exit from SaveSetting";
    //qDebug("with par[27][4]=%f",par[27][4]);
    setGuiMute(false);
}


void ksemawc::msgErrLoad(QString where,QString fnERR){
    fflush(stdout);
    QMessageBox msgBox;
    msgBox.setText(where+": error loading file:\n"+fnERR);
    msgBox.setStandardButtons(QMessageBox::Ok);
    msgBox.exec();
}


void ksemawc::LoadFilenk(){
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
    qDebug()<<"LoadFilenk-> load file nk-solutions= "<<fnFnk;
    QFile file(fnFnk);
    if (!file.open (QIODevice::ReadOnly | QIODevice::Text)){
        qWarning()<<"\t ERROR!!!";
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
    qDebug()<<"\tInfo: "<<line.simplified();
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
    qDebug()<<"\tN. (L,n,k)= "<<N<<" a2nm="<<a2nm;
    rxy[16][3]=1000.;
    rxy[16][4]=-1000.;
    rxy[17][3]=1000.;
    rxy[17][4]=-1000.;
    Sol.clear();
    Sol.resize(N);
    //load data
    for(int L=0;L<N;L++){
        line=stream.readLine();
        line=line.simplified();
        QStringList List;
        List =line.split(" ");
        int nV=List.count();
        if(nV!=5){
            List=line.split("\t");//retry with tab separator
            nV=List.count();
            if(nV!=5){
                qDebug()<<"nV="<<nV<<" line="<<line;
                qWarning()<<"LoadFilenk-> ERROR reading file-nk= "<<fnam<<": separator char invalid";
                QMessageBox msgBox;
                msgBox.setText("LoadFilenk-> ERROR reading file-nk="+fnam+"\n: separator char invalid or ERRn and ERRk missing");
                msgBox.setStandardButtons(QMessageBox::Ok);
                msgBox.exec();
                return;
            }
        }
        for(int kk=1;kk<=nV;kk++)
            Sol[L][kk]=List.at(kk-1).toDouble();
        Sol[L][0]=1.;//enable data
        Sol[L][1]=Sol[L][1]*a2nm;//wavelength in nm
        //rxy[20][3]=min(rxy[20][3],Sol[L][1]);
        //rxy[20][4]=max(rxy[20][4],Sol[L][1]);
        rxy[16][3]=min(rxy[16][3],Sol[L][2]-Sol[L][4]);
        rxy[16][4]=max(rxy[16][4],Sol[L][2]+Sol[L][4]);
        rxy[17][3]=min(rxy[17][3],Sol[L][3]-Sol[L][5]);
        rxy[17][4]=max(rxy[17][4],Sol[L][3]+Sol[L][5]);
    }
    file.close();
    ui->DP_RXY_16_3->setText(QString::number(rxy[16][3]));
    ui->DP_RXY_16_4->setText(QString::number(rxy[16][4]));
    ui->DP_RXY_17_3->setText(QString::number(rxy[17][3]));
    ui->DP_RXY_17_4->setText(QString::number(rxy[17][4]));
    int NORD=2;
    std::vector<std::vector<double>> NKNEW(N, std::vector<double>(6, 0.0));
    for(int I=0;I<N;I++){
        double WWM=1.E30;
        for(int H=0;H<N;H++){
            if(Sol[H][1]>=.0 && Sol[H][1]<WWM){
                WWM=Sol[H][1];
                NORD=H;
            }
        }
        for(int k=1;k<=5;k++){
            NKNEW[I][k]=Sol[NORD][k];
        }
        Sol[NORD][1]=-1.;
    }
    int i=0;
    int inew=0;
    while(i<N){
        if(NKNEW[i][4]<=.0)
            NKNEW[i][4]=0.001*NKNEW[i][2];//set minimal n-error
        if(NKNEW[i][5]<=0.){//set minimal k-error  (was Sol[i][5]: unsorted row; must test the sorted row like NKNEW[i][4] above)
            if(NKNEW[i][3]>0.)
                NKNEW[i][5]=0.001*NKNEW[i][3];
            else
                NKNEW[i][5]=0.001*(rxy[17][4]-rxy[17][3]);
        }
        for(int k=1;k<=5;k++)
            Sol[inew][k]=NKNEW[i][k];
        Sol[inew][0]=1.;//enabled for fit
        inew++;
        i++;
    }
    rxy[20][3]=min(rxy[20][3],Sol[0][1]);
    rxy[20][4]=max(rxy[20][4],Sol[N-1][1]);
}


void ksemawc::ClrFnk(){
    ui->lineEdit_Fnk -> setText("mate/aa999.9");
    ui->lineEdit_infoFnk ->setText("");
    NANK[9]="mate/aa999.9";
    Sol.clear();
}


void ksemawc::SaveFnk(){
    qDebug()<<"-> SaveFnk with lastAction="<<lastAction;
    QString line;
    int iCase=0;
    int N;
    if(lastAction.contains("FitExpMeas") || lastAction.contains("CalcMis") || lastAction.contains("IbridCurrent"))
        N=NeV;
    else
        N=Sol.size();
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
    pippo.setOkButtonText("OK");
    //pippo.setCancelButtonText("No, thanks!");
    bool ok=pippo.exec();
    QString text =pippo.textValue();
    if (ok && !text.isEmpty()){
        line=text;
    }
    QFile file(fnFnk);
    if(!file.open(QIODevice::WriteOnly | QIODevice::Text)){
        msgErrLoad("SaveFnk",fnFnk);
        qWarning()<<"IONK-> ERROR opening file= "<<fnFnk;
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
        if(lastAction.contains("FitExpMeas") || lastAction.contains("CalcMis") || lastAction.contains("IbridCurrent")){
            for(int L=0;L<N;L++)
                stream<<ms.lambda[L]<<"\t"<<nkSol.n[L]<<"\t"<<nkSol.k[L]<<"\t"<<nkSol.n[L]*0.001<<"\t"<<nkSol.k[L]*0.001<<"\n";
        }
        else{
            for(int L=0;L<N;L++)
                stream<<Sol[L][1]<<"\t"<<Sol[L][2]<<"\t"<<Sol[L][3]<<"\t"<<Sol[L][4]<<"\t"<<Sol[L][5]<<"\n";
        }
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
    qDebug()<<"->Setnk("<<ifile<<" with fnk["<<ifile<<"]="<<fnk[ifile];
    double val1,val2;
    QTextStream stream ( &file );
    QString line,line2;
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
    if(val1<val2){
        idToLineEdit["WLminNK"+QString::number(ifile)] -> setText(QString::number(nint(val1)));
        idToLineEdit["WLmaxNK"+QString::number(ifile)] -> setText(QString::number(nint(val2)));
    }
    else{
        idToLineEdit["WLminNK"+QString::number(ifile)] -> setText(QString::number(nint(val2)));
        idToLineEdit["WLmaxNK"+QString::number(ifile)] -> setText(QString::number(nint(val1)));
    }
    file.close();
    MCRange();
    line2=idToLineEdit["lineEdit"+QString::number(ifile)] -> text();
    updateMatenk(7+ifile,line2);
}


void ksemawc::mDwUp(int iLayer, int Dw1UpM1){
    if(Dw1UpM1 == 1)
        stack.moveDown(iLayer - 1);  // 0-based index
    else
        stack.moveUp(iLayer - 1);
    stackToPm(stack);
    SetModel(nlayer);
}


void ksemawc::Clrfn(){
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
        idToLineEdit["WLmin"+QString::number(imis)] -> clear();
        idToLineEdit["WLmax"+QString::number(imis)] -> clear();
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
    ui->cBteE1 -> clear();
    ui->cBteE2 -> clear();
    ui->cBteE3 -> clear();
    ui->cBteE4 -> clear();
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


void ksemawc::listMeas(){
    if(fnSample.isEmpty() || fnSample.contains("mate/aa999"))
        return;
    qDebug()<<"->listMeas  with fnSample= "<<fnSample;
    QString lab,estens[7];
    estens[1]=".tn";
    estens[2]=".tp";
    estens[3]=".rn";
    estens[4]=".rp";
    estens[5]=".r1";
    estens[6]=".an";

    ui->lineEdit_Tn->clear();
    ui->lineEdit_Tp->clear();
    ui->lineEdit_Rn->clear();
    ui->lineEdit_Rp->clear();
    ui->lineEdit_R1->clear();
    ui->lineEdit_Apds->clear();
    for(int imis=1;imis<=6;imis++){
        //idToCheckBox["checkB_mis"+QString::number(imis)+"_1"] ->setCheckState(Qt::Unchecked);
        idToComboBox["cBmis"+QString::number(imis)] -> clear();
        idToLineEdit["WLmin"+QString::number(imis)] -> clear();
        idToLineEdit["WLmax"+QString::number(imis)] -> clear();
        for(int j=0;j<2;j++){
            lab="v";
            if(j==1) lab="i";
            for(int i=0;i<10;i++){
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
    for(int i=0;i<10;i++){
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
        // After a fresh manual data load all four ELI rows default to the first
        // measurement file (listMeas selects index 0 in every cBmis7/9/11/13).
        // Spread the angles (theta) of that file across consecutive rows
        // (ELI-1 -> 1st theta, ELI-2 -> 2nd theta, ...) and DESELECT the rows that
        // have no theta left. Otherwise every row reloads the same file at theta 0,
        // duplicating the data into ELI-3/ELI-4. The number of active rows must
        // equal the number of angles in the file. (Rows pointing to a different
        // measurement are left untouched: the user selected those manually.)
        QString labRef=ui->cBmis7->currentText();
        int nThetaRef=ui->cBteE1->count();
        int iThetaPre=ui->cBteE1->currentIndex();
        for(int j=2;j<=4;j++){
            QString lab=idToComboBox["cBmis"+QString::number(5+2*j)]->currentText();
            if(labRef.isEmpty() || lab!=labRef)
                continue;//different (manually chosen) measurement: keep as is
            int iThetaNext=iThetaPre+1;
            if(iThetaNext<nThetaRef){
                idToComboBox["cBteE"+QString::number(j)]->setCurrentIndex(iThetaNext);
                iThetaPre=iThetaNext;
            }
            else
                idToComboBox["cBmis"+QString::number(5+2*j)]->setCurrentIndex(-1);
        }
        iNewSample=0;
    }
}


void ksemawc::pwTn(){
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


void ksemawc::pwTp(){
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
    int s1p2=nint(par[27][2]);
    ui->cB_PAR_27_2 -> setCurrentIndex(s1p2-1);
    MCRange();
}


void ksemawc::pwRn(){
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


void ksemawc::pwRp(){
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


void ksemawc::pwR1(){
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


void ksemawc::pwApds(){
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
            qDebug("pwEj(%d) fnEli=%s",j,fnEli.toStdString().c_str());
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
                qDebug()<<"->pwEj("<<j<<") file="<<fnEli<<" Ntheta="<<ntel<<" Ndat="<<ndat<<" format="<<info;
                EliTab[j][0][0]=ntel;
                for(i=0;i<ntel;i++){
                    stream >> theta >> ndat >> wmin >> wmax;
                    EliTab[j][i][0]=theta;//per-theta row is 0-based [i] (was constant [ntel]); matches data rows below and the VASE branch
                    EliTab[j][i][1]=ndat;
                    EliTab[j][i][2]=wmin;
                    EliTab[j][i][3]=wmax;
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
    int ntel=static_cast<int>(EliTab[j][0][0]+0.5);
    if(ntel>0){
        int index=idToComboBox["cBteE"+QString::number(j)]->currentIndex();//ui->cBteE1 -> currentIndex();
        if(index<0)
            return;
        qDebug()<<"\t->pwSubEj("<<j<<") ntel="<<ntel<<" index="<<index;
        double theta=idToComboBox["cBteE"+QString::number(j)] -> currentText().toDouble();
        idToDoubleSpinBox["dSB_PM_"+QString::number(85+j)+"_1"]->setValue(theta);
        pmTe[j][1]=theta;
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
    if(isGuiMuted())
        return;
    qDebug("->MCRange");
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
    qDebug()<<"\t Lmin="<<Lmin<<" Lmax="<<Lmax;
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
    // Renaming items can clear the popup view's row-hidden flags, so re-assert
    // the nk-file filter to keep only the loaded nk files visible.
    updateNkRows(ui->sB_nnk->value());
}


void previewFile(QString filename, QString lab,QString& info,double& wmin,double& wmax){
    if(lab.isEmpty())
        return;
    qDebug()<<"->previewFile "<<filename<<" lab="<<lab;
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
        if(info.contains("PE UV", Qt::CaseInsensitive)){
            // file Perkin Elmer L900 - L950
            for(i=0;i<6;i++) line = stream.readLine();
            info = stream.readLine();
            for(i=0;i<5;i++) line = stream.readLine();
            if(line.contains("PerkinElmer UV WinLab 7.1", Qt::CaseInsensitive))
                ilinrim=69;
            else if(line.contains("PerkinElmer UV WinLab 7.4", Qt::CaseInsensitive))
                ilinrim=73;
            else if(line.contains("PerkinElmer UV WinLab 5", Qt::CaseInsensitive))
                ilinrim=65;
            else
                ilinrim=69;
            for(i=0;i<ilinrim;i++) line = stream.readLine();
            stream >> wmax;
            stream >> ridelta;
            stream >> ndati;
            wmin=wmax+ridelta*(ndati-1);
        }
        else if(line0.contains("#####SCALED") || line0.contains("#####SCALEDA") ||
                line0.contains("#####SCALED%")){
            // SCALED file
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
}


void ksemawc::SetModel(const int &){
    int nlayerOld=nlayer;
    nlayer=ui->sB_PAR_51_2 -> value();
    int Dnl=nlayer-nlayerOld;
    qDebug()<<"-> SetModel with Dnl="<<Dnl;
    setGuiMute(true);
    if(nlayer != stack.nLayers()){
        if(Dnl > 0)
            stack.layers.insert(stack.layers.begin(), Dnl, Layer{});
        else if(Dnl < 0)
            stack.layers.erase(stack.layers.begin(), stack.layers.begin() + (-Dnl));
        stackToPm(stack);
    }
    // Enable/disable Up/Down buttons based on nlayer
    for(int i=1; i<=N_LAYER_MAX; i++){
        idToPushButton["pBmDw"+QString::number(i)]->setEnabled(false);
        idToPushButton["pBmUp"+QString::number(i)]->setEnabled(false);
    }
    idToPushButton["pBmDw1"]->setEnabled(nlayer > 0);
    for(int i=2; i<=N_LAYER_MAX; i++){
        idToPushButton["pBmUp"+QString::number(i)]->setEnabled(nlayer >= i);
        idToPushButton["pBmDw"+QString::number(i)]->setEnabled(nlayer > i);
    }
    {
        int nL = stack.nLayers();
        ensureLayerRow(nL);
        int nLoop = std::max(m_maxLayerRow, nL);
    for(int i=1;i<=nLoop;i++){
        bool vis=(i<=nL);
        idToComboBox[wCB3(i)]->setVisible(vis);
        idToComboBox[wCB1(i)]->setVisible(vis);
        idToDoubleSpinBox[wD(i)]->setVisible(vis);
        idToDoubleSpinBox[wDn(i)]->setVisible(vis);
        idToDoubleSpinBox[wNc(i)]->setVisible(vis);
        idToDoubleSpinBox[wDk(i)]->setVisible(vis);
        idToDoubleSpinBox[wKc(i)]->setVisible(vis);
        idToDoubleSpinBox[wRg(i)]->setVisible(vis);
        idToDoubleSpinBox[wSd(i)]->setVisible(vis);
        idToDoubleSpinBox[wNu(i)]->setVisible(vis);
        idToPushButton["pBmUp"+QString::number(i)]->setVisible(vis);
        idToPushButton["pBmDw"+QString::number(i)]->setVisible(vis);
        if(m_layerGrid)
            if(auto* it = m_layerGrid->itemAtPosition(i, 0))
                if(auto* lbl = qobject_cast<QLabel*>(it->widget()))
                    lbl->setVisible(vis);
        if(i <= nL){
            Layer& l = stack.layers[i-1];
            if(l.type < 0 || l.type > 4) l.type = 0;
            if(l.materialIndex < 1 || l.materialIndex > N_CNK_MAX) l.materialIndex = 1;
            idToComboBox[wCB1(i)]->setEnabled(true);
            idToComboBox[wCB3(i)]->setEnabled(true);
            idToComboBox[wCB1(i)]->setCurrentIndex(l.materialIndex - 1);
            idToComboBox[wCB3(i)]->setCurrentIndex(l.type);
            if(l.type == 0){
                idToDoubleSpinBox[wD(i)]  ->setEnabled(true);
                idToDoubleSpinBox[wD(i)]  ->setValue(l.thickness.value * 1.E-6);
                idToDoubleSpinBox[wDn(i)] ->setEnabled(false);
                idToDoubleSpinBox[wNc(i)] ->setEnabled(false);
                idToDoubleSpinBox[wDk(i)] ->setEnabled(false);
                idToDoubleSpinBox[wKc(i)] ->setEnabled(false);
                idToDoubleSpinBox[wSd(i)] ->setEnabled(false);
                idToDoubleSpinBox[wRg(i)] ->setEnabled(true);
                idToDoubleSpinBox[wRg(i)] ->setValue(l.roughness.value);
                idToDoubleSpinBox[wNu(i)] ->setEnabled(false);
            } else if(l.type == 1){
                idToDoubleSpinBox[wD(i)]  ->setEnabled(true);
                idToDoubleSpinBox[wD(i)]  ->setValue(l.thickness.value);
                idToDoubleSpinBox[wDn(i)] ->setEnabled(false);
                idToDoubleSpinBox[wNc(i)] ->setEnabled(false);
                idToDoubleSpinBox[wDk(i)] ->setEnabled(false);
                idToDoubleSpinBox[wKc(i)] ->setEnabled(false);
                idToDoubleSpinBox[wSd(i)] ->setEnabled(false);
                idToDoubleSpinBox[wRg(i)] ->setEnabled(true);
                idToDoubleSpinBox[wRg(i)] ->setValue(l.roughness.value);
                idToDoubleSpinBox[wNu(i)] ->setEnabled(true);
                idToDoubleSpinBox[wNu(i)] ->setValue(l.nonUniformity.value);
            } else if(l.type == 2){
                idToDoubleSpinBox[wD(i)]  ->setEnabled(true);
                idToDoubleSpinBox[wDn(i)] ->setEnabled(true);
                idToDoubleSpinBox[wNc(i)] ->setEnabled(true);
                idToDoubleSpinBox[wDk(i)] ->setEnabled(true);
                idToDoubleSpinBox[wKc(i)] ->setEnabled(true);
                idToDoubleSpinBox[wRg(i)] ->setEnabled(true);
                idToDoubleSpinBox[wSd(i)] ->setEnabled(true);
                idToDoubleSpinBox[wNu(i)] ->setEnabled(true);
                idToDoubleSpinBox[wD(i)]  ->setValue(l.thickness.value);
                idToDoubleSpinBox[wDn(i)] ->setValue(l.nGrad.value);
                idToDoubleSpinBox[wNc(i)] ->setValue(l.nCurv.value);
                idToDoubleSpinBox[wDk(i)] ->setValue(l.kGrad.value);
                idToDoubleSpinBox[wKc(i)] ->setValue(l.kCurv.value);
                idToDoubleSpinBox[wRg(i)] ->setValue(l.roughness.value);
                idToDoubleSpinBox[wSd(i)] ->setValue(l.slopeNGrad.value);
                idToDoubleSpinBox[wNu(i)] ->setValue(l.nonUniformity.value);
            }
        } else{
            idToComboBox[wCB1(i)]->setEnabled(false);
            idToComboBox[wCB3(i)]->setEnabled(false);
            idToDoubleSpinBox[wD(i)]  ->setEnabled(false);
            idToDoubleSpinBox[wDn(i)] ->setEnabled(false);
            idToDoubleSpinBox[wNc(i)] ->setEnabled(false);
            idToDoubleSpinBox[wDk(i)] ->setEnabled(false);
            idToDoubleSpinBox[wKc(i)] ->setEnabled(false);
            idToDoubleSpinBox[wRg(i)] ->setEnabled(false);
            idToDoubleSpinBox[wSd(i)] ->setEnabled(false);
            idToDoubleSpinBox[wNu(i)] ->setEnabled(false);
        }
    }
    } // end nLoop block
    //refresh EMA fraction — all user materials + special media
    for(int i=1;i<=N_CNK_MAX;i++)
        idToDoubleSpinBox[emaKey(i)] -> setValue(pmFe[i][1]);
    setGuiMute(false);
}


void ksemawc::updateNkRows(int n){
    for(int i=1; i<=N_CNK_USER_MAX; i++){
        bool vis=(i<=n);
        idToPushButton["pBnk"+QString::number(i)]->setVisible(vis);
        idToPushButton["pBclearnk"+QString::number(i)]->setVisible(vis);
        idToLineEdit["lineEdit"+QString::number(i)]->setVisible(vis);
        idToLineEdit["lineEdit_infoNK_"+QString::number(i)]->setVisible(vis);
        idToLineEdit["WLminNK"+QString::number(i)]->setVisible(vis);
        idToLineEdit["WLmaxNK"+QString::number(i)]->setVisible(vis);
    }
    // Filter NK-file items in all material source combos (cB_cnk{i}a and cB_cnk{i}b).
    // Items 0-7 = "constant nk" + "Fit #1-7": always visible.
    // Item k >= 8 = "nk-(k-7)": visible only if k-7 <= n (file actually loaded).
    const QString suffixes[2] = {"a", "b"};
    for(int i = 1; i <= N_CNK_MAX; i++){
        QString si = QString::number(i);
        for(const QString& suf : suffixes){
            QComboBox* cb = idToComboBox.value("cB_cnk"+si+suf, nullptr);
            if(!cb) continue;
            for(int k = 8; k < cb->count(); k++){
                // "Material#1" is a trailing non-nk source option (index NSRC_MATREF1):
                // always keep it visible, only the unloaded nk-file slots get hidden.
                bool hide = (k - 7 > n) && cb->itemText(k) != "Material#1";
                if(auto* lv = qobject_cast<QListView*>(cb->view()))
                    lv->setRowHidden(k, hide);
            }
        }
    }
}


void ksemawc::updateMatRows(int n){
    if(!m_cnkGrid) return;
    // Row 0 of m_cnkGrid = column-header row (always visible).
    // Rows 1..N_CNK_USER_MAX-1 correspond to materials 2..N_CNK_USER_MAX.
    for(int r = 1; r <= N_CNK_USER_MAX-1; r++){
        bool vis = (r < n);   // visible when material (r+1) <= n
        for(int c = 0; c <= 8; c++){
            auto* item = m_cnkGrid->itemAtPosition(r, c);
            if(item && item->widget()) item->widget()->setVisible(vis);
        }
    }
    // Filter layer material combos: show only items 0..n-1 (Materials #1..#n)
    for(int i = 1; i <= m_maxLayerRow; i++){
        QComboBox* cb = idToComboBox.value(wCB1(i), nullptr);
        if(!cb) continue;
        for(int j = 0; j < cb->count(); j++)
            if(auto* lv=qobject_cast<QListView*>(cb->view())) lv->setRowHidden(j, j >= n);
        if(cb->isEnabled() && cb->currentIndex() >= n)
            cb->setCurrentIndex(n - 1);
    }
}


void ksemawc::RefreshModel(){
    if(isGuiMuted())
        return;
    qDebug()<<"-> RefreshModel";
    for(int i=1;i<=std::max(m_maxLayerRow, stack.nLayers());i++){
        if(i <= stack.nLayers()){
            Layer& l        = stack.layers[i-1];
            l.type          = idToComboBox[wCB3(i)]->currentIndex();
            l.materialIndex = idToComboBox[wCB1(i)]->currentIndex() + 1;
            par[50+i][1]    = l.materialIndex;
            par[50+i][3]    = l.type + 1;
            if(l.type == 0){
                idToDoubleSpinBox[wD(i)]  ->setEnabled(true);
                idToDoubleSpinBox[wDn(i)] ->setEnabled(false);
                idToDoubleSpinBox[wNc(i)] ->setEnabled(false);
                idToDoubleSpinBox[wDk(i)] ->setEnabled(false);
                idToDoubleSpinBox[wKc(i)] ->setEnabled(false);
                idToDoubleSpinBox[wSd(i)] ->setEnabled(false);
                idToDoubleSpinBox[wRg(i)] ->setEnabled(true);
                idToDoubleSpinBox[wNu(i)] ->setEnabled(false);
            } else if(l.type == 1){
                idToDoubleSpinBox[wD(i)]  ->setEnabled(true);
                idToDoubleSpinBox[wDn(i)] ->setEnabled(false);
                idToDoubleSpinBox[wNc(i)] ->setEnabled(false);
                idToDoubleSpinBox[wDk(i)] ->setEnabled(false);
                idToDoubleSpinBox[wKc(i)] ->setEnabled(false);
                idToDoubleSpinBox[wSd(i)] ->setEnabled(false);
                idToDoubleSpinBox[wRg(i)] ->setEnabled(true);
                idToDoubleSpinBox[wRg(i)] ->setValue(l.roughness.value);
                idToDoubleSpinBox[wNu(i)] ->setEnabled(true);
                idToDoubleSpinBox[wNu(i)] ->setValue(l.nonUniformity.value);
            } else if(l.type == 2){
                idToDoubleSpinBox[wD(i)]  ->setEnabled(true);
                idToDoubleSpinBox[wDn(i)] ->setEnabled(true);
                idToDoubleSpinBox[wNc(i)] ->setEnabled(true);
                idToDoubleSpinBox[wDk(i)] ->setEnabled(true);
                idToDoubleSpinBox[wKc(i)] ->setEnabled(true);
                idToDoubleSpinBox[wRg(i)] ->setEnabled(true);
                idToDoubleSpinBox[wSd(i)] ->setEnabled(true);
                idToDoubleSpinBox[wNu(i)] ->setEnabled(true);
                idToDoubleSpinBox[wDn(i)] ->setValue(l.nGrad.value);
                idToDoubleSpinBox[wNc(i)] ->setValue(l.nCurv.value);
                idToDoubleSpinBox[wDk(i)] ->setValue(l.kGrad.value);
                idToDoubleSpinBox[wKc(i)] ->setValue(l.kCurv.value);
                idToDoubleSpinBox[wRg(i)] ->setValue(l.roughness.value);
                idToDoubleSpinBox[wSd(i)] ->setValue(l.slopeNGrad.value);
                idToDoubleSpinBox[wNu(i)] ->setValue(l.nonUniformity.value);
            }
        } else{
            idToComboBox[wCB1(i)]->setEnabled(false);
            idToComboBox[wCB3(i)]->setEnabled(false);
            idToDoubleSpinBox[wD(i)]  ->setEnabled(false);
            idToDoubleSpinBox[wDn(i)] ->setEnabled(false);
            idToDoubleSpinBox[wNc(i)] ->setEnabled(false);
            idToDoubleSpinBox[wDk(i)] ->setEnabled(false);
            idToDoubleSpinBox[wKc(i)] ->setEnabled(false);
            idToDoubleSpinBox[wRg(i)] ->setEnabled(false);
            idToDoubleSpinBox[wSd(i)] ->setEnabled(false);
            idToDoubleSpinBox[wNu(i)] ->setEnabled(false);
        }
    }
    //refresh EMA fraction — all user materials + special media
    for(int i=1;i<=N_CNK_MAX;i++)
        idToDoubleSpinBox[emaKey(i)] -> setValue(pmFe[i][1]);
}


void ksemawc::setMat(int nM){
    qDebug()<<"->setMat called by cnk "<<nM;
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
    qDebug()<<"->setEMA called by checkBox "<<nM;
    idToComboBox["cB_cnk"+QString::number(nM)+"b"] -> setEnabled(state==Qt::Checked);
    idToDoubleSpinBox[emaKey(nM)] -> setEnabled(state==Qt::Checked);
    if(state!=Qt::Checked)
        idToComboBox["cB_cnk"+QString::number(nM)+"b"] ->setCurrentIndex(-1);
    int curI=idToComboBox["cB_cnk"+QString::number(nM)+"b"] ->currentIndex();
    idToLineEdit["LEcnk"+QString::number(nM)+"_5"] -> setEnabled(curI==0 && state==Qt::Checked);
    idToLineEdit["LEcnk"+QString::number(nM)+"_6"] -> setEnabled(curI==0 && state==Qt::Checked);
}


void ksemawc::setKindOsc(int N) {
    if(isGuiMuted())
        return;
    int pmIdx = 101 + 5 * (N - 1);   // legacy index: only used to build the widget name
    QString name = QString("cBpm_%1_1").arg(pmIdx);
    QComboBox* w = findChild<QComboBox*>(name);
    // Oscillator kind storage is pmOs[1+5*(N-1)] (what listOsc reads); pmAt(pmIdx)
    // under the uniform ip encoding would resolve to a layer gradient, not here.
    if (w) pmOs[1 + 5 * (N - 1)][1] = w->currentIndex() + 1;
    listOsc();
}


void ksemawc::setOscN(int k){
    if(isGuiMuted())
        return;
    setTabSim();//save current setting
    int ioptFit=ui->cB_cnk1a -> currentIndex();//Fit#
    int Nosc=nint(pf[ioptFit][1]);//N. Osc
    Qt::CheckState state;
    state=idToCheckBox["cBosc_"+QString::number(k)]->checkState();
    qDebug()<<"->setOscN ioptFit="<<ioptFit<<" N.Osc="<<Nosc<<" call-by-oscN="<<k;
    if(state==Qt::Checked){
        pf[ioptFit][Nosc+1+1]=k;
        pf[ioptFit][1]=Nosc+1;
    }
    else{
        Nosc++;//in the call of setTabSim Nosc has been already decreased
        int jj=1;
        while(nint(pf[ioptFit][jj+1])!=k && jj<Nosc){
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
    if(isGuiMuted())
        return;
    int iok;
    int ioptFit=ui->cB_cnk1a -> currentIndex();
    qDebug()<<"->ListOsc: ioptFit="<<ioptFit;
    setGuiMute(true);
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
                if(nint(pf[ioptFit][i])==k)iok=1;
            }
            if(iok==1){
                int ifu=nint(pmOs[1+(k-1)*5][1]);
                idToCheckBox["cBosc_"+QString::number(k)]-> setCheckState ( Qt::Checked );
                idToComboBox["cBpm_"+QString::number(100+1+(k-1)*5)+"_1"] -> setEnabled(true);
                idToComboBox["cBpm_"+QString::number(100+1+(k-1)*5)+"_1"] -> setCurrentIndex(nint(pmOs[1+(k-1)*5][1])-1);
                //idToComboBox["cBpm_"+QString::number(100+1+(k-1)*5)+"_1"] -> setSizeAdjustPolicy(QComboBox::AdjustToContents);
                for(int iC=2;iC<=5;iC++){
                    idToLineEdit["LEpm_"+QString::number(100+iC+(k-1)*5)+"_1"] -> setEnabled(true);
                    idToLineEdit["LEpm_"+QString::number(100+iC+(k-1)*5)+"_1"] -> setText("");
                }
                if(ifu<4){//Lorentz, Q-homo, Q-inhomo
                    idToLineEdit["LEpm_"+QString::number(100+5+(k-1)*5)+"_1"] -> setEnabled(false);
                    idToLineEdit["LEpm_"+QString::number(100+5+(k-1)*5)+"_1"] -> setText("unused");
                }
                else if(ifu==4){//Flat
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
                else if(ifu==19){
                    idToLineEdit["LEpm_"+QString::number(100+3+(k-1)*5)+"_1"] -> setEnabled(false);
                    idToLineEdit["LEpm_"+QString::number(100+3+(k-1)*5)+"_1"] -> setText("unused");
                }
                for(int iC=2;iC<=5;iC++){
                    QString content=idToLineEdit["LEpm_"+QString::number(100+iC+(k-1)*5)+"_1"] -> text();
                    if(content!="unused")
                       idToLineEdit["LEpm_"+QString::number(100+iC+(k-1)*5)+"_1"] -> setText(QString::number(std::abs(pmOs[iC+(k-1)*5][1])));
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
    setGuiMute(false);
}


void ksemawc::saveNKmodel(){
    int iok2;
    int ioptFit=ui->cB_cnk1a -> currentIndex();
    qDebug()<<"->saveOscModel ioptFit="<<ioptFit;
    if(ioptFit>0 && ioptFit<8){
        QString fn= QFileDialog::getSaveFileName(
                    this,
                    "Filename to save",
                    pathroot,
                    "nkMod file (*.nkMod)");
        if(fn.isEmpty())
            return;
        if(!fn.contains(".nkMod"))
            fn=fn+".nkMod";
        QInputDialog pippo;
        pippo.setLabelText("Please set a comment!");
        pippo.setInputMode(QInputDialog::TextInput);
        //pippo.setTextValue(QString::number(expTime));
        pippo.setOkButtonText("OK");
        //pippo.setCancelButtonText("No, thanks!");
        bool ok=pippo.exec();
        QString line;
        QString text =pippo.textValue();
        if (ok && !text.isEmpty()){
            line=text;
        }
        QFile file(fn);
        if(!file.open(QIODevice::WriteOnly | QIODevice::Text)){
            msgErrLoad("Save Osc Model",fn);
            qWarning()<<"IONK-> ERROR opening file= "<<fn;
            return;
        }
        {// (overwrite confirmation removed: it was dead code — iCase was always 0; getSaveFileName already confirms overwrite)
            QTextStream stream (&file);
            stream<<line<<"\n";
            for(int k=1;k<=20;k++){
                iok2=0;
                for(int i=2;i<=21;i++) {
                    if(nint(pf[ioptFit][i])==k) iok2=1;
                }
                if(iok2==1){
                    stream<<k<<"\t";//#Osc
                    for(int J=1;J<=5;J++)
                        stream<<pmOs[J+(k-1)*5][1]<<"\t";//Osc kind and parameter
                    stream<<"\n";
                }
            }
        }
        file.close();
    }
}


void ksemawc::loadNKmodel(){
    qDebug()<<"-> LoadNKmodel";
    int ioptFit=ui->cB_cnk1a -> currentIndex();
    if(ioptFit==0 || ioptFit>7){
        QMessageBox msgBox;
        msgBox.setText("Please set a valid Fit option");
        msgBox.setStandardButtons(QMessageBox::Ok);
        msgBox.exec();
        return;
    }
    QString fn = QFileDialog::getOpenFileName(
                this,
                "Choose a NK model", //window title
                pathroot,                 //initial directory
                "nkMod file (*.nkMod)"); //file extension
    if(fn.isEmpty())
        return;
    QFile file(fn);
    if (!file.open (QIODevice::ReadOnly | QIODevice::Text))
        return;
    setGuiMute(true);
    // Loading a model REPLACES the current one: clear every oscillator first, so that
    // oscillators absent from the file do not survive into the loaded model.
    for(int k=1;k<=20;k++){
        idToCheckBox["cBosc_"+QString::number(k)]-> setCheckState ( Qt::Unchecked );
        idToComboBox["cBpm_"+QString::number(100+1+(k-1)*5)+"_1"] -> setEnabled(false);
        idToLineEdit["LEpm_"+QString::number(100+2+(k-1)*5)+"_1"] -> setEnabled(false);
        idToLineEdit["LEpm_"+QString::number(100+3+(k-1)*5)+"_1"] -> setEnabled(false);
        idToLineEdit["LEpm_"+QString::number(100+4+(k-1)*5)+"_1"] -> setEnabled(false);
        idToLineEdit["LEpm_"+QString::number(100+5+(k-1)*5)+"_1"] -> setEnabled(false);
    }
    qDebug()<<"\tfile: "<<fn;
    QTextStream stream (&file);
    QString line;
    line=stream.readLine();
    qDebug()<<"\tcomment: "<<line;
    do{
        line=stream.readLine();
        line=line.simplified();
        if(line.isEmpty())
            continue;//tolerate blank lines (e.g. a trailing one)
        QStringList List;
        List =line.split(" ");
        int nV=List.count();
        if(nV!=6){
            qDebug("Error: expectet 6 values at line, but found %d",nV);
            break;//give up on the malformed line, but keep what was read and restore the GUI
        }
        int OscN=nint(List.at(0).toDouble());
        for(int J=1;J<=5;J++)
            pmOs[J+(OscN-1)*5][1]=List.at(J).toDouble();
        int ifu=nint(pmOs[1+(OscN-1)*5][1]);
        int k=OscN;
        idToCheckBox["cBosc_"+QString::number(k)]-> setCheckState ( Qt::Checked );
        idToComboBox["cBpm_"+QString::number(100+1+(k-1)*5)+"_1"] -> setEnabled(true);
        idToComboBox["cBpm_"+QString::number(100+1+(k-1)*5)+"_1"] -> setCurrentIndex(nint(pmOs[1+(k-1)*5][1])-1);
        //idToComboBox["cBpm_"+QString::number(100+1+(k-1)*5)+"_1"] -> setSizeAdjustPolicy(QComboBox::AdjustToContents);
        for(int iC=2;iC<=5;iC++){
            idToLineEdit["LEpm_"+QString::number(100+iC+(k-1)*5)+"_1"] -> setEnabled(true);
            idToLineEdit["LEpm_"+QString::number(100+iC+(k-1)*5)+"_1"] -> setText("");
        }
        if(ifu<4){//Lorentz, Q-homo, Q-inhomo
            idToLineEdit["LEpm_"+QString::number(100+5+(k-1)*5)+"_1"] -> setEnabled(false);
            idToLineEdit["LEpm_"+QString::number(100+5+(k-1)*5)+"_1"] -> setText("unused");
        }
        else if(ifu==4){//Flat
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
        else if(ifu==19){
            idToLineEdit["LEpm_"+QString::number(100+3+(k-1)*5)+"_1"] -> setEnabled(false);
            idToLineEdit["LEpm_"+QString::number(100+3+(k-1)*5)+"_1"] -> setText("unused");
        }
        for(int iC=2;iC<=5;iC++){
            QString content=idToLineEdit["LEpm_"+QString::number(100+iC+(k-1)*5)+"_1"] -> text();
            if(content!="unused")
               idToLineEdit["LEpm_"+QString::number(100+iC+(k-1)*5)+"_1"] -> setText(QString::number(std::abs(pmOs[iC+(k-1)*5][1])));
        }
    }while(!stream.atEnd());
    file.close();
    setGuiMute(false);
    // The checkboxes were ticked with the GUI muted, so setOscN() bailed out and
    // pf[ioptFit][] still lists the oscillators of the PREVIOUS model. Rebuild it from
    // the widgets, otherwise saveNKmodel() — which iterates over pf — writes back only
    // the oscillators that happened to be in the old list.
    setTabSim();
    listOsc();
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
    qDebug()<<"->saveNKsim with lastAction="<<lastAction;
    int ierr,iok=1;
    QString sample,fn2s;
    sample=ui->lineEdit_sample-> text();
    if(sample.contains("mate/aa999.9"))
        sample="";
    sample=sample+"_"+lastAction+".nk";
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
            break;
        }
    }
    if(iok==1){
        qDebug()<< "fNKsim= "<<fNKsim;
        qDebug()<< "fn2s= "<<fn2s;
        // QFile::copy (no shell) — safe with '+'/spaces/long names in fn2s.
        ierr = copyFileOverwrite(fNKsim, fn2s) ? 0 : 1;
        if(ierr != 0)
            qWarning()<<"Error copying NKsim!!!";
        else
            qDebug()<<"NKsim saved as "<<fn2s;
    }
}


void ksemawc::saveSim(){
    int ierr,iok=1;
    QString fn2s;
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
            break;
        }
    }
    if(iok==1){
        qDebug()<< "fMisSim= "<<fMisSim;
        qDebug()<< "fn2s= "<<fn2s;
        // QFile::copy (no shell) — safe with '+'/spaces/long names in fn2s.
        ierr = copyFileOverwrite(fMisSim, fn2s) ? 0 : 1;
        if(ierr != 0){
            qWarning()<<"Error copying MisSim.dat!!!";
            msgErrLoad("saveSim","Error copying MisSim.dat!!!");
        }
        else
            qDebug()<<"MisSim saved as "<<fn2s;
    }
}


void ksemawc::tabChanged(){
    QMessageBox msgBox;
    int itab=ui->tabWidget -> currentIndex();
    qDebug()<<"-> tabChanged:  itab="<<itab<<" lastTab="<<lastTab;
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
                    break;
                }
            }
        }
        readRangeEli();
        SPADA();//load measurements, file-nk and solution-nk
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
    qDebug()<<"-> setTabSim";
    int ioptFit=ui->cB_cnk1a -> currentIndex();
    if(ioptFit>0 && ioptFit<8){
        for(int i=1;i<=21;i++) pf[ioptFit][i]=0.;
        int j=2;
        Qt::CheckState state;
        for(int i=1;i<=20;i++){
            state=idToCheckBox["cBosc_"+QString::number(i)]-> checkState();
            if( state == Qt::Checked ) {
                pf[ioptFit][j]=i;
                pmOs[1+(i-1)*5][1]=idToComboBox["cBpm_"+QString::number(100+1+(i-1)*5)+"_1"] -> currentIndex();
                pmOs[1+(i-1)*5][1]++;
                for(int ii=2;ii<6;ii++){
                    QString content=idToLineEdit["LEpm_"+QString::number(100+ii+(i-1)*5)+"_1"] -> text();
                    if(content!="unused")
                        pmOs[ii+(i-1)*5][1]=content.toDouble();
                }
                if(pmOs[1+(i-1)*5][1]==15 || pmOs[1+(i-1)*5][1]==16){
                    if(pmOs[4+(i-1)*5][1]>pmOs[5+(i-1)*5][1]*0.45){
                        pmOs[4+(i-1)*5][1]=pmOs[5+(i-1)*5][1]*0.45;
                        idToLineEdit["LEpm_"+QString::number(100+4+(i-1)*5)+"_1"] -> setText(QString::number(pmOs[4+(i-1)*5][1]));
                    }
                }
                j++;
            }
        }
        pf[ioptFit][1]=j-2;
        //qDebug("\tNosc=%d",nint(pf[ioptFit][1]));
    }
}


void ksemawc::refreshFitPar(){
    int icoherent,klim,ivp,n_fresh,iosc,Jcombo,ip;
    setGuiMute(true);
    qDebug()<<"-> refreshFitPar";
    fflush(stdout);
    npp=ui->sB_PAR_34_5 -> value();
    n_fresh=0;
    QString Lj,Lparametro;
    //setting the parameter choice
    for(int j=1;j<=m_maxFitRow;j++){
        idToComboBox["cBParFit_"+QString::number(j)] -> clear();
        idToComboBox["cBParFit_"+QString::number(j)] -> addItem("none");
    }
    for(int i=1;i<=nlayer;i++){//layers
        icoherent=idToComboBox[wCB3(i)] -> currentIndex();
        Lj="L"+QString::number(i);
        if(icoherent==0)
            klim=-1;
        else if(icoherent==1)
            klim=1;
        else if(icoherent==2)
            klim=7;
        for(int j=1;j<=npp;j++){
            idToComboBox["cBParFit_"+QString::number(j)]-> addItem(Lj+ParFitLab[6]);//roughness is always included
            for(int k=0;k<=klim;k++){
                if(k!=6)
                    idToComboBox["cBParFit_"+QString::number(j)]-> addItem(Lj+ParFitLab[k]);
            }
        }
    }
    for(int i=1;i<=N_CNK_USER_MAX;i++){//EMA
        Qt::CheckState state=idToCheckBox["cB_EMA_"+QString::number(i)]->checkState();
        if(state==Qt::Checked){
            for(int j=1;j<=npp;j++)
                idToComboBox["cBParFit_"+QString::number(j)]-> addItem(ParFitLab[12]+QString::number(i));
        }
    }
    for(int i=1;i<=4;i++){//Thetas
        Qt::CheckState state=idToCheckBox["checkB_mis"+QString::number(5+2*i)+"_1"]->checkState();
        if(state==Qt::Checked){
            for(int j=1;j<=npp;j++)
                idToComboBox["cBParFit_"+QString::number(j)]-> addItem(ParFitLab[13]+QString::number(i));
        }
    }
    for(int i=1;i<=20;i++){//oscillator
        Lj="O"+QString::number(i);
        Qt::CheckState state=idToCheckBox["cBosc_"+QString::number(i)]-> checkState();
        if( state == Qt::Checked ) {
            int ifu=idToComboBox["cBpm_"+QString::number(100+1+(i-1)*5)+"_1"]->currentIndex();
            ifu++;
            for(int j=1;j<=npp;j++){
                idToComboBox["cBParFit_"+QString::number(j)]-> addItem(Lj+ParFitLab[8]);//"C"
                if(ifu!=4 && ifu!=5 && ifu!=17 && ifu!=19)
                    idToComboBox["cBParFit_"+QString::number(j)]-> addItem(Lj+ParFitLab[9]);//"E0"
                if(ifu!=4){
                    idToComboBox["cBParFit_"+QString::number(j)]-> addItem(Lj+ParFitLab[10]);//"D"
                }
                if(ifu>5 && ifu!=17)
                    idToComboBox["cBParFit_"+QString::number(j)]-> addItem(Lj+ParFitLab[11]);//"W"
            }
        }
    }

    //setting the parameter value
    for(int j=1;j<=m_maxFitRow;j++){
        if(j<=npp){
            //if(lastTab != itab){
            ivp=nint(pmAt(nint(ppm[j]))[2]);
            ip=ppm[j];
            if(ip >= 1 && ip <= 792){              // layer parameter
                int prop = (ip-1)/99;
                int n1   = (ip-1)%99 + 1;
                int n2   = (prop==7) ? 0 : prop+1;
                Lparametro = "L"+QString::number(n1)+ParFitLab[n2];
            }
            else if(ip >= 793 && ip <= 891){       // EMA mat 1..99
                Lparametro = ParFitLab[12]+QString::number(ip-792);
            }
            else if(ip >= 892 && ip <= 895){       // ThetaEli 1..4
                Lparametro = ParFitLab[13]+QString::number(ip-891);
            }
            else if(ip >= 897 && ip <= 996){       // oscillator params
                int param0 = ip-897;
                iosc = param0/5 + 1;
                Lparametro = "O"+QString::number(iosc)+ParFitLab[7+param0%5];
            }
            Jcombo=idToComboBox["cBParFit_"+QString::number(j)] -> findText(Lparametro);
            if(Jcombo>=0) idToComboBox["cBParFit_"+QString::number(j)] -> setCurrentIndex(Jcombo);
            idToLineEdit["DPparFitV_"+QString::number(j)] -> setText(QString::number(pmAt(ip)[1]));
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
                idToLineEdit["DPparFitErr_"+QString::number(j)] -> setText(QString::number(pmAt(ip)[4]));
                idToLineEdit["DPparFitGC_"+QString::number(j)] -> setText(QString::number(pmAt(ip)[5]));
            }
        }
        else{
            idToCheckBox["chBeParFit_"+QString::number(j)] -> setCheckState ( Qt::Unchecked );
            idToCheckBox["chBeParFit_"+QString::number(j)] -> setEnabled(false);
        }
    }
    nPar=n_fresh;
    ui->sB_PAR_35_5 -> setValue(nPar);
    setGuiMute(false);
    SaveSetting(-1);
}


void ksemawc::PanFitEnable(){
    int n_fresh,n1,n2,ip;
    int itab=ui->tabWidget -> currentIndex();

    if(isGuiMuted())
        return;
    qDebug()<<"-> PanFitEnable with itab="<<itab;
    fflush(stdout);
    setGuiMute(true);
    npp=ui->sB_PAR_34_5 -> value();
    n_fresh=0;
    Qt::CheckState state;
    //state = Qt::Checked;
    QString Lj,Valore;
    if(itab==3 && itab==lastTab){
        fflush(stdout);
        for(int j=1;j<=m_maxFitRow;j++){
            if(j<=npp){
                idToLineEdit["DPparFitV_"+QString::number(j)] -> setEnabled(true);
                idToCheckBox["chBeParFit_"+QString::number(j)] -> setEnabled(true);
                Lj=idToComboBox["cBParFit_"+QString::number(j)] -> currentText();
                if(Lj=="none")// unassigned row: cannot be a fit parameter (mirror PanFitChoice) — must uncheck before n_fresh
                    idToCheckBox["chBeParFit_"+QString::number(j)] -> setCheckState( Qt::Unchecked );
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
                ip=0;
                if(Lj=="none"){
                    ppm[j]=0;// nothing to write: leave parameter arrays untouched (avoids corrupting ThetaEli via bogus ip=893)
                }
                else{
                    if(Lj.at(0)=='L'){
                        n1=(Lj.at(1)).digitValue();
                        if(Lj.size()>2 && Lj.at(2).isDigit()){
                            n1=n1*10+Lj.at(2).digitValue();
                        }
                        n2=-1;
                        do{
                            n2++;
                        } while(!Lj.contains(ParFitLab[n2],Qt::CaseSensitive) && n2<7);
                        ip=ipLayer(n1,n2);
                    }
                    else if(Lj.at(0)=='f'){
                        ip = 792 + Lj.mid(5).toInt();  // "fEMA_" = 5 chars
                    }
                    else if(Lj.at(0)=='T')
                        ip = 891 + Lj.at(6).digitValue();  // "Theta_" = 6 chars
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
                        ip = 896 + (n1-1)*5 + (n2-7) + 1;
                    }
                    ppm[j]=ip;
                    Valore=idToLineEdit["DPparFitV_"+QString::number(j)] -> text();
                    pmAt(ip)[1]=Valore.toDouble();
                    if(ip>=897 && ip<=996)//oscillator params (uniform ip encoding v5.2)
                        pmAt(ip)[1]=std::abs(pmAt(ip)[1]);
                    if( state == Qt::Unchecked )
                        pmAt(ip)[2]=0;
                    else{
                        pmAt(ip)[2]=j;
                        pmAt(n_fresh)[3]=ip;
                    }
                }
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
    }
    setGuiMute(false);
}


void ksemawc::PanFitPar(){
    if(isGuiMuted())
        return;
    qDebug()<<"-> PanFitPar";
    int icoherent,klim,nppNew;
    int itab=ui->tabWidget -> currentIndex();
    setGuiMute(true);
    nppNew=ui->sB_PAR_34_5 -> value();
    Qt::CheckState state;
    //state = Qt::Checked;
    QString Lj;
    if(itab==3 && itab==lastTab){
        ensureFitRow(nppNew);   // create widget rows if nppNew > m_maxFitRow
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
                icoherent=idToComboBox[wCB3(i)] -> currentIndex();
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
                        if(ifu!=4 && ifu!=5 && ifu!=17 && ifu!=19)
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
        npp=nppNew;
        PanFitEnable();
    }
    setGuiMute(false);
}


void ksemawc::PanFitChoice(){
    if(isGuiMuted())
        return;
    int n1,n2,ip;
    int itab=ui->tabWidget -> currentIndex();
    QString Lj;
    Qt::CheckState state;
    qDebug()<<"-> PanFitChoice";
    setGuiMute(true);
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
                if(n2==0 || (1<=n2 && n2<=7))
                    ip=ipLayer(n1,n2);
                else if(8<=n2 && n2<=11)
                    ip = 896 + (n1-1)*5 + (n2-7) + 1;
                else if(n2==12)
                    ip = 792 + Lj.mid(5).toInt();
                else if(n2==13)
                    ip = 891 + Lj.at(6).digitValue();
                if(1<=ip && ip<=996){
                    idToLineEdit["DPparFitV_"+QString::number(j)] -> setText(QString::number(pmAt(ip)[1]));
                    state=idToCheckBox["chBeParFit_"+QString::number(j)]-> checkState();
                    if( state == Qt::Checked ){
                        idToLineEdit["DPparFitErr_"+QString::number(j)] -> setText(QString::number(pmAt(ip)[4]));
                        idToLineEdit["DPparFitGC_"+QString::number(j)] -> setText(QString::number(pmAt(ip)[5]));
                    }
                    else{
                        idToLineEdit["DPparFitErr_"+QString::number(j)] -> clear();
                        idToLineEdit["DPparFitGC_"+QString::number(j)]  -> clear();
                    }
                }
            }
        }
    }
    setGuiMute(false);
}


void ksemawc::PlotMENK(){
    qDebug()<<"-> PlotMENK()";
    iColor=0;//reset to black color
    //SPADA();
    PlotME();
    PlotNK(1);
}


void ksemawc::PlotME(){
    int Ndata=NeV;
    int iRD,iPSI=0,iDELTA=0;
    std::vector <double> Xp(Ndata),Yp(Ndata),ErrXp(Ndata),ErrYp(Ndata);
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
            qDebug()<<"->PlotMe: DATO["<<i<<"]="<<DATO[i]<<" ifirstWarning="<<ifirstWarning;
        if(DATO[i]>0){
            double X,Y,errY;
            for(int L=1;L<=Ndata;L++){
                X=ms.lambda[L-1];
                if(i<=6){
                    Y=ms.measures[i-1].value[L-1]*100.;
                    errY=ms.measures[i-1].error[L-1]*100.;
                }
                else if(i>=7){
                    Y=ms.measures[i-1].value[L-1];
                    errY=ms.measures[i-1].error[L-1];
                }
                Xp[L-1]=X;
                Yp[L-1]=Y;
                ErrYp[L-1]=errY;
            }
            iColor=0;
            PLOTline1bar2(2,iRD,0,iGraph,Ndata,Xp.data(),Yp.data(),ErrXp.data(),ErrYp.data());
        }
        else if(DATO[i]!=0)
            PLOTline1bar2(2,iRD,0,iGraph,0,Xp.data(),Yp.data(),ErrXp.data(),ErrYp.data());//erase without data plot
    }

    iRD=1;
    if(DATO[1]>0 && DATO[3]>0){
        int iPlotA=0;
        rxy[18][3]=1000.;
        rxy[18][4]=-1000.;
        for(int L=1;L<=Ndata;L++){
            Xp[L-1]=ms.lambda[L-1];
            Yp[L-1]=(1.-ms.measures[0].value[L-1]-ms.measures[2].value[L-1])*100.;
            ErrYp[L-1]=(ms.measures[0].error[L-1]+ms.measures[2].error[L-1])*100.;
            if(Yp[L-1]<0.) iPlotA++;
            rxy[18][3]=min(rxy[18][3],Yp[L-1]);
            rxy[18][4]=max(rxy[18][4],Yp[L-1]);
        }
        ui->DP_RXY_18_3 -> setText(QString::number(rxy[18][3]));
        ui->DP_RXY_18_4 -> setText(QString::number(rxy[18][4]));
        labelQwt="A_front";
        PLOTline1bar2(2,iRD,0,9,Ndata,Xp.data(),Yp.data(),ErrXp.data(),ErrYp.data());
        if(ifirstWarning==0 && iPlotA>0){
            QMessageBox msgBox;
            msgBox.setText("ATTENTION: A = 1 - Tn - Rn < 0 !!!\nPlease check your measurements.");
            msgBox.setStandardButtons(QMessageBox::Ok);
            msgBox.exec();
            ifirstWarning=-1;
        }
    }
    else if(DATO[1]!=0 && DATO[3]!=0)
        PLOTline1bar2(2,iRD,0,9,0,Xp.data(),Yp.data(),ErrXp.data(),ErrYp.data());

    if(DATO[1]>0 && DATO[5]>0){
        int iPlotA=0;
        rxy[15][3]=1000.;
        rxy[15][4]=-1000.;
        for(int L=1;L<=Ndata;L++){
            Xp[L-1]=ms.lambda[L-1];
            Yp[L-1]=(1.-ms.measures[0].value[L-1]-ms.measures[4].value[L-1])*100.;
            ErrYp[L-1]=(ms.measures[0].error[L-1]+ms.measures[4].error[L-1])*100.;
            if(Yp[L-1]<0.) iPlotA++;
            rxy[15][3]=min(rxy[15][3],Yp[L-1]);
            rxy[15][4]=max(rxy[15][4],Yp[L-1]);
        }
        ui->DP_RXY_15_3 -> setText(QString::number(rxy[15][3]));
        ui->DP_RXY_15_4 -> setText(QString::number(rxy[15][4]));
        labelQwt="A_back";
        PLOTline1bar2(2,iRD,0,16,Ndata,Xp.data(),Yp.data(),ErrXp.data(),ErrYp.data());
        if(ifirstWarning<=0 && iPlotA>0){
            QMessageBox msgBox;
            msgBox.setText("ATTENTION: A= 1 - Tn - R1 <0 !!!\nPlease check your measurements.");
            msgBox.setStandardButtons(QMessageBox::Ok);
            msgBox.exec();
        }
        ifirstWarning=1;
    }
    else if(DATO[1]!=0 && DATO[5]!=0)
        PLOTline1bar2(2,iRD,0,16,0,Xp.data(),Yp.data(),ErrXp.data(),ErrYp.data());


    // save to HD re-sampled experimental measurements and loaded NK
    QFile file0(pathroot+"expo/MeasExp.dat");
    if (!file0.open(QIODevice::WriteOnly | QIODevice::Text)){
        msgErrLoad("Save file0","Error opening expo/MisSFexp.dat");
    }
    QTextStream stream0 (&file0);

    QFile file1(pathroot+"expo/ErrMeas.dat");
    if (!file1.open(QIODevice::WriteOnly | QIODevice::Text)){
        msgErrLoad("Save file1","Error opening expo/ErrMeas.dat");
    }
    QTextStream stream1 (&file1);

    QFile file2(pathroot+"expo/loadedNK.dat");
    if (!file2.open(QIODevice::WriteOnly | QIODevice::Text)){
        msgErrLoad("Save file2","Error opening expo/loadedNK.dat");
    }
    QTextStream stream2 (&file2);

    stream0<<"wl(nm)\tTn\tTp\tRn\tRp\tR1\tAPDS\tDEL_1\tPSI_1\tDEL_2\tPSI_2\tDEL_3\tPSI_3\tDEL_4\tPSI_4"<<"\n";
    stream1<<"wl(nm)\teTn\teTp\teRn\teRp\teR1\teAPDS\teDEL_1\tePSI_1\teDEL_2\tePSI_2\teDEL_3\tePSI_3\teDEL_4\tePSI_4"<<"\n";
    // Export every loaded nk material (no fixed 8-pair cap: a project may use up to N_CNK_USER_MAX).
    int nNkExport=0;
    for(int j=0;j<N_CNK_USER_MAX;j++)
        if(nkMaterials[j].isLoaded()) nNkExport=j+1;//highest loaded slot; keeps column j <-> material j+1
    stream2<<"wl(nm)";
    for(int j=1;j<=nNkExport;j++) stream2<<"\tn"<<j<<"\tk"<<j;
    stream2<<"\n";
    for(int i=0;i<NeV;i++){
        stream0<<QString::number(ms.lambda[i])<<"\t";
        for(int j=0;j<14;j++) stream0<<QString::number(ms.measures[j].value[i],'f',7)<<"\t";
        stream1<<QString::number(ms.lambda[i])<<"\t";
        for(int j=0;j<14;j++) stream1<<QString::number(ms.measures[j].error[i],'f',7)<<"\t";
        stream2<<QString::number(ms.lambda[i])<<"\t";
        for(int j=0;j<nNkExport;j++) stream2<<QString::number(nkMaterials[j].n[i],'f',7)<<"\t"<<QString::number(nkMaterials[j].k[i],'f',7)<<"\t";
        stream0<<"\n";
        stream1<<"\n";
        stream2<<"\n";
    }
    file0.close();
    file1.close();
    file2.close();
}


void ksemawc::PlotNK(int iRD){
    int Ndata=Sol.size();
    //if(Ndata==0)
    //    return;
    std::vector <double> Xp(Ndata),Yp(Ndata),ErrXp(Ndata),ErrYp(Ndata);
    // n-data plot
    for(int i=0;i<Ndata;i++){
        Xp[i]=Sol[i][1];
        Yp[i]=Sol[i][2];
        ErrYp[i]=Sol[i][4];
    }
    PLOTline1bar2(2,iRD,iColor,12,Ndata,Xp.data(),Yp.data(),ErrXp.data(),ErrYp.data());
    // k-data plot
    for(int i=0;i<Ndata;i++){
        Xp[i]=Sol[i][1];
        Yp[i]=Sol[i][3];
        ErrYp[i]=Sol[i][5];
    }
    PLOTline1bar2(2,iRD,iColor,13,Ndata,Xp.data(),Yp.data(),ErrXp.data(),ErrYp.data());
    //eps1 eps2 plot
    if(nint(par[10][1])==1){
        //eps1-data plot
        double v1,v2,v3,v4,yd,yu;
        rxy[26][3]=1000.;
        rxy[26][4]=-1000.;
        for(int i=0;i<Ndata;i++){
            v1=pow(Sol[i][2]-Sol[i][4],2.)-pow(Sol[i][3]-Sol[i][5],2.);
            v2=pow(Sol[i][2]-Sol[i][4],2.)-pow(Sol[i][3]+Sol[i][5],2.);
            v3=pow(Sol[i][2]+Sol[i][4],2.)-pow(Sol[i][3]-Sol[i][5],2.);
            v4=pow(Sol[i][2]+Sol[i][4],2.)-pow(Sol[i][3]+Sol[i][5],2.);
            yd=min(v1,v2);
            yd=min(yd,v3);
            yd=min(yd,v4);
            yu=max(v1,v2);
            yu=max(yu,v3);
            yu=max(yu,v4);
            rxy[26][3]=min(rxy[26][3],yd);
            rxy[26][4]=max(rxy[26][4],yu);
            Xp[i]=Sol[i][1];
            Yp[i]=(yu+yd)/2.;
            ErrYp[i]=(yu-yd)/2.;
        }
        PLOTline1bar2(2,iRD,iColor,14,Ndata,Xp.data(),Yp.data(),ErrXp.data(),ErrYp.data());
        ui->DP_RXY_26_3->setText(QString::number(rxy[26][3]));
        ui->DP_RXY_26_4->setText(QString::number(rxy[26][4]));
        //eps2-data plot
        rxy[27][3]=1000.;
        rxy[27][4]=-1000.;
        for(int i=0;i<Ndata;i++){
            v1=2.*(Sol[i][2]-Sol[i][4])*(Sol[i][3]-Sol[i][5]);
            v2=2.*(Sol[i][2]-Sol[i][4])*(Sol[i][3]+Sol[i][5]);
            v3=2.*(Sol[i][2]+Sol[i][4])*(Sol[i][3]-Sol[i][5]);
            v4=2.*(Sol[i][2]+Sol[i][4])*(Sol[i][3]+Sol[i][5]);
            yd=min(v1,v2);
            yd=min(yd,v3);
            yd=min(yd,v4);
            yu=max(v1,v2);
            yu=max(yu,v3);
            yu=max(yu,v4);
            rxy[27][3]=min(rxy[27][3],yd);
            rxy[27][4]=max(rxy[27][4],yu);
            Xp[i]=Sol[i][1];
            Yp[i]=(yu+yd)/2.;
            ErrYp[i]=(yu-yd)/2.;
        }
        PLOTline1bar2(2,iRD,iColor,15,Ndata,Xp.data(),Yp.data(),ErrXp.data(),ErrYp.data());
        ui->DP_RXY_27_3->setText(QString::number(rxy[27][3]));
        ui->DP_RXY_27_4->setText(QString::number(rxy[27][4]));
    }
}


void ksemawc::Simula(){
    qDebug()<<"->Simula";
    lastAction="Simula";
    AdjRoughMax();
    SaveSetting(-1);
    std::vector<std::vector<double>> mc(15, std::vector<double>(NeV+1, 0.0));
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
        double INSTR[21];
        for(int j=1;j<=20;j++)
            INSTR[j]=par[j][3];
        for(int ic=1;ic<=14;ic++){
            if(DATO[ic]!=0){
                DATO[ic]=1;
                par[ic][4]=1.;
                idToCheckBox["checkB_mis"+QString::number(ic)+"_3"] -> setEnabled(true);
                for(int i=1;i<=NeV;i++){
                    ms.measures[ic-1].value[i-1]=mc[ic][i];//value/error are 0-based (size NeV); mc is 1-based (size NeV+1)
                    if(ic<=6){
                        double DL=ms.driftBaseline[i-1];
                        if(ic==1 || ic==2)
                            ms.measures[ic-1].error[i-1]=ms.measures[ic-1].value[i-1]*INSTR[3]+DL;
                        if(ic==3 || ic==4 || ic==5)
                            ms.measures[ic-1].error[i-1]=ms.measures[ic-1].value[i-1]*(INSTR[3]+INSTR[4])+DL;
                    }
                    else if(ic>=7 && ic<=14)
                        ms.measures[ic-1].error[i-1]=0.01;
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
        return;
    }
    else{
        int ipr;
        QString line,pezzo;
        QTextStream stream (&file);
        line = stream.readLine();//info line
        qDebug() << line;
        line = stream.readLine();//Ndat & "nm" for new file
        line=line.simplified();
        QStringList List0;
        List0 =line.split(" ");
        pezzo=List0.at(0).toLocal8Bit().constData();
        ipr=pezzo.toInt();
        double a2nm=0.1;
        if(line.contains("nm"))
            a2nm=1.;
        qDebug("Ndat= %d a2nm=%f",ipr,a2nm);
        std::vector<std::vector<double>> wp(ipr+1, std::vector<double>(3, 0.0));
        for(int i=1; i<=ipr; i++){
            double tmp1, tmp2;
            stream >> tmp1 >> tmp2;
            wp[i][1] = tmp1 * a2nm;
            wp[i][2] = tmp2;
            //qDebug()<<i<<" "<<wp[i][1]<<" "<<wp[i][2];
        }
        file.close();
        //weight
        double speso=0.;
        int jjmin=1;
        while(wp[jjmin][1] < ms.lambda[0] && jjmin<ipr)
            jjmin++;
        int jj=jjmin;
        int jjmax=jj;
        while(jj<=ipr && wp[jj][1] >= ms.lambda[0] && wp[jj][1] <= ms.lambda[NeV-1] ){
            speso=speso+wp[jj][2];
            jjmax=jj;
            jj++;
            //qDebug("while jj=%d",jj);
        }
        // reset adders
        double av[7],avs[7],val[7][3];
        for( int i=1;i<=6;i++){
            av[i]=.0;
            avs[i]=.0;
        }
        // mean @theta
        for(int ii=jjmin;ii<=jjmax;ii++){
            int ij=1;
            while(wp[ii][1]>ms.lambda[ij-1] && ij<=nint(par[28][2]-1))
                ij++;
            if(ij>1) ij--;
            double x=(wp[ii][1]-ms.lambda[ij-1])/(ms.lambda[ij]-ms.lambda[ij-1]);
            for(int j=1;j<=2;j++){
                int ik=ij+j-1;
                if(j==1){
                    for(int iv1=1;iv1<=6;iv1++){
                        val[iv1][1]=ms.measures[iv1-1].value[ik-1]*(1.-x);
                        val[iv1][2]=mc[iv1][ik]*(1.-x);
                    }
                }
                else{
                    for(int iv1=1;iv1<=6;iv1++){
                        val[iv1][1]=val[iv1][1]+ms.measures[iv1-1].value[ik-1]*x;
                        val[iv1][2]=val[iv1][2]+mc[iv1][ik]*x;
                    }
                }
            }
            for(int j=1;j<=6;j++){
                av[j]=av[j]+val[j][1]*wp[ii][2]/speso;
                avs[j]=avs[j]+val[j][2]*wp[ii][2]/speso;
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


void ksemawc::CalcMis(std::vector<std::vector<double>>& mc){
    qDebug("->CalcMis with cnk[1].forceMode=%d",nint(cnk[1].forceMode));
    nextColor();
    int s1p2u3=ui->cB_PAR_35_2 ->currentIndex();
    s1p2u3++;
    Qt::CheckState state;
    state=ui->checkB_PAR_54_2 -> checkState ();
    if(state==Qt::Unchecked)
        par[54][2]=0;
    else
        par[54][2]=1;
    int iSpec0Hemi1=par[54][2];
    int nwl=NeV;
    double teta[6],wl,VNK[N_CNK_HARD_MAX+1][3],vot[6][3],delOld[4],del=0.;
    std::vector <double> Xp(nwl),Yp(nwl),ErrXp(nwl),ErrYp(nwl),nn(nwl),kk(nwl),e1(nwl),e2(nwl);
    fill_n(ErrXp.data(), nwl, 0.);
    fill_n(ErrYp.data(), nwl, 0.);
    teta[1]=ui->dSB_PAR_6_1->value();
    teta[2]=ui->dSB_PM_86_1->value();
    teta[3]=ui->dSB_PM_87_1->value();
    teta[4]=ui->dSB_PM_88_1->value();
    teta[5]=ui->dSB_PM_89_1->value();
    state=ui->checkBox_deltaConnect->checkState();
    int iDelCon=0;
    if(state==Qt::Checked)
        iDelCon=1;
    for(int k=1;k<6;k++){
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
        wl=ms.lambda[i-1];
        Xp[i-1]=wl;
        COSVNK(VNK,i);
        nn[i-1]=VNK[1][1];
        kk[i-1]=VNK[1][2];
        e1[i-1]=VNK[1][1]*VNK[1][1]-VNK[1][2]*VNK[1][2];
        e2[i-1]=2.*VNK[1][1]*VNK[1][2];
        nkSol.n[i-1]=VNK[1][1];
        nkSol.k[i-1]=VNK[1][2];
        streamNK<<wl<<"\t"<<VNK[1][1]<<"\t"<<VNK[1][2]<<"\t"<<VNK[1][1]*0.001<<"\t"<<VNK[1][2]*0.001<<"\n";
        if(DATO[1]!=0 || DATO[3]!=0 || DATO[5]!=0){
            if(iSpec0Hemi1==1)
                iShemispherical=true;
            else
                iShemispherical=false;
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
            par[54][2]=0;//Tp and Rp are always direct/specular
            iShemispherical=false;
            ASSEMBLER(i,1,teta[1],vot);
            par[54][2]=iSpec0Hemi1;
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
            if(iSpec0Hemi1==1)
                iShemispherical=true;
            else
                iShemispherical=false;
            ASSEMBLER(i,2,0.,vot);
            mc[6][i]=vot[4][1];
        }
        else
            mc[6][i]=.0;
        for(int j=1;j<=4;j++){
            iShemispherical=false;
            if(DATO[7+2*(j-1)]!=0 || DATO[8+2*(j-1)]!=0){
                ASSEMBLER(i,2+j,teta[j+1],vot);
                del=vot[5][2];//Delta
                if(iDelCon==1){
                    if(i==1){
                        delOld[j-1]=del;
                    }
                    del=del-Dperiod*nint((del-delOld[j-1])/Dperiod);
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
                    Diff=Diff+(mc[ic][i]-ms.measures[ic-6+5].value[i-1])/NeV;
                }
                Nperiod=nint(Diff/Dperiod);
            }
            for(int i=1;i<=NeV;i++){
                Yp[i-1]=mc[ic][i];
                if(ic<=5)
                    Yp[i-1]=Yp[i-1]*100.;
                if(ic<=6)
                    FM=FM+pow((mc[ic][i]-ms.measures[ic-1].value[i-1])/ms.measures[ic-1].error[i-1],2.)/NeV;
                else{
                    FM=FM+pow((mc[ic][i]-ms.measures[ic-1].value[i-1]-Nperiod*Dperiod)/ms.measures[ic-1].error[i-1],2.)/NeV;
                    Yp[i-1]=Yp[i-1]-Nperiod*Dperiod;
                }
            }
            par[35+ic][1]=sqrt(FM);
            if(ic<=6)
                PLOTline1bar2(1,0,iColor,ic,NeV,Xp.data(),Yp.data(),ErrXp.data(),ErrYp.data());
            else if(ic==7 || ic==9 || ic==11 || ic==13)
                PLOTline1bar2(1,0,iColor,7,NeV,Xp.data(),Yp.data(),ErrXp.data(),ErrYp.data());
            else if(ic==8 || ic==10 || ic==12 || ic==14)
                PLOTline1bar2(1,0,iColor,8,NeV,Xp.data(),Yp.data(),ErrXp.data(),ErrYp.data());
            SIMchi2=SIMchi2+FM;
        }
    }
    ui->lineEdit_SIMchi2->setText("chi2= "+QString::number(nMIS>0. ? SIMchi2/nMIS : 0.));//guard nMIS==0 (no enabled measurement -> NaN)
    if(DATO[1]!=0 && DATO[3]!=0){
        for(int i=1;i<=NeV;i++){
            Yp[i-1]=(1.-mc[1][i]-mc[3][i])*100.;
        }
        labelQwt="A_front";
        PLOTline1bar2(1,0,iColor,9,NeV,Xp.data(),Yp.data(),ErrXp.data(),ErrYp.data());
    }
    if(DATO[1]!=0 && DATO[5]!=0){
        for(int i=1;i<=NeV;i++){
            Yp[i-1]=(1.-mc[1][i]-mc[5][i])*100.;
        }
        labelQwt="A_back";
        PLOTline1bar2(1,0,iColor,16,NeV,Xp.data(),Yp.data(),ErrXp.data(),ErrYp.data());
    }
    PLOTline1bar2(1,0,iColor,12,NeV,Xp.data(),nn.data(),ErrXp.data(),ErrYp.data());
    PLOTline1bar2(1,0,iColor,13,NeV,Xp.data(),kk.data(),ErrXp.data(),ErrYp.data());
    if(nint(par[10][1])==1){//plot permittivity
        PLOTline1bar2(1,0,iColor,14,NeV,Xp.data(),e1.data(),ErrXp.data(),ErrYp.data());
        PLOTline1bar2(1,0,iColor,15,NeV,Xp.data(),e2.data(),ErrXp.data(),ErrYp.data());
    }
    QFile file(fMisSim);
    if (!file.open(QIODevice::WriteOnly | QIODevice::Text)){
        msgErrLoad("CalcMis-fMisSim",fMisSim);
        return;
    }
    QTextStream stream(&file);
    stream<<"wl\tTn\tTp\tRn\tRp\tR1\tApds\tDELTA\tPSI\n";
    for(int i=1;i<=NeV;i++){
        stream<<ms.lambda[i-1];
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
    fill_n(ErrYp, 92, 0.);
    SaveSetting(-1);
    int s1p2u3=ui->cB_PAR_35_2 ->currentIndex();
    s1p2u3++;
    int iSpec0Hemi1=nint(par[54][2]);
    // kind of average
    QString fname=pathroot+NANK[10].simplified();
    QFile file(fname);
    if (!file.open(QIODevice::ReadOnly | QIODevice::Text)){
        msgErrLoad("PlotAve fname",fname);
        return;
    }
    QTextStream stream(&file);
    line=stream.readLine();//info line
    qDebug()<<line;
    line=stream.readLine();
    line=line.simplified();
    QStringList List;
    List=line.split(" ");
    pezzo=List.at(0).toLocal8Bit().constData();
    ipr=pezzo.toInt();
    double wl2nm=1.;
    if(line.contains("nm"))
        wl2nm=0.1;
    std::vector<std::vector<double>> wp(ipr+1, std::vector<double>(3, 0.0));
    double speso=0.;
    for(int i=1;i<=ipr;i++){
        double tmp1, tmp2;
        stream >> tmp1 >> tmp2;
        wp[i][1]=tmp1;
        wp[i][2]=tmp2;
        wp[i][1]=wl2nm*wp[i][1];
        qDebug()<< wp[i][1] <<"\t"<< wp[i][2];
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
    if(iSpec0Hemi1==1)
        iShemispherical=true;
    else
        iShemispherical=false;
    for(int I=1;I<=NT;I++){
        teta=rxy[21][1]+(rxy[21][2]-rxy[21][1])*(I-1)/static_cast<double>(NT-1);
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
            if(wp[ii][1]<=ms.lambda[0]){
                ij=1;
                x=0.;
                if(wp[ii][1]<ms.lambda[0] && I==1)
                    qWarning()<<"Attention wl_mean="<<wp[ii][1]<<" <WLmin!!!";
            }
            else if(wp[ii][1]>ms.lambda[0] && wp[ii][1]<ms.lambda[NeV-1]){
                ij=1;
                while(wp[ii][1]>ms.lambda[ij-1])
                    ij++;
                if(ij>1)
                    ij--;
                x=(wp[ii][1]-ms.lambda[ij-1])/(ms.lambda[ij]-ms.lambda[ij-1]);
            }
            else if(wp[ii][1]>=ms.lambda[NeV-1]){
                ij=NeV-1;
                x=1.;
                if(wp[ii][1]>ms.lambda[NeV-1] && I==1)
                    qWarning()<<"Attention wl_mean="<<wp[ii][1]<<" >WLmax!!!!";
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
        ui->DP_RXY_23_1 -> setText(QString::number(static_cast<int>(ymin*100.)));
        ui->DP_RXY_23_2 -> setText(QString::number(nint(ymax*100.)));
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
        return;
    }
    QTextStream out(&file2);
    out<<" Teta(deg)  Tau  Rho  Rho1"<<"\n";
    out<<NT;
    for(int i=0;i<NT;i++)//0-based: tet/tau/... are filled tet[I-1] for I=1..NT (indices 0..NT-1)
        out<< tet[i]<<"\t"<<tau[i]<<"\t"<<rho[i]<<"\t"<<rho1[i]<<"\n";
    file2.close();
    qInfo()<<"Data vs theta saved in: expo/ThetaRhoVStheta.dat";
}


void ksemawc::PlotAbsEL(){
    std::vector <int> LastMediumAtWl(NeV+1);
    double vosi[6][3],denom;
    std::vector <double> wwl(NeV+1),x(NeV+1),y(NeV+1),ErrXp(NeV+1),ErrYp(NeV+1),Trasm1(NeV+1);
    std::vector<std::vector<double>> ABS(NeV+1, std::vector<double>(10, 0.0));
    fill_n(ErrXp.data(), NeV+1, 0.);
    fill_n(ErrYp.data(), NeV+1, 0.);
    fill_n(Trasm1.data(), NeV+1, 1.);
    std::vector<std::vector<long double>> Ps(NeV+1, std::vector<long double>(11, 0.0));
    std::vector<std::vector<long double>> Pp(NeV+1, std::vector<long double>(11, 0.0));
    std::vector<std::vector<std::complex<long double>>> Bs(NeV+1, std::vector<std::complex<long double>>(11, 0.0));
    std::vector<std::vector<std::complex<long double>>> Cs(NeV+1, std::vector<std::complex<long double>>(11, 0.0));
    std::vector<std::vector<std::complex<long double>>> Bp(NeV+1, std::vector<std::complex<long double>>(11, 0.0));
    std::vector<std::vector<std::complex<long double>>> Cp(NeV+1, std::vector<std::complex<long double>>(11, 0.0));
    complex<double> irup,irdw,irup2,irdw2,mups,mupp,mdws,mdwp,pq,ir[10], NQ, MUS, DELTA;
    int icol,ivnkdw,ivnkup,s1p2u3, LastMedium,ivnk;
    double VNK[N_CNK_HARD_MAX+1][3];
    SaveSetting(-1);

    s1p2u3=nint(par[35][2]);
    double teta=par[6][1]*deg2rad;
    int iFirstCoe=1;
    int Nlayers=nint(par[51][2]);
    ivnkdw=CNK_OUT_SF;//output medium for SF measurements
    ivnkup=CNK_IN_SF; //refractive index of input medium for SF measurements
    int iRD=1;//redraw the Absorptance plot

    qDebug() << "plot Abs at each layer of " << QString::number(nint(par[51][2]))<<" at theta="<<par[6][1]<<" with polarization "<<s1p2u3<<" Nlayers="<<Nlayers<<" ivnkdw="<<ivnkdw;

    //set wwl array
    for(int i=1;i<=NeV;i++){
        wwl[i]=ms.lambda[i-1];
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
            if(nint(par[50+j][3])==1 && j>1)
                iWorn=1;
            ivnk=static_cast<int>(par[50+j][1]);//index of the material to be assigned to the ii layer
            ir[j]=complex<double>(VNK[ivnk][1],-VNK[ivnk][2]);
            NQ=ir[j]*ir[j];
            MUS=sqrt(NQ-pq);
            DELTA=2.*PIG*pmAt(j)[1]*MUS/wwl[i];
            if(std::abs(imag(DELTA))>ABSmax) {
                //if(i==iwl2print)
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
        if(i==iwl2print) qDebug()<<"LastMediumAtWl["<<i<<"]="<<LastMediumAtWl[i];
    }

    //computing at each wavelength
    for(int i=1;i<=NeV;i++){
        //refractive indices @ wl
        COSVNK(VNK,i);

        //set last medium
        ivnkdw=CNK_OUT_SF;
        LastMedium=LastMediumAtWl[i];
        if(LastMedium<Nlayers+1){
            ivnkdw=nint(par[50+LastMedium][1]);
            for(int m=LastMedium+1;m<=Nlayers+1;m++){
                Ps[i][m]=0.;
                Pp[i][m]=0.;
            }
        }
        //loop on the layers
        for(int iLayer=1;iLayer<=LastMedium;iLayer++){
            if(i==iwl2print){
                qDebug()<<"\n|PlotAbsEL@wl="<<ms.lambda[iwl2print-1]<<"> for-loop with iLayer="<<iLayer<<"<="<<LastMedium;
            }
            if(iLayer==1 && (nint(par[51][3])==1 || LastMedium==1)){//the first layer is inchoerent or completely absorbing
                irup=complex<double>(VNK[ivnkup][1],-VNK[ivnkup][2]);
                pq=pow(irup*sin(teta),2.);// (n_input*sin(teta))**2.
                if(i==iwl2print)
                    qDebug()<<"\t Incoherent ivnkup="<<ivnkup<<" nup="<<VNK[ivnkup][1]<<" kup="<<-VNK[ivnkup][2];
                ivnkdw=nint(par[51][1]);//refractive index of first layer
                irup=complex<double>(VNK[ivnkup][1],-VNK[ivnkup][2]);
                irup2=irup*irup;
                pq=pow(irup*sin(teta),2.);// (n_input*sin(teta))**2.
                irdw=complex<double>(VNK[ivnkdw][1],-VNK[ivnkdw][2]);
                irdw2=irdw*irdw;
                mups=sqrt(irup2-pq);
                mupp=irup2/mups;
                mdws=sqrt(irdw2-pq);
                mdwp=irdw2/mdws;
                double atten=exp(-4.*PIG*pmD[1][1]/wwl[i]*imag(-mdws));
                double RaS=pow(std::abs((mups-mdws)/(mups+mdws)),2.);
                double RaP=pow(std::abs((mupp-mdwp)/(mupp+mdwp)),2.);
                //double TaS=1.-RaS;//T=1-R because A=0 at the interface
                //double TaP=1.-RaP;
                //double R1aS=RaS;
                //double R1aP=RaP;
                if(i==iwl2print){
                    qDebug()<<"The first layer is inchoerent with first interface:";
                    qDebug()<<"nUP="<<VNK[ivnkup][1]<<" kUP="<<VNK[ivnkup][2]<<" nDW="<<VNK[ivnkdw][1]<<" kDW="<<VNK[ivnkdw][2]<<" d="<<pmD[1][1];
                    qDebug()<<"atten="<<atten<<" RaS="<<RaS<<" RaP="<<RaP;
                    qDebug()<<"2nd interface: call BUILDER: iFirstLayer=2 ncoe="<<Nlayers-1;
                }
                BUILDER(i,1,2,Nlayers-1,pq,vosi);
                double RbS=vosi[2][1];
                double RbP=vosi[2][2];
                if(i==iwl2print) qDebug()<<"RbS="<<RbS<<" RbP="<<RbP;
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
            }
            else{//thin layer
                if(iLayer<LastMedium){
                    BUILDER(i,1,iLayer,LastMedium-iLayer,pq,vosi);
                    Bs[i][iLayer]=freCoeff[4][0];
                    Cs[i][iLayer]=freCoeff[4][1];
                    Bp[i][iLayer]=freCoeff[5][0];
                    Cp[i][iLayer]=freCoeff[5][1];
                    Ps[i][iLayer]=real(freCoeff[0][0]);
                    Pp[i][iLayer]=real(freCoeff[0][1]);
                    if(iLayer==iFirstCoe){
                        Bs[i][iFirstCoe-1]=complex<double>(vosi[2][1],0.);//Rs full stack
                        Bp[i][iFirstCoe-1]=complex<double>(vosi[2][2],0.);
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
            }
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
            PLOTline1bar2(1,iRD,icol,9,Np,x.data(),y.data(),ErrXp.data(),ErrYp.data());
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
    int iSpec0Hemi1=nint(par[54][2]);
    if(iSpec0Hemi1==1)
        iShemispherical=true;
    else
        iShemispherical=false;
    ui->lineEdit_texture->setText("theta1="+QString::number(theta1)+" theta2="+QString::number(theta2));
    qDebug()<<"->Plot texturized";
    AdjRoughMax();
    SaveSetting(-1);
    double wl,VNK[N_CNK_HARD_MAX+1][3],vot[6][3];
    std::vector <double> Xp(NeV),Atext(NeV),Rtext(NeV),Ttext(NeV),ErrXp(NeV),ErrYp(NeV);
    fill_n(ErrXp.data(), NeV, 0.);
    fill_n(ErrYp.data(), NeV, 0.);
    for(int i=1;i<=NeV;i++){
        wl=ms.lambda[i-1];;
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
    }
    nextColor();
    PLOTline1bar2(1,0,iColor,1,NeV,Xp.data(),Ttext.data(),ErrXp.data(),ErrYp.data());
    PLOTline1bar2(1,0,iColor,3,NeV,Xp.data(),Rtext.data(),ErrXp.data(),ErrYp.data());
    PLOTline1bar2(1,0,iColor,9,NeV,Xp.data(),Atext.data(),ErrXp.data(),ErrYp.data());
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
    double WL;
    std::vector <double> Xp(Ndata),Yp(Ndata),ErrXp(Ndata),ErrYp(Ndata);
    int iChoiceWL=ui->comboBox_searchNK -> currentIndex();
    if(iChoiceWL==0)
        WL=ui->dSB_WLsearchNK -> value();
    else if(iChoiceWL==1)
        WL=ui->dSB_PAR_4_1 -> value();
    else
        WL=ui->dSB_PAR_4_2 -> value();
    int iWL=1;
    while(ms.lambda[iWL-1]<WL && iWL<NeV)
        iWL++;
    WL=ms.lambda[iWL-1];
    qDebug()<<"->searchNK @ WL= "<<WL;
    SaveSetting(-1);
    PLOTline1bar2(1,1,iCol,11,0,Xp.data(),Yp.data(),ErrXp.data(),ErrYp.data());//erase plot
    cnk[1].forceMode=3;//nk-uncknown will be set to forced values
    for(int imis=1;imis<=14;imis++){
        if(DATO[imis]==2){
            if(imis<=6){
                iCol=imis;
                qDebug("imis=%d: %f +-%f\n ",imis,ms.measures[imis-1].value[iWL-1],ms.measures[imis-1].error[iWL-1]);
                Ndat=SOLVE(imis,iWL,Xp.data(),Yp.data(),ErrYp.data());
                PLOTline1bar2(2,0,iCol,11,Ndat,Xp.data(),Yp.data(),ErrXp.data(),ErrYp.data());
            }
            else{
                int I=imis-6;
                if(I==1 || I==3 || I==5 || I==7) iCol=1;
                if(I==2 || I==4 || I==6 || I==8) iCol=4;
                qDebug("imis=%d: %f +-%f\n ",imis,ms.measures[I+5].value[iWL-1],ms.measures[I+5].error[iWL]);
                Ndat=SOLVE(imis,iWL,Xp.data(),Yp.data(),ErrYp.data());
                PLOTline1bar2(2,0,iCol,11,Ndat,Xp.data(),Yp.data(),ErrXp.data(),ErrYp.data());
            }
        }
    }
    cnk[1].forceMode=0;//nk-uncknown by the chosen option
}


void ksemawc::RefTrackG(){
    SaveSetting(-1);
    //SPADA();
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
    int NmisEnab=nint(par[22][2]);
    qDebug("Exhaustive numerical n-search in wl-n and wl-k spaces with N.meas=%d sampled on %d points",NmisEnab,NeV);
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
        qDebug()<<"results="<<result;
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
    qDebug()<<"calN1K2NK3= "<<calN1K2NK3;
    double n,k,kMin=100,kMax=-100.;
    double FM,FMold,nMin,dFMold=0.;
    par[1][1]=rxy[16][1];//nMin
    par[1][2]=rxy[16][2];//nMax
    par[2][1]=rxy[17][1];//kMin
    par[2][2]=rxy[17][2];//kMax
    par[6][1]=ui->dSB_PAR_6_1->value();
    pmTe[1][1]=ui->dSB_PM_86_1->value();
    pmTe[2][1]=ui->dSB_PM_87_1->value();
    pmTe[3][1]=ui->dSB_PM_88_1->value();
    pmTe[4][1]=ui->dSB_PM_89_1->value();
    int IL1=1;
    int IL2=NeV;
    Sol.clear();//reset Sol
    Sol.reserve(2*NeV);
    double Dn=(par[1][2]-par[1][1])/100.;
    k=par[2][1];
    for(int IL=IL1;IL<=IL2;IL++){
        ui->progressBar_ENS->setValue(static_cast<int>(static_cast<double>(IL-IL1)/static_cast<double>(IL2-IL1)*100.));
        ui->progressBar_ENS->update();
        par[24][2]=IL;
        double lam=ms.lambda[IL-1];
        double Trasm=1.,dt,dlim,tol;
        if(nint(par[50+nint(par[53][2])][3])==1){
            if(DATO[1]>0) Trasm=ms.measures[0].value[IL-1];
            if(DATO[2]>0) Trasm=ms.measures[1].value[IL-1];
            double dinco=pmD[nint(par[53][2])][1];
            if(iw==3) qDebug()<<"Trasm = "<<Trasm;
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
                cnk[1].forceMode=3;//nk-uncknown will be set to forced values
                cnk[1].nForced=n;
                cnk[1].kForced=k;
                double dt=0.001*std::abs(rxy[17][2]-rxy[17][1]);
                FM=FindRoot([this](double k){ return FMER(k); },k,dt,dlim,tol);
                k=cnk[1].kForced;
            }
            else if(calN1K2NK3==1){//n calculation with k=k_simulation
                double VNK[N_CNK_HARD_MAX+1][3];
                cnk[1].forceMode=1;//n-unknown will be set to forced value
                COSVNK(VNK,IL);
                k=VNK[1][2];
                cnk[1].nForced=n;
                FM=FMER(k);
            }
            else if(calN1K2NK3==2){//k calculation with n=n_simulation
                double VNK[N_CNK_HARD_MAX+1][3];
                cnk[1].forceMode=2;//k-unknown will be set to forced value
                COSVNK(VNK,IL);
                n=VNK[1][1];
                cnk[1].kForced=k;
                double dt=0.001*std::abs(rxy[17][2]-rxy[17][1]);
                FM=FindRoot([this](double k){ return FMER(k); },k,dt,dlim,tol);
                if(IL==iwl2print)
                    qDebug("IL=%d FM=%f n=%f k_fin=%f",IL,FM,n,cnk[1].kConst);
                k=cnk[1].kForced;
                Sol.emplace_back(std::array<double, 6>{1., lam, n, k, n*0.001, 0.});
                double Dk=k/100.;//error evaluation
                do{
                    k=k+Dk;
                    FMold=FMER(k*1.01);
                }while(FMold<1.);
                int iSol=Sol.size();
                iSol--;
                Sol[iSol][5]=std::abs(k-Sol[iSol][3]);
                i=101;
                //continue;
            }
            if(calN1K2NK3!=2){
                if((FMold>1.||i==0) && FM<=1.){
                    nMin=n;
                    kMin=k;
                    kMax=k;
                }
                else if(i>0 && FMold<=1. && FM<=1.){
                    kMin=min(kMin,k);
                    kMax=max(kMax,k);
                }
                else if(i>0 && FMold<=1. && FM>1.){
                    Sol.emplace_back(std::array<double, 6>{1., lam, (nMin+n)/2., (kMin+kMax)/2., (n-nMin)/2., (kMax-kMin)/2.});
                }
                else if(irm==1&& dFMold<0. && (FM-FMold)>0. && FM>1. && FMold>1.){//found a minimum
                    Sol.emplace_back(std::array<double, 6>{1., lam, n, k, n*0.0001, k*0.0001});
                }
            }
            dFMold=FM-FMold;
            FMold=FM;
            QCoreApplication::processEvents();
            if(iStop){
                qInfo()<<"The exhaustive computing was stopped as requested by the user!";
                return;
            }
        }//end for n
    }//end for wl
    int iSol=Sol.size();
    qDebug()<<"\nFound "<<iSol<<" nk-solutions";
    if(iSol==0){
        QMessageBox msgBox;
        msgBox.setText("No solution has been found!\nPlease change the n-range of the search");
        msgBox.setStandardButtons(QMessageBox::Ok);
        msgBox.exec();
    }
    nextColor();
    PlotNK(0);
    cnk[1].forceMode=0;//nk-unknown by the chosen option
}


void ksemawc::selectNsol(){
    if(m_graphs[12] == nullptr){
        QMessageBox msgBox;
        msgBox.setText("Please open the n(λ) graph first by clicking 'Plot Exp. Measures!'");
        msgBox.setStandardButtons(QMessageBox::Ok);
        msgBox.exec();
        return;
    }

    if(d_picker){
        delete d_picker;
        d_picker = nullptr;
        picker_m = nullptr;
    }
    polygonF.clear();

    QMessageBox msgBox;
    msgBox.setText("Please delimit with a polygon the n-solutions to be deleted!\n"
                   "Left click to set a point; click on the first point to terminate!");
    msgBox.setStandardButtons(QMessageBox::Ok);
    msgBox.exec();

    picker_m = new QwtPickerClickPointMachine();
    d_picker = new QwtPlotPicker(
        QwtAxis::XBottom, QwtAxis::YLeft,
        QwtPlotPicker::CrossRubberBand, QwtPicker::AlwaysOn,
        m_graphs[12]->canvas());   // ← sostituito G12_wn
    d_picker->setStateMachine(picker_m);
    connect(d_picker, qOverload<const QPointF &>(&QwtPlotPicker::selected),
            this, &ksemawc::drawPolygon);
}


void ksemawc::drawPolygon(QPointF pos){
    double x1,y1;
    x1=pos.x();
    y1=pos.y();
    if(L1E2==2)
        x1=1240./x1;
    polygonF<<QPointF(x1,y1);
    int Np=polygonF.count();
    if(Np<=1)
        return;
    std::vector <double> x(Np),y(Np),ErrXp(Np),ErrYp(Np);
    for(int i=0;i<Np;i++){
        x[i]=polygonF.at(i).x();
        y[i]=polygonF.at(i).y();
    }
    PLOTline1bar2(1,0,1,12,Np,x.data(),y.data(),ErrXp.data(),ErrYp.data());//n plot
    double delta=pow((x[Np-1]-x[0])/(rxy[20][2]-rxy[20][1]),2.)+
                   pow((y[Np-1]-y[0])/(rxy[16][2]-rxy[16][1]),2.);
    if(sqrt(delta)>0.01)
        return;
    QPointF ptF;
    int nSol=Sol.size();
    for(int i=0;i<nSol;i++){
        ptF.setX(Sol[i][1]);
        ptF.setY(Sol[i][2]);
        if(polygonF.containsPoint(ptF,Qt::OddEvenFill)){
            qDebug("i=%d wl=%f n=%f is contained! It will be canceled\n",i,x1,y1);
            Sol[i][0]=0.;//to be canceled
        }
        else
            Sol[i][0]=1.;
    }
    fflush(stdout);
    //cancel selected nk-solutions
    int newN=-1;
    std::vector<std::vector<double>> NKNEW(nSol, std::vector<double>(6, 0.0));
    for(int I=0;I<nSol;I++){
        if(Sol[I][0]>0.5){
            newN++;
            for(int k=1;k<=5;k++)
                NKNEW[newN][k]=Sol[I][k];
        }
    }
    Sol.clear();//reset Sol
    Sol.resize(newN+1);//newN is the last filled index (= #kept-1); size must be newN+1 (and 0 when nothing kept)
    for(int I=0;I<=newN;I++){
        for(int k=1;k<=5;k++)
            Sol[I][k]=NKNEW[I][k];
        Sol[I][0]=1;
    }
    qDebug("now the nk-solutions are %d",newN);
    if (d_picker) {
        d_picker->deleteLater();
        d_picker = nullptr;
        picker_m = nullptr; // È gestito dal picker, quindi basta azzerarlo
    }
    PlotNK(1);//plot nk
}


void ksemawc::IbridPlotFit(){
    lastAction="IbridCurrent";
    SaveSetting(-1);
    IbridKernel("g");
}


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
    state=ui->cB_EMA_1->checkState();
    if(state==Qt::Checked){
        QMessageBox msgBox;
        msgBox.setText("ATTENTION: EMA will be disabled to best-fit n, k , epsi1 and epsi2");
        msgBox.setInformativeText("Do you agree?");
        msgBox.setStandardButtons(QMessageBox::Yes | QMessageBox::No);
        msgBox.setDefaultButton(QMessageBox::No);
        int ret = msgBox.exec();
        switch (ret) {
        case QMessageBox::Yes:
            ui->cB_EMA_1->setCheckState(Qt::Unchecked);
            break;
        case QMessageBox::No:
            return;
            break;
        }
    }
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
    Sol.clear();
}


void ksemawc::setStop(){
    iStop=true;
    qDebug()<<"-> setStop: iStop= true";
    fflush(stdout);
}


void ksemawc::IbridKernel(QString rc){
    PanFitEnable();
    // data number
    qDebug()<<"-> IbridKernel with rc="<<rc;
    iDelCon=(ui->checkBox_deltaConnect->checkState()==Qt::Checked)?1:0;//refresh global from GUI (was stale from last Read/SaveSetting), cf. CalcMis
    ui->lineEdit_chi2Best->setStyleSheet("background-color: white;");
    int mwl=NeV;//SF and ELI
    int mwla=NeV;//enabled wavelengths
    int md=Sol.size();//nk-data
    int mda=0;//enabled nk-data
    int ial=0;
    constexpr int nparmax=17;
    int ibridok,m=0,ifit;
    int Info=0;
    Qt::CheckState state;
    state=ui->checkB_PAR_54_2 -> checkState ();
    if(state==Qt::Unchecked)
        par[54][2]=0;
    else
        par[54][2]=1;
    int iSpec0Hemi1=par[54][2];
    double pst[nparmax+1][3],vot[6][3];
    std::vector<std::vector<double>> miss(6, std::vector<double>(NeV+1, 0.0));
    std::vector<std::vector<double>> elis(9, std::vector<double>(NeV+1, 0.0));
    std::vector<std::vector<double>> PsiMat(mwl, std::vector<double>(4, 0.0));
    std::vector<std::vector<double>> DelMat(mwl, std::vector<double>(4, 0.0));
    std::vector <double> wwl(mwl),kk(mwl),nn(mwl),e1(mwl),e2(mwl),tn(mwl),tp(mwl),rn(mwl),rp(mwl),r1(mwl),ErrXp(mwl),ErrYp(mwl),
            ncent(mwl),ndev(mwl),kcent(mwl),kdev(mwl),psi(mwl),delta(mwl),ps(std::max(NeV,997));
    double chi2ini=1.E+09;
    par[38][1]=1.e+20;
    par[38][2]=0.;
    for(int i=0;i<md;i++){
        if(nint(Sol[i][0])>0){
            if(Sol[i][1]>=par[4][1] && Sol[i][1]<=par[4][2]){
                mda=mda+1;
                par[38][1]=min(par[38][1],Sol[i][1]);
                par[38][2]=max(par[38][2],Sol[i][1]);
            }
            else
                Sol[i][0]=0.;
        }
        if(Sol[i][4]<=0.){
            if(ial==0){
                qDebug("ERRn = 0 @%f  => ERRn = Dn =%f",Sol[i][1],par[21][3]);
                ial=1;
            }
            Sol[i][4]=std::abs(par[21][3]);
        }
    }
    qDebug("WLmin=%f WLmax=%f mwl=%d md=%d mda=%d",par[4][1],par[4][2],mwl,md,mda);
    par[38][3]=1.e+20;
    par[38][4]=0.;
    for(int i=1;i<=NeV;i++){
            par[38][3]=min(par[38][3],ms.lambda[i-1]);
            par[38][4]=max(par[38][4],ms.lambda[i-1]);
    }
    //initialization par for fit
    int n=ui->sB_PAR_35_5 -> value();//nint(par[35][5]);
    std::vector <double> p(n);
    //double* p=nullptr;
    //p = new double[n];
    for(int i=1;i<=n;i++){
        int ipm=nint(pmAt(i)[3]);
        p[i-1]=pmAt(ipm)[1];
    }
    if(ifirstcall==0){// IbridKernel has been called
        ifirstcall++;
    }
    // initialization
    double chi2fin=par[55][2];
    double fredeg=par[56][2];
    // minimal initialization
    double te=0.;//incident angle
    int s1p2=nint(par[27][2]);//polarization
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
                        ms.measures[i-1].value[ii-1]=miss[i][ii];//restore measured values
                    for(int i=1;i<=8;i++)
                        ms.measures[i+5].value[ii-1]=elis[i][ii];//idem
                }
                for(int i=1;i<=n;i++){
                    int ip=nint(pmAt(i)[3]);
                    p[i-1]=pst[i-1][1];//central value (jie=0) — pst is stored 0-based (see jie==0 / accumulo branches)
                    pmAt(ip)[1]=p[i-1];
                    if(pst[i-1][2]>.0)
                        pmAt(ip)[4]=sqrt(pst[i-1][2]);//error = deviation rms
                    else
                        pmAt(ip)[4]=0.;
                }
                for(int i=0;i<NeV;i++){
                    //store central value of solution
                    nn[i]=ncent[i];
                    kk[i]=kcent[i];
                    nkSol.n[i]=nn[i];
                    nkSol.k[i]=kk[i];
                    if(ndev[i]>0.)
                        ndev[i]=sqrt(ndev[i]);
                    else
                        ndev[i]=0.;
                    if(kdev[i]>0.)
                        kdev[i]=sqrt(kdev[i]);
                    else
                        kdev[i]=0.;
                }
                qDebug()<<"Saving nk-solutions in temp";
                Sol.resize(mwl);
                for(int i=0;i<mwl;i++){
                    Sol[i][1]=wwl[i];
                    Sol[i][2]=ncent[i];
                    Sol[i][3]=kcent[i];
                    Sol[i][4]=ndev[i];
                    Sol[i][5]=kdev[i];
                    Sol[i][0]=1.;//data enabled
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
                        miss[i][ii]=ms.measures[i-1].value[ii-1];//store measured values
                    for(int i=1;i<=8;i++)
                        elis[i][ii]=ms.measures[i+5].value[ii-1];
                }
            }
            else
                r=" ";//stop error computing
        }

        if(r=="i" && ieon==1 && jie>0){
            // measurement perturbation with the error
            int imisura=1;
            qDebug("..... IbridOne with error computing: step %d of %d",jie,jiemax);
            for(int i=1;i<=5;i++){
                if(DATO[i]==2){
                    int iespo=static_cast<int>((jie+1)/imisura);
                    qDebug("  measure= %d (-1)^%d= %d",i,iespo,static_cast<int>(pow(-1,iespo)));
                    for(int ii=1;ii<=NeV;ii++)
                        ms.measures[i-1].value[ii-1]=miss[i][ii]+ms.measures[i-1].error[ii-1]*pow(-1.,iespo);
                    imisura++;
                }
            }
            for(int i=8;i<=14;i=i+2){
                if(DATO[i]==2){
                    int iespo=static_cast<int>((jie+1)/imisura);
                    qDebug("  measure= %d (-1)^%d= %d",i,iespo,static_cast<int>(pow(-1,iespo)));
                    for(int ii=1;ii<=NeV;ii++)
                        ms.measures[i-7+5].value[ii-1]=elis[i-1][ii]+ms.measures[i-7+5].error[ii-1]*pow(-1.,iespo);
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
                    qDebug("-> launch IbridOne: compute k from Tn");
                    ibridok=1;
                }
                else if(DATO[2]==2){
                    qDebug("-> launch IbridOne: compute k from Tp teta= %f deg  S1P2= %d",par[6][1],s1p2);
                    ibridok=1;
                }
                else if(DATO[1]!=2 && DATO[2]!=2){
                    qDebug("-> ABORT IbridOne because transmittance is not enabled");
                    ibridok=-100;
                    r=" ";//abort
                    QMessageBox msgBox;
                    msgBox.setText("IbridOne aborted because transmittance is not enabled");
                    msgBox.setStandardButtons(QMessageBox::Ok);
                    msgBox.exec();
                    return;
                }
                if(DATO[3]==2){
                    qDebug("\tfit of Rn");
                    ibridok++;
                }
                if(DATO[4]==2){
                    qDebug("\tfit of Rp  teta= %f deg  S1P2= %d",par[6][1],s1p2);
                    ibridok++;
                }
                if(DATO[5]==2){
                    qDebug("\tfit of R1");
                    ibridok++;
                }
                if(DATO[8]==2){
                    qDebug("\tfit of PSI_1");
                    ibridok++;
                }
                if(DATO[10]==2){
                    qDebug("\tfit of PSI_2");
                    ibridok++;
                }
                if(DATO[12]==2){
                    qDebug("\tfit of PSI_3");
                    ibridok++;
                }
                if(DATO[14]==2){
                    qDebug("\tfit of PSI_4");
                    ibridok++;
                }
                if(ibridok<2){
                    qDebug("IbridOne aborted because needs T&R or T&DELTA-PSI!");
                    r=" ";//abort
                    QMessageBox msgBox;
                    msgBox.setText("IbridOne aborted because needs T&R or T&DELTA-PSI!");
                    msgBox.setStandardButtons(QMessageBox::Ok);
                    msgBox.exec();
                    return;
                }
                else{
                    m=(ibridok-1)*mwla;
                    qDebug("     ... on %d data and %d parameters",m,n);
                }
            }
            else if(r=="f"){// suppress p_fit ><  SELMQ
                if(nint(par[32][5])==0)
                    m=mda;
                else
                    m=2*mda;
                int npfit=n;
                for(int ipf=npfit;ipf>=1;ipf--){
                    int ipm=nint(pmAt(ipf)[3]);
                    qDebug("ipf = %d ipm = %d",ipf,ipm);
                    if(ipm<100){
                        pmAt(ipm)[2]=.0;
                        for(int iii=ipf;iii<=n;iii++){
                            pmAt(iii)[3]=pmAt(iii+1)[3];// move 1 step down
                            pmAt(iii)[4]=pmAt(iii+1)[4];
                            pmAt(iii)[5]=pmAt(iii+1)[5];
                        }
                        n=n-1;
                        qDebug("n = %d",n);
                    }
                    par[35][5]=n;
                }
                QString fitwhat;
                if(nint(par[32][5])==0)
                    fitwhat="n";
                else if(nint(par[32][5])==1)
                    fitwhat="n & k";
                else if(nint(par[32][5])==2)
                    fitwhat="epsi1 & epsi2";
                qDebug("launch fit of %s",fitwhat.toStdString().c_str());
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
                qDebug("launch fit of Selected Experimental Measurements: Nsem=%d Ndata=%d",iSEM,m);
            }

            //Fit, compute jac & cov matrices, err. e corr.
            fredeg=static_cast<double>(m>n?m-n:1);//guard m==n (avoid /0)
            par[56][2]=fredeg;
            ifit=1;
            QFile fchi2(filechi2);
            fchi2.resize(0);
            fchi2.close();
            int lwa=m*n+5*n+m;
            std::vector <int> iwa(n);
            std::vector <double> wa(lwa),fvec(m),FJAC(m*n);
            double tol=sqrt(dpmpar(1));
            qDebug("n(param)= %d m(Ndata)= %d lwa= %d fredeg=%f tol=%e s1p2=%d",
                   n,m,lwa,fredeg,tol,s1p2);
            if(n==0){
                QMessageBox msgBox;
                msgBox.setText("ATTENTION: please enable at least one parameter for best-fit!!!");
                msgBox.setStandardButtons(QMessageBox::Ok);
                msgBox.exec();
                return;
            }
            if(m==0){
                QMessageBox msgBox;
                msgBox.setText("ATTENTION: please select at least one measurement!!!");
                msgBox.setStandardButtons(QMessageBox::Ok);
                msgBox.exec();
                return;
            }
            iRecChi2=1;//save chi2 values at the end of called functions
            if(r=="f" && n>0){
                Info=lmdif1(FSQ, this,m,n,p.data(),fvec.data(),tol,iwa.data(),wa.data(),lwa);
            }
            else if(r=="fsem" && n>0 && m>0){
                Info=lmdif1(FSEM, this,m,n,p.data(),fvec.data(), tol, iwa.data(), wa.data(), lwa);
            }
            else if(r=="i"){
                Info=lmdif1(FRCK, this,m,n,p.data(), fvec.data(), tol, iwa.data(), wa.data(), lwa);
            }
            iRecChi2=0;

            //analysis chi2
            QFile fchi2b(filechi2);
            if(fchi2b.open(QIODevice::ReadOnly | QIODevice::Text)){
                QTextStream inp(&fchi2b);
                QString linea,st1;
                int iLine=0;
                double chi2;
                do{
                    linea=inp.readLine();
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
                    iLine++;
                }while(!inp.atEnd());
                fchi2b.close();
                par[27][1]=chi2ini;
                par[55][2]=chi2fin;
                ui->DP_PAR_27_1->setText(QString::number(chi2ini));
                ui->DP_PAR_55_2->setText(QString::number(chi2fin));
            }

            // set BF parameters
            //fnorm=enorm(m,fvec);
            qDebug("BF: chi2= %.9g",chi2fin);
            for(int i=1;i<=n;i++){
                int ip=nint(pmAt(i)[3]);
                pmAt(ip)[1]=p[i-1];
                //check value
                if((ip>=1 && ip<=99) || (ip>=496 && ip<=594) || (ip>=897 && ip<=996))//D, Rg, oscillators
                    if(pmAt(ip)[1]<.0)
                        pmAt(ip)[1]=0.;
                if(ip>=694 && ip<=891){//Nu (694-792) + EMA fractions (793-891)
                    if(pmAt(ip)[1]<.0)
                        pmAt(ip)[1]=0.;
                    if(pmAt(ip)[1]>1.)
                        pmAt(ip)[1]=1.;
                }
                pmAt(ip)[4]=.0;
                qDebug("\tp[%d]=pm[%d][1]=%.9g",i-1,ip,p[i-1]);
            }
            if(ifit==1){//write Info fit
                qDebug("********************INFO FIT***********************");
                qDebug("Enabled point m=%d",m);
                qDebug("N parameters n=%d => degree of freedom=%d",n,nint(fredeg));
                //qDebug("Final L2 norm of the residual = %f",fnorm);
                if(Info<0)
                    qDebug("excetution terminated by user");
                else if(Info==0)
                    qDebug("Info=0: improper input parameters");
                else if(Info==1)
                    qDebug("Info=1: algorithm estimates that the relative error in the sum of squares is at most tol");
                else if(Info==2)
                    qDebug("Info=2: algorithm estimates that the relative error between x and the solution is at most tol");
                else if(Info==3)
                    qDebug("Info=3: conditions for Info = 1 and Info = 2 both hold");
                else if(Info==4)
                    qDebug("Info=4:  fvec is orthogonal to the columns of the jacobian to machine precision");
                else if(Info==5)
                    qDebug("Info=5:  number of calls to fcn has reached or exceeded %d",200*(n+1));
                else if(Info==6)
                    qDebug("Info=6: tol is too small. no further reduction in the sum of squares is possible");
                else if(Info==7)
                    qDebug("Info=7: tol is too small. no further improvement in the approximate solution x is possible");
                qDebug("*************************************************");
                ifit=0;
            }

            if(m>0 && ieon==0){
                //computing jabobiana fjac
                for(int i=1;i<=996;i++)//save actual status
                    ps[i]=pmAt(i)[1];
                //int iflag=2;
                double epsfcn=100.*tol;//0.;
                qDebug("epsfcn=%e",epsfcn);
                if(r=="f" && n>0)
                    fdjac2(FSQ,this,m,n,p.data(),fvec.data(),FJAC.data(),m,epsfcn,wa.data());
                else if(r=="fsem" && n>0)
                    fdjac2(FSEM,this,m,n,p.data(),fvec.data(),FJAC.data(),m,epsfcn,wa.data());
                else if(r=="i" && n>0)
                    fdjac2(FRCK,this,m,n,p.data(),fvec.data(),FJAC.data(),m,epsfcn,wa.data());
                for(int i=1;i<=996;i++)//return to best fit
                    pmAt(i)[1]=ps[i];
                for(int i=1;i<=n;i++){
                    p[i-1]=pmAt(nint(pmAt(i)[3]))[1];
                }

                //computing matrix inv covar (fjac_trasposta)*(fjac)
                std::vector<std::vector<double>> cinv(n, std::vector<double>(n, 0.0));
                std::vector<std::vector<double>> cov(n, std::vector<double>(n, 0.0));
                for(int i=0;i<n;i++){
                    for(int j=0;j<n;j++){
                        cinv[i][j]=.0;
                        for(int i2=0;i2<m;i2++)
                            cinv[i][j]=cinv[i][j]+FJAC[i2+i*m]*FJAC[i2+j*m];//fjac[i + j * ldfjac]
                    }
                }

                //compunting matrix varianza-covarianza
                qDebug("computing MATINV...");
                // costruisci array di puntatori per MATINV
                std::vector<double*> cinv_ptr(n), cov_ptr(n);
                for(int i=0; i<n; i++){
                    cinv_ptr[i] = cinv[i].data();
                    cov_ptr[i]  = cov[i].data();
                }
                bool berr=MATINV(n,n,cinv_ptr.data(),cov_ptr.data());
                if(berr){
                    qDebug("done!");
                    //computing errors at 99% of confidency
                    for(int i=1;i<=n;i++){
                        if(cov[i-1][i-1]<1.e3){
                            pmAt(nint(pmAt(i)[3]))[4]=sqrt(std::abs(6.635*cov[i-1][i-1]));
                        }
                    }
                    fflush(stdout);

                    //computing global correlation
                    if(n>1){
                        for(int i=1;i<=n;i++){
                            if(1.-1./(cov[i-1][i-1]*cinv[i-1][i-1])>0.)
                                pmAt(nint(pmAt(i)[3]))[5]=sqrt(1.-1./(cov[i-1][i-1]*cinv[i-1][i-1]));
                            else
                                pmAt(nint(pmAt(i)[3]))[5]=1.;
                            if(pmAt(nint(pmAt(i)[3]))[5]>1.)
                                pmAt(nint(pmAt(i)[3]))[5]=1.;
                        }
                    }
                }
                //refresh PanFit
                Qt::CheckState state;
                for(int i=1;i<=npp;i++){
                    state=idToCheckBox["chBeParFit_"+QString::number(i)]-> checkState();
                    int ip=nint(ppm[i]);
                    //check value (uniform ip encoding v5.2)
                    if((ip>=1 && ip<=99) || (ip>=496 && ip<=594) || (ip>=897 && ip<=996))//D, Rg, oscillators
                        if(pmAt(ip)[1]<.0)
                            pmAt(ip)[1]=0.;
                    if(ip>=694 && ip<=891){//Nu (694-792) + EMA fractions (793-891)
                        if(pmAt(ip)[1]<.0)
                            pmAt(ip)[1]=0.;
                        if(pmAt(ip)[1]>1.)
                            pmAt(ip)[1]=1.;
                    }
                    idToLineEdit["DPparFitV_"+QString::number(i)] -> setText(QString::number(pmAt(ip)[1]));
                    if( state == Qt::Checked ){
                        idToLineEdit["DPparFitErr_"+QString::number(i)] -> setText(QString::number(pmAt(ip)[4]));
                        idToLineEdit["DPparFitGC_"+QString::number(i)] -> setText(QString::number(pmAt(ip)[5]));
                    }
                    else{
                        idToLineEdit["DPparFitErr_"+QString::number(i)] -> setText(QString::number(0));
                        idToLineEdit["DPparFitGC_"+QString::number(i)] -> setText(QString::number(0));
                    }
                }
            }
        }

        if(r=="f" || r=="fsem" || r=="i" || r=="g"){
            // The best-fit values were written back into the pm* arrays (and the
            // jacobian pass restored them from ps[]); propagate them into stack
            // so the recomputed curves use the fitted layer parameters.
            pmValuesToStack(stack);
            //refresh plots
            nextColor();
            double ymin=1.e6;
            double ymax=-1.e6;
            double chi2=.0;
            double delOld[4],del,delDiff[4]={0},VNK[N_CNK_HARD_MAX+1][3];
            qDebug()<<"r="<<r;
            for(int i=1;i<=mwl;i++){
                par[24][2]=i;//iWL
                double wl=ms.lambda[i-1];
                cnk[1].forceMode=0;//nk-unknown by the chosen option
                COSVNK(VNK,i);
                cnk[1].nForced=VNK[1][1];
                cnk[1].kForced=VNK[1][2];
                nkSol.n[i-1]=VNK[1][1];
                if(r=="i" || r=="g"){//IbridOne, compute k with best-fit parameters
                    cnk[1].forceMode=2;//k-unknown will be set to forced value
                    par[24][2]=i;
                    nkSol.n[i-1]=VNK[1][1];
                    double Trasm=1.,dt,dlim,tol;
                    if(nint(par[50+nint(par[53][2])][3])==1){
                        if(DATO[1]>0) Trasm=ms.measures[0].value[i-1];
                        if(DATO[2]>0) Trasm=ms.measures[1].value[i-1];
                        double dinco=pmD[nint(par[53][2])][1];
                        dt=-0.1*wl/4./3.14/dinco*log(Trasm);
                        dlim=1.e-10;
                        if(iw>0)
                            qDebug("Trasm = %f dt=%f",Trasm,dt);
                    }
                    else{
                        dt=0.001;//DELTA_k for thin film
                        dlim=1.e-5;
                    }
                    tol=dt/1000.;
                    cnk[1].forceMode=2;//k-unknown will bet set to forced values
                    double k=cnk[1].kForced;
                    FindRoot([this](double k){ return DELTAT(k); },k,dt,dlim,tol);
                    nkSol.k[i-1]=cnk[1].kForced;
                    cnk[1].forceMode=0;
                }
                else if(r=="f" || r=="fsem"){
                    nkSol.n[i-1]=VNK[1][1];
                    nkSol.k[i-1]=VNK[1][2];
                }
                wwl[i-1]=wl;
                nn[i-1]=nkSol.n[i-1];
                kk[i-1]=nkSol.k[i-1];
                e1[i-1]=nkSol.n[i-1]*nkSol.n[i-1]-nkSol.k[i-1]*nkSol.k[i-1];
                e2[i-1]=2.*nkSol.n[i-1]*nkSol.k[i-1];
                ymax=max(kk[i-1],ymax);
                ymin=min(kk[i-1],ymin);
                ErrXp[i-1]=0.;
                ErrYp[i-1]=0.;
                if(r=="f" || r=="fsem"){
                    cnk[1].forceMode=0;//nk-unknown by the chosen option
                }
                else if(r=="i" || r=="g"){
                    cnk[1].forceMode=2;//k-unknown set to forced value
                }

                //computing Tn Rn R1
                if(DATO[1]==2 || DATO[3]==2 || DATO[5]==2){
                    if(iSpec0Hemi1==1)
                        iShemispherical=true;
                    else
                        iShemispherical=false;
                    ASSEMBLER(i,1,0.,vot);
                    tn[i-1]=vot[1][s1p2]*100.;
                    rn[i-1]=vot[2][s1p2]*100.;
                    r1[i-1]=vot[3][s1p2]*100.;
                    if(DATO[1]==2)
                        chi2=chi2+pow((ms.measures[0].value[i-1]-tn[i-1]/100.)/ms.measures[0].error[i-1],2.)/fredeg;
                    if(DATO[3]==2)
                        chi2=chi2+pow((ms.measures[2].value[i-1]-rn[i-1]/100.)/ms.measures[2].error[i-1],2.)/fredeg;
                    if(DATO[5]==2)
                        chi2=chi2+pow((ms.measures[4].value[i-1]-r1[i-1]/100.)/ms.measures[4].error[i-1],2.)/fredeg;
                }
                if(DATO[2]==2 || DATO[4]==2){
                    te=par[6][1]*deg2rad;
                    par[54][2]=0;//Tp and Rp are always direct/specular
                    iShemispherical=false;
                    ASSEMBLER(i,1,te,vot);
                    par[54][2]=iSpec0Hemi1;
                    tp[i-1]=vot[1][s1p2]*100.;
                    rp[i-1]=vot[2][s1p2]*100.;
                    if(DATO[2]==2)
                        chi2=chi2+pow((ms.measures[1].value[i-1]-tp[i-1]/100.)/ms.measures[1].error[i-1],2.)/fredeg;
                    if(DATO[4]==2)
                        chi2=chi2+pow((ms.measures[3].value[i-1]-rp[i-1]/100.)/ms.measures[3].error[i-1],2.)/fredeg;
                }
                for(int ieli=0;ieli<4;ieli++){
                    if(DATO[7+2*ieli]==2 ||DATO[8+2*ieli]==2){
                        te=pmTe[1+ieli][1]*deg2rad;;//par[14+ieli][1]*deg2rad;
                        iShemispherical=false;
                        ASSEMBLER(i,1,te,vot);
                        del=vot[5][2];//Delta
                        if(iDelCon==1){
                            if(i==1){
                                delOld[ieli]=del;
                            }
                            del=del-Dperiod*nint((del-delOld[ieli])/Dperiod);
                            delOld[ieli]=del;
                        }
                        PsiMat[i-1][ieli]=vot[5][1];
                        DelMat[i-1][ieli]=del;
                        if(DATO[7+2*ieli]>0)
                            delDiff[ieli]=delDiff[ieli]+(del-ms.measures[6+2*ieli].value[i-1])/mwl;
                    }
                }
                cnk[1].forceMode=0;
            }
            for(int ieli=0;ieli<4;ieli++){
                if(DATO[7+2*ieli]>0)
                    qDebug("|IbridKernel> delDiff[%d]=%f",ieli,delDiff[ieli]);
            }
            if(ieon==1){// save central value and RMS deviation
                if(jie==0){
                    qDebug("***setting ncent & kcent and initialize ndev kdev****");
                    for(int i=0;i<NeV;i++){
                        ncent[i]=nn[i];
                        ndev[i]=0.;
                        kcent[i]=kk[i];
                        kdev[i]=0.;
                        qDebug("wl=%f ncent=%f kcent=%f",wwl[i],ncent[i],kcent[i]);
                    }
                    for(int i=0;i<n;i++){
                        pst[i][1]=p[i];
                        pst[i][2]=0.;
                    }
                }
                else{
                    qDebug("***computing ndev & kdev****");
                    for(int i=0;i<NeV;i++){
                        ndev[i]=ndev[i]+pow(ncent[i]-nn[i],2.)/jiemax;
                        kdev[i]=kdev[i]+pow(kcent[i]-kk[i],2.)/jiemax;
                    }
                    for(int i=0;i<n;i++)
                        pst[i][2]=pst[i][2]+pow(pst[i][1]-p[i],2.)/jiemax;
                }
            }

            //refresh plot
            PLOTline1bar2(1,0,iColor,12,mwl,wwl.data(),nn.data(),ErrXp.data(),ErrYp.data());//n plot
            PLOTline1bar2(1,0,iColor,13,mwl,wwl.data(),kk.data(),ErrXp.data(),ErrYp.data());//k plot
            if(nint(par[10][1])==1){
                PLOTline1bar2(1,0,iColor,14,mwl,wwl.data(),e1.data(),ErrXp.data(),ErrYp.data());//e1 plot
                PLOTline1bar2(1,0,iColor,15,mwl,wwl.data(),e2.data(),ErrXp.data(),ErrYp.data());//e2 plot
            }
            if(DATO[1]==2)
                PLOTline1bar2(1,0,iColor,1,mwl,wwl.data(),tn.data(),ErrXp.data(),ErrYp.data());//Tn plot
            if(DATO[2]==2)
                PLOTline1bar2(1,0,iColor,2,mwl,wwl.data(),tp.data(),ErrXp.data(),ErrYp.data());//Tp plot
            if(DATO[3]==2)
                PLOTline1bar2(1,0,iColor,3,mwl,wwl.data(),rn.data(),ErrXp.data(),ErrYp.data());//Rn plot
            if(DATO[4]==2)
                PLOTline1bar2(1,0,iColor,4,mwl,wwl.data(),rp.data(),ErrXp.data(),ErrYp.data());//Rp plot
            if(DATO[5]==2)
                PLOTline1bar2(1,0,iColor,5,mwl,wwl.data(),r1.data(),ErrXp.data(),ErrYp.data());//R1 plot
            for(int ieli=0;ieli<4;ieli++){
                if(DATO[7+2*ieli]==2){
                    int Nperiod=nint(delDiff[ieli]/Dperiod);
                    for(int i=0;i<NeV;i++)
                        delta[i]=DelMat[i][ieli]-Nperiod*Dperiod;
                    PLOTline1bar2(1,0,iColor,7,mwl,wwl.data(),delta.data(),ErrXp.data(),ErrYp.data());//DELTA plot
                }
                if(DATO[8+2*ieli]==2){
                    for(int i=0;i<NeV;i++)
                        psi[i]=PsiMat[i][ieli];
                    PLOTline1bar2(1,0,iColor,8,mwl,wwl.data(),psi.data(),ErrXp.data(),ErrYp.data());//PSI plot
                }
            }
        }

        //save Job# in /temp and set njobBest
        if(ieon==2 || (ieon==0 && (r=="f" || r=="i" || r=="fsem"))){
            StoreFitSet();
            double chi2Best= ui->lineEdit_chi2Best->text().toDouble();
            if(chi2fin<chi2Best){
                ui->lineEdit_chi2Best->setStyleSheet("background-color: yellow;");
                ui->lineEdit_chi2Best->setText(QString::number(chi2fin));
                ui->spinBox_jobBest->setValue(jobtot);
            }
        }

        if(ieon==0)
            r="X";
    }
}


void ksemawc::StoreFitSet(){
    jobtot++;
    QString fileStorecurrent=fileStore;
    fileStore=pathroot+"temp/semaw"+QString::number(jobtot)+".Spj";
    SaveSetting(3);
    fileStore=fileStorecurrent;
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

void ksemawc::reset(){
    ui->sB_PAR_8_1 -> setValue(0);
    ui->sB_PAR_8_2 -> setValue(0);
    jobtot=0;
    ifirstcall=0;
    ui->spinBox_jobBest -> setValue(0);
    ui->lineEdit_chi2Best->setText("1.E+99");
}


void ksemawc::SPADA(){
    qDebug("-> SPADA ....");

    // Fase 1 — resize strutture nuove
    ms.resize(NeV);
    nkMaterials.resize(N_CNK_USER_MAX);
    for(auto& nk : nkMaterials)
        nk.resize(NeV);
    nkSol.resize(NeV);

    // Tipo misura — popolato qui una volta sola
    static const MeasureType sfTypes[6] = {
        MeasureType::TransmittanceNormal,
        MeasureType::TransmittancePolarised,
        MeasureType::ReflectanceNormal,
        MeasureType::ReflectancePolarised,
        MeasureType::ReflectanceBack,
        MeasureType::Absorptance
    };
    for(int i = 0; i < 6; ++i)
        ms.measures[i].type = (DATO[i+1] >= 1.) ? sfTypes[i] : MeasureType::None;

    int N,MR,Ndati,ilinrim,NANG;
    double instr[21],Xtmp,Ytmp,wmin,wmax,dmin,dmax,Div,rstep,DLAM;
    std::vector<double> X,Y,Z;
    X.reserve(1000);  // initial dimension to avoid frequent re-allocation
    Y.reserve(1000);
    Z.reserve(1000);
    std::vector<double> Delta,Psi,ErrD,ErrP;
    Delta.reserve(1000);  // initial dimension to avoid frequent re-allocation
    Psi.reserve(1000);
    ErrD.reserve(1000);
    ErrP.reserve(1000);
    QString ST1,ST2,DIS[7],EST[8],ER,SPA,SAL,fnam,line,line2,line3,pezzo;
    //reset Min Max
    for(int I=1;I<=22;I++){
        rxy[I][3]=1.e+36;
        rxy[I][4]=-1.e+36;
    }
    int IUVIR=nint(par[24][1]);
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
            ST1=QString::number(nint(par[20+i][4]));
            EST[i]=ST2+ST1+DIS[i];
        }
    }
    if(!NANK[9].contains("mate/aa999")){
        qDebug()<<"\tSPADA-> call LoadFilenk with file-solutions "<<NANK[9];
        LoadFilenk();
    }
    for(int i=1;i<=20;i++)
        instr[i]=par[i][3];
    if(instr[6]<.5)
        ER="y";
    else
        ER=" ";
    SPA="L";
    if(nint(rxy[25][4])==2)
        SPA="E";
    SAL="n";
    MR=nint(par[22][1]);
    int nrif=nint(par[10][3]);
    dmin=par[4][1];
    dmax=par[4][2];
    int STEP=1;
    // set WL and baseline error
    for(int L=1;L<=NeV;L++){
        double xwe;
        if(SPA=="L")
            xwe=dmin+static_cast<double>(L-1)/(NeV-1)*(dmax-dmin);
        else
            xwe=1./(1./dmin+static_cast<double>(L-1)/(NeV-1)*(1./dmax-1./dmin));
        ms.lambda[L-1]=xwe;
        ms.measures[5].value[L-1]=1.;//Rreference stored in temporary allocation
        double errBL=instr[5];
        ms.driftBaseline[L-1]=instr[5];
        if(ER!="y"){
            if(IUVIR==1){
                if(xwe>184.9 && xwe<860.81)
                    errBL=.0005+.838951/(1+pow((xwe-183.5)/0.7,2.))+.0026/(1.+pow((xwe-870.)/25.,4.));
                if(xwe>3100. && xwe<3200.1)
                    errBL=-4.66E-04+.00476/(1+pow((xwe-3500.)/200.,2.));
            }
            if(errBL<instr[5]) errBL=instr[5];
            ms.driftBaseline[L-1]=errBL;
        }
    }
    //file-nk load (materials 1..N_CNK_USER_MAX). fnk[] holds the full path (with
    // ".nk") for EVERY nk slot, while the legacy NANK[] only stores slots 1..8.
    // The engine supports N_CNK_USER_MAX nk materials (SETVNK reads io=8..7+N_CNK_USER_MAX,
    // nkMaterials sized N_CNK_USER_MAX), so iterate over fnk[] up to N_CNK_USER_MAX —
    // otherwise files 9.. (e.g. one assigned to a Material) are never loaded and the
    // simulation reads an empty nkMaterials slot -> absurd results.
    for(int I=1;I<=N_CNK_USER_MAX;I++){
        fnam=fnk[I];
        if(!fnam.contains("mate/aa999.9") && !fnam.isEmpty()){
            qDebug()<<"\tSPADA-> load file-nk "<<fnam;
            QFile file(fnam);
            if(!file.open(QIODevice::ReadOnly | QIODevice::Text)){
                msgErrLoad("SPADA-fnam",fnam);
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
            qDebug()<<"\tNdata= "<<Ndati;
            double A2nm=0.1;
            if(line.contains("nm"))
                A2nm=1;
            for(int J=1;J<=Ndati;J++){
                line=stream.readLine();
                line=line.simplified();
                QStringList List;
                List =line.split(" ");
                int nV=List.count();
                if(nV!=5){
                    List =line.split("\t");
                    nV=List.count();
                    if(nV!=5){
                        qDebug()<<"nV="<<nV<<" line="<<line;
                        QMessageBox msgBox;
                        msgBox.setText("LoadFilenk-> ERROR reading file-nk="+fnam+"\n: separator char invalid or ERRn and ERRk missing");
                        msgBox.setStandardButtons(QMessageBox::Ok);
                        msgBox.exec();
                        continue;
                    }
                }
                X.push_back(List.at(0).toDouble()*A2nm);//wl
                Y.push_back(List.at(1).toDouble());//n
                Z.push_back(List.at(2).toDouble());//k
            }
            file.close();
            if(X.front()<X.back()){
                N=1;
                STEP=1;
            }
            else{
                N=Ndati;
                STEP=-1;
            }
            nkMaterials[I-1].resize(NeV);
            CONVER(X.data(),Y.data(),Ndati,N,STEP,7+I,1);
            CONVER(X.data(),Z.data(),Ndati,N,STEP,7+I,2);
            X.clear();
            Y.clear();
            Z.clear();
        }
    }
    //file-SF load
    if(nint(par[23][1])==1){
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
                continue;
            }
            else{
                qDebug()<<"\tSPADA-> i="<<i<<" loading fileSF "<<fnam;
                QTextStream stream (&file);
                line = stream.readLine();
                line2 = stream.readLine();
                if(line.contains("PE UV")){// file Perkin Elmer L900 - L950
                    qDebug()<<"\tfile Perkin Elmer L900 - L950";
                    Div=100.;
                    for(int ir=1;ir<=11;ir++)
                        line = stream.readLine();
                    line3 = stream.readLine();
                    if(line3.contains("PerkinElmer UV WinLab 7.1", Qt::CaseInsensitive))
                        ilinrim=70;
                    else if(line3.contains("PerkinElmer UV WinLab 7.4", Qt::CaseInsensitive))
                        ilinrim=74;
                    else if(line3.contains("PerkinElmer UV WinLab 5"))
                        ilinrim=66;
                    else
                        ilinrim=70;
                    for(int ir=1;ir<=ilinrim;ir++)
                        line = stream.readLine();
                    stream >> rstep;
                    STEP=nint(rstep);
                    stream >> Ndati;
                    N=Ndati;
                    qDebug()<<"Ndati="<<Ndati;
                    do{
                        line = stream.readLine();
                    }while(!line.contains("#DATA"));
                    for(int j=1;j<=Ndati;j++){
                        stream >> Xtmp >> Ytmp;
                        X.push_back(Xtmp);
                        Y.push_back(Ytmp);
                    }
                    fflush(stdout);
                }
                else if(line2.contains("#####SCALED")){//SCALED file
                    qDebug()<<"\tSCALED file";
                    double xin,xfi;
                    stream >> xin;
                    stream >> xfi;
                    stream >> Ndati;
                    if(xin<xfi){
                        STEP=1;
                        N=1;
                    }
                    else{
                        STEP=-1;
                        N=Ndati;
                    }
                    if(line2.contains("#####SCALED ")){
                        qDebug() << "\ttype: #####SCALED";

                        bool ok;
                        int ncolo = QInputDialog::getInt(this,
                            "Which column have to be loaded?",
                            "Column number (2..6):",
                            2,    // default
                            2,    // min
                            6,    // max
                            1,    // step
                            &ok);
                        if(!ok) return;  // utente ha premuto Cancel

                        for(int j=1; j<=Ndati; j++){
                            double Xtmp, Ytmp, zz;
                            if(ncolo==2) stream >> Xtmp >> Ytmp;
                            if(ncolo==3) stream >> Xtmp >> zz >> Ytmp;
                            if(ncolo==4) stream >> Xtmp >> zz >> zz >> Ytmp;
                            if(ncolo==5) stream >> Xtmp >> zz >> zz >> zz >> Ytmp;
                            if(ncolo==6) stream >> Xtmp >> zz >> zz >> zz >> zz >> Ytmp;
                            X.push_back(Xtmp);
                            Y.push_back(Ytmp);
                        }

                        int unita = QInputDialog::getInt(this,
                            "Absolute or %?",
                            "Absolute (0) o percent (1):",
                            0,    // valore default
                            0,    // minimo
                            1,    // massimo
                            1,    // step
                            &ok);
                        if(!ok) return;
                        Div = (unita == 0) ? 1.0 : 100.0;
                    }
                    else if(line2.contains("#####SCALEDA")){
                        qDebug("\ttype: #####SCALEDA");
                        for(int j=1;j<=Ndati;j++){
                            stream >> Xtmp >> Ytmp;
                            X.push_back(Xtmp);
                            Y.push_back(Ytmp);
                        }
                        Div=1.;
                    }
                    else if(line2.contains("#####SCALED%")){
                        qDebug("\ttype: #####SCALED%%");
                        for(int j=1;j<=Ndati;j++){
                            stream >> Xtmp >> Ytmp;
                            X.push_back(Xtmp);
                            Y.push_back(Ytmp);
                            //cout<< j<<"\t"<<X[j]<<"\t"<<Y[j]<<"\n";
                        }
                        Div=100.;
                    }
                }
                else if(line2.contains("RTmethod[")){
                    // file Transmittance by VASE
                    qDebug("\ttype VASE transmittance");
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
                        Xtmp=pezzo.toDouble();
                        if(iL0E1==1)
                            Xtmp=1240./Xtmp;
                        X.push_back(Xtmp);
                        pezzo=List.at(3).toLocal8Bit().constData();
                        Y.push_back(pezzo.toDouble());
                    }while(!stream.atEnd());
                    Ndati=j;
                    qDebug("\tNdat=%d",Ndati);
                    Div=1.;
                    STEP=1;
                    N=1;
                }
                else if(line2.contains("##XYUNITS= W")){//Lamba9 or old FTIR
                    qDebug("\ttype: Lambda9 or old FTIR");
                    line = stream.readLine();
                    line = stream.readLine();
                    stream >> line >> DLAM >> ST1>> Div;
                    qDebug("DLAM= %f Div= %f",DLAM,Div);
                    Div=Div*100.;
                    stream >> line >> Ndati>>ST1;
                    qDebug("Ndati= %d",Ndati);
                    do{
                        line = stream.readLine();
                        line=line.simplified();
                    }while(line.isEmpty());
                    int NRIGHE=static_cast<int>(Ndati/14.);
                    double DE=14.*((Ndati/14.)-NRIGHE);
                    if(DE>.5) NRIGHE++;
                    qDebug()<<"NRIGHE="<<NRIGHE<<" DE="<<DE;
                    double X1;
                    for(int J=1;J<=NRIGHE;J++){
                        line=stream.readLine();
                        line=line.simplified();
                        QStringList List;
                        List =line.split(" ");
                        int nV=List.count();
                        //qDebug()<<"nV="<<nV;
                        QString pezzo;
                        for(int iv=0;iv<nV;iv++){
                            pezzo=List.at(iv).toLocal8Bit().constData();
                            if(iv==0){
                                X1=pezzo.toDouble();
                            }
                            else{
                                X.push_back(X1+(iv-1)*DLAM);
                                Y.push_back(pezzo.toDouble());
                                //qDebug()<<"X="<<X.back()<<" Y="<<Y.back();
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
                    double xf,yf;
                    if(line3.contains("##DATA TYPE= UL")){//new file Lamba9)
                        qDebug("\ttype: new Lambda9");
                        for(int J=1;J<=6;J++)
                            line = stream.readLine();
                        stream >> line >> xf;
                        stream >> line >> yf;
                        stream >> line >> wmin;
                        stream >> line >> wmax;
                        stream >> line >> Ndati;
                        line = stream.readLine();
                        line = stream.readLine();
                        line = stream.readLine();
                        DLAM=(wmax-wmin)/static_cast<double>(Ndati-1);
                        Div=1./yf;
                        int NRIGHE=static_cast<int>(Ndati/10.);
                        double DE=10.*((Ndati/10.)-static_cast<int>(Ndati/10.));
                        if(DE>.5) NRIGHE++;
                        //int K;
                        double X1;
                        for(int J=1;J<=NRIGHE;J++){
                            K=(J-1)*10;
                            line=stream.readLine();
                            line=line.simplified();
                            //cout<<line.toStdString()<<"\n";
                            QStringList List;
                            List =line.split(" ");
                            int nV=List.count();
                            QString pezzo;
                            for(int iv=0;iv<nV;iv++){
                                pezzo=List.at(iv).toLocal8Bit().constData();
                                if(iv==0){
                                    X1=pezzo.toDouble();
                                }
                                else{
                                    X.push_back(X1*xf+(iv-1)*DLAM);
                                    Y.push_back(pezzo.toDouble());
                                }
                            }
                        }
                        N=1;
                        STEP=1;
                    }
                    else if(line3.contains("##DATA TYPE= UV")){//file Lambda19
                        qDebug("\ttype: Lambda19");
                        for(int J=1;J<=9;J++)
                            line = stream.readLine();
                        stream>>line>>xf;
                        stream>>line>>yf;
                        stream>>line>>wmax;
                        stream>>line>>wmin;
                        stream>>line>>Ndati;
                        for(int i1=1;i1<=4;i1++)
                            line = stream.readLine();
                        DLAM=(wmax-wmin)/(Ndati-1);
                        Div=100./yf;
                        int NRIGHE=static_cast<int>(Ndati/5.);
                        double DE=5.*((Ndati/5.)-static_cast<int>(Ndati/5.));
                        if(DE>.5) NRIGHE++;
                        //int K;
                        double X1;
                        for(int J=1;J<=NRIGHE;J++){
                            K=(J-1)*5;
                            line=stream.readLine();
                            line=line.simplified();
                            QStringList List;
                            List =line.split(" ");
                            int nV=List.count();
                            QString pezzo;
                            for(int iv=0;iv<nV;iv++){
                                pezzo=List.at(iv).toLocal8Bit().constData();
                                if(iv==0){
                                    X1=pezzo.toDouble();
                                }
                                else{
                                    X.push_back(X1*xf-(iv-1)*DLAM);
                                    Y.push_back(pezzo.toDouble());
                                }
                            }
                        }
                        N=Ndati;
                        STEP=-1;
                    }
                    else if(line3.contains("##DATA TYPE= IN")){//new IR file
                        qDebug("\ttype: new IR");
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
                                        X.push_back(pezzo.toDouble()*xf);
                                    else{
                                        X.push_back(X.back() + DLAM);
                                        Y.push_back(pezzo.toDouble()*yf);
                                    }
                                }
                            }
                            else{//FTIR SPECTRUM GS
                                if(K==0) Div=Div*100.;
                                double Xt1,Yt[4];
                                stream>>Xt1>>Yt[0]>>Yt[1]>>Yt[2]>>Yt[3];
                                //stream>>X1>>Y[K+1]>>Y[K+2]>>Y[K+3]>>Y[K+4];
                                for(int iiii=1;iiii<=4;iiii++){
                                    X.push_back(Xt1*xf+(iiii-1)*DLAM);
                                    Y.push_back(Yt[iiii-1]);
                                }
                                K=K+4;
                            }
                        }
                        N=1;
                        STEP=1;
                    }
                }
                file.close();
                if(Ndati<=0){
                    QMessageBox msgBox;
                    msgBox.setText("ATTENTION: the file "+fnam+" has no data!");
                    msgBox.setStandardButtons(QMessageBox::Ok);
                    msgBox.exec();
                    continue;
                }
                for(int L=0;L<Ndati;L++){
                    if(IUVIR==2) X[L]=1./X[L]*1.E7;//cm^-1 -> nm
                    Y[L]=Y[L]/Div;
                }
                if(i==0){//reference mirror
                    CONVER(X.data(),Y.data(),Ndati,N,STEP,6,IUVIR);
                    rxy[6][3]=1.e+36;
                    rxy[6][4]=-1.e+36;
                }
                else{
                    CONVER(X.data(),Y.data(),Ndati,N,STEP,i,IUVIR);
                    idToLineEdit["DP_RXY_"+QString::number(i)+"_"+QString::number(3)]-> setText(QString::number(rxy[i][3]));
                    idToLineEdit["DP_RXY_"+QString::number(i)+"_"+QString::number(4)]-> setText(QString::number(rxy[i][4]));
                }
                X.clear();
                Y.clear();
            }
        }
    }
    //file ELI load
    if(nint(par[23][2])==1){
        //int DP0cosDtanP1=ui->comboBox_DeltaPsiScale->currentIndex();
        for(int J=1;J<=4;J++){
            if(DATO[5+2*J]>0 && !NANK[11+J].contains("mate/aa999.9")){
                qDebug()<<NANK[11+J];
                fnam=pathroot+NANK[11+J].simplified()+".el";
                qDebug()<<"\tSPADA-> J="<<J<<" loading fileELI "<<fnam;
                QFile file(fnam);
                if(!file.open(QIODevice::ReadOnly | QIODevice::Text)){
                    msgErrLoad("SPADA-fnam tris",fnam);
                    continue;
                }
                else{
                    int iL0E1=0;
                    QString pezzo;
                    QTextStream stream (&file);
                    line = stream.readLine();
                    qDebug()<<line;
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
                            qDebug()<<line;
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
                    qDebug()<<"Nang= "<<NANG<<" NdatiEL= "<<Ndati<<" st4= "<<line2;
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
                        QMessageBox msgBox;
                        msgBox.setText("ATTENTION: set the format of file ELI!!!");
                        msgBox.setStandardButtons(QMessageBox::Ok);
                        msgBox.exec();
                        return;
                    }
                    qDebug("iformat=%d par[13+J][1]=%f",iformat,par[13+J][1]);
                    double theta2load=idToComboBox["cBteE"+QString::number(J)]->currentText().toDouble();
                    double yOld=0.;
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
                            PS=atan(PS)/deg2rad;
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
                        //if(std::abs(TE-par[13+J][1])<=1.){//0.001
                        if(std::abs(TE-theta2load)<=1.){
                            //cout<<TE<<"\t"<<LAM<<"\t"<<DE<<"\t"<<PS<<"\t"<<DDE<<"\t"<<DPS<<"\n";
                            if(iDelCon==1){
                                if(I==1)
                                    yOld=DE;
                                //qDebug("eV=%f DE=%f DE-yOld=%f (DE-yOld)/Dperiod=%f nint((DE-yOld)/Dperiod)=%d DEfinal=%f",1240./LAM,DE,yOld,(DE-yOld)/Dperiod,nint((DE-yOld)/Dperiod),DE-Dperiod*nint((DE-yOld)/Dperiod));
                                if(std::abs((DE-yOld)/Dperiod)>1.5)
                                    DE=DE-Dperiod*nint((DE-yOld)/Dperiod);
                                yOld=DE;
                            }
                            NDEL=NDEL+1;
                            X.push_back(LAM);
                            Delta.push_back(DE);
                            ErrD.push_back(DDE);
                            Psi.push_back(PS);
                            ErrP.push_back(DPS);
                            // X[NDEL]=LAM;
                            // EPR[1][NDEL][1]=DE;
                            // EPR[1][NDEL][2]=DDE;
                            // EPR[2][NDEL][1]=PS;
                            // EPR[2][NDEL][2]=DPS;
                        }
                    }
                    file.close();
                    qDebug("Theta= %f  NdatEl= %d",theta2load,NDEL);
                    if(NDEL==0)
                        return;
                    int Nstart=1;
                    STEP=1;
                    CONVER(X.data(),Delta.data(),NDEL,Nstart,STEP,N_CONV_ELI_BASE+2*(J-1)+0,1);//Delta
                    CONVER(X.data(),ErrD.data() ,NDEL,Nstart,STEP,N_CONV_ELI_BASE+2*(J-1)+0,2);//ErrDelta
                    CONVER(X.data(),Psi.data()  ,NDEL,Nstart,STEP,N_CONV_ELI_BASE+2*(J-1)+1,1);//Psi
                    CONVER(X.data(),ErrP.data() ,NDEL,Nstart,STEP,N_CONV_ELI_BASE+2*(J-1)+1,2);//ErrPsi
                    X.clear();
                    Delta.clear();
                    ErrD.clear();
                    Psi.clear();
                    ErrP.clear();
                    idToLineEdit["DP_RXY_"+QString::number(7+2*(J-1))+"_"+QString::number(3)]-> setText(QString::number(rxy[7+2*(J-1)][3]));
                    idToLineEdit["DP_RXY_"+QString::number(7+2*(J-1))+"_"+QString::number(4)]-> setText(QString::number(rxy[7+2*(J-1)][4]));
                    idToLineEdit["DP_RXY_"+QString::number(8+2*(J-1))+"_"+QString::number(3)]-> setText(QString::number(rxy[8+2*(J-1)][3]));
                    idToLineEdit["DP_RXY_"+QString::number(8+2*(J-1))+"_"+QString::number(4)]-> setText(QString::number(rxy[8+2*(J-1)][4]));
                }
            }
        }
        setRangeEli();
    }
    par[38][3]=par[4][1];
    par[38][4]=par[4][2];
    //par[21][1]=NeV;
    par[28][2]=NeV;
}


void ksemawc::CONVER(double *X,double *Y,int NDATI,int N1,int STEP,int I,int IUVIR){
    //int iw=1;//verbose (0 silent, 1 verbose)
    iw=nint(par[18][1]);
    qDebug("CONVER: NDATI=%d N1=%d STEP=%d I=%d IUVIR=%d",NDATI,N1,STEP,I,IUVIR);
    int ALARM=0;
    double cosDelMin=10.;
    double cosDelMax=-10.;
    double tanPsiMin=1.e+36;
    double tanPsiMax=-1.;
    double INSTR[21];
    double WL,WLa,WLb;

    for(int j=1;j<=20;j++){
        INSTR[j]=par[j][3];
    }
    for(int H=1;H<=NeV;H++){
        std::vector<double> xs,ys;
        xs.reserve(10);
        ys.reserve(10);
        WL=ms.lambda[H-1];
        double FACO=1.;
        if((I==3 || I==4 || I==5) && IUVIR==1){
            if(nint(INSTR[10])>0)
                FACO=FACO*ms.measures[5].value[H-1];
        }

        // points to average
        if(H>1)
            WLa=WL-(ms.lambda[H-1]-ms.lambda[H-2])/2.;
        else
            WLa=WL-(ms.lambda[H]-ms.lambda[H-1])/2.;
        if(H<NeV)
            WLb=WL+(ms.lambda[H]-ms.lambda[H-1])/2.;
        else
            WLb=WL+(ms.lambda[H-1]-ms.lambda[H-2])/2.;
        int N=N1;
        while(X[N-1]<WLa && (N+STEP)<=NDATI && (N+STEP)>=1)
            N=N+STEP;
        if(X[N-1]>WLb && N-STEP>=1 && N-STEP<=NDATI)
            N=N-STEP;
        if(H==iwl2print || iw==1){
            qDebug("\n>>> eV=%f WL= %f WLa= %f  WLb= %f",1240./WL,WL,WLa,WLb);
            qDebug("starting-N= %d X[%d]= %f",N,N-1,X[N-1]);
            qDebug("Data in [WLa , WLb]:");
        }

        //counting data in [WLa,WLb]
        int NPM=0;
        int J=0;
        do {
            NPM=NPM+1;
            xs.push_back(X[N-1]);
            ys.push_back(Y[N-1]);
            if(H==iwl2print ||iw==1)
                qDebug()<<N<<xs[J]<<"\t"<<ys[J];
            J++;
            N=N+STEP;
        }while(N<=NDATI && N>=1 && X[N-1]>=WLa && X[N-1]<=WLb);
        if(N<NDATI && N>1)
            N=N-STEP;
        else if(N<=1 && STEP==-1)
            N=2;
        else if(N>=NDATI && STEP ==1)
            N=NDATI-1;

        // check
        if(NPM<=2){//target: X[N-1] <= WL <= X[N-1+STEP]
            if(H==iwl2print ||iw==1)
                qDebug("... N=%d with NPM= %d -> set NPM=2",N,NPM);
            if(N>1 && N<NDATI){
                bool adj=true;
                while(adj && N>=1 && N<=NDATI && N+STEP>=1 && N+STEP<=NDATI){//keep X[N-1] and X[N+STEP-1] in bounds (was reading X[NDATI] past end)
                    if(X[N-1]>WL)
                        N=N-STEP;
                    if(X[N+STEP-1]<WL)
                        N=N+STEP;
                    if(X[N-1]<=WL && X[N-1+STEP]>=WL)
                        adj=false;
                    if(H==iwl2print ||iw==1)
                        qDebug("end while adj loop: N=%d X[N-1]=%f WL=%f X[N-1+STEP]=%f adj=%b",N,X[N-1],WL,X[N-1+STEP],adj);
                }
            }
            NPM=2;
            //ns.resize(2);
            xs.resize(2);ys.resize(2);
            if(N<=1){
                if(STEP==1)
                    N=1;
                else
                    N=2;
            }else if(N>=NDATI){
                if(STEP==1)
                    N=NDATI-1;
                else
                    N=NDATI;
            }
            for(int j=0; j<NPM; j++){
                //ns[j]=N;
                xs[j]=X[N-1];
                ys[j]=Y[N-1];
                if(H==iwl2print ||iw==1)
                    qDebug("N=%d xs[%d]=%f ys[%d]=%f",N,j,xs[j],j,ys[j]);
                N = N+STEP;
            }
            if((H==iwl2print ||iw==1) && X[N-1]<WL)
                qDebug("ATTENTION: X[N-1]<WL !!!");
        }

        if(NPM>3 && I>=N_CONV_ELI_BASE && nint(par[10][2])==0){ //ellissometric measurement and 3 points condition
            if(H==iwl2print ||iw==1)
                qDebug()<<"... ELLI meas -> set NPM=3!";
            do{
                N=N-STEP;
            }while(X[N-1]>WL);
            N=N-STEP;
            NPM=3;
            for(int j=0;j<NPM;j++){
                //ns[j]=N;
                xs[j]=X[N-1];
                ys[j]=Y[N-1];
                if(H==iwl2print ||iw==1)
                    qDebug()<<xs[j]<<"\t"<<ys[j];
                N=N+STEP;
            }
        }

        double a=0.,b=0.,c=0.;
        //linear extrapolation
        if(NPM==2){
            a=0;
            b=(ys[1]-ys[0])/(xs[1]-xs[0]);
            c=-b*xs[0]+ys[0];
            //YFIN=b*(WL-xs[0])+ys[0]
        }
        else if(NPM==3){//interpolation on 3 points
            if(std::abs(xs[0]-xs[1])<.01 || std::abs(xs[1]-xs[2])<.01){
                qDebug("ATTENTION: 2 data with same abscissa WL=%f xs[0]=%f xs[1]=%f xs[2]=%f\n!!!",
                       WL,xs[0],xs[1],xs[2]);
            }
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
            if(H==iwl2print ||iw==1 ){
                qDebug("Coeff parabola at 3 points:");
                qDebug("X: %f %f %f",XD1,XD2,XD3);
                qDebug("Y: %f %f %f",YD1,YD2,YD3);
                qDebug("a= %f  b= %f  c= %f",a,b,c);
                qDebug("FACO=%f ", FACO);
            }
        }

        // fit/interpolation
        else{ //if(NPM>3)
            int mFit=NPM;
            int nFit=3;
            int lwa=mFit*nFit+5*nFit+mFit;
            std::vector <int> iwa(nFit);
            std::vector <double> p(nFit),fvec(mFit),wa(lwa);
            p[2]=0.;
            pTF2.resize(NPM);
            for(int j=0;j<NPM;j++){
                pTF2[j].x = xs[j];
                pTF2[j].y = ys[j];
                p[2]=p[2]+ys[j];
            }
            p[0]=0.;
            p[1]=0.;
            p[2]=p[2]/static_cast<double>(NPM);
            double tol=sqrt(dpmpar(1));
            int Info=lmdif1(FPAR, pTF2.data(), mFit, nFit, p.data(),fvec.data(), tol, iwa.data(), wa.data(), lwa);
            if(Info<=0)
                qWarning("improper input parameters");
            c=p[2];
            b=p[1];
            a=p[0];
            if(H==iwl2print ||iw==1){
                qDebug("coeff of BF parabolic:\na= %f  b= %f  c= %f",a,b,c);
                for(int iii=0;iii<NPM;iii++)
                    qDebug()<< xs[iii]<<"\t"<<ys[iii];
            }
        }

        double YFIN=(a*WL*WL+b*WL+c)*FACO;
        if(H==iwl2print ||iw==1)
            qDebug("Final: WL=%f YFIN=%f",WL,YFIN);
        if(I>=8 && I<=7+N_CNK_USER_MAX){//nk file
            // file number = I-7. NANK[] only names slots 1..8 (and indices 9.. hold
            // unrelated entries / are out of bounds), so use it only for 1..8.
            QString nkName = (I-7>=1 && I-7<=8) ? NANK[I-7].simplified()
                                                : QString("nk-")+QString::number(I-7);
            if(YFIN<.0 && ALARM==0 && nint(par[9][1])!=1){
                if(IUVIR==1){
                    QMessageBox msgBox;
                    msgBox.setText("n-file "+nkName+" is <0 !!!\n Do you want impose n>=0?");
                    msgBox.setStandardButtons(QMessageBox::Yes | QMessageBox::No);
                    msgBox.setDefaultButton(QMessageBox::Yes);
                    int ret = msgBox.exec();
                    switch (ret) {
                    case QMessageBox::Yes:
                        par[9][1]=1;
                        break;
                    }
                    qDebug()<<"n-file "<<nkName<<" is <0 !!!";
                }
                else if(IUVIR==2){
                    QMessageBox msgBox;
                    msgBox.setText("k-file "+nkName+" is <0 !!!\n Do you want impose k>=0?");
                    msgBox.setStandardButtons(QMessageBox::Yes | QMessageBox::No);
                    msgBox.setDefaultButton(QMessageBox::Yes);
                    int ret = msgBox.exec();
                    switch (ret) {
                    case QMessageBox::Yes:
                        par[9][1]=1;
                        break;
                    }
                    qDebug()<<"k-file "<<nkName<<" is <0 !!!";
                }
                ALARM=1;
            }
            if(nint(par[9][1])==1 && YFIN<.0)
                YFIN=.0;
            int jM = I - 8;  // index in nkMaterials: 0..N_CNK_USER_MAX-1
            if(IUVIR == 1)
                nkMaterials[jM].n[H-1] = YFIN;
            else
                nkMaterials[jM].k[H-1] = YFIN;
        }
        if(I<=6){//SF measurements
            ms.measures[I-1].value[H-1] = YFIN;
            rxy[I][3]=min(rxy[I][3],YFIN*100.);
            rxy[I][4]=max(rxy[I][4],YFIN*100.);
            double DL=ms.driftBaseline[H-1];
            double errY=0;
            if(I==1 || I==2)
                errY=YFIN*INSTR[3]+DL;
            if(I==3 || I==4 || I==5)
                errY=YFIN*(INSTR[3]+INSTR[4])+DL;
            if(iw==1)
                qDebug("ms.measures[%d].value[%d]= %f",I-1,H-1,YFIN);
            ms.measures[I-1].error[H-1] = errY;
        }
        if(I>=N_CONV_ELI_BASE && I<=N_CONV_ELI_BASE+7){//eli measurements
            int eOff = I - N_CONV_ELI_BASE;  // 0..7: even=Delta, odd=Psi
            if(eOff%2==0){
                cosDelMin=min(cosDelMin,cos(YFIN*deg2rad));
                cosDelMax=max(cosDelMax,cos(YFIN*deg2rad));
            }
            else{
                tanPsiMin=min(tanPsiMin,tan(YFIN*deg2rad));
                tanPsiMax=max(tanPsiMax,tan(YFIN*deg2rad));
            }
            if(IUVIR==1){
                ms.measures[eOff+6].value[H-1] = YFIN;
                rxy[7+eOff][3]=min(rxy[7+eOff][3],YFIN);
                rxy[7+eOff][4]=max(rxy[7+eOff][4],YFIN);
            }
            else
                ms.measures[eOff+6].error[H-1] = YFIN;
        }
        if(iw==1){
            QMessageBox msgBox;
            msgBox.setText("Do you want to continue step by step?");
            msgBox.setStandardButtons(QMessageBox::Yes | QMessageBox::No);
            msgBox.setDefaultButton(QMessageBox::Yes);
            int ret = msgBox.exec();
            switch (ret) {
            case QMessageBox::No:
                iw=0;
                par[18][1]=0;
            }
        }
    }
    if(I>=N_CONV_ELI_BASE && I<=N_CONV_ELI_BASE+7 && IUVIR==1){
        if((I-N_CONV_ELI_BASE)%2==0){
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
    PointToFit2 *punti = static_cast<PointToFit2*>(p);
    for(int i=0;i<m;i++){
        double xi = punti[i].x;
        double yi = punti[i].y;
        fvec[i] = yi - (x[0] * xi * xi + x[1] * xi + x[2]);
    }
    return(0);
}


void ksemawc::ensureGraph(int ic, int iRD){
    if(ic < 1 || ic > N_GRAPHS)
        return;
    //qDebug("ensureGraph: ic=%d iRD=%d",ic,iRD);
    int wGhx = static_cast<int>(rxy[24][1]);
    int wGhy = static_cast<int>(rxy[24][2]);
    int graphDrift = 40;

    if(m_graphs[ic] == nullptr){
        m_graphs[ic] = new QwtPlot();
        m_graphs[ic]->setAttribute(Qt::WA_DeleteOnClose);
        m_graphs[ic]->setGeometry(nOpenGraph * graphDrift,
                                  nOpenGraph * graphDrift,
                                  wGhx, wGhy);
        nOpenGraph++;
        // La connect usa la copia di ic per cattura — corretta per tutti i 16
        connect(m_graphs[ic], &QObject::destroyed, this, [this, ic]() {
            m_graphs[ic]  = nullptr;
            m_pickers[ic] = nullptr;
        });

        m_pickers[ic] = new QwtPlotPicker(
            QwtAxis::XBottom, QwtAxis::YLeft,
            QwtPlotPicker::CrossRubberBand, QwtPicker::AlwaysOn,
            m_graphs[ic]->canvas());
        m_pickers[ic]->setStateMachine(new QwtPickerDragPointMachine());
        m_pickers[ic]->setRubberBandPen(QPen(Qt::red));
        m_pickers[ic]->setTrackerPen(QPen(Qt::black));
    }
    else {
        if(iRD == 1){
            m_graphs[ic]->detachItems(QwtPlotItem::Rtti_PlotCurve, true);
            m_graphs[ic]->detachItems(QwtPlotItem::Rtti_PlotIntervalCurve, true);
        }
        m_graphs[ic]->raise();
        m_graphs[ic]->activateWindow();
    }
}


void ksemawc::PLOTline1bar2(int iL1B2, int iRD, int iCol, int ic,int Ndata, double *Xp, double *Yp,
                                    double *ErrXp, double *ErrYp){
    qDebug("PLOT: iL1B2=%d iRD=%d iCol=%d ic=%d Ndata=%d L1E2=%d",iL1B2,iRD,iCol,ic,Ndata,L1E2);

    // --- Trasform Y  ---
    double Ymin = rxy[17][1], Ymax = rxy[17][2];
    std::vector<double> Px(Ndata), Py(Ndata);
    for(int i = 0; i < Ndata; i++){
        Px[i] = Xp[i];
        Py[i] = Yp[i];
    }
    if(ic == 13 && nint(par[31][5]) == 1){
        if(Ymin > 0.) Ymin = log10(Ymin); else Ymin = log10(Ymax / 1.E+04);
        Ymax = log10(Ymax);
        for(int i = 0; i < Ndata; i++){
            double Yup = Yp[i] + ErrYp[i];
            double Ydw = Yp[i] - ErrYp[i];
            if(Yup < Ymin) Yup = Ymin;
            if(Ydw < Ymin) Ydw = Ymin;
            Py[i]     = 0.5 * (log10(Yup) + log10(Ydw));
            ErrYp[i]  = 0.5 * (log10(Yup) - log10(Ydw));
        }
    }
    if(ic == 7 && nint(par[33][5]) == 1){
        for(int i = 0; i < Ndata; i++){
            double Yup = cos((Yp[i] + ErrYp[i]) * deg2rad);
            double Ydw = cos((Yp[i] - ErrYp[i]) * deg2rad);
            Py[i]     = 0.5 * (Yup + Ydw);
            ErrYp[i]  = 0.5 * std::abs(Yup - Ydw);
        }
    }
    if(ic == 8 && nint(par[33][5]) == 1){
        for(int i = 0; i < Ndata; i++){
            double Yup = tan((Yp[i] + ErrYp[i]) * deg2rad);
            double Ydw = tan((Yp[i] - ErrYp[i]) * deg2rad);
            Py[i]     = 0.5 * (Yup + Ydw);
            ErrYp[i]  = 0.5 * std::abs(Yup - Ydw);
        }
    }

    //check Graph number
    if(ic < 1 || ic > N_GRAPHS)
        return;

    // open Graph
    ensureGraph(ic, iRD);

    // --- Trasform X ---
    double WEmin = rxy[20][1], WEmax = rxy[20][2];
    if(L1E2 == 2 && ic != 10 && ic != 11){
        WEmin = 1240. / WEmin;
        WEmax = 1240. / WEmax;
        for(int i = 0; i < Ndata; i++)
            Px[i] = 1240. / Px[i];
    }

    // --- Lines and error bars ---
    QwtPlot* graph = m_graphs[ic];
    QwtPlotCurve *dataPlot = new QwtPlotCurve(labelQwt);
    dataPlot->setSamples(Px.data(), Py.data(), Ndata);
    int wT = static_cast<int>(rxy[24][3]);
    if(iL1B2 == 1){
        dataPlot->setPen(QPen(myColor[iCol], wT, Qt::DotLine));
    }
    else{
        QwtPlotIntervalCurve *range_plot = new QwtPlotIntervalCurve("range");
        dataPlot->setSymbol(new QwtSymbol(QwtSymbol::Ellipse,
                                          Qt::NoBrush, QPen(myColor[iCol]), QSize(2, 2)));
        dataPlot->setStyle(QwtPlotCurve::NoCurve);
        dataPlot->setRenderHint(QwtPlotItem::RenderAntialiased);
        QVector<QwtIntervalSample> range(Ndata);
        for(int i = 0; i < Ndata; i++)
            range[i] = QwtIntervalSample(Px[i], Py[i] - ErrYp[i], Py[i] + ErrYp[i]);
        QwtIntervalSymbol *errorbar = new QwtIntervalSymbol(QwtIntervalSymbol::Bar);
        errorbar->setPen(QPen(myColor[iCol], wT));//was iColor (global): use the curve's iCol like the symbol above
        errorbar->setWidth(wT);
        range_plot->setSamples(new QwtIntervalSeriesData(range));
        range_plot->setSymbol(errorbar);
        range_plot->setStyle(QwtPlotIntervalCurve::NoCurve);
        range_plot->attach(graph);
    }

    // --- Title axes and scale ---
    // Tabella: yTitle, rxy_yMin_idx, rxy_yMax_idx
    // (indici in rxy per i casi standard; i casi speciali sono gestiti sotto)
    struct GraphMeta {
        const char* yTitle;
        int         rxy_y;   // indice riga di rxy per scala Y
        bool        hasWE;   // true se usa asse WEmin/WEmax, false se usa asse proprio
    };
    static const GraphMeta meta[N_GRAPHS + 1] = {
//      Y_label                     index of rxy[n][4] to set Y_range
        {nullptr,                   0,  false}, // [0] non usato
        {"Tnormal (%)",             1,  true }, // ic=1 Tn
        {"Tpolarised (%)",          2,  true }, // ic=2 Tp
        {"Rnormal (%)",             3,  true }, // ic=3 Rn
        {"Rpolarized (%)",          4,  true }, // ic=4 Rp
        {"R1normal (%)",            5,  true }, // ic=5 R1
        {"Apds",                    6,  true }, // ic=6 Apds
        {nullptr,                   7,  true }, // ic=7 Delta (dynamic title)
        {nullptr,                   8,  true }, // ic=8 Psi   (idem )
        {"Absorptance_front (%)",  18,  true }, // ic=9 Afront
        {"Tau Rho Rho1 (%)",       23,  false}, // ic=10 Tau,Rho, Rho1 (X axis = Theta)
        {"k",                      17,  false}, // ic=11 nk space
        {"n",                      16,  true }, // ic=12 n
        {nullptr,                  17,  true }, // ic=13 k or log(k) (dynamic title)
        {"eps1",                   26,  true }, // ic=14 epsi1
        {"eps2",                   27,  true }, // ic=15 epsi2
        {"Absorptance_back (%)",   15,  true }, // Aback
    };

    const GraphMeta& m = meta[ic];

    // Titolo asse Y
    if(ic == 7){
        graph->setAxisTitle(0, nint(par[33][5]) == 0 ? "Delta (deg)" : "cos(Delta)");
    } else if(ic == 8){
        graph->setAxisTitle(0, nint(par[33][5]) == 0 ? "Psi (deg)" : "tan(Psi)");
    } else if(ic == 13){
        graph->setAxisTitle(0, nint(par[31][5]) != 1 ? "k" : "log10(k)");
    } else {
        graph->setAxisTitle(0, m.yTitle);
    }

    // Title X axis
    if(m.hasWE){
        graph->setAxisTitle(2, L1E2 == 1 ? "wavelength (nm)" : "photon energy (eV)");
    } else if(ic == 10){
        graph->setAxisTitle(2, "Theta (deg)");
    } else if(ic == 11){
        graph->setAxisTitle(2, "n");
    }

    // Scale Y axis
    if(ic == 7 && nint(par[33][5]) != 0)
        graph->setAxisScale(0, rxy[19][1], rxy[19][2], 0);
    else if(ic == 8 && nint(par[33][5]) != 0)
        graph->setAxisScale(0, rxy[22][1], rxy[22][2], 0);
    else if(ic == 13)
        graph->setAxisScale(0, Ymin, Ymax, 0);
    else
        graph->setAxisScale(0, rxy[m.rxy_y][1], rxy[m.rxy_y][2], 0);

    // Scale X axis
    if(ic == 10)
        graph->setAxisScale(2, rxy[21][1], rxy[21][2], 0);
    else if(ic == 11)
        graph->setAxisScale(2, rxy[16][1], rxy[16][2], 0);
    else
        graph->setAxisScale(2, WEmin, WEmax, 0);

    // Attach curve
    dataPlot->attach(graph);

    // Legend (only ic 9 e 16)
    if((ic == 9 || ic == 16) && labelQwt.contains("layer"))
        graph->insertLegend(new QwtLegend(), QwtPlot::RightLegend);

    graph->setAutoReplot();
    graph->show();
}


void ksemawc::COSVNK(double VNK[][3],int L){
    // L = index of wavelenght
    if(L==iwl2print)
        qDebug("COSVNK @ iWL=%d WL=%f",L,ms.lambda[L-1]);
    double f2;
    for(int J=1;J<=N_CNK_MAX;J++){
        int iForce=nint(cnk[J].forceMode);
        int i1=nint(cnk[J].nSrc);
        SETVNK(i1,J,VNK,L,2);
        int i2=nint(cnk[J].emaSrc);
        if(i2>=0){
            f2=pmFe[J][1];
            double n1=VNK[J][1];
            double k1=VNK[J][2];
            SETVNK(i2,J,VNK,L,5);
            double n2=VNK[J][1];
            double k2=VNK[J][2];
            auto [nEff, kEff] = EMA(n1,k1,n2,k2,f2);
            VNK[J][1]=nEff;
            VNK[J][2]=kEff;
        }
        if(iForce==1)
            VNK[J][1]=cnk[J].nForced;
        else if(iForce==2)
            VNK[J][2]=cnk[J].kForced;
        else if(iForce==3){
            VNK[J][1]=cnk[J].nForced;
            VNK[J][2]=cnk[J].kForced;
        }
        if(L==iwl2print)
            qDebug("COSVNK: cnk[%d].forceMode=%d i1=%d i2=%d VNK[%d][1]=%f VNK[%d]=%f",J,nint(cnk[J].forceMode),i1,i2,J,VNK[J][1],J,VNK[J][2]);
    }
}


void ksemawc::SETVNK(int io,int J,double VNK[][3],int L,int Kcte){
    // L = index of wavelenght
    // set VNK(J,1) and VNK(J,2) given CNK(J,k)
    if(io==0){
        VNK[J][1]=cnk[J].col(Kcte);
        VNK[J][2]=cnk[J].col(Kcte+1);
    }
    else if(io>=1 && io<=7){ // 1-7 Fit options
        double WL=ms.lambda[L-1];
        double ev=1240./WL;
        auto fd = FDISP(io,ev);
        VNK[J][1]=fd.sqn;
        VNK[J][2]=fd.sqk;
    }
    else if(io>=8 && io<=7+N_CNK_USER_MAX){
        VNK[J][1]=nkMaterials[io-8].n[L-1];
        VNK[J][2]=nkMaterials[io-8].k[L-1];
    }
    else if(io==NSRC_MATREF1){// "Material#1": reuse material #1's n,k (COSVNK fills VNK[1] first)
        VNK[J][1]=VNK[1][1];
        VNK[J][2]=VNK[1][2];
    }
}


EmaResult EMA(double N,double K,double NA,double KA,double FA){
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
    return{NE, KE};
}


FdispResult FDISP(int iopt,double eV){
    // nk computing by oscillators
    //qDebug("->FDISP iopt=%d @ E=%f",iopt,eV);
    double cln2=0.693147181;
    double sqrln2=0.832554611;
    double kVolSca=0;
    long double E3, Ex, Ey, W, KTi, r, Eb, Eh, El, derexp, derparab, TUx, TUy;

    // epsilon and nk by sum of oscillators
    int noscFT=nint(pf[iopt][1]);
    long double E=eV;
    long double epr=0.;
    long double epi=0.;
    long double eprj,epij;
    for(int io=2;io<=noscFT+1;io++){
        int i=nint(pf[iopt][io]);
        int ifu=nint(pmOs[1+(i-1)*5][1]);
        long double C=std::abs(pmOs[2+(i-1)*5][1]);
        long double E0=std::abs(pmOs[3+(i-1)*5][1]);
        long double D=std::abs(pmOs[4+(i-1)*5][1]);
        long double K=std::abs(pmOs[5+(i-1)*5][1]);
        //qDebug("i=%d ifu=%d C=%Lf E0=%Lf D=%Lf K=%Lf",i,ifu,C,E0,D,K);
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
            err=gsl_sf_complex_psi_e(static_cast<double>(real(xi)), static_cast<double>(imag(xi)), &re, &im );
            psi=complex<double>(re.val, im.val);
            cmpxOut=due*(log(xi)-pigc/tan(pigc*xi)-(complex<long double>)psi)-uno/xi;
            xi=sqrt(K/(E0+E+aim*D));
            err=gsl_sf_complex_psi_e(static_cast<double>(real(xi)), static_cast<double>(imag(xi)), &re, &im );
            psi=complex<double>(re.val, im.val);
            cmpxOut=cmpxOut+due*(log(xi)-pigc/tan(pigc*xi)-(complex<long double>)psi)-uno/xi;
            xi=sqrt(K/E0);
            err=gsl_sf_complex_psi_e(static_cast<double>(real(xi)), static_cast<double>(imag(xi)), &re, &im );
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
        else if(ifu==17){// Drude ionized
            long double DOM=D;
            long double Etr=sqrt(2.*C/PIG)*0.75;
            if(E>Etr) DOM=D*pow(Etr/E,1.5);
            eprj=-2.*C/PIG/(DOM*DOM+E*E);
            epij= 2.*C/PIG*DOM/(pow(E,3.)+E*DOM*DOM);
        }
        else if(ifu==18){// Lorentz-Dirac oscillator
            long double den=pow(E0*E0-E*E,2.)+pow(2.*E*D+E*E*E/K,2.);
            eprj=2.*C/PIG*(E0*E0-E*E)/den;
            epij=2.*C/PIG*(2.*E*D+E*E*E/K)/den;
// Expression as given in Prokopidis_2014
//            long double den=pow(E0*E0-E*E,2.)+pow(E*(2.*D+E0*E0/K),2.);
//            eprj=2.*C/PIG*((E0*E0-E*E)+E*E/K*(2.*D+E0*E0/K))/den;
//            epij=2.*C/PIG*(2.*E*D+E*E*E/K)/den;
        }
        else if(ifu==19){// k-oscillator for volume scattering
            if(iShemispherical)
                kVolSca=C*pow(E,K);
            else
                kVolSca=D*pow(E,K);
            eprj=0.;
            epij=0.;
        }
        epr=epr+eprj;
        epi=epi+epij;
        if(isnan(eprj) || isnan(epij)){
            qWarning("->FDISP-ERROR: iopt=%d OscNumber=%d ifu=%d E(eV)=%f eprj=%Lf epij=%Lf C=%Lf",iopt,io-1,ifu,eV,eprj,epij,C);
            iStop=true;
        }
    }
    double r_sqn = sqrt(sqrt(epr*epr + epi*epi) + epr) / sqrt(2.);
    double r_sqk = par[11][1] * (sqrt(sqrt(epr*epr + epi*epi) - epr) / sqrt(2.)+kVolSca);
    return {static_cast<double>(epr),
            static_cast<double>(epi),
            r_sqn,
            r_sqk
    };//par[11][1] is the attenuation coeff of k by Fit#N
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
            cDAW[i]=exp(-pow((2.*static_cast<double>(i)-1.)*H,2.));
    }
    if(std::abs(x)<0.2){// Use series expansion
        x2=x*x;
        value=x*(1.-A1*x2*(1.-A2*x2*(1.-A3*x2)));
    }
    else{// Use sampling theorem representation
        xx=std::abs(x);
        n0=2*nint(0.5*xx/H);
        xp=xx-static_cast<double>(n0)*H;
        e1=exp(2.*xp*H);
        e2=e1*e1;
        d1=static_cast<double>(n0+1);
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
            n1=2*static_cast<int>((Emax-Emin)/de/2.); //n1 must be even
            n2=0.;
        }
        else{
            de=(Emax-Emin)/nstep;
            if((Emax-E)<de*2.){
                Emax=E+de*2.;
                n1=2;
            }
            else{
                n1=2*static_cast<int>((Emax-E)/de/2.); //n1 must be even
            }
            de=(Emax-E)/n1;
            n2=2*static_cast<int>((E-Emin)/de/2.); //n2 must be even
            Emin=E-de*n2;
        }
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
    }while(std::abs((chi1-chi1old)/chi1)>0.001 && nstep<801);
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

    chi1=0.;
    do{
        chi1old=chi1;
        chi1=0.;
        Emax=Emaxin;
        Emin=Eminin;
        if (E>=Emax || E<=Emin){
            de=(Emax-Emin)/nstep;
            n1=2*static_cast<int>((Emax-Emin)/de/2.); //n1 must be even
            n2=0.;
        }
        else{
            de=(Emax-Emin)/nstep;
            if((Emax-E)<de*2.){
                Emax=E+de*2.;
                n1=2;
            }
            else{
                n1=2*static_cast<int>((Emax-E)/de/2.); //n1 must be even
            }
            de=(Emax-E)/n1;
            n2=2*static_cast<int>((E-Emin)/de/2.); //n2 must be even
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
    }while(std::abs((chi1-chi1old)/chi1)>0.001 && nstep<801);
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
            n1=2*static_cast<int>((Emax-Emin)/de/2.); //n1 must be even
            n2=0.;
        }
        else{
            de=(Emax-Emin)/nstep;
            if((Emax-E)<de*2.){
                Emax=E+de*2.;
                n1=2;
            }
            else{
                n1=2*static_cast<int>((Emax-E)/de/2.); //n1 must be even
            }
            de=(Emax-E)/n1;
            n2=2*static_cast<int>((E-Emin)/de/2.); //n2 must be even
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
    }while(std::abs((chi1-chi1old)/chi1)>0.001 && nstep<801);
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
            n1=2*static_cast<int>((Emax-Emin)/de/2.); //n1 must be even
            n2=0.;
        }
        else{
            de=(Emax-Emin)/nstep;
            if((Emax-E)<de*2.){
                Emax=E+de*2.;
                n1=2;
            }
            else{
                n1=2*static_cast<int>((Emax-E)/de/2.); //n1 must be even
            }
            de=(Emax-E)/n1;
            n2=2*static_cast<int>((E-Emin)/de/2.); //n2 must be even
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
    }while(std::abs((chi1-chi1old)/chi1)>0.001 && nstep<801);
    return(chi1);
}


long double ImCodyM1M2(long double E,long double C,long double E0,long double E3,long double D){
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
    else if (std::abs(E-EM)>=D) {
        chi2=chi2-K1*sqrt(std::abs(E-EM)-D);
    }
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
            n1=2*static_cast<int>((Emax-Emin)/de/2.); //n1 must be even
            n2=0.;
        }
        else{
            de=(Emax-Emin)/nstep;
            if((Emax-E)<de*2.){
                Emax=E+de*2.;
                n1=2;
            }
            else{
                n1=2*static_cast<int>((Emax-E)/de/2.); //n1 must be even
            }
            de=(Emax-E)/n1;
            n2=2*static_cast<int>((E-Emin)/de/2.); //n2 must be even
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
        chi1=chi1*2.*de*2./PIG;
        nstep=nstep*2;
    }while(std::abs((chi1-chi1old)/chi1)>0.001 && nstep<801);
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
    else if (std::abs(E-EM)>=D) {
        chi2=chi2-K1*sqrt(std::abs(E-EM)-D);
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
            n1=2*static_cast<int>((Emax-Emin)/de/2.); //n1 must be even
            n2=0.;
        }
        else{
            de=(Emax-Emin)/nstep;
            if((Emax-E)<de*2.){
                Emax=E+de*2.;
                n1=2;
            }
            else{
                n1=2*static_cast<int>((Emax-E)/de/2.); //n1 must be even
            }
            de=(Emax-E)/n1;
            n2=2*static_cast<int>((E-Emin)/de/2.); //n2 must be even
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
    }while(std::abs((chi1-chi1old)/chi1)>0.001 && nstep<801);
    return(chi1);
}


void ksemawc::ASSEMBLER(int iwl, int ikind, double teta, double vot[6][3]){
    int ivnkdw,ivnkup;
    double VNK[N_CNK_HARD_MAX+1][3],vosi[6][3];
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
    double wl=ms.lambda[iwl-1];//wavelength

    //vot inizialization
    vot[1][1]=1.;
    vot[1][2]=1.;
    for(int j=2;j<=5;j++){
        vot[j][1]=0.;
        vot[j][2]=0.;
    }

    //complex refractive indices and initialization entrance medium
    COSVNK(VNK,iwl);
    // all ELI measurements share a single input medium (CNK_IN_ELI)
    ivnkup=cnkInputMedium(ikind);
    irup=complex<double>(VNK[ivnkup][1],-VNK[ivnkup][2]);
    pq=pow(irup*sin(teta),2.);// (N_entrance*sin(teta))**2.
    //*** multilayer parameters
    int Nlayer=stack.nLayers();

    //*** loop on the interfaces = Nlayer+1
    int i=1;
    while(i<=Nlayer+1){
        int ncoe=0;
        //*** computing R T at i interface
        // Guard i<=Nlayer REQUIRED here: the last iteration (i==Nlayer+1)
        // represents the exit medium interface and has no corresponding layer.
        if(i <= Nlayer && (stack.layers[i-1].type > 0 || stack.layers[i-1].roughness.value > 0.)){//film or roughness
            int ifst=i;
            //thin film
            if(stack.layers[i-1].type > 0){//film!
                ncoe=1;
                while(ifst+ncoe<=Nlayer && stack.layers[ifst+ncoe-1].type > 0){//multilayer
                    ncoe=ncoe+1;
                }
                BUILDER(iwl,ikind,ifst,ncoe,pq,vosi);
                //save the physical parameters affering to the single interface
                vot[4][1]=vosi[4][1];// save Apds
                vot[5][1]=vosi[5][1];// save PSI
                vot[5][2]=vosi[5][2];// salva DELTA
            }

            //bulk with roughness
            else if(stack.layers[i-1].roughness.value > 0.){ //rough bulk!
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
                ivnkup=cnkInputMedium(ikind);// ELI all share CNK_IN_ELI
            else
                ivnkup=stack.layers[i-2].materialIndex;
            if(i<=Nlayer)
                ivnkdw=stack.layers[i-1].materialIndex;
            else //the next medium is the exit one for SF
                ivnkdw=CNK_OUT_SF;
            irup=complex<double>(VNK[ivnkup][1],-VNK[ivnkup][2]);
            irup2=irup*irup;
            irdw=complex<double>(VNK[ivnkdw][1],-VNK[ivnkdw][2]);
            irdw2=irdw*irdw;
            mups=sqrt(irup2-pq);
            mupp=irup2/mups;
            mdws=sqrt(irdw2-pq);
            mdwp=irdw2/mdws;
            //symmetric coating on back surface
            if(i==Nlayer+1 && stack.layers[i-2].type == 0 && nint(par[52][2])==1){ //symmetric multistrate on the rear face last layer
                for(int ii=1;ii<=2;ii++){
                    vosi[1][ii]=vot[1][ii];
                    vosi[2][ii]=vot[3][ii];
                    vosi[3][ii]=vot[2][ii];
                }
            }
            else{
                //vosi[1][1]=4.*real(mups)*real(mdws)/pow(std::abs(mups+mdws),2.);//requires transparent in out media
                //vosi[1][2]=4.*real(mupp)*real(mdwp)/pow(std::abs(mupp+mdwp),2.);
                vosi[2][1]=pow(std::abs((mups-mdws)/(mups+mdws)),2.);
                vosi[2][2]=pow(std::abs((mupp-mdwp)/(mupp+mdwp)),2.);
                vosi[1][1]=1.-vosi[2][1];//T=1-R because A=0 at the interface
                vosi[1][2]=1.-vosi[2][2];
                vosi[3][1]=vosi[2][1];
                vosi[3][2]=vosi[2][2];
                if(i==1 && Nlayer==1){ //BARE SUBSTRATE
                    //*** Apds
                    vot[4][1]=4.*real(irup*(conj(mdws)-irdw))/pow(std::abs(irup+mdws),2.);
                    //*** PSI e DELTA
                    rr=((mupp-mdwp)/(mupp+mdwp))/((mups-mdws)/(mups+mdws));
                    vot[5][1]=atan(std::abs(rr))/deg2rad;//PSI
                    vot[5][2]=arg(rr)/deg2rad;//DELTA
                }
            }
        }


        // combination with R T of the previous interface
        if(i==1)
            alfa=1.;
        else{
            ivnkup=stack.layers[i-2].materialIndex;
            irup=complex<double>(VNK[ivnkup][1],-VNK[ivnkup][2]);
            irup2=irup*irup;
            mups=sqrt(irup2-pq);
            alfa=exp(-4.*PIG*stack.layers[i-2].thickness.value/wl*imag(-mups));
        }
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
        i=i+ncoe+1;
    }
}


void ksemawc::BUILDER(int iwl,int ikind,int ifst,int ncoe, std::complex<double> pq,double vosi[6][3]){
    // Sized on N_LAYER_HARD_MAX (stack.nLayers() can reach it) +2: the layer loop
    // below also addresses index ifst+ncoe == imax+1. Zero-initialised: only rows
    // 1..imax are filled from the stack, the extra row must not hold garbage.
    constexpr int NLB=N_LAYER_HARD_MAX+2;
    int ivnk,NFA,iv[NLB]={},irougfa[NLB]={},irougms[NLB]={},iv2[NLB]={};
    int nu[NLB][3]={};
    //           [ i][0]=Nlayer with thickness nonuniform
    //           [ i][1]=NFAmin
    //           [ i][1]=NFAmax
    double ds[NLB]={},gn[NLB]={},cn[NLB]={},gk[NLB]={},
            ck[NLB]={},ru[NLB]={},nmedio,kmedio,al[NLB]={},VNK[N_CNK_HARD_MAX+1][3],
            dz,BpN,zz,ApN,ApK,BpK,CpK,CpN,Rs,R1s,Rp,R1p,Ts,Tp,As,Ap,ra,Ps,Pp,
            dsnu[NLB]={},DELTA,PSI;
    std::vector<double> d;
    complex<double> NQ,mus,mup,rhoS,rhopS,tauS,rhoP,rhopP,tauP,nq1,mus1,mup1,out[10][3],Bs,Cs,Bp,Cp,denoS,denoP;
    std::vector<std::complex<double>> ir;

    // Subroutine for build-up the indicated coherent part of the multilayer and/or roughness
    double wl=ms.lambda[iwl-1];//wavelength

    //layer parameters
    int imax=stack.nLayers();
    int nino=nint(par[29][1]);//discretization inhomogeneity
    int N=nint(par[28][1]);//   discretization roughness
    for(int i=1;i<=imax;i++){
        const Layer& l = stack.layers[i-1];
        ds[i]=l.thickness.value;
        al[i]=l.type + 1;                            // 0-based type → 1-based legacy al
        gn[i]=l.nGrad.value + 1240./wl*l.slopeNGrad.value;
        cn[i]=l.nCurv.value;
        gk[i]=l.kGrad.value;
        ck[i]=l.kCurv.value;
        ru[i]=l.roughness.value;
        dsnu[i]=l.nonUniformity.value;
        if(ru[i]>0.){//reduce the thickness of the layers which interface is rough for 3*sigma
            ds[i]=ds[i]-3.*ru[i];
            if(i>1)
                ds[i-1]=ds[i-1]-3.*ru[i];
        }
    }

    //refractive indices
    COSVNK(VNK,iwl);
    int ivnk1=cnkInputMedium(ikind);
    ir.push_back({0.,0.});//to create ir[0] which is not used
    ir.push_back(complex<double>(VNK[ivnk1][1],-VNK[ivnk1][2]));//this is ir[1]
    d.push_back(0.);//to create d[0] which is not used
    d.push_back(0.);//d[1]
    if(ifst>1){
        ivnk1=stack.layers[ifst-2].materialIndex;
        ir[1]=complex<double>(VNK[ivnk1][1],-VNK[ivnk1][2]);
    }
    NFA=1;

    // set index vector and refractive indices
    int ims=ifst;
    int nroug=0;
    int Nnouni=0;
    int iDWroug=0;
    // The trailing (ims==ifst+ncoe) pass adds the rough interface shared with the
    // NEXT layer down; ims<=imax is REQUIRED because that layer only exists when
    // the coherent group does not end on the last one. Without it, a rough bulk as
    // last layer (ncoe=1, ifst=imax) reads ru[imax+1] and then indexes
    // stack.layers[imax] out of bounds.
    while(ims<=(ifst+ncoe-1) || (ims==ifst+ncoe && ims<=imax && ru[ims]>0. && iDWroug==0)){
        if(ds[ims]<0.00001){
            if(iwl==iwl2print) qDebug("Skip layer %d because has zero thickness",ims);
            ims++;
            continue;
        }
        // roughness
        if(ru[ims]>0.01){// add double layer for roughness
            nroug=nroug+1;//one more rough-interface
            NFA++;
            ir.push_back(ir[NFA-1]);// same refractive index of the previous fase
            d.push_back(3.*ru[ims]);//first half rough layer
            NFA++;
            irougfa[nroug]=NFA;//NFA assigned to the new layer
            irougms[nroug]=ims;//index of the concerning model-layer
            ivnk=stack.layers[ims-1].materialIndex;
            ir.push_back(complex<double>(VNK[ivnk][1],-VNK[ivnk][2]));
            d.push_back(3.*ru[ims]);//second half rough thickness
            if(ims==ifst+ncoe){
                iDWroug++;
                continue;
            }
        }
        if(ds[ims]>0.00001){
            if(nint(al[ims])==3){ //inomogeneous layer
                ivnk=stack.layers[ims-1].materialIndex;
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
                    ir.push_back(complex<double>(ApN*zz*zz+BpN*zz+CpN,-(ApK*zz*zz+BpK*zz+CpK)));
                    d.push_back(dz*ds[ims]);
                }
                if(ru[ims]>0.)//roughness
                    ir[NFA-nino]=ir[NFA-nino+1];//nk for roughness
                if(dsnu[ims]>0.){//nonuniform thickness
                    Nnouni++;
                    nu[Nnouni][0]=ims;//N layer
                    nu[Nnouni][2]=NFA;//NFA_max
                    nu[Nnouni][1]=NFA-nino+1;//NFA_min
                }
            }
            else if(nint(al[ims])==2){// homogeneous layer
                NFA=NFA+1;
                ivnk=stack.layers[ims-1].materialIndex;
                ir.push_back(complex<double>(VNK[ivnk][1],-VNK[ivnk][2]));
                d.push_back(ds[ims]);
                if(dsnu[ims]>0.){//nonuniform thickness
                    Nnouni++;
                    nu[Nnouni][0]=ims;//N layer
                    nu[Nnouni][2]=NFA;//NFA_max
                    nu[Nnouni][1]=NFA;//NFA_min
                }
            }
            // if(iwl==iwl2print)
            //     qDebug("BUILDER: wl= %f ds[%d]=%f ir[%d]=%f+i*%f",wl,ims,ds[ims],NFA,real(ir[NFA]),imag(ir[NFA]));
        }
        ims++;
    }
    // substrate
    NFA=NFA+1;
    if(ims<=imax)
        ivnk=stack.layers[ims-1].materialIndex;
    else if(imax==1 && nint(al[1])==1)
        ivnk=stack.layers[imax-1].materialIndex;
    else{ //without substrate: the exit medium closes the stack
        ivnk=CNK_OUT_SF;
    }
    ir.push_back(complex<double>(VNK[ivnk][1],-VNK[ivnk][2]));
    d.push_back(0.);
    // if(iwl==iwl2print)
    //     qDebug("BUILDER: wl=%f ds[%d]=%f ivnk=%d NFA=%d ir[%d][1]=%f+i*%f",wl,ims,ds[ims],ivnk,NFA,NFA,real(ir[NFA]),imag(ir[NFA]));

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
    double Dit2=static_cast<double>(Nit2);
    for(int it2=0;it2<Nit2;it2++){//iteration over combinations of d-nonuniform
        for(int Jnu=0;Jnu<Nnouni;Jnu++){
            if(Jnu+1<Nnouni)
                iv2[Jnu]=static_cast<int>(it2/pow(M,Nnouni-Jnu-1));
            else{
                if(it2>0)
                    iv2[Jnu]++;
            }
            if(iv2[Jnu]+1>M)
                iv2[Jnu]=iv2[Jnu]-M*static_cast<int>(iv2[Jnu]/M);//was iv[Jnu] (roughness array): wrap must use iv2, cf. roughness loop iv[j] below
            int NFAmin=nu[Jnu+1][1];
            int NFAmax=nu[Jnu+1][2];
            double Dd=1./static_cast<double>(NFAmax-NFAmin+1);
            for(int NFAi=NFAmin;NFAi<=NFAmax;NFAi++){
                d[NFAi]=ds[nu[Jnu+1][0]]*Dd*(1.+dsnu[nu[Jnu+1][0]]*(iv2[Jnu]-N)/N);
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
                    iv[j]=static_cast<int>(it/pow(M,nroug-j-1));
                else{
                    if(it>0)
                        iv[j]++;
                }
                if(iv[j]+1>M)
                    iv[j]=iv[j]-M*static_cast<int>(iv[j]/M);
                ra=pow(-1.,iv[j])*dx[N][static_cast<int>((iv[j]+1)/2.+1)];//normalized thickness variation
                wtot=wtot*w[N][static_cast<int>((iv[j]+1)/2.+1)];         //cumulative weight
                d[irougfa[j+1]-1]=(3.+ra)*ru[irougms[j+1]];   //thickness first half rough layer
                d[irougfa[j+1]]=(3.-ra)*ru[irougms[j+1]];//thickness seconf half rough layer
            }
            CALFRE(NFA,wl,pq,ir,d,out);
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
            PSI=PSI+atan(std::abs(rhoP/rhoS))/deg2rad*wtot;
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
        if(nint(par[54][2])==0){ //specular SF measurements
            NQ=ir[NFA]*ir[NFA];
            mus=sqrt(NQ-pq);
            mup=NQ/mus;
            nq1=ir[1]*ir[1];
            mus1=sqrt(nq1-pq);
            mup1=nq1/mus1;

            vosi[1][1]=vosi[1][1]+real(mus)/real(mus1)*pow(std::abs(tauS),2.)/Dit2;
            vosi[1][2]=vosi[1][2]+real(mup)/real(mup1)*pow(std::abs(tauP),2.)/Dit2;
            vosi[2][1]=vosi[2][1]+pow(std::abs(rhoS),2.)/Dit2;
            vosi[2][2]=vosi[2][2]+pow(std::abs(rhoP),2.)/Dit2;
            vosi[3][1]=vosi[3][1]+pow(std::abs(rhopS),2.)/Dit2;
            vosi[3][2]=vosi[3][2]+pow(std::abs(rhopP),2.)/Dit2;
        }
        else{ // hemispherical SF measurements
            vosi[1][1]=vosi[1][1]+Ts/Dit2;
            //vosi[1][2]=vosi[1][2]+Tp/Dit2;//Tp always direct
            vosi[2][1]=vosi[2][1]+Rs/Dit2;
            //vosi[2][2]=vosi[2][2]+Rp/Dit2;//Rp always specular
            vosi[3][1]=vosi[3][1]+R1s/Dit2;
            vosi[3][2]=vosi[3][2]+R1p/Dit2;
        }
        vosi[4][1]=vosi[4][1]+As/Dit2;
        vosi[4][2]=vosi[4][2]+Ap/Dit2;
        // vosi[5][1]=vosi[5][1]+atan(std::abs(rhoP/rhoS))/deg2rad/Dit2;//PSI
        // vosi[5][2]=vosi[5][2]+arg(rhoP/rhoS)/deg2rad/Dit2;//DELTA
        vosi[5][1]=vosi[5][1]+PSI/Dit2;//PSI
        vosi[5][2]=vosi[5][2]+DELTA/Dit2;//DELTA
    }
}


void CALFRE(int NFAin,double wl,complex<double> pq,const std::vector<std::complex<double>>& IR,std::vector<double> d,complex<double> out[10][3]){
    complex<double> SC[3][3],PC[3][3],MUS,MUP,NQ,CI,DELTA,C,S,B,BP,CP,deno,nq1,mus1,mup1,
            rhoS,rhopS,tauS,rhoP,rhopP,tauP;
    complex<double> a11,a12,a21,a22,b11,b12,b21,b22;
    CI=complex<double>(0.,1.);

//check if a layer absorbs completely the light. In that case this layer becomes the exit medium
    int NFA=NFAin;
    for(int I=2;I<=NFAin-1;I++){
        NQ=IR[I]*IR[I];
        MUS=sqrt(NQ-pq);
        MUP=NQ/MUS;
        DELTA=2.*PIG*d[I]*MUS/wl;
        if(std::abs(imag(DELTA))>ABSmax) {
            NFA=I;
            break;
        }
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
        //if(nint(wl)==2500)
        //    qDebug("CALFRE: wl=%f IR[%d]=%f-i*%f d[%d]=%f",wl,I,real(IR[I]),imag(IR[I]),I,d[I]);

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
    }
    NQ=IR[NFA]*IR[NFA];
    MUS=sqrt(NQ-pq);
    MUP=NQ/MUS;
    nq1=IR[1]*IR[1];
    mus1=sqrt(nq1-pq);
    mup1=nq1/mus1;
    //if(nint(wl)==2500)
    //    qDebug("CALFRE: wl=%f IR[%d]=%f+i*%f mus1=%f+i%f mup1=%f+i*%f",wl,NFA,real(IR[NFA]),imag(IR[NFA]),real(mus1),imag(mus1),real(mup1),imag(mup1));

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
    double Rs=pow(std::abs(rhoS),2.);
    double R1s=pow(std::abs(rhopS),2.);
    double Ts=real(MUS)/real(mus1)*pow(std::abs(tauS),2.);
    double As=4.*real(mus1)*real(B*conj(C)-MUS)/pow(std::abs(deno),2.);
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
    double Rp=pow(std::abs(rhoP),2.);
    double R1p=pow(std::abs(rhopP),2.);
    double Tp=real(MUP)/real(mup1)*pow(std::abs(tauP),2.);
    double Ap=4.*real(mup1)*real(B*conj(C)-MUP)/pow(std::abs(deno),2.);
    out[2][0]=B;
    out[3][0]=C;
    out[0][2]=real(B*conj(C));//Poynting vector p-pol

    // PSI & DELTA
    double PSI=atan(std::abs(rhoP/rhoS))/deg2rad;
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
}


int ksemawc::SOLVE(int imis,int iWL,double *Xp,double *Yp,double *ErrYp){
    int Ndat=0;
    int ikind=0,ivot=0,IPD=0;
    double tetar=0.,KCL,vot[6][3];
    int iSpec0Hemi1=nint(par[54][2]);
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
        tetar=pmTe[1][1];
        ivot=5;
        IPD=imis-static_cast<int>(imis/2.)*2;
    }
    else if(imis==9 || imis==10){
        ikind=4;
        tetar=pmTe[2][1];
        ivot=5;
        IPD=imis-static_cast<int>(imis/2.)*2;
    }
    else if(imis==11 ||imis==12){
        ikind=5;
        tetar=pmTe[3][1];
        ivot=5;
        IPD=imis-static_cast<int>(imis/2.)*2;
    }
    else if(imis==13||imis==14){
        ikind=6;
        tetar=pmTe[4][1];
        ivot=5;
        IPD=imis-static_cast<int>(imis/2.)*2;
    }
    iShemispherical=false;
    if(imis==1 || imis==3 ||imis==5)
        if(iSpec0Hemi1==1)
            iShemispherical=true;
    if(imis==2 || imis==4)
        par[54][2]=0.;
    int s1p2=nint(par[27][2]);
    qDebug("-> SOLVE: imis=%d ikind=%d iWL=%d s1p2=%d teta=%f",imis,ikind,iWL,s1p2,tetar);
    tetar=tetar*deg2rad;
    double Dn=(par[1][2]-par[1][1])/100.;//n-step
    double Dk=(par[2][2]-par[2][1])/500.;//k-step
    double k1,FM,FMold=1.E+10;
    for(int iN=0;iN<101;iN++){
        cnk[1].nForced=par[1][1]+iN*Dn;
        FMold=1.E+10;
        for(int iK=0;iK<501;iK++){
            KCL=par[2][1]+iK*Dk;
            cnk[1].kForced=KCL;
            ASSEMBLER(iWL,ikind,tetar,vot);
            if(imis<=6)
                FM=std::abs(ms.measures[imis-1].value[iWL-1]-vot[ivot][s1p2])/ms.measures[imis-1].error[iWL-1];
            else
                FM=std::abs(sin(deg2rad*(ms.measures[imis-1].value[iWL-1]-vot[ivot][IPD+1])/2.))/(deg2rad*ms.measures[imis-1].error[iWL-1]);
            if(FM<=1. && FMold>1.){
                k1=KCL;
            }
            else if(FM>1. && FMold<=1.){
                Xp[Ndat]=cnk[1].nForced;
                Yp[Ndat]=0.5*(k1+KCL);
                ErrYp[Ndat]=0.5*(KCL-k1);
                if(ErrYp[Ndat]<Dk/4.)
                    ErrYp[Ndat]=Dk/4.;
                Ndat++;
            }
            FMold=FM;
        }
    }
    par[54][2]=iSpec0Hemi1;
    return(Ndat);
}


double ksemawc::FMER(double k){
    int iWL=nint(par[24][2]);
    int s1p2=nint(par[27][2]);
    int iSpec0Hemi1=nint(par[54][2]);
    double teta[6];
    cnk[1].kForced=k;
    teta[1]=par[6][1];
    teta[2]=pmTe[1][1];
    teta[3]=pmTe[2][1];
    teta[4]=pmTe[3][1];
    teta[5]=pmTe[4][1];
    for(int i=1;i<=5;i++)
        teta[i]=teta[i]*deg2rad;
    double vot[6][3];
    double FM=.0;
    int NMIS=0;
    if(DATO[1]==2 || DATO[3]==2 ||DATO[5]==2){
        if(iSpec0Hemi1==1)
            iShemispherical=true;
        else
            iShemispherical=false;
        ASSEMBLER(iWL,1,0.,vot);
        for(int i=1;i<=3;i++){
            int j=2*i-1;
            if(DATO[j]==2){
                FM=FM+pow((vot[i][s1p2]-ms.measures[j-1].value[iWL-1])/ms.measures[j-1].error[iWL-1],2.);
                NMIS++;
            }
        }
    }
    if(DATO[2]==2 || DATO[4]==2){
        par[54][2]=0;//Tp and Rp are always direct/specular
        iShemispherical=false;
        ASSEMBLER(iWL,1,teta[1],vot);
        par[54][2]=iSpec0Hemi1;//restore (uses the function-scope iSpec0Hemi1; inner shadow removed)
        for(int i=1;i<=2;i++){
            int j=2*i;
            if(DATO[j]==2){
                FM=FM+pow((vot[i][s1p2]-ms.measures[j-1].value[iWL-1])/ms.measures[j-1].error[iWL-1],2.);
                NMIS++;
            }
        }
    }
    if(DATO[6]==2){
        if(iSpec0Hemi1==1)
            iShemispherical=true;
        else
            iShemispherical=false;
        ASSEMBLER(iWL,2,0.,vot);
        FM=FM+pow(vot[4][1]-ms.measures[5].value[iWL-1],2.);
        NMIS++;
    }
    for(int I=1;I<=4;I++){
        iShemispherical=false;
        if(DATO[5+2*I]==2 || DATO[6+2*I]==2){
            ASSEMBLER(iWL,I+2,teta[1+I],vot);
            double DE=ms.measures[6+2*(I-1)].value[iWL-1];
            double DC=vot[5][2];//Delta;
            double DDE=ms.measures[6+2*(I-1)].error[iWL-1]*deg2rad;
            double PE=ms.measures[7+2*(I-1)].value[iWL-1];
            double PC=vot[5][1];//PSI
            double DP=ms.measures[7+2*(I-1)].error[iWL-1];
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
    //qDebug("FM=%f",FM);
    return(FM);
}


int ksemawc::FSQ(void *p, int m, int n, const double *x, double *fvec, int iflag){
    /* calculate the functions at x and return the values in fvec[0] through fvec[m-1] */
    //struct pointToFit2 *pTF2 = (struct pointToFit2 *)p;
    int ioptf=nint(pmOs[0][1]);
    // addressing fit-parameters to PM
    for(int i=0;i<n;i++){
        int ip=nint(pmAt(i+1)[3]);
        pmAt(ip)[1]=x[i];
    }

    // computing fvec for n fit
    int j=-1;
    double chi2=0.;
    double fredeg=real(m>n?m-n:1);//guard m==n (avoid /0)
    double erynMIN=(rxy[16][2]-rxy[16][1])/1.e4;
    double erykMIN=(rxy[17][2]-rxy[17][1])/1.e4;
    double ere1MIN=(rxy[26][2]-rxy[26][1])/1.e4;
    double ere2MIN=(rxy[27][2]-rxy[27][1])/1.e4;
    double wl,eV,y,ery,nVal,kVal,epsR,epsI;
    int md=Sol.size();
    for(int i=0;i<md;i++){
        if(nint(Sol[i][0])==1){//the data is enabled
            wl=Sol[i][1];
            eV=1240./wl;
            auto fd=FDISP(ioptf,eV);
            nVal=fd.sqn;
            kVal=fd.sqk;
            if(nint(par[32][5])<=1){
                j++;
                y=Sol[i][2];
                ery=Sol[i][4];
                if(ery<erynMIN) ery=erynMIN;
                fvec[j]=(y-nVal)/ery;
                chi2=chi2+pow(fvec[j],2.)/fredeg;
            }
            if(nint(par[32][5])==1){
                j++;
                y=Sol[i][3];
                ery=Sol[i][5];
                if(ery<erykMIN) ery=erykMIN;
                fvec[j]=(y-kVal)/ery;
                chi2=chi2+pow(fvec[j],2.)/fredeg;
            }
            if(nint(par[32][5])==2){
                epsR=nVal*nVal-kVal*kVal;
                epsI=2.*nVal*kVal;
                j++;
                ery=2.*(Sol[i][2]*Sol[i][4]+Sol[i][3]*Sol[i][5]);
                if(ery<ere1MIN) ery=ere1MIN;
                fvec[j]=((Sol[i][2]*Sol[i][2]-Sol[i][3]*Sol[i][3])-epsR)/ery;
                chi2=chi2+pow(fvec[j],2.)/fredeg;
                j++;
                ery=2.*(Sol[i][4]*Sol[i][3]+Sol[i][2]*Sol[i][5]);
                if(ery<ere2MIN) ery=ere2MIN;
                fvec[j]=(2.*Sol[i][2]*Sol[i][3]-epsI)/ery;;
                chi2=chi2+pow(fvec[j],2.)/fredeg;
            }
        }
    }
    // monitor values
    qDebug("|FSQ> chi2=%.9g\t",chi2);
//    for(int i=0;i<n;i++)
//        qDebug("p[%d]=%.9g\t",i,x[i]);
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


int ksemawc::FSEM(void *p, int m, int n, const double *x, double *fvec, int iflag){
    //***** Function to fit selected experimental measurements
    ksemawc* self = static_cast<ksemawc*>(p);

    cnk[1].forceMode=0;//n-unknown by the chosen option
    int mwl=NeV;
    int s1p2=nint(par[27][2]);
    int iSpec0Hemi1=nint(par[54][2]);
    double vot[6][3];
    double te=0.;
    // addressing fit parameter onto PM
    int ip;
    for(int i=1;i<=n;i++){
        ip=nint(pmAt(i)[3]);
        pmAt(ip)[1]=x[i-1];
        if(iw>0) qDebug("|FSEM> p[%d]=%f",i,x[i-1]);
    }
    // Sync stack so ASSEMBLER sees the updated layer parameters
    pmValuesToStack(self->stack);
    double chi2=.0;
    double fredeg=static_cast<double>(m>n?m-n:1);//guard m==n (avoid /0)
    //computing fvec
    int ivec=0;
    for(int i=1;i<=mwl;i++){
        par[24][2]=i;
        //computing Tn Rn R1
        if(DATO[1]==2 || DATO[3]==2 || DATO[5]==2){
            if(iSpec0Hemi1==1)
                iShemispherical=true;
            else
                iShemispherical=false;
            self->ASSEMBLER(i,1,0.,vot);
            if(DATO[1]==2){
                fvec[ivec]=(self->ms.measures[0].value[i-1]-vot[1][s1p2])/self->ms.measures[0].error[i-1];
                chi2=chi2+pow(fvec[ivec],2.)/fredeg;
                ivec++;
            }
            if(DATO[3]==2){
                fvec[ivec]=(self->ms.measures[2].value[i-1]-vot[2][s1p2])/self->ms.measures[2].error[i-1];
                chi2=chi2+pow(fvec[ivec],2.)/fredeg;
                ivec++;
            }
            if(DATO[5]==2){
                fvec[ivec]=(self->ms.measures[4].value[i-1]-vot[3][s1p2])/self->ms.measures[4].error[i-1];
                chi2=chi2+pow(fvec[ivec],2.)/fredeg;
                ivec++;
            }
        }
        //computing Tp Rp
        if(DATO[2]==2 || DATO[4]==2){
            te=par[6][1]*deg2rad;
            par[54][2]=0;//Tp and Rp are always direct/specular
            iShemispherical=false;
            self->ASSEMBLER(i,1,te,vot);
            par[54][2]=iSpec0Hemi1;
            if(DATO[2]==2){
                fvec[ivec]=(self->ms.measures[1].value[i-1]-vot[1][s1p2])/self->ms.measures[1].error[i-1];
                chi2=chi2+pow(fvec[ivec],2.)/fredeg;
                ivec++;
            }
            if(DATO[4]==2){
                fvec[ivec]=(self->ms.measures[3].value[i-1]-vot[2][s1p2])/self->ms.measures[3].error[i-1];
                chi2=chi2+pow(fvec[ivec],2.)/fredeg;
                ivec++;
            }
        }
        //computing PSI DELTA
        iShemispherical=false;
        for(int j=1;j<=4;j++){
            if(DATO[7+2*(j-1)]==2 || DATO[8+2*(j-1)]==2){
                te=pmTe[j][1]*deg2rad;
                self->ASSEMBLER(i,2+j,te,vot);
                if(DATO[7+2*(j-1)]==2){//DELTA
                    double offset=nint((self->ms.measures[6+2*(j-1)].value[i-1]-vot[5][2])/Dperiod)*Dperiod;
                    fvec[ivec]=(self->ms.measures[6+2*(j-1)].value[i-1]-vot[5][2]-offset)/self->ms.measures[6+2*(j-1)].error[i-1];
                    chi2=chi2+pow(fvec[ivec],2.)/fredeg;
                    ivec++;
                }
                if(DATO[8+2*(j-1)]==2){//PSI
                    fvec[ivec]=(self->ms.measures[7+2*(j-1)].value[i-1]-vot[5][1])/self->ms.measures[7+2*(j-1)].error[i-1];
                    chi2=chi2+pow(fvec[ivec],2.)/fredeg;
                    ivec++;
                }
            }
        }
    }
    // monitor values
    qDebug("|FSEM> chi2=%.9g\t",chi2);
//    for(int i=0;i<n;i++)
//        qDebug("p[%d]=%.9g\t",i,x[i]);
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


int ksemawc::FRCK(void *p, int m, int n, const double *x, double *fvec, int iflag){
    //***** Subroutine of IbridOne: computing k from T given n
    ksemawc* self = static_cast<ksemawc*>(p);


    cnk[1].forceMode=2;//k-unknown will bet set to forced values
    int mwl=NeV;
    int s1p2=nint(par[27][2]);
    int iSpec0Hemi1=nint(par[54][2]);
    double vot[6][3],VNK[N_CNK_HARD_MAX+1][3];
    double te=0.;
    if(DATO[4]==2)
        te=par[6][1]*deg2rad;
    else if(DATO[8]==2)
        te=pmTe[1][1]*deg2rad;
    else if(DATO[10]==2)
        te=pmTe[2][1]*deg2rad;
    else if(DATO[12]==2)
        te=pmTe[3][1]*deg2rad;
    else if(DATO[14]==2)
        te=pmTe[4][1]*deg2rad;
    // addressing fit parameter onto PM
    int ip;
    for(int i=1;i<=n;i++){
        ip=nint(pmAt(i)[3]);
        pmAt(ip)[1]=x[i-1];
        if(iw>0) qDebug("|FRCK> p[%d]=%f",i,x[i-1]);
    }
    // Sync stack so ASSEMBLER sees the updated layer parameters
    pmValuesToStack(self->stack);
    double chi2=.0;
    double fredeg=static_cast<double>(m>n?m-n:1);//guard m==n (avoid /0)
    //computing fvec
    int ivec=0;
    double wl;
    for(int i=1;i<=mwl;i++){
        wl=self->ms.lambda[i-1];
        par[24][2]=i;
        self->COSVNK(VNK,i);
        self-> nkSol.n[i-1]=VNK[1][1];
        cnk[1].nForced=VNK[1][1];
        cnk[1].kForced=self->nkSol.k[i-1];
        // computing k from T
        //DELTA_k
        double Trasm=1.,dt,dlim,tol;;
        if(nint(par[50+nint(par[53][2])][3])==1){//incoherent layer
            if(DATO[1]>0) Trasm=self->ms.measures[0].value[i-1];
            if(DATO[2]>0) Trasm=self->ms.measures[1].value[i-1];
            double dinco=pmD[nint(par[53][2])][1];
            if(iw==1) qDebug("Trasm = %f",Trasm);
            dt=-0.1*wl/4./3.14/dinco*log(Trasm);
            dlim=1.e-10;//dlim
        }
        else{
            dt=0.001;
            dlim=1.e-5;//dlim
        }
        tol=dt/1000.;//tol
        double k=cnk[1].kForced;
        //FindRoot([this](double k){ return DELTAT(k); },k,dt,dlim,tol);
        self->FindRoot([self](double k) {
                return self->DELTAT(k);
            }, k, dt, dlim, tol);//root is delivered via cnk[1].kForced side effect; return value unused here
        self->nkSol.k[i-1]=cnk[1].kForced;
        //computing R
        if(DATO[3]==2 || DATO[5]==2){
            if(iSpec0Hemi1==1)
                iShemispherical=true;
            else
                iShemispherical=false;
            self->ASSEMBLER(i,1,0.,vot);
            if(DATO[3]==2){
                fvec[ivec]=(self->ms.measures[2].value[i-1]-vot[2][s1p2])/self->ms.measures[2].error[i-1];
                chi2=chi2+pow(fvec[ivec],2.)/fredeg;
                ivec++;
            }
            if(DATO[5]==2){
                fvec[ivec]=(self->ms.measures[4].value[i-1]-vot[3][s1p2])/self->ms.measures[4].error[i-1];
                chi2=chi2+pow(fvec[ivec],2.)/fredeg;
                ivec++;
            }
        }
        par[54][2]=0;//Tp and Rp and ELI are always direct/specular
        iShemispherical=false;
        if(DATO[4]==2){
            self->ASSEMBLER(i,1,te,vot);
            fvec[ivec]=(self->ms.measures[3].value[i-1]-vot[2][s1p2])/self->ms.measures[3].error[i-1];
            chi2=chi2+pow(fvec[ivec],2.)/fredeg;
            ivec++;
        }
        if(DATO[8]==2){
            self->ASSEMBLER(i,1,te,vot);
            fvec[ivec]=(self->ms.measures[7].value[i-1]-vot[5][1])/self->ms.measures[7].error[i-1];
            chi2=chi2+pow(fvec[ivec],2.)/fredeg;
            ivec++;
        }
        if(DATO[10]==2){
            self->ASSEMBLER(i,1,te,vot);
            fvec[ivec]=(self->ms.measures[9].value[i-1]-vot[5][1])/self->ms.measures[9].error[i-1];
            chi2=chi2+pow(fvec[ivec],2.)/fredeg;
            ivec++;
        }
        if(DATO[12]==2){
            self->ASSEMBLER(i,1,te,vot);
            fvec[ivec]=(self->ms.measures[11].value[i-1]-vot[5][1])/self->ms.measures[11].error[i-1];
            chi2=chi2+pow(fvec[ivec],2.)/fredeg;
            ivec++;
        }
        if(DATO[14]==2){
            self->ASSEMBLER(i,1,te,vot);
            fvec[ivec]=(self->ms.measures[13].value[i-1]-vot[5][1])/self->ms.measures[13].error[i-1];
            chi2=chi2+pow(fvec[ivec],2.)/fredeg;
            ivec++;
        }
        par[54][2]=iSpec0Hemi1;
    }
    cnk[1].forceMode=0;//nk-unknown by the chosen option
    // monitor values
    // qDebug("|FRCK> chi2=%g\t",chi2);
    // for(int i=0;i<n;i++)
    //     qDebug("p[%d]=%g\t",i,x[i]);
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


double ksemawc::DELTAT(double k){
    //subroutine for computing the difference Texp and Tcalc
    int s1p2=1;
    int iSpec0Hemi1=par[54][2];
    double texp=0.,tcalc,te=0.,vot[6][3];
    int iwl=nint(par[24][2]);
    if(DATO[1]==2){
        texp=ms.measures[0].value[iwl-1];
        te=0.;
        s1p2=1;
        if(iSpec0Hemi1==1)
            iShemispherical=true;
        else
            iShemispherical=false;
    }
    else if(DATO[2]==2){
        texp=ms.measures[1].value[iwl-1];
        te=par[6][1]*deg2rad;
        s1p2=nint(par[27][2]);
        par[54][2]=0;//Tp and Rp are always direct/specular
        iShemispherical=false;
    }
    //cnk[1].nConst=sqn;
    cnk[1].kForced=k;
    ASSEMBLER(iwl,1,te,vot);
    par[54][2]=iSpec0Hemi1;
    tcalc=vot[1][s1p2];
    double dT=tcalc-texp;
    return(dT);
}


double ksemawc::FindRoot(std::function<double(double)> FT,double t,double dt,double dlim,double tol){
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
        qDebug("\n FindRoot: t=%e dt=%e dlim=%e",t,dt,dlim);
    if(std::abs(d1) <= dlim) istop=1;
    while(istop == 0){
        t=t+dt;
        d1=FT(t);
        if(iWarning==1){
            qDebug("t: %e -> %e dt=%e\td: %e -> %e DELTAd=%e",t0,t,dt,d0,d1,d1-d0);
            //           cout << "To continue press a key & enter: ";
            //           cin  >> ch;
        }
        if(std::abs(d1) <= dlim){
            istop=1;
            //    qDebug("  STOP per raggiunto limite precisione!\n");
        }
        else if(d0*d1 <.0 && std::abs(d1) >= (10.*std::abs(d0))){
            t=t0;
            d1=d0;
            dt=dt/10.;
            //    qDebug("dt decimato = %e\n",dt);
            if(std::abs(dt) <= dlim){
                istop=1;
                //      qDebug(" STOP per raggiunto limite precisione!\n");
            }
        }
        else if(d0*d1 <.0 && std::abs(d1) < (10.*std::abs(d0))) {// Brent
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
                if(std::abs(fc) < std::abs(fb)){
                    a=b;
                    b=c;
                    c=a;
                    fa=fb;
                    fb=fc;
                    fc=fa;
                }
                //      qDebug("\nBrent: N_iteration=%d \n",iter);
                //   "a,b,c=",a,b,c
                //   "fa,fb,fc=",fa,fb,fc
                tol1=2.*eps*std::abs(b)+0.5*tol; // Convergence check.
                xm=.5*(c-b);
                if(std::abs(xm)<=tol1 || fb==0. || iter>=itmax){
                    //	qDebug("eps=%e tol=%e tol1=%e \n xm=%e\n",eps,tol,tol1,xm);
                    t=b;
                    d1=fb;
                    istop=1;
                    //        if(iter<itmax)
                    //          qDebug("STOP per raggiunto limite di precisione\n");
                    //        else
                    //          dDebug("STOP per raggiunto N.max iterazioni\n");
                }
                if(istop==0){
                    if(std::abs(e)>=tol1 && std::abs(fa)>std::abs(fb)){
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
                        pp=std::abs(pp);
                        if(2.*pp<min(3.*xm*q-std::abs(tol1*q),std::abs(e*q))){
                            e=d; //Accept interpolation.
                            d=pp/q;
                            //            qDebug("...quadratic interpolation...\n");
                        }
                        else{
                            //            qDebug("...bisection..\n");
                            d=xm; //Interpolation failed, use bisection.
                            e=d;
                        }
                    }
                    else{ //Bounds decreasing too slowly, use bisection.
                        //          qDebug("...bisection to speed up..\n");
                        d=xm;
                        e=d;
                    }
                    a=b; //Move last best guess to a.
                    fa=fb;
                    if(std::abs(d)>tol1) //Evaluate new trial root.
                        b=b+d;
                    else{
                        if(xm>0)
                            b=b+tol1;
                        else
                            b=b-tol1;
                    }
                    t=b;
                    fb=FT(t);
                    //        qDebug("'b=%e,fb=%e \n",b,fb);
                }
            }
        }
        else if(std::abs(d1)<std::abs(d0)) {
            //      qDebug("  -> better\n");
            t0=t;
            d0=d1;
            irif=0;
            imi=imi+1;
        }
        else if(std::abs(d1)>std::abs(d0) && d0*d1>0.){
            //      qDebug("  -> worst\n");
            t=t0;
            dt=-dt;
            irif=irif+1;
            imi=0;
            if(irif>2){
                if(std::abs(d0-d1)>0.3*dlim){
                    dt=dt/2.;
                    //          qDebug("step->step/2 !\n");
                }
                else
                    irif=ilim+10;
            }
            d1=d0;
        }
        else
            inul=inul+1;
        if(irif>ilim || imi>ilim || inul>ilim){
            //      qDebug("C0: irif o imi o inul > 1000! -> stop\n");
            istop=1;
        }
    }
    d1=FT(t);
    if(iWarning==1)
        qDebug("Root -> t=%e d1=%e",t,d1);
    //   std::cout << "To continue press a key & enter: ";
    //   std::cin  >> ch;
    return(d1);
}


void nextColor(){
    iColor++;
    if(iColor>=7)
        iColor=1;
}


bool MATINV(int n, int np, double **as, double **b){
    /* subroutine for computing the inverse matrix by subroutines given by Numerical Recipe
      n = dimension of the matrix to be inverted
      np = max dimension of matrices used in called subroutine
      as(n,n) = matrix to be inverted
      b(n,n) = inverted matrix
    */
    //qDebug("-> MATINV\n");
    std::vector<int> indx(n);
    std::vector<std::vector<double>> y(np, std::vector<double>(np, 0.0));
    std::vector<std::vector<double>> a(np, std::vector<double>(np, 0.0));
    std::vector<double*> a_ptr(np);
    for(int i=0; i<np; i++) a_ptr[i] = a[i].data();

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
    bool berr=ludcmp(np,a_ptr.data(),n,indx.data(),d);// Decompose the matrix just once.
    if(!berr)
        return(false);
    std::vector<double> yCol(n);
    for(int j=0;j<n;j++){//  Find inverse by columns.
        for(int j2=0;j2<n;j2++)
            yCol[j2]=y[j2][j];
        lubksb(np,a_ptr.data(),n,indx.data(),yCol.data());
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
    return(true);
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
    //qDebug("-> lubksb\n");
    int i,ii,j,ll;
    double sum;
    ii=-1;// sentinel "no nonvanishing element of b seen yet" (0-based: -1, NOT 0, which is a valid index)
    /* When ii is set to a positive value, it will become the index
    of the first nonvanishing element of b. We now do
    the forward substitution, equation (2.3.6). The only new
    wrinkle is to unscramble the permutation as we go.
    */
    for(i=0;i<n;i++){
        ll=indx[i];
        sum=b[ll];
        b[ll]=b[i];
        if(ii>=0){
            for(j=ii;j<i;j++)// FORTRAN do j=ii,i-1 -> 0-based j<i
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


bool ludcmp(int np,double **a,int n,int *indx,double& d){
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
    //qDebug("-> ludcmp\n");
    int i,imax,j,k;
    double aamax,dum,sum;// vv stores the implicit scaling of each row.
    std::vector<double> vv(n);
    d=1.;// No row interchanges yet.
    for(i=0;i<n;i++){// Loop over rows to get the implicit scaling information.
        aamax=0.;
        for(j=0;j<n;j++)
            if(std::abs(a[i][j])>aamax) aamax=std::abs(a[i][j]);
        if(aamax==0.){
            qWarning("singular matrix in ludcmp!!!!");
            return(false);
        }
        vv[i]=1./aamax;//Save the scaling.
    }
    for(j=0;j<n;j++){// This is the loop over columns of Crout's method.
        for(i=0;i<j;i++){// This is equation (2.3.12) except for i=j.  FORTRAN do i=1,j-1 -> 0-based i<j
            sum=a[i][j];
            for(k=0;k<i;k++)// FORTRAN do k=1,i-1 -> 0-based k<i
                sum=sum-a[i][k]*a[k][j];
            a[i][j]=sum;
        }
        aamax=0.;// Initialize for the search for largest pivot element.
        for(i=j;i<n;i++){// This is i=j of equation (2.3.12) and i=j+1...N of equation (2.3.13).
            sum=a[i][j];
            for(k=0;k<j;k++)// FORTRAN do k=1,j-1 -> 0-based k<j
                sum=sum-a[i][k]*a[k][j];
            a[i][j]=sum;
            dum=vv[i]*std::abs(sum);// Figure of merit for the pivot.
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
        if(j!=n-1){//         Now,finally, divide by the pivot element.  FORTRAN if(j.ne.n) -> 0-based j!=n-1
            dum=1./a[j][j];
            for(i=j+1;i<n;i++)
                a[i][j]=a[i][j]*dum;
        }
    }//! Go back for the next column in the reduction.
    return(true);
}
