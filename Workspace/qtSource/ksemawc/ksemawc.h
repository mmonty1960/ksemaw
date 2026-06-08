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

#include <complex>
#include <qwt_plot.h>
#include <qwt_plot_curve.h>
#include <qwt_series_data.h>
#include <gsl/gsl_math.h>
#include <qwt_plot_picker.h>
#include <qwt_picker_machine.h>

#include "ui_ksemawc.h"
class QwtPlotZoomer;
class QwtPlotPicker;
class QwtPlotPanner;
class Plot;
class QPolygon;
using namespace std;

struct PointToFit2 {
    double x;
    double y;
};

// ─── Contants ───────────────────────────────────────────────────────────────
static constexpr int N_MEAS   = 20;
//static constexpr int N_GRAPHS = 16;

// ─── MeasureType ────────────────────────────────────────────────────────────
enum class MeasureType {
    None = 0,
    TransmittanceNormal,
    TransmittancePolarised,
    ReflectanceNormal,
    ReflectancePolarised,
    ReflectanceBack,
    Absorptance,
    Delta1, Psi1,
    Delta2, Psi2,
    Delta3, Psi3,
    Delta4, Psi4
};

// ─── MeasureData ────────────────────────────────────────────────────────────
struct MeasureData {
    MeasureType         type     = MeasureType::None;
    QString             filename;
    double              angle    = 0.0;
    std::vector<double> value;
    std::vector<double> error;

    bool isActive() const { return type != MeasureType::None; }
    bool isValid()  const { return isActive() && !filename.isEmpty() && !value.empty(); }
};

// ─── MeasureSet ─────────────────────────────────────────────────────────────
struct MeasureSet {
    std::vector<double>             lambda;
    std::vector<double>             driftBaseline;
    std::array<MeasureData, N_MEAS> measures;

    int  nPoints() const { return static_cast<int>(lambda.size()); }

    void resize(int Nwl){
        lambda.assign(Nwl, 0.0);
        driftBaseline.assign(Nwl, 0.0);
        for(auto& m : measures){
            m.value.assign(Nwl, 0.0);
            m.error.assign(Nwl, 0.0);
        }
    }
};

// ─── NkData ─────────────────────────────────────────────────────────────────
struct NkData {
    QString             filename;
    std::vector<double> n;
    std::vector<double> k;

    bool isLoaded() const { return !filename.isEmpty() && !n.empty(); }

    void resize(int Nwl){
        n.assign(Nwl, 0.0);
        k.assign(Nwl, 0.0);
    }
};

// ─── NkReference ────────────────────────────────────────────────────────────
struct NkReference {
    QString             filename;
    std::vector<double> lambda;
    std::vector<double> n;
    std::vector<double> k;
    std::vector<double> err_n;
    std::vector<double> err_k;

    bool isLoaded() const { return !filename.isEmpty() && !n.empty(); }
};

// ─── NkSolution ─────────────────────────────────────────────────────────────
struct NkSolution {
    std::vector<double> n;
    std::vector<double> k;

    void resize(int Nwl){
        n.assign(Nwl, 0.0);
        k.assign(Nwl, 0.0);
    }
};

// Results of Effective Medium Approximation
struct EmaResult {
    double n;   // refrative index effective
    double k;   // extinction coefficinet effective
};

// Results oscillators
struct FdispResult {
    double epsiR;  // real part of electrical permittivity
    double epsiI;  // imaginary part of electrical permittivity
    double sqn;    // n computed from complex electrical permittivity
    double sqk;    // k computed from electrical permittivity
};

namespace Ui{
class ksemawc;
}

class ksemawc : public QWidget//, private Ui::ksemawcDLG
{
    Q_OBJECT

private:
    MeasureSet              ms;
    std::vector<NkData>     nkMaterials;
    NkReference             nkRef;
    NkSolution              nkSol;
    std::vector<PointToFit2> pTF2;

    bool m_guiIsMuted = false; // setGuiMute is not active
    void setGuiMute(bool mute);
    void Setnk(int ifile);
    void Clrnk(int ifile);
    void updateMatenk(int, QString);
    void ReadSetting(QString);

    Ui::ksemawc *ui;

    static constexpr int N_GRAPHS = 16;
    QwtPlot*        m_graphs[N_GRAPHS + 1]={};   // index 1..16, [0] not used
    QwtPlotPicker*  m_pickers[N_GRAPHS + 1]={};  // idem

    void ensureGraph(int ic, int iRD);

    QwtPlotPicker* m_picker= nullptr;
    QwtPlotPicker* d_picker= nullptr;
    QwtPickerClickPointMachine* picker_m;

public:
    ksemawc(QWidget *parent=0);    
    QMap<QString, QComboBox*> idToComboBox;
    QMap<QString, QCheckBox*> idToCheckBox;
    QMap<QString, QLineEdit*> idToLineEdit;
    QMap<QString, QDoubleSpinBox*> idToDoubleSpinBox;
    bool isGuiMuted() const { return m_guiIsMuted; }// To know the state of setGuiMute

public Q_SLOTS:
    void PLOTline1bar2(int iL1B2,int iRD,int iCol,int ic,int Ndata,double *Xp,double *Yp,double *ErrXp,double *ErrYp);
    void about();
    void LoadProject();
    void SaveProject();
    void LoadFilenk();
    void ClrFnk();
    void Clrfn();
    void SaveFnk();
    void SaveSetting(int iCall);
    void mDwUp(int iLayer, int Dw1UpM1);
    void callSetSample();
    void setSample();
    void listMeas();
    void pwTn ();
    void pwTp ();
    void pwRn ();
    void pwRp ();
    void pwR1 ();
    void pwApds ();
    void pwEj(int j);
    void pwSubEj(int j);
    void MCRange();
    void listOsc();
    void loadNKmodel();
    void saveNKmodel();
    void tabChanged();
    void setTabSim();
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
    void CalcMis(std::vector<std::vector<double>>& mc);
    void PlotAve();
    void PlotTexturized();
    void searchNK();
    void RefTrackG();
    void NumericalSearch();
    void FitN();
    void FitNK();
    void FitE1E2();
    void FitSelExpMeas();
    void IbridPlotFit();
    void IbridFit();
    void IbridOne();
    void setStop();
    void IbridOneStore();
    void IbridKernel(QString r);
    void ClearTempIbri();
    void StoreFitSet();
    void GoBest();
    void GoPrevious();
    void GoNext();
    void rangeWL();
    void refreshFitPar();
    void saveNKsim();
    void saveSim();
    void PlotAbsEL();
    void setOscN(int k);
    void manageLEwl();
    void setMat(int nM);
    void setEMA(int nM);
    void setRifMir();
    void SPADA();
    void CONVER(double *X,double *Y,int NDATI,int N1,int STEP,int I,int IUVIR);
    void selectNsol();
    void drawPolygon(QPointF pos);
    void setKindOsc(int N);
    void setFontDia();
    void msgErrLoad(QString where, QString fnERR);
    void setRangeEli();
    void readRangeEli();
    void reset();
    void SETVNK(int io,int J,double VNK[17][3],int L,int Kcte);
    void COSVNK(double VNK[17][3],int L);
    void ASSEMBLER(int iwl, int ikind, double teta, double vot[6][3]);
    void BUILDER(int iwl,int ikind,int ifst,int ncoe, std::complex<double> pq,double vosi[6][3]);
    int SOLVE(int imis,int iWL,double *Xp,double *Yp,double *ErrYp);
    double FMER(double k);
    static int FSEM(void *p, int m, int n, const double *x, double *fvec, int iflag);
    static int FRCK(void *p, int m, int n, const double *x, double *fvec, int iflag);
    double DELTAT(double k);
    double FindRoot(std::function<double(double)> f,double t,double dt,double dlim,double tol);
    static int FSQ(void *p, int m, int n, const double *x, double *fvec, int iflag);
private Q_SLOTS:

};

#endif