#ifndef CHACKJKDATA_H
#define CHACKJKDATA_H
#include <string>
#include "TFile.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TGraphAsymmErrors.h"
#include "TMultiGraph.h"
#include "TTree.h"
#include "TColor.h"
#include "TLegend.h"
#include "TF1.h"
#include "TH2D.h"
#include "TCanvas.h"
#include "color_linux.h"

using namespace std;

const std::vector<std::vector<int>> rgb_Petroff = {
    {87, 144, 252},    // #5790fc
    {248, 156, 32},   // #f89c20
    {228, 37, 54},    // #e42536
    {150, 74, 139},   // #964a8b
    {156, 156, 161},  // #9c9ca1
    {122, 33, 221},   // #7a21dd
    //kP6
    {24, 69, 251},
    {255, 94, 2},
    {201, 31, 22 },
    {200, 73, 169 },
    {173, 173, 125 },
    {134, 200, 221 },
    {87, 141, 255},
    {101, 99, 100 },
    //kP8
    {63, 144, 218 },
    {255, 169, 14 },
    {189, 31, 1},
    {148, 164, 162 },
    {131, 45, 182 },
    {169, 107, 89 },
    {231, 99, 0 },
    {185, 172, 112 },
    {113, 117, 129},
    {146, 218, 221} 
    //kP10
  };

const std::vector<std::vector<int>> rgb_chat = {
    {51, 102, 229},  // 蓝色
    {51, 179, 76},   // 绿色
    {204, 102, 26}   // 琥珀
};

int getBranchSize(TTree *t, std::string tree_name);

TMultiGraph *checkjkdata(std::string filename, std::string tree_name);

TGraph *checkjkdata(std::string filename, std::string tree_name, int ien);

TGraphErrors *getAve(string filename, string gname);

double LegendreP(unsigned int n, double x);

double G(int n, double Delta, double R, double r);

// P''_n(x)
double D2P(int n, double x);

// 1 / r * d^2 (r * G_n^{R, Delta}(r) ) / dr^2
double D2rG_r(int n, double Delta, double R, double r);

TGraphAsymmErrors *getFitErrband(TF1 *f, int N = 100);

TList *getGraphList(TFile *f);

void AddTwoSet(double a, string filename1, double b, string filename2, string outputname);

TGraphErrors *DivideG(TGraphErrors *g1, TGraphErrors *g2);

void CleanOld(string name);

class ColorsSchemes
{
  public:
  int kCstart = 103;
  int kCend = 103;
  int kP6Blue = 103;
  int kP8Blue = 109;
  int kP10Blue = 117;
  std::vector<TColor *> kcs;

  ColorsSchemes(std::vector<std::vector<int>> rgb, std::string tit);
  ColorsSchemes(ColorsSchemes &cs);

  ColorsSchemes();

  ~ColorsSchemes();

  void register_custom_colors(std::vector<std::vector<int>> rgb, std::string tit);

  void unregister_custom_colors();

};

class graphdata
{
public:
  TMultiGraph *mg = nullptr;
  int ngraphs = 0;
  TGraphErrors *aveGraph = nullptr;
  TLegend *leg = new TLegend(0.6, 0.4, 0.85, 0.65);

  double a_x = 1;
  double a_y = 1;
  TF1 *f_x = new TF1("f_x", "x", 0, 100);

  int Color = kAzure;
  int Pallete = EColorPalette::kLightTemperature;
  double Alpha = 0.3;
  int MarkerStyle = 20;
  double MarkerSize = 0.5;
  string drawoption = "AP";

  graphdata();
  graphdata(double a_x, double a_y);
  graphdata(double a_x, double a_y, string drawoption);
  graphdata(const graphdata &other);
  ~graphdata();

  TLegend *resetLeg(vector<double> Locate={0.6, 0.4, 0.85, 0.65}, double TextSize = 0.05, int TextFont = 42, string suffix = "");

  void setAxAy(double ax, double ay);

  void setFormula(string formula);

  void setColor(int color);

  void setAlpha(double alpha);

  void setPallete(int pallete);

  void setMarkerStyle(int style);

  void setMarkerSize(double size);

  void setDrawOption(string drawopt);

  void makeUP(TGraph *g, int Ng, int ien) const;

  void makeUP(TGraph *g, int ien = -1) const;

  TGraph *getiGraph(int ien) const;

  TMultiGraph *getAllGraph() const;

  TGraphErrors *getAveGraph() const;
 };

class jkdata : public graphdata
{
public:
  jkdata();

  jkdata(string filename, string plotname);
  jkdata(string filename, string plotname, double a_x, double a_y);
  jkdata(const jkdata &other);

  void init(string filename, string plotname);

  ~jkdata();

  jkdata Clone() const;

  void GetPoint(int i, double &x, vector<double> &y, bool ifunited = true) const;

  void AddPoint(double x, vector<double> y, bool ifunited = true) const;

  void RemovePoint(int i);

  TGraphErrors *getiGraphErrors(int ien);

  jkdata operators(const jkdata &other, int oper) const;

  jkdata operators(const double oy, int oper) const;

  jkdata operators(const TGraph *oG, int oper) const;

  jkdata operators(TF1 *f, int oper) const;

  jkdata operator=(const jkdata &other);

  jkdata operator+(const jkdata &other) const;

  jkdata operator-(const jkdata &other) const;

  jkdata operator*(const jkdata &other) const;

  jkdata operator/(const jkdata &other) const;

  jkdata operator+(const double other) const;

  jkdata operator-(const double other) const;

  jkdata operator*(const double other) const;

  jkdata operator/(const double other) const;

  jkdata operator+(const TGraph *other) const;

  jkdata operator-(const TGraph *other) const;

  jkdata operator*(const TGraph *other) const;

  jkdata operator/(const TGraph *other) const;

  jkdata operator+(TF1 *f) const;

  jkdata operator-(TF1 *f) const;

  jkdata operator*(TF1 *f) const;

  jkdata operator/(TF1 *f) const;

  void print(string filename) const;
  
};

class timedata : public graphdata
{
public:
  timedata(string filename_model, int time_start, int time_end, string plotname, string drawopt = "AP");
  timedata(string filename_model, int time_start, int time_end, string plotname, double a_x, double a_y, string drawopt = "AP");
  timedata(string filename_model, vector<int> t_list, string plotname, string drawopt = "AP");
  timedata(string filename_model, vector<int> t_list, string plotname, double a_x, double a_y, string drawopt = "AP");
  timedata(vector<jkdata> &datalist, vector<int> t_list, string plotname, string drawopt = "AP");
  timedata(vector<jkdata> &datalist, vector<int> t_list, string plotname, double a_x, double a_y, string drawopt = "AP");
  ~timedata();

  void init(string filename_model, string plotname, int time_start, int time_end);
  void init(string filename_model, string plotname, vector<int> t_list);
  void init(const vector<jkdata> &datalist, string plotname, vector<int> t_list);

  TGraphErrors *getiGraph(int ien);

  TGraphErrors *getSlide(int ix);

  void GetPoint(int i, double &x, vector<double> &y, vector<double> &ey, bool ifunited = true);

  void AddPoint(double x, vector<double> y, vector<double> ey, bool ifunited = true);

  void RemovePoint(int i);
};

void chackjkdata();

class jkfitdata : public jkdata
{
private:
  TGraph *actfitfunc = nullptr;
  TMultiGraph *actmf = nullptr;
  double act_a_x;
  double act_a_y;
  double calfitfunc(double x, int ien);
  TGraph *FitGraph(int ien, int N = 100);
  
public:
  // TList *fitlist = new TList();
  TF1 *fitfunc;
  vector<double> ave_chi2dof;
  vector<double> chi2dof_jk;
  TGraphAsymmErrors *fitgraph = nullptr;
  TGraphErrors *mf = nullptr;
  string fitdrawopt = "R";
  vector<vector<double>> cov_par;
  vector<vector<double>> cov_data;
  vector<int> pick_ip;

  jkfitdata(string tar, string plotname, string func_name = "fit_func", string opt = "R");

  jkfitdata(string tar, string plotname, double a_x, double a_y, string func_name = "fit_func", string opt = "R");
  jkfitdata(const jkfitdata &other);

  jkfitdata(const jkdata &other);

  ~jkfitdata();

  void getFitFunc(string filename, string func_name);

  TMultiGraph *GetFitGraph(int ien = -1, int color = kRed);

  TF1 *GetiFit(int ien = -1) const;

  vector<double> GetChi2Dof(int ien = -1) const;

  TGraphAsymmErrors *GetFitErrorBand(bool ori = true);

  void makeUPf(TF1 *f, int ien = -1);

  void setAxAy(double ax, double ay);

  void setFormula(string formula);

  jkfitdata operator=(const jkfitdata &other);

  jkfitdata operator=(const jkdata &other);

  jkfitdata Clone() const;

  jkfitdata operators(TF1 *f, int oper);

  void AddPoint(double x, vector<double> y, bool ifunited = true);

  void RemovePoint(int i);

  jkfitdata operators(const double other, int oper);

  jkfitdata operator+(TF1 *f);
  jkfitdata operator-(TF1 *f);
  jkfitdata operator*(TF1 *f);
  jkfitdata operator/(TF1 *f);

  jkfitdata operator+(const double other);
  jkfitdata operator-(const double other);
  jkfitdata operator*(const double other);
  jkfitdata operator/(const double other);

  static jkfitdata FitS(jkdata &data, TF1 *f, string opt = "");
  void Fit(TF1 *f);
  /**
   * @brief Perform a correlated (or optionally uncorrelated) fit on jackknife data.
   *
   * This function fits the provided TF1 function to the averaged jackknife data,
   * taking into account the correlations among data points through the covariance matrix.
   * It supports several options for fitting range, minimizer choice, and fit mode.
   *
   * The fit can be performed in three modes:
   * - Correlated fit (default): uses the full covariance matrix.
   * - Uncorrelated fit (opt::U): uses only the diagonal elements of the covariance matrix.
   * - Ignore-error fit (opt::W): treats all errors as equal (identity covariance matrix).
   *
   * Additional options can control fitting behavior:
   * - "R": draw range with the specified range of the function not for fit.
   * - "Q": quiet mode, suppressing printout during the fit.
   * - "S[n]": select Nstep fitting points within the range (default is nentries/2).
   * - "M[i]": choose minimizer type:
   *      - 1 = Minuit
   *      - 2 = Minuit2 (default)
   *      - 3 = GSLMultiMin
   *      - 4 = GSLMultiFit
   *      - 5 = GSLSimAn
   *
   * @param data Reference to the jackknife data container (jkdata) to fit.
   * @param f Pointer to the TF1 function to be fitted.
   * @param opt Optional string containing fit options as described above.
   *
   * @return A jkfitdata object containing the fit results, including:
   * - Fitted parameters and their uncertainties
   * - Chi2/dof of the fit
   * - Cloned function and graph objects for further inspection or plotting
   *
   * @note Currently, this function provides an overview of the fit result.
   *       It can be extended to estimate phase shift errors via Monte Carlo sampling
   *       based on the covariance matrix.
   *
   * @warning Ensure that the number of measurements (nentries) is larger than
   *          the number of data points to avoid ill-conditioned covariance matrices.
   */
  static jkfitdata CorrelatedFitS(jkdata &data, TF1 *f, string opt = "");
  void CorrelatedFit(TF1 *f);

  vector<vector<double>> GetDataConvMat(bool is_norm = false);
  void SetDataConvMat(vector<vector<double>> cov, bool from_norm = false, vector<int> pick = {});
  vector<vector<double>> GetFitConvMat();

  vector<vector<double>> GetjkFitParameters() const;

  static jkfitdata Construct(jkdata &data, TF1 *f, vector<vector<double>> par, string opt);

  static jkfitdata Construct(jkdata &data, TF1 *f, jkfitdata &fitdata, string opt = "R");

  void setFitDrawOpt(string opt);

  void print(string filename) const;
  void print() const;
};

struct CorrMetrics
{
  double mean_absR;          // 平均绝对相关系数
  double max_absR;           // 最大绝对相关系数
  double frob_offdiag_ratio; // 非对角项占比
  double cond_number;        // 条件数 (SVD)
  double effective_rank;     // 有效秩
};

vector<vector<double>> CovToCorr(const vector<vector<double>> &Co);
CorrMetrics ComputeCorrMetrics(const vector<vector<double>> &C);
vector<vector<double>> GetCov(const std::string &filename, bool is_norm);
vector<vector<double>> trimCov(const vector<vector<double>> & cov_o, vector<int> pick_ip);

void CompareCorrelation(const std::string &filename1,
                        const std::string &filename2,
                        bool is_norm);
void CompareCorrelation(const std::vector<std::vector<double>> &C_ori,
                        const std::vector<std::vector<double>> &yCsub);
TH2D *DirectComp(const std::string &filename1,
                 const std::string &filename2,
                 bool is_norm, vector<int> pick_ip = {});
TH2D *DirectComp(const std::vector<std::vector<double>> &C1,
                const std::vector<std::vector<double>> &C2, vector<int> pick_ip = {});
// name is for give an unique identity for TH2D also a simple control for data
// like name =="xxx abs" for draw a abs value for matrix
TH2D *DrawCorrHeatmap(const vector<vector<double>> &Co, const std::string &name = "corr", vector<int> pick_ip = {});

#endif