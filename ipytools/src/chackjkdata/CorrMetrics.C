#include "chackjkdata.h"
#include "TROOT.h"
#include <TMatrixDSym.h>
#include <Eigen/Dense>

using namespace std;

vector<vector<double>> CovToCorr(const vector<vector<double>> &Co) {
  int N = Co.size();
  vector<double> diag_sqrt(N);
  vector<vector<double>> R(N, vector<double>(N));

  for (int i = 0; i < N; i++)
    diag_sqrt[i] = sqrt(Co[i][i]);

  for (int i = 0; i < N; i++)
    for (int j = 0; j < N; j++)
    {  
      R[i][j] = (diag_sqrt[i] * diag_sqrt[j] != 0) ? Co[i][j] / (diag_sqrt[i] * diag_sqrt[j]) : 0;
    }
  return R;
}

Eigen::MatrixXd CovToCorr(const vector<vector<double>> &Co, Eigen::MatrixXd &C) {
  int N = Co.size();
  Eigen::VectorXd diag_sqrt(N);
  Eigen::MatrixXd R(N, N);
  C.resize(N, N);

  for (int i = 0; i < N; i++)
    diag_sqrt(i) = sqrt(Co[i][i]);

  for (int i = 0; i < N; i++)
    for (int j = 0; j < N; j++)
    {  
      R(i, j) = (diag_sqrt[i] * diag_sqrt[j] != 0) ? Co[i][j] / (diag_sqrt[i] * diag_sqrt[j]) : 0;
      C(i, j) = Co[i][j];
    }
  return R;
}

CorrMetrics ComputeCorrMetrics(const vector<vector<double>> &Co) {
    CorrMetrics M;
    int N = Co.size();
    Eigen::MatrixXd C;
    Eigen::MatrixXd R = CovToCorr(Co, C);

    double sum_abs = 0;
    double max_abs = 0;
    int count = 0;
    for (int i=0;i<N;i++)
        for (int j=0;j<N;j++) {
            if (i==j) continue;
            double a = std::abs(R(i,j));
            sum_abs += a;
            if (a > max_abs) max_abs = a;
            count++;
        }
    M.mean_absR = (count>0) ? sum_abs/count : 0;
    M.max_absR  = max_abs;

    // Frobenius off-diagonal ratio
    double frob_total = C.norm();
    double frob_diag = C.diagonal().norm();
    double frob_off = std::sqrt(std::max(0.0, frob_total*frob_total - frob_diag*frob_diag));
    M.frob_offdiag_ratio = (frob_total != 0) ? frob_off/frob_total : 0;

    // SVD
    Eigen::JacobiSVD<Eigen::MatrixXd> svd(C, Eigen::ComputeThinU | Eigen::ComputeThinV);
    Eigen::VectorXd s = svd.singularValues();
    double smax = s.maxCoeff();
    double smin = s.minCoeff();
    M.cond_number = (smin>0) ? smax/smin : std::numeric_limits<double>::infinity();

    double sum_s = s.sum();
    double sum_s2 = s.array().square().sum();
    M.effective_rank = (sum_s2>0) ? (sum_s*sum_s)/sum_s2 : 0;

    return M;
}

vector<vector<double>> GetCov(const std::string &filename, bool is_norm)
{
  TFile *f = new TFile(filename.c_str());
  TMatrixDSym *cov = (TMatrixDSym *)f->Get("cov_data");
  if (!cov)
  {
    cERR << "There is no cov matrix in" << filename << endl;
    return vector<vector<double>>();
  }

  vector<vector<double>> cov_data(cov->GetNcols(), vector<double>(cov->GetNrows()));
  vector<double> D(cov->GetNcols());
  for (int i = 0; i < cov->GetNcols(); i++)
  {
    D[i] = sqrt((*cov)(i, i));
  }
  for (int i = 0; i < cov->GetNcols(); i++)
  {
    for (int j = 0; j < cov->GetNrows(); j++)
    {
      if (is_norm)
        cov_data[i][j] = (*cov)(i, j) / D[i] / D[j];
      else
        cov_data[i][j] = (*cov)(i, j);
    }
  }

  f->Close();
  delete f;
  return cov_data;
}

vector<vector<double>> trimCov(const vector<vector<double>> & cov_o, vector<int> pick_ip)
{
  int N = pick_ip.size();
  vector<vector<double>> cov_t(N, vector<double>(N));
  for (int i = 0; i < N; i++)
  {
    int ip = pick_ip[i];
    for (int j = 0; j < N; j++)
    {
      int jp = pick_ip[j];
      cov_t[i][j] = cov_o[ip][jp];
    }
  }
  return cov_t;
}

void CompareCorrelation(const std::string &filename1,
                        const std::string &filename2,
                        bool is_norm)
{
  vector<vector<double>> cov1 = GetCov(filename1, is_norm);
  vector<vector<double>> cov2 = GetCov(filename2, is_norm);
  CompareCorrelation(cov1, cov2);
}

void CompareCorrelation(const std::vector<std::vector<double>> &C_ori,
                        const std::vector<std::vector<double>> &C_sub)
{
  auto metrics_ori = ComputeCorrMetrics(C_ori);
  auto metrics_sub = ComputeCorrMetrics(C_sub);

  std::cout << "======= Original Data =======\n";
  std::cout << "mean|R|=" << metrics_ori.mean_absR
            << "  max|R|=" << metrics_ori.max_absR
            << "  frob_off/total=" << metrics_ori.frob_offdiag_ratio
            << "  cond=" << metrics_ori.cond_number
            << "  eff_rank=" << metrics_ori.effective_rank
            << std::endl;

  std::cout << "======= Subset Data =======\n";
  std::cout << "mean|R|=" << metrics_sub.mean_absR
            << "  max|R|=" << metrics_sub.max_absR
            << "  frob_off/total=" << metrics_sub.frob_offdiag_ratio
            << "  cond=" << metrics_sub.cond_number
            << "  eff_rank=" << metrics_sub.effective_rank
            << std::endl;
}

TH2D *DirectComp(const std::string &filename1,
                 const std::string &filename2,
                 bool is_norm, vector<int> pick_ip)
{
  vector<vector<double>> cov1 = GetCov(filename1, is_norm);
  vector<vector<double>> cov2 = GetCov(filename2, is_norm);
  return DirectComp(cov1, cov2, pick_ip);
}

TH2D *DirectComp(const std::vector<std::vector<double>> &C1,
                 const std::vector<std::vector<double>> &C2, vector<int> pick_ip)
{
  int N = C1.size();
  if (N != C2.size())
  {
    cERR << "the two covariance have differnet size." << endl;
    return nullptr;
  }
  Eigen::MatrixXd C_1;
  Eigen::MatrixXd R_1 = CovToCorr(C1, C_1);
  Eigen::MatrixXd C_2;
  Eigen::MatrixXd R_2 = CovToCorr(C2, C_2);
  auto x = gROOT->FindObject("compare_cov");
  if (x)
    delete x;

  bool use_pick = (pick_ip.size() != 0);
  int Ndata = use_pick ? pick_ip.size() : N;
  TH2D *h = new TH2D("compare_cov", "compare correlated matrix: |R1-R2|", Ndata, 0, Ndata, Ndata, 0, Ndata);
  for (int i = 0; i < Ndata; i++)
  {
    int ip = use_pick ? pick_ip[i] : i;
    for (int j = 0; j < Ndata; j++)
    {
      int jp = use_pick ? pick_ip[j] : j;
      h->SetBinContent(i + 1, j + 1, abs((R_1(ip, jp) - R_2(ip, jp)) / (R_1(ip, jp) + R_2(ip, jp))));
    }
  }
  h->SetDrawOption("COLZ");
  h->SetMaximum(1);
  h->SetMinimum(0);
  return h;
}

TH2D *DrawCorrHeatmap(const vector<vector<double>> &Co, const std::string &name, vector<int> pick_ip)
{
  Eigen::MatrixXd C;
  Eigen::MatrixXd R = CovToCorr(Co, C);
  int N = R.rows();
  auto x = gROOT->FindObject(name.c_str());
  if (x)
    delete x;
  bool use_pick = (pick_ip.size() != 0);
  int Ndata = use_pick ? pick_ip.size() : N;
  TH2D *h = new TH2D(name.c_str(), "Correlation matrix; i; j", Ndata, 0, Ndata, Ndata, 0, Ndata);
  for (int i = 0; i < Ndata; i++)
  {
    int ip = use_pick ? pick_ip[i] : i;
    for (int j = 0; j < Ndata; j++)
    {
      int jp = use_pick ? pick_ip[j] : j;
      if (name.find("abs") != string::npos)
        h->SetBinContent(i + 1, j + 1, abs(R(ip, jp)));
      else
        h->SetBinContent(i + 1, j + 1, R(ip, jp));
    }
  }
  h->SetDrawOption("COLZ");
  return h;
}