#include "chackjkdata.h"
#include "TKey.h"
#include "TROOT.h"

using namespace std;

int getBranchSize(TTree *t, string tree_name)
{
  int datasize;
  TBranch *b = t->GetBranch(tree_name.c_str());
  if (!b)
  {
    printf("No branch %s found in tree %s\n", tree_name.c_str(), t->GetName());
    return -1;
  }
  string b_tit = b->GetTitle();
  size_t start = b_tit.find('[');
  size_t end = b_tit.find(']');

  if (start != std::string::npos && end != std::string::npos && end > start)
  {
    // 提取方括号中的内容
    std::string numberStr = b_tit.substr(start + 1, end - start - 1);
    datasize = std::stoi(numberStr);
  }
  else
  {
    std::cout << "No array size found" << std::endl;
  }
  return datasize;
}

TMultiGraph *checkjkdata(std::string filename, std::string tree_name)
{
  int L = 96;
  int ReIm = 2; // 0 for real part, 1 for imaginary part
  if(filename.find("observables") != std::string::npos)
    ReIm = 1; // only real
  TFile *f = new TFile(filename.c_str());
  TTree *t = (TTree *)f->Get("jackknifed_data");
  if (!t)
  {
    printf("No 'jackknifed_data' tree found in file %s\nTry to find '%s_tree'\n", filename.c_str(), tree_name.c_str());
    t = (TTree *)f->Get((tree_name + "_tree").c_str());
    if (!t)
    {
      printf("Still No '%s_tree' tree found in file %s\n", tree_name.c_str(), filename.c_str());
      return nullptr;
    }
    printf("Found '%s_tree' tree in file %s!\n", tree_name.c_str(), filename.c_str());
  }
  TGraphErrors *gr = (TGraphErrors *)(f->Get(tree_name.c_str()));
  if (!gr)
  {
    printf("No graph %s found in file %s\n", tree_name.c_str(), filename.c_str());
    return nullptr;
  }

  int datasize = gr->GetN();
  int flg_mis = 0;
  if (datasize != (48 + 44 * L + 12 * L * L + L * L * L) / 48)
    flg_mis = 1;
  double data[datasize * ReIm];
  t->SetBranchAddress(tree_name.c_str(), data);
  int Nentry = t->GetEntries();
  TMultiGraph *mg = new TMultiGraph();
  // for (int i = 0; i < 1; i++)
  for (int i = 0; i < Nentry; i++)
  {
    t->GetEntry(i);

    TGraph *g = new TGraph(datasize);
    if (flg_mis == 1)
      for (int id = 0; id < datasize; id++)
      {
        double r = gr->GetPointX(id);
        g->SetPoint(id, r, data[ReIm * id]);
      }
    else
      for (int iz = 0; iz < L / 2 + 1; iz++)
        for (int iy = iz; iy < L / 2 + 1; iy++)
          for (int ix = iy; ix < L / 2 + 1; ix++)
          {
            int id = (6 * ix - 3 * iy * iy +
                      iy * (-3 + 6 * (L / 2 + 1)) +
                      iz * (-1 + iz * iz - 3 * iz * (L / 2 + 1) +
                            3 * pow(L / 2 + 1, 2))) /
                     6;
            double r = sqrt(ix * ix + iy * iy + iz * iz);
            g->SetPoint(id, r, data[2 * id]);
          }
    g->SetMarkerColor(i - 9 + 860);
    g->SetMarkerStyle(20);
    g->SetMarkerSize(0.5);
    mg->Add(g);
  }
  f->Close();
  delete f;
  return mg;
}

TGraph *checkjkdata(std::string filename, std::string tree_name, int ien)
{
  int L = 96;
  TFile *f = new TFile(filename.c_str());
  TTree *t = (TTree *)f->Get("jackknifed_data");
  if (!t)
  {
    printf("No jackknifed_data tree found in file %s\n", filename.c_str());
    return nullptr;
  }
  TGraphErrors *gr = (TGraphErrors *)(f->Get(tree_name.c_str()));
  if (!gr)
  {
    printf("No graph %s found in file %s\n", tree_name.c_str(), filename.c_str());
    return nullptr;
  }

  int datasize = gr->GetN();
  int flg_mis = 0;
  if (datasize != (48 + 44 * L + 12 * L * L + L * L * L) / 48)
    flg_mis = 1;

  double data[datasize * 2];
  t->SetBranchAddress(tree_name.c_str(), data);
  int Nentry = t->GetEntries();
  if (ien >= Nentry || ien < 0)
  {
    printf("Invalid entry number %d, should be between -1 and %d\n", ien, Nentry - 1);
    return nullptr;
  }
  TGraph *g = new TGraph(datasize);
  // for (int i = 0; i < 1; i++)
  {
    t->GetEntry(ien);

    if (flg_mis == 1)
      for (int id = 0; id < datasize; id++)
      {
        double r = gr->GetPointX(id);
        g->SetPoint(id, r, data[2 * id]);
      }
    else
      for (int iz = 0; iz < L / 2 + 1; iz++)
        for (int iy = iz; iy < L / 2 + 1; iy++)
          for (int ix = iy; ix < L / 2 + 1; ix++)
          {
            int id = (6 * ix - 3 * iy * iy +
                      iy * (-3 + 6 * (L / 2 + 1)) +
                      iz * (-1 + iz * iz - 3 * iz * (L / 2 + 1) +
                            3 * pow(L / 2 + 1, 2))) /
                     6;
            double r = sqrt(ix * ix + iy * iy + iz * iz);
            g->SetPoint(id, r, data[2 * id]);
          }
    g->SetMarkerColor(ien - 9 + 860);
    g->SetMarkerStyle(20);
    g->SetMarkerSize(0.5);
  }
  f->Close();
  delete f;
  return g;
}

TGraphErrors *getAve(string filename, string gname)
{
  TFile *f = new TFile(filename.c_str());
  f->Get(gname.c_str())->Draw("AP");
  TGraphErrors *g = new TGraphErrors(*(TGraphErrors *)(f->Get(gname.c_str())));
  g->SetMarkerColor(1);
  g->SetLineColor(1);
  g->Draw("same *");
  if (!g)
  {
    printf("No graph %s found in file %s\n", gname.c_str(), filename.c_str());
    return nullptr;
  }
  g->SetMarkerColor(860);
  g->SetLineColor(860);
  g->SetMarkerStyle(20);
  g->SetMarkerSize(0.5);
  delete f;
  return g;
}

double LegendreP(unsigned int n, double x)
{
  if (n == 0)
    return 1.0;
  if (n == 1)
    return x;

  double P0 = 1;
  double P1 = x;

  for (int i = 2; i <= n; ++i)
  {
    double Pn = ((2.0 * i - 1.0) * x * P1 - (i - 1.0) * P0) / i; // 递推关系
    P0 = P1;
    P1 = Pn;
  }

  return P1;
}

double G(int n, double Delta, double R, double r)
{
  double x = (r - R) / Delta;
  // if(abs(x)>1)
  //	return 0;
  // return 1. / r * sqrt((2. * n + 1) / Delta) * ROOT::Math::legendre(n, x);
  return 1. / r * sqrt((2. * n + 1) / 2. / Delta) * LegendreP(n, x);
}

// P''_n(x)
double D2P(int n, double x)
{
  if (abs(abs(x) - 1) < 1.e-60)
    return pow(x, n) * n * (n + 1.) * (n * n + n - 2.);
  return (1. + n) * ((2 * x * x + n * (x * x - 1.)) * LegendreP(n, x) - 2 * x * LegendreP(n + 1, x)) / pow(x * x - 1., 2);
}

// 1 / r * d^2 (r * G_n^{R, Delta}(r) ) / dr^2
double D2rG_r(int n, double Delta, double R, double r)
{
  return 1. / (r * Delta * Delta) * sqrt((2. * n + 1) / Delta / 2.) * D2P(n, (r - R) / Delta);
}

TGraphAsymmErrors *getFitErrband(TF1 *f, int N)
{
  TGraphAsymmErrors *graph = new TGraphAsymmErrors(N);
  double dx = (f->GetXmax() - f->GetXmin()) / (N - 1);
  int Npar = f->GetNpar();
  double paras[Npar];
  f->GetParameters(paras);
  const double *perr = f->GetParErrors();
  for (int i = 0; i < N; i++)
  {
    TF1 *f_copy = new TF1(*f);
    double x = f->GetXmin() + i * dx;
    double y = f->Eval(x);

    double y_up = y, y_down = y;
    for (int j = 0; j < pow(2, Npar); j++)
    {
      vector<int> sign(Npar, -1);
      for (int k = 0; k < Npar; k++)
        sign[k] = sign[k] + 2 * ((j >> k) & 1);
      for (int k = 0; k < Npar; k++)
        f_copy->SetParameter(k, paras[k] + sign[k] * perr[k]);
      double y_now = f_copy->Eval(x);
      if (y_now > y_up)
        y_up = y_now;
      if (y_now < y_down)
        y_down = y_now;
    }
    delete f_copy;

    graph->SetPoint(i, x, y);
    graph->SetPointError(i, 0, 0, (y - y_down), (y_up - y));
  }
  graph->SetTitle(f->GetTitle());
  graph->SetName(Form("%s_ErrBand", f->GetName()));
  graph->SetFillColorAlpha(6, 0.3);
  // graph->SetFillStyle(3005);
  graph->SetDrawOption("A4");
  return graph;
}

TList *getGraphList(TFile *f)
{
  TList *glist = new TList();
  TList *list = f->GetListOfKeys();
  for (int i = 0; i < list->GetSize(); i++)
  {
    TKey *key = (TKey *)list->At(i);
    if (string(key->ReadObj()->ClassName()) == "TGraphErrors")
    {
      TGraph *g = (TGraph *)key->ReadObj();
      glist->Add(g);
    }
  }
  return glist;
}

void AddTwoSet(double a, string filename1, double b, string filename2, string outputname)
{
  int L = 96;
  TFile *f1 = new TFile(filename1.c_str());
  TFile *f2 = new TFile(filename2.c_str());
  TTree *t1 = (TTree *)f1->Get("jackknifed_data");
  TTree *t2 = (TTree *)f2->Get("jackknifed_data");
  TList *glist1 = getGraphList(f1);
  TList *glist2 = getGraphList(f2);
  int ng = glist1->GetSize();
  int ndata = ((TGraph *)(glist1->At(0)))->GetN();
  int flg_mis = 0;
  if (ndata != (48 + 44 * L + 12 * L * L + L * L * L) / 48)
    flg_mis = 1;
  int Njack = t1->GetEntries();
  if (ng != glist2->GetSize())
  {
    cERR << "Error: two files have different number of graphs!" << endl;
    return;
  }
  if (Njack != t2->GetEntries())
  {
    cERR << "Error: two files have different number of entries!" << endl;
    return;
  }
  if (ndata != ((TGraphErrors *)(glist2->At(0)))->GetN())
  {
    cERR << "Error: two files have different number of data points!" << endl;
    return;
  }

  TFile *fout = new TFile(outputname.c_str(), "recreate");
  TTree *tout = new TTree("jackknifed_data", "jackknifed_data");
  tout->SetEntries(Njack);
  for(int ig = 0; ig < ng; ++ig)
  {
    TGraphErrors *g1 = (TGraphErrors *)(glist1->At(ig));
    TGraphErrors *g2 = (TGraphErrors *)(glist2->At(ig));
    TGraphErrors *gout = new TGraphErrors(ndata);
    gout->SetNameTitle(g1->GetName(), g1->GetTitle());

    double jkdata1[ndata * 2];
    double jkdata2[ndata * 2];
    double jkdataout[ndata * 2];
    t1->SetBranchAddress(g1->GetName(), jkdata1);
    t2->SetBranchAddress(g1->GetName(), jkdata2);
    TBranch *bn = tout->Branch(g1->GetName(), jkdataout, Form("ijkdata[%d]/D", 2 * ndata));
    vector<double> ave(ndata * 2, 0);
    vector<double> err(ndata * 2, 0);
    for (int i = 0; i < Njack; i++)
    {
      t1->GetEntry(i);
      t2->GetEntry(i);
      for (int j = 0; j < ndata * 2; j++)
      {
        double pave = (jkdata1[j] *a + jkdata2[j]*b);
        jkdataout[j] = pave;
        ave[j] += pave / Njack;
        err[j] += pave * pave / Njack;
      }
      bn->Fill();
    }
    if (flg_mis == 1)
    {
      for (int j = 0; j < ndata; j++)
      {
        err[2 * j] = sqrt((err[2 * j] - ave[2 * j] * ave[2 * j]) * (Njack - 1));
        gout->SetPoint(j, g1->GetX()[j], ave[2 * j]);
        gout->SetPointError(j, 0, err[2 * j]);
      }
    }
    else
    {
      for (int iz = 0; iz < L / 2 + 1; iz++)
        for (int iy = iz; iy < L / 2 + 1; iy++)
          for (int ix = iy; ix < L / 2 + 1; ix++)
          {
            int j = (6 * ix - 3 * iy * iy +
                      iy * (-3 + 6 * (L / 2 + 1)) +
                      iz * (-1 + iz * iz - 3 * iz * (L / 2 + 1) +
                            3 * pow(L / 2 + 1, 2))) /
                     6;
            double r = sqrt(ix * ix + iy * iy + iz * iz);
            err[2 * j] = sqrt((err[2 * j] - ave[2 * j] * ave[2 * j]) * (Njack - 1));
            gout->SetPoint(j, r, ave[2 * j]);
            gout->SetPointError(j, 0, err[2 * j]);
          }
      }
      gout->Write();
  }
  tout->Write();
  fout->Close();
  delete f1;
  delete f2;
  delete fout;
  delete glist1;
  delete glist2;
}

TGraphErrors *DivideG(TGraphErrors *g1, TGraphErrors *g2)
{
  if (g1->GetN() != g2->GetN())
  {
    cERR << "Error: two graphs have different number of data points!" << endl;
    return nullptr;
  }
  int ndata = g1->GetN();
  TGraphErrors *gout = new TGraphErrors(ndata);
  gout->SetName("divide");
  for (int i = 0; i < ndata; i++)
  {
    double x = g1->GetX()[i];
    double y1 = g1->GetY()[i];
    double ey1 = g1->GetEY()[i];
    int i2 = -1;
    for (int j = 0; j < ndata; j++)
    {
      if (abs(g2->GetX()[j] - x) < 1.e-3)
      {
        i2 = j;
        break;
      }
      if (j == ndata - 1)
      {
        cERR << "Error: two graphs have different x values! for x = " << x << endl;
        delete gout;
        return nullptr;
      }
    }
    double y2 = g2->GetY()[i2];
    double ey2 = g2->GetEY()[i2];
    if (y2 == 0)
    {
      cWARN << "Divide by zero at x = " << x << endl;
      y2 = 1.e-30;
    }
    double y = y1 / y2;
    double ey = sqrt(pow(ey1 / y2, 2) + pow(ey2 * y1 / y2 / y2, 2));
    gout->SetPoint(i, x, y);
    gout->SetPointError(i, 0, ey);
  }
  gout->SetMarkerSize(0.5);
  gout->SetMarkerStyle(20);
  gout->SetMarkerColor(kTeal);
  gout->SetLineColor(kViolet);
  return gout;
}

void CleanOld(string name)
{
  TObject * old = gROOT->FindObject(name.c_str());
  if(old)
    delete old;
}

void chackjkdata()
{

  TCanvas *c = new TCanvas("c", "c", 800, 600);
  // jkdata *jk = new jkdata("output/set1/pot_3s1_cen_FconfNLmdc_t008_1600conf_080bin.root", "Rcorr_t07");
  timedata *jk = new timedata("output/set1/pot_3s1_cen_FconfNOmgccc_t%03d_1600conf_080bin.root", 8, 15, "pot");
  TMultiGraph *mg = jk->getAllGraph();
  if (!mg)
  {
    cerr << "No data here!" << endl;
    return;
  }
  // mg->GetXaxis()->SetRangeUser(-1, 10);
  mg->Draw("AP");
}