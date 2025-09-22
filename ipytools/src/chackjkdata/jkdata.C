#include "chackjkdata.h"
#include <fstream>

using namespace std;
//class jkdata : public graphdata
  jkdata::jkdata() : graphdata()
  {
  }
  jkdata::jkdata(string filename, string plotname) : graphdata()
  {
    init(filename, plotname);
  }
  jkdata::jkdata(string filename, string plotname, double a_x, double a_y) : graphdata(a_x, a_y)
  {
    init(filename, plotname);
  }
  jkdata::jkdata(const jkdata &other) : graphdata(other)
  {
  }

  void jkdata::init(string filename, string plotname)
  {
    aveGraph = getAve(filename, plotname);
    for (int i = 0; i < aveGraph->GetN(); ++i)
    {
      double x = aveGraph->GetX()[i] * a_x;
      x = f_x->Eval(x);
      aveGraph->SetPoint(i, x, aveGraph->GetY()[i] * a_y);
      aveGraph->SetPointError(i, 0, aveGraph->GetErrorY(i) * a_y);
    }
    makeUP(aveGraph);
    aveGraph->SetName(Form("%s_ave", plotname.c_str()));

    mg = checkjkdata(filename, plotname);
    mg->SetName(plotname.c_str());
    ngraphs = mg->GetListOfGraphs()->GetSize();
    TList *lg = mg->GetListOfGraphs();

    TIter next(lg);
    TObject *obj;
    int ig = 0;
    while ((obj = next()))
    {
      TGraph *graph = dynamic_cast<TGraph *>(obj);
      graph->SetName(Form("%s_jk%d", plotname.c_str(), ig));
      graph->SetTitle(Form("jk-%d", ig));
      if (graph)
      {
        // 获取图的点数
        int n = graph->GetN();
        if (n == 0)
          continue;

        for (int i = 0; i < n; ++i)
        {
          double x = graph->GetX()[i] * a_x;
          x = f_x->Eval(x);
          graph->SetPoint(i, x, graph->GetY()[i] * a_y);
        }
        makeUP(graph, ig);
        leg->AddEntry(graph, Form("jk-%d", ig), "lp");
        ig++;
      }
    }
    mg->SetDrawOption(drawoption.c_str());
  }

  jkdata::~jkdata()
  {
  }

  jkdata jkdata::Clone() const
  {
    jkdata res;
    res.aveGraph = (TGraphErrors *)(aveGraph->Clone());
    res.mg = (TMultiGraph *)(mg->Clone());
    res.ngraphs = ngraphs;
    res.f_x = (TF1 *)(f_x->Clone());
    res.leg = (TLegend *)(leg->Clone());
    res.a_x = a_x;
    res.a_y = a_y;
    res.Color = Color;
    res.Pallete = Pallete;
    res.Alpha = Alpha;
    res.MarkerStyle = MarkerStyle;
    res.MarkerSize = MarkerSize;
    res.drawoption = drawoption;
    return res;
  }

  void jkdata::GetPoint(int i, double &x, vector<double> &y, bool ifunited) const
  {
    x=aveGraph->GetX()[i];
    if(!ifunited)
      x/=a_x;
    y.resize(ngraphs);
    for (int j = 0; j < ngraphs; j++)
    {
      TGraph *g = (TGraph *)mg->GetListOfGraphs()->At(j);
      y[j] = g->GetY()[i];
      if(!ifunited)
        y[j] /= a_y;
    }
  }

  void jkdata::AddPoint(double x, vector<double> y, bool ifunited) const
  {
    double ax = 1;
    double ay = 1;
    if (!ifunited)
    {
      ax = a_x;
      ay = a_y;
    }
    x *= ax;
    double avey = 0;
    double erry = 0;
    if(y.size() != ngraphs)
    {
      cERR << "import data size not equal to ngraphs!" << endl;
      return;
    }
    for (int i = 0; i < y.size(); i++)
    {
      TGraph *g = (TGraph *)mg->GetListOfGraphs()->At(i);
      g->AddPoint(x, y[i] * ay);
      avey += y[i];
      erry += y[i] * y[i];
    }
    avey /= y.size();
    erry /= y.size();
    erry = sqrt((erry - avey * avey) * (y.size() - 1));
    aveGraph->AddPoint(x, avey * ay);
    aveGraph->SetPointError(aveGraph->GetN() - 1, 0, erry * ay);
  }

  void jkdata::RemovePoint(int i)
  {
    aveGraph->RemovePoint(i);
    for (int j = 0; j < ngraphs; j++)
    {
      TGraph *g_j = getiGraph(j);
      g_j->RemovePoint(i);
    }
  }

  TGraphErrors *jkdata::getiGraphErrors(int ien)
  {
    TGraph *ig = (TGraph *)(mg->GetListOfGraphs()->At(ien));
    TGraphErrors *g = new TGraphErrors(ig->GetN());
    for (int i = 0; i < g->GetN(); i++)
    {
      g->SetPoint(i, ig->GetX()[i], ig->GetY()[i]);
      g->SetPointError(i, 0, aveGraph->GetErrorY(i));
    }
    makeUP(g, ien);

    return g;
  }

  jkdata jkdata::operators(const jkdata &other, int oper) const
  {
    jkdata res = Clone(); // proprety copy as same as rhs
    if (ngraphs != other.ngraphs)
    {
      cERR << "ngraphs not equal! Can't calculate them." << endl;
      return jkdata();
    }
    int datasize = aveGraph->GetN();
    if (datasize != other.aveGraph->GetN())
    {
      cERR << "datasize aveGraph->GetN() not equal! Can't calculate them." << endl;
      return jkdata();
    }
    if (a_x != other.a_x || a_y != other.a_y)
      cWARN << "a_x or a_y not equal! Using rhs's a_x and a_y." << endl;

    delete res.mg;
    delete res.aveGraph;
    res.mg = new TMultiGraph();
    res.aveGraph = new TGraphErrors(datasize);
    vector<double> ave(datasize, 0), ave2(datasize, 0);
    for (int i = -1; i < ngraphs; i++)
    {
      TGraph *g;
      TGraph *rg;
      TGraph *og;
      if (i == -1)
      {
        g = new TGraphErrors(datasize);
        rg = aveGraph;
        og = other.aveGraph;
      }
      else
      {
        g = new TGraph(datasize);
        rg = getiGraph(i);
        og = other.getiGraph(i);
      }
      for (int j = 0; j < datasize; j++)
      {
        double x = rg->GetX()[j] / a_x;
        double ox = og->GetX()[j] / other.a_x;
        if (ox != x)
        {
          cERR << "x not equal! Can't calculate them." << endl;
          return jkdata();
        }
        double ry = rg->GetY()[j] / a_y;
        double oy = og->GetY()[j] / other.a_y;

        double y;
        switch (oper)
        {
        case 0: // add
          y = (ry + oy) * res.a_y;
          break;
        case 1: // sub
          y = (ry - oy) * res.a_y;
          break;
        case 2: // mul
          y = (ry * oy) * res.a_y * res.a_y;
          break;
        case 3: // div
          y = ry / oy;
          break;

        default:
          y = 0. / 0; // NaN
          break;
        }
        g->SetPoint(j, x * res.a_x, y);
        if (i != -1)
        {
          ave[j] += y;
          ave2[j] += y * y;
        }
      }
      makeUP(g, i);
      if (i == -1)
        res.aveGraph = (TGraphErrors *)g;
      else
        res.mg->Add(g);
    }

    for (int j = 0; j < datasize; j++)
    {
      double avej = ave[j] / ngraphs;
      double ave2j = ave2[j] / ngraphs;
      double errj = sqrt((ave2j - avej * avej) * (ngraphs - 1));
      res.aveGraph->SetPointError(j, 0, errj);
    }
    return res;
  }

  jkdata jkdata::operators(const double oy, int oper) const
  {
    jkdata res = Clone(); // proprety copy as same as rhs
    int datasize = aveGraph->GetN();

    delete res.mg;
    delete res.aveGraph;
    res.mg = new TMultiGraph();
    res.aveGraph = new TGraphErrors(datasize);
    vector<double> ave(datasize, 0), ave2(datasize, 0);
    for (int i = -1; i < ngraphs; i++)
    {
      TGraph *g;
      TGraph *rg;
      TGraph *og;
      if (i == -1)
      {
        g = new TGraphErrors(datasize);
        rg = aveGraph;
      }
      else
      {
        g = new TGraph(datasize);
        rg = getiGraph(i);
      }
      for (int j = 0; j < datasize; j++)
      {
        double x = rg->GetX()[j] / a_x;
        double ry = rg->GetY()[j] / a_y;

        double y;
        switch (oper)
        {
        case 0: // add
          y = (ry + oy);
          break;
        case 1: // sub
          y = (ry - oy);
          break;
        case 2: // mul
          y = ry * oy;
          break;
        case 3: // div
          y = ry / oy;
          break;

        default:
          y = 0. / 0; // NaN
          break;
        }
        g->SetPoint(j, x * res.a_x, y *res.a_y);
        if (i != -1)
        {
          ave[j] += y;
          ave2[j] += y * y;
        }
      }
      makeUP(g, i);
      if (i == -1)
        res.aveGraph = (TGraphErrors *)g;
      else
        res.mg->Add(g);
    }

    for (int j = 0; j < datasize; j++)
    {
      double avej = ave[j] / ngraphs;
      double ave2j = ave2[j] / ngraphs;
      double errj = sqrt((ave2j - avej * avej) * (ngraphs - 1));
      res.aveGraph->SetPointError(j, 0, errj * res.a_y);
    }
    return res;
  }

  jkdata jkdata::operators(const TGraph* oG, int oper) const
  {
    jkdata res = Clone(); // proprety copy as same as rhs
    int datasize = aveGraph->GetN();

    delete res.mg;
    delete res.aveGraph;
    res.mg = new TMultiGraph();
    res.aveGraph = new TGraphErrors(datasize);
    vector<double> ave(datasize, 0), ave2(datasize, 0);
    for (int i = -1; i < ngraphs; i++)
    {
      TGraph *g;
      TGraph *rg;
      TGraph *og;
      if (i == -1)
      {
        g = new TGraphErrors(datasize);
        rg = aveGraph;
      }
      else
      {
        g = new TGraph(datasize);
        rg = getiGraph(i);
      }
      for (int j = 0; j < datasize; j++)
      {
        double x = rg->GetX()[j];
        double ry = rg->GetY()[j];
        double oy = oG->Eval(x);

        double y;
        switch (oper)
        {
        case 0: // add
          y = (ry + oy);
          break;
        case 1: // sub
          y = (ry - oy);
          break;
        case 2: // mul
          y = (ry * oy);
          break;
        case 3: // div
          y = ry / oy;
          break;

        default:
          y = 0. / 0; // NaN
          break;
        }
        g->SetPoint(j, x, y);
        if (i != -1)
        {
          ave[j] += y;
          ave2[j] += y * y;
        }
      }
      makeUP(g, i);
      if (i == -1)
        res.aveGraph = (TGraphErrors *)g;
      else
        res.mg->Add(g);
    }

    for (int j = 0; j < datasize; j++)
    {
      double avej = ave[j] / ngraphs;
      double ave2j = ave2[j] / ngraphs;
      double errj = sqrt((ave2j - avej * avej) * (ngraphs - 1));
      res.aveGraph->SetPointError(j, 0, errj);
    }
    return res;
  }

  jkdata jkdata::operators(TF1 *f, int oper) const
  {
    jkdata res = Clone();
    int datasize = aveGraph->GetN();

    delete res.mg;
    delete res.aveGraph;
    res.mg = new TMultiGraph();
    res.aveGraph = new TGraphErrors(datasize);
    vector<double> ave(datasize, 0), ave2(datasize, 0);
    for (int i = -1; i < ngraphs; i++)
    {
      TGraph *g;
      TGraph *rg;
      if (i == -1)
      {
        g = new TGraphErrors(datasize);
        rg = aveGraph;
      }
      else
      {
        g = new TGraph(datasize);
        rg = getiGraph(i);
      }
      for (int j = 0; j < datasize; j++)
      {
        double x = rg->GetX()[j] / a_x;
        double ry = rg->GetY()[j] / a_y;

        double y;
        switch (oper)
        {
        case 0: // add
          y = (ry + f->Eval(x));
          break;
        case 1: // sub
          y = (ry - f->Eval(x));
          break;
        case 2: // mul
          y = (ry * f->Eval(x));
          break;
        case 3: // div
          y = ry / f->Eval(x);
          break;

        default:
          y = 0. / 0; // NaN
          break;
        }
        g->SetPoint(j, x * res.a_x, y * res.a_y);
        if (i != -1)
        {
          ave[j] += y;
          ave2[j] += y * y;
        }
      }
      makeUP(g, i);
      if (i == -1)
        res.aveGraph = (TGraphErrors *)g;
      else
        res.mg->Add(g);
    }

    for (int j = 0; j < datasize; j++)
    {
      double avej = ave[j] / ngraphs;
      double ave2j = ave2[j] / ngraphs;
      double errj = sqrt((ave2j - avej * avej) * (ngraphs - 1));
      res.aveGraph->SetPointError(j, 0, errj * res.a_y);
    }
    return res;
  }

  jkdata jkdata::operator=(const jkdata &other)
  {
    if (this != &other)
    {
      if (mg != nullptr)
        delete mg;
      if (aveGraph != nullptr)
        delete aveGraph;
      if (f_x != nullptr)
        delete f_x;
      if (leg != nullptr)
        delete leg;
      mg = (TMultiGraph *)(other.mg->Clone());
      aveGraph = (TGraphErrors *)(other.aveGraph->Clone());
      f_x = (TF1 *)(other.f_x->Clone());
      leg = (TLegend *)(other.leg->Clone());
      ngraphs = other.ngraphs;
      a_x = other.a_x;
      a_y = other.a_y;
      Color = other.Color;
      Pallete = other.Pallete;
      Alpha = other.Alpha;
      MarkerStyle = other.MarkerStyle;
      MarkerSize = other.MarkerSize;
      drawoption = other.drawoption;
    }
    return *this;
  }

  jkdata jkdata::operator+(const jkdata &other) const
  {
    return operators(other, 0);
  }

  jkdata jkdata::operator-(const jkdata &other) const
  {
    return operators(other, 1);
  }

  jkdata jkdata::operator*(const jkdata &other) const
  {
    return operators(other, 2);
  }

  jkdata jkdata::operator/(const jkdata &other) const
  {
    return operators(other, 3);
  }

  jkdata jkdata::operator+(const double other) const
  {
    return operators(other, 0);
  }

  jkdata jkdata::operator-(const double other) const
  {
    return operators(other, 1);
  }

  jkdata jkdata::operator*(const double other) const
  {
    return operators(other, 2);
  }

  jkdata jkdata::operator/(const double other) const
  {
    return operators(other, 3);
  }

  jkdata jkdata::operator+(const TGraph* other) const
  {
    return operators(other, 0);
  }

  jkdata jkdata::operator-(const TGraph* other) const
  {
    return operators(other, 1);
  }

  jkdata jkdata::operator*(const TGraph* other) const
  {
    return operators(other, 2);
  }

  jkdata jkdata::operator/(const TGraph* other) const
  {
    return operators(other, 3);
  }

  jkdata jkdata::operator+(TF1 *f) const
  {
    return operators(f, 0);
  }

  jkdata jkdata::operator-(TF1 *f) const 
  {
    return operators(f, 1);
  }

  jkdata jkdata::operator*(TF1 *f) const 
  {
    return operators(f, 2);
  }

  jkdata jkdata::operator/(TF1 *f) const 
  {
    return operators(f, 3);
  }

  void jkdata::print(string filename) const 
  {
    ofstream out(filename.c_str());
    if (out.fail())
    {
      cERR << "Can't open file " << filename << endl;
      return;
    }
    out << "datasize: " << aveGraph->GetN() << endl;
    out << "jksize: " << ngraphs << endl;
    out << "aveData:" << endl;
    for (int i = 0; i < aveGraph->GetN(); i++)
      out << aveGraph->GetX()[i] << " " << aveGraph->GetY()[i] << " " << aveGraph->GetErrorY(i) << endl;
    out << "jkdata:" << endl;
    for (int i = 0; i < ngraphs; i++)
    {
      TGraph *g = getiGraph(i);
      for (int j = 0; j < g->GetN(); j++)
        out << g->GetX()[j] << " " << g->GetY()[j] << " " << g->GetErrorY(j) << endl;
      out << endl;
    }
    out.close();
  }