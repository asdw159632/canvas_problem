#include "chackjkdata.h"
#include "TLegendEntry.h"

using namespace std;
//class timedata : public graphdata
  timedata::timedata(string filename_model, int time_start, int time_end, string plotname, string drawopt) : graphdata()
  {
    drawoption = drawopt;
    init(filename_model, plotname, time_start, time_end);
  }
  timedata::timedata(string filename_model, int time_start, int time_end, string plotname, double a_x, double a_y, string drawopt) : graphdata(a_x, a_y, drawopt)
  {
    init(filename_model, plotname, time_start, time_end);
  }
  timedata::timedata(string filename_model, vector<int> t_list, string plotname, string drawopt) : graphdata()
  {
    drawoption = drawopt;
    init(filename_model, plotname, t_list);
  }
  timedata::timedata(string filename_model, vector<int> t_list, string plotname, double a_x, double a_y, string drawopt) : graphdata(a_x, a_y, drawopt)
  {
    init(filename_model, plotname, t_list);
  }
  timedata::timedata(vector<jkdata> &datalist, vector<int> t_list, string plotname, string drawopt) : graphdata()
  {
    drawoption = drawopt;
    init(datalist, plotname, t_list);
  }
  timedata::timedata(vector<jkdata> &datalist, vector<int> t_list, string plotname, double a_x, double a_y, string drawopt) : graphdata(a_x, a_y, drawopt)
  {
    init(datalist, plotname, t_list);
  }
  timedata::~timedata()
  {
  }

  void timedata::init(string filename_model, string plotname, int time_start, int time_end)
  {
    if (time_end < time_start)
    {
      cERR << "Invalid time range!" << endl;
      exit(1);
    }
    vector<int> t_list(0);
    for (int t = time_end; t >=time_start; t--)//倒序是因为一般时间大误差大，可以当作背景先画上
      t_list.push_back(t);
    init(filename_model, plotname, t_list);
  }

  void timedata::init(string filename_model, string plotname, vector<int> t_list)
  {
    vector<string> filenames(t_list.size());

    ngraphs = t_list.size();
    mg = new TMultiGraph();
    mg->SetName(Form("%s_%s_mg_timedata", filename_model.c_str(), plotname.c_str()));

    if (plotname == "Rcorr")//Rcorr一定是按时间排序,倒序是因为一般时间大误差大，可以当作背景先画上
    {
      sort(t_list.begin(), t_list.end(), std::greater<int>());

      ngraphs += 2;
    }
    for (int i = 0; i < t_list.size(); i++)
      filenames[i] = Form(filename_model.c_str(), t_list[i]);
    for (int i = 0; i < ngraphs; i++)
    {
      TGraphErrors *g;
      string filename;
      int time_index;
      if (plotname != "Rcorr")
      {
        filename = filenames[i];
        g = getAve(filename, plotname);
        time_index = t_list[i];
      }
      else
      {
        if (i == 0)
          filename = filenames[0];
        else if (i > ngraphs - 2)
          filename = filenames[i - 2];
        else
          filename = filenames[i - 1];
        time_index = t_list[0] - i + 1;
        g = getAve(filename, Form("Rcorr_t%03d", time_index));
      }
      if (g)
      {
        for (int j = 0; j < g->GetN(); j++)
        {
          double x = g->GetX()[j] * a_x;
          x = f_x->Eval(x);
          g->SetPoint(j, x, g->GetY()[j] * a_y);
          g->SetPointError(j, 0, g->GetErrorY(j) * a_y);
        }
        g->SetTitle(Form("t/a = %d", time_index));
        if (plotname != "Rcorr")
          g->SetName(Form("%s_timedata", filename.c_str()));
        else
          g->SetName(Form("%s_timedata", Form("Rcorr_t0%d", time_index)));
        makeUP(g, i);
        if (drawoption.find("C4") != string::npos || drawoption.find("C3") != string::npos)
          leg->AddEntry(g, Form("t/a = %d", time_index), "lf");
        else
          leg->AddEntry(g, Form("t/a = %d", time_index), "ep");
        mg->Add(g);
      }
    }
    mg->SetDrawOption(drawoption.c_str());
  }

  void timedata::init(const vector<jkdata>& datalist, string plotname, vector<int> t_list)
  {
    if (t_list.size() != datalist.size())
    {
      cerr << "Invalid time range! t_list size is not equal to datalist" << endl;
      exit(1);
    }

    ngraphs = datalist.size();
    mg = new TMultiGraph();
    mg->SetName(Form("%s_mg_timedata", plotname.c_str()));
    for (int i = 0; i < ngraphs; i++)
    {
      TGraphErrors *g = (TGraphErrors*)datalist[i].aveGraph->Clone();
      int time_index = t_list[i];
      if (g)
      {
        for (int j = 0; j < g->GetN(); j++)
        {
          double x = g->GetX()[j] * a_x;
          x = f_x->Eval(x);
          g->SetPoint(j, x, g->GetY()[j] * a_y);
          g->SetPointError(j, 0, g->GetErrorY(j) * a_y);
        }
        g->SetTitle(Form("t/a = %d", time_index));
        g->SetName(Form("%s_t%03d_timedata", plotname.c_str(), time_index));
        makeUP(g, i);
        if (drawoption.find("C4") != string::npos || drawoption.find("C3") != string::npos)
          leg->AddEntry(g, Form("%s t/a = %d", plotname.c_str(), time_index), "lf");
        else
          leg->AddEntry(g, Form("%s t/a = %d", plotname.c_str(), time_index), "ep");
        mg->Add(g);
      }
    }
    mg->SetDrawOption(drawoption.c_str());
  }


  TGraphErrors *timedata::getiGraph(int ien)
  {
    TGraphErrors *g = (TGraphErrors *)(mg->GetListOfGraphs()->At(ien));

    return g;
  }

  TGraphErrors *timedata::getSlide(int ix)
  {
    TGraphErrors *g = new TGraphErrors(ngraphs);
    double x = getiGraph(0)->GetX()[ix];
    int Nx = getiGraph(0)->GetN();
    for (int i = 0; i < ngraphs; i++)
    {
      TGraphErrors *g_i = getiGraph(i);
      TList *leg_lf = leg->GetListOfPrimitives();
      TLegendEntry *le = (TLegendEntry *)(leg->GetListOfPrimitives()->At(i));
      const char *enTit = le->GetLabel();
      int time_i;
      sscanf(enTit, "t/a = %d", &time_i);
      g->SetPoint(i, (double)time_i, g_i->GetY()[ix]);
      g->SetPointError(i, 0, g_i->GetErrorY(ix));
    }
    g->SetName(Form("slide_r = %f", x));
    g->SetTitle(Form("slide_r = %f", x));
    makeUP(g, Nx, ix);
    return g;
  }

  void timedata::GetPoint(int i, double &x, vector<double> &y, vector<double> &ey, bool ifunited)
  {
    x = getiGraph(0)->GetX()[i];
    if (!ifunited)
      x /= a_x;
    y.resize(ngraphs);
    ey.resize(ngraphs);
    for (int j = 0; j < ngraphs; j++)
    {
      TGraphErrors *g_j = getiGraph(j);
      y[j] = g_j->GetY()[i];
      ey[j] = g_j->GetErrorY(i);
      if (!ifunited)
      {
        y[j] /= a_y;
        ey[j] /= a_y;
      }
    }
  }

  void timedata::AddPoint(double x, vector<double> y, vector<double> ey, bool ifunited)
  {
    double ay = 1;
    if(!ifunited)
    {
      x *= a_x;
      ay = a_y;
    }
    if (y.size() != ngraphs)
    {
      cERR << "import data size not equal to ngraphs!" << endl;
      return;
    }
    if (y.size() != ey.size())
    {
      cERR << "import data size not equal to its error size!" << endl;
      return;
    }
    for (int i = 0; i < ngraphs; i++)
    {
      TGraphErrors *g_i = getiGraph(i);
      g_i->AddPoint(x, y[i] * ay);
      g_i->SetPointError(g_i->GetN() - 1, 0, ey[i] * ay);
    }
  }

  void timedata::RemovePoint(int i)
  {
    for (int j = 0; j < ngraphs; j++)
    {
      TGraphErrors *g_j = getiGraph(j);
      g_j->RemovePoint(i);
    }
  }