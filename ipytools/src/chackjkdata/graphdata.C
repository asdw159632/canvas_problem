#include "chackjkdata.h"
#include "TROOT.h"
#include "TLegendEntry.h"
#include "TStyle.h"
#include <fstream>

using namespace std;
//class graphdata

  graphdata::graphdata()
  {
  }
  graphdata::graphdata(double a_x, double a_y) : a_x(a_x), a_y(a_y)
  {
  }
  graphdata::graphdata(double a_x, double a_y, string drawoption) : a_x(a_x), a_y(a_y), drawoption(drawoption)
  {
  }
  graphdata::graphdata(const graphdata &other) : aveGraph((TGraphErrors *)(other.aveGraph->Clone())),
                                      mg((TMultiGraph *)(other.mg->Clone())),
                                      ngraphs(other.ngraphs),
                                      leg((TLegend *)(other.leg->Clone())),
                                      a_x(other.a_x),
                                      a_y(other.a_y),
                                      f_x((TF1 *)(other.f_x->Clone())),
                                      Color(other.Color),
                                      Pallete(other.Pallete),
                                      Alpha(other.Alpha),
                                      MarkerStyle(other.MarkerStyle),
                                      MarkerSize(other.MarkerSize),
                                      drawoption(other.drawoption)
  {
  }
  graphdata::~graphdata()
  {
    if (mg)
    {
      delete mg;
      mg = nullptr;
    }
    if (aveGraph)
    {
      delete aveGraph;
      aveGraph = nullptr;
    }
    delete f_x;
    if (leg)
    {
      delete leg;
      leg = nullptr;
    }
  }

  TLegend* graphdata::resetLeg(vector<double> Locate, double TextSize, int TextFont, string suffix)
  {
    TObject *obj = (TObject *)gROOT->FindObject("rleg");
    if (obj)
      delete obj;
    TLegend *rleg = new TLegend(Locate[0], Locate[1], Locate[2], Locate[3]);
    rleg->SetName("rleg");
    TList *entries = leg->GetListOfPrimitives();

    // 遍历所有 TLegendEntry 并修改
    for (Int_t i = 0; i < entries->GetSize(); i++)
    {
      TLegendEntry *entry = (TLegendEntry *)entries->At(i);
      rleg->AddEntry(entry->GetObject(), entry->GetLabel(), entry->GetOption());
    }
    return rleg;
  }

  void graphdata::setAxAy(double ax, double ay)
  {
    a_x = ax;
    a_y = ay;

    if (aveGraph != nullptr)
    {
      for (int i = 0; i < aveGraph->GetN(); ++i)
      {
        double x = aveGraph->GetX()[i];
        double y = aveGraph->GetY()[i];
        double ey = aveGraph->GetEY()[i];
        x *= a_x;
        y *= a_y;
        ey *= a_y;
        aveGraph->SetPoint(i, x, y);
        aveGraph->SetPointError(i, 0, ey);
      }
      f_x->SetRange(aveGraph->GetXaxis()->GetXmin(), aveGraph->GetXaxis()->GetXmax());
    }
    if (mg != nullptr)
      for (int ien = 0; ien < ngraphs; ++ien)
      {
        TGraph *graph = getiGraph(ien);
        if (graph)
        {
          // 获取图的点数
          int n = graph->GetN();
          if (n == 0)
            continue;

          for (int i = 0; i < n; ++i)
          {
            double x = graph->GetX()[i];
            double y = graph->GetY()[i];
            x *= a_x;
            y *= a_y;
            graph->SetPoint(i, x, y);
          }
        }
      }
  }

  void graphdata::setFormula(string formula)
  {
    delete f_x;
    f_x = new TF1("f_x", formula.c_str(), 0, 100);

    if (aveGraph != nullptr)
    {
      for (int i = 0; i < aveGraph->GetN(); ++i)
      {
        double x = aveGraph->GetX()[i];
        x = f_x->Eval(x);
        aveGraph->SetPointX(i, x);
      }
      f_x->SetRange(aveGraph->GetXaxis()->GetXmin(), aveGraph->GetXaxis()->GetXmax());
    }
    if (mg != nullptr)
      for (int ien = 0; ien < ngraphs; ++ien)
      {
        TGraph *graph = getiGraph(ien);
        if (graph)
        {
          // 获取图的点数
          int n = graph->GetN();
          if (n == 0)
            continue;

          for (int i = 0; i < n; ++i)
          {
            double x = graph->GetX()[i];
            x = f_x->Eval(x);
            graph->SetPointX(i, x);
          }
        }
      }
  }

  void graphdata::setColor(int color)
  {
    Color = color;
    if (aveGraph != nullptr)
    {
      aveGraph->SetLineColor(color);
      aveGraph->SetFillColorAlpha(color, Alpha);
      aveGraph->SetMarkerColor(color);
    }
  }

  void graphdata::setAlpha(double alpha)
  {
    Alpha = alpha;
    if (aveGraph != nullptr)
      aveGraph->SetFillColorAlpha(Color, Alpha);
    if (mg != nullptr)
      for (int ien = 0; ien < ngraphs; ++ien)
        getiGraph(ien)->SetFillColorAlpha(Color, Alpha);
  }

  void graphdata::setPallete(int pallete)
  {
    Pallete = pallete;
    gStyle->SetPalette(Pallete);
    if (mg != nullptr)
      for (int ien = 0; ien < ngraphs; ++ien)
      {
        TGraph *g = getiGraph(ien);
        int color = Color;
        if (ien >= 0)
        {
          int Ncolor = gStyle->GetNumberOfColors();
          int icolor = gStyle->GetColorPalette(Ncolor / ngraphs * ien);
          color = gROOT->GetColor(icolor)->GetNumber();
        }
        g->SetLineColor(color);
        g->SetFillColorAlpha(color, Alpha);
        g->SetMarkerColor(color);
      }
  }

  void graphdata::setMarkerStyle(int style)
  {
    MarkerStyle = style;
    if (aveGraph != nullptr)
      aveGraph->SetMarkerStyle(style);
    if (mg != nullptr)
      for (int ien = 0; ien < ngraphs; ++ien)
        getiGraph(ien)->SetMarkerStyle(style);
  }

  void graphdata::setMarkerSize(double size)
  {
    MarkerSize = size;
    if (aveGraph != nullptr)
      aveGraph->SetMarkerSize(size);
    if (mg != nullptr)
      for (int ien = 0; ien < ngraphs; ++ien)
        getiGraph(ien)->SetMarkerSize(size);
  }

  void graphdata::setDrawOption(string drawopt)
  {
    drawoption = drawopt;

    int linewid = 1;
    if (drawopt.find("C4") != string::npos)
      linewid = 4;
    if (aveGraph != nullptr)
    {
      aveGraph->SetLineWidth(linewid);
      aveGraph->SetDrawOption(drawopt.c_str());
    }
    if (mg != nullptr)
      for (int ien = 0; ien < ngraphs; ++ien)
      {
        TGraph *g = getiGraph(ien);
        g->SetLineWidth(linewid);
        g->SetDrawOption(drawopt.c_str());
      }
  }

  void graphdata::makeUP(TGraph *g, int Ng, int ien) const
  {
    g->SetMarkerStyle(MarkerStyle);
    g->SetMarkerSize(MarkerSize);
    gStyle->SetPalette(Pallete);
    int color = Color;
    if (ien >= 0)
    {
      int Ncolor = gStyle->GetNumberOfColors();
      int icolor = gStyle->GetColorPalette(Ncolor / Ng * ien);
      color = gROOT->GetColor(icolor)->GetNumber();
    }
    g->SetLineColor(color);
    if (drawoption.find("C") != string::npos)
    {
      g->SetLineWidth(4);
    }
    g->SetFillColorAlpha(color, Alpha);
    g->SetMarkerColor(color);
  }

  void graphdata::makeUP(TGraph *g, int ien) const 
  {
    makeUP(g, ngraphs, ien);
  }

  TGraph *graphdata::getiGraph(int ien) const 
  {
    TGraph *g = (TGraph *)(mg->GetListOfGraphs()->At(ien));

    return g;
  }

  TMultiGraph *graphdata::getAllGraph() const 
  {
    return mg;
  }

  TGraphErrors *graphdata::getAveGraph() const 
  {
    return aveGraph;
  }

