#include <iomanip>
#include "chackjkdata.h"
#include "TROOT.h"
#include <fstream>

using namespace std;

//class ColorsSchemes

  ColorsSchemes::ColorsSchemes(){}
  ColorsSchemes::ColorsSchemes(ColorsSchemes &cs): kCstart(cs.kCstart),
                                    kCend(cs.kCend),
                                    kP6Blue(cs.kP6Blue),
                                    kP8Blue(cs.kP8Blue),
                                    kP10Blue(cs.kP10Blue)
  {
    int nColor = kCend - kCstart + 1;
    for (int i = 0; i < nColor; i++)
      kcs.push_back((TColor*)cs.kcs[i]->Clone());
  }
  ColorsSchemes::ColorsSchemes(std::vector<std::vector<int>> rgb, std::string tit)
  {
    register_custom_colors(rgb, tit);
  }
  ColorsSchemes::~ColorsSchemes()
  {
    unregister_custom_colors();
  }

  void ColorsSchemes::register_custom_colors(std::vector<std::vector<int>> rgb, std::string tit)
  {
    int rows = rgb.size();
    for (int i = 0; i < rows; ++i)
    {
      int r = rgb[i][0];
      int g = rgb[i][1];
      int b = rgb[i][2];

      kCend = kCstart + i;
      TColor *old_color = gROOT->GetColor(kCend);
      if (old_color != nullptr)
      {
        delete old_color; // 显式删除旧颜色
        old_color = nullptr;
      }
      string name = Form("%s_%d", tit.c_str(), i);
      CleanOld(name);
      TColor *color = new TColor(kCend, r / 255.0, g / 255.0, b / 255.0, name.c_str());
      kcs.push_back(color);
    }
  }

  void ColorsSchemes::unregister_custom_colors()
  {
    for (int i = 0; i < kcs.size(); ++i)
      if (kcs[i])
      {
        delete kcs[i];
        kcs[i] = nullptr;
      }
  }
