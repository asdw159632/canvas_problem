# core is python 3.6.8
import ROOT
import math
from ctypes import c_double
import numpy as np
import pandas as pd

# 初始化标记（防止重复执行）
_INITIALIZED = False

def initialize(force=False):
    """
    初始化函数，可重复调用但只执行一次
    Args:
        force (bool): 强制重新初始化
    """
    global _INITIALIZED, c1, c2, c3, zero, a_fm, a_m1, mpi, color, chatGPTAdColor, m_h1,m_h2,delta,m_red
    
    if _INITIALIZED and not force:
        return
    
    # 1. initial Canvases
    c1 = ROOT.TCanvas("c", "c", 800, 600)
    c3 = ROOT.TCanvas("c3", "c3", 800, 600)
    c2 = ROOT.TCanvas("c2", "c2", 1600, 600)
    c2.Divide(2, 1)

    # 质量获取
    m_h1=4796.8
    m_h2=938.9
    m_red=785.2
    delta=0.672612
    m_red= 785.2

    
    # 设置画布边距
    for pad in [c2.cd(1), c2.cd(2), c1.cd()]:
        ROOT.gPad.SetLeftMargin(0.13)
        ROOT.gPad.SetBottomMargin(0.13)
        ROOT.gPad.SetTopMargin(0.08)
        ROOT.gPad.SetRightMargin(0.04)
    
    # 2. some needed staff
    zero = ROOT.TF1("0", "0", 0, 100)
    zero.SetLineWidth(4)
    zero.SetLineColor(ROOT.kBlack)
    
    # 3. lattice parameters
    a_fm = 0.084372  # fm
    a_m1 = 2338.8    # MeV
    mpi = 138        # MeV
    
    # 4. 加载C++头文件
    ROOT.gSystem.Load("./ipytools/lib/libjkdata.so")
    ROOT.gInterpreter.ProcessLine('#include "./ipytools/inc/chackjkdata.h"')


    # color scheme
    chatGPTAdColor=[
        [51, 179, 76],   #绿色
        [51, 102, 229],  #蓝色
        [255, 15, 24],   #红色
        [128, 0, 128],   #紫色
        [221, 160, 221], #淡紫
        [139, 69, 19],   #棕色
        [222, 184, 135], #米棕
        [204, 102, 26]   #琥珀   
    ]

    #color=ROOT.ColorsSchemes(ROOT.rgb_chat,"chat")
    color=ROOT.ColorsSchemes(chatGPTAdColor,"chat")
    
    _INITIALIZED = True

# 首次自动初始化
initialize()

# some function is widely used in this notebook
def printfitfunc(fitcheck, ien=-1, ifprint=True):
    """打印拟合函数参数"""
    fitfunc = fitcheck.GetiFit(ien)
    colnames = ["parvalue", "parerror"]
    rows = []
    data = []
    
    for i in range(fitfunc.GetNpar()):
        if ien >= 0:
            data.append([fitfunc.GetParameter(i), "/"])
        else:
            data.append([fitfunc.GetParameter(i), fitfunc.GetParError(i)])
        rows.append(fitfunc.GetParName(i))
    
    chi2dof = fitcheck.GetChi2Dof(ien)
    if ien >= 0:
        data.append([chi2dof[0], "/"])
    else:
        data.append([chi2dof[0], chi2dof[1]])
    rows.append("chi2/DOF")
    
    df = pd.DataFrame(data, index=rows, columns=colnames)
    if ifprint:
        fitfunc.Print()
        print(df)
    return df

def printfitpars(fitcheck):
 colnames=["parvalue","parerror"]
 for ien in range(fitcheck.ngraphs):
  colnames.append("jk_"+str(ien))
 rows=[]
 data=[]
 for ip in range(fitcheck.fitfunc.GetNpar()):
  rows.append(fitcheck.fitfunc.GetParName(ip))
  datai=[]
  for ien in range(-1,fitcheck.ngraphs):
   fitfunc=fitcheck.GetiFit(ien)
   if ien>=0:
    datai.append(fitfunc.GetParameter(ip))
   else:
    datai.append(fitfunc.GetParameter(ip))
    datai.append(fitfunc.GetParError(ip))
  data.append(datai)
 rows.append("chi2/DOF")
 datai=[]
 for ien in range(-1,fitcheck.ngraphs):
  chi2dof=fitcheck.GetChi2Dof(ien)
  if ien>=0:
    datai.append(chi2dof[0])
  else:
   datai.append(chi2dof[0])
   datai.append(chi2dof[1])
 data.append(datai) 
 df=pd.DataFrame(data,index=rows,columns=colnames)
 print(df)
 return df
 
#Add a table print of the time-dependent fit results
def show_fitres(trange,tarfit="fit_func",spin=1,Set="aveset",sys="FconfNOmgccc",info="B01_1600conf_080bin_misnered"):
 colnames=[]
 rows=[]
 data=[]
 suffix=""
 graph="pot"
 if tarfit!="fit_func":
  suffix="_observables"
  graph="phase_shift_kcotd_nocoul"
 for time in range(trange[0],trange[1]+1):
  rows.append(f"t={time}")
  jkfit=ROOT.jkfitdata(f"output/{Set}/pot_{2*spin+1}s{spin}_cen_{sys}_t{time:03d}_{info}{suffix}",graph,tarfit)
  func=jkfit.fitfunc
  if time==trange[0]:
   colnames=[func.GetParName(i) for i in range(func.GetNpar())]
   colnames.append("Chi2/Ndf")
  rowdata=[f"{func.GetParameter(i):.5f} +- {func.GetParError(i):.5f}" for i in range(func.GetNpar())]
  rowdata.append(f"{jkfit.ave_chi2dof[0]:.5f} +- {jkfit.ave_chi2dof[1]:.5f}")
  data.append(rowdata)
 df=pd.DataFrame(data,columns=colnames,index=rows)
 return df,func

def df_trans_to_pave(df,pave, func, x="r", unit="fm"):
    """将DataFrame转换为TPaveText"""
    table_text = df.to_string()
    pave.SetTextAlign(12)
    pave.SetTextSize(0.022)
    pave.AddText(func.GetExpFormula().Data())
    
    fitrange = np.array([func.GetXmin(), func.GetXmax()])
    pave.AddText(f"fit range: {x}#in[{fitrange[0]:.2f},{fitrange[1]:.2f}] {unit}")
    
    for line in table_text.split("\n"):
        pave.AddText(line)
    return pave

def setcolors(smg, color=color, Alpha=0.3):
  ncolor=color.kCend-color.kCstart+1
  for i in range(smg.GetListOfGraphs().GetEntries()):
   #using chatGPT advised
   ithcolor=i%ncolor
   smg.GetListOfGraphs().At(i).SetMarkerColor(color.kCstart + ithcolor)
   smg.GetListOfGraphs().At(i).SetLineColor(color.kCstart + ithcolor)
   smg.GetListOfGraphs().At(i).SetFillColorAlpha(color.kCstart + ithcolor, Alpha)

def showColorLeg(showCS=chatGPTAdColor):
  legc=ROOT.TLegend(0.5,0.3,0.9,0.6)
  legc.SetName("leg-color")
  m = [ROOT.TMarker(0,0,21) for i in range(len(showCS))]
  ncolor=color.kCend-color.kCstart+1
  for i in range(len(showCS)):
    ithcolor=i%ncolor
    m[i].SetMarkerColor(color.kCstart+ithcolor)
    m[i].SetMarkerSize(2)
    legc.AddEntry(m[i],f"{i}", "p")
  return legc,m
