# ipytools/__init__.py
from .head import (  # 从 head.py 导入所有需要暴露的内容
    initialize,      # 初始化函数
    c1, c2, c3, zero,    # 全局变量（Canvas等）
    a_fm, a_m1, mpi, # 其他配置变量
    color, chatGPTAdColor,cell2_runtime,
    # 其他需要暴露的函数...
    printfitfunc, printfitpars, show_fitres, df_trans_to_pave, setcolors, showColorLeg,m_h1,m_h2,delta,m_red
)

# 可选：声明允许通过 `from ipytools import *` 导入的内容
__all__ = [
    "initialize", "c1", "c2", "c3", "zero", "m_h1","m_h2","delta","m_red","cell2_runtime",
    "a_fm", "a_m1", "mpi", "color", "chatGPTAdColor",
    "printfitfunc", "printfitpars", "show_fitres", "df_trans_to_pave", "setcolors", "showColorLeg"
]