#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jun 12 14:12:37 2017

@author: malfatti
"""
import numpy as np

from datetime import datetime


## Level 0
def Set(Backend='TkAgg', AxesObj=(), FigObj=(), FigTitle='', Params=False, 
            Plot=False, Axes=False):
    if Params:
#            print('Set plot parameters...')
        Params = {'backend': Backend,
#                  'text.usetex': True, 'text.latex.unicode': True,
#                  'text.latex.preamble': '\\usepackage{siunitx}',
#                  'font.family': 'serif', 'font.serif': 'Computer Modern Roman',
                  'axes.titlesize': 'medium', 'axes.labelsize': 'medium',
                  'xtick.labelsize': 'small', 'xtick.direction': 'out',
                  'ytick.labelsize': 'small', 'ytick.direction': 'out',
                  'legend.fontsize': 'small', 'legend.labelspacing': 0.4,
                  'figure.titlesize': 'large', 'figure.titleweight': 'normal',
                  
                  'image.cmap': 'cubehelix', 'savefig.transparent': True,
                  'svg.fonttype': 'none'}
        return(Params)
    
    elif Plot:
#            print('Set plot figure...')
        FigObj.suptitle(FigTitle); FigObj.tight_layout(); 
        FigObj.subplots_adjust(top=0.925)
    
    elif Axes:
#            print('Set plot axes...')
        AxesObj.spines['right'].set_visible(False)
        AxesObj.spines['top'].set_visible(False)
        AxesObj.yaxis.set_ticks_position('left')
        AxesObj.xaxis.set_ticks_position('bottom')
        AxesObj.locator_params(tight=True)
    
    else: print("'Params', 'Plot' or 'Axes' must be True.")
    
    return(None)


def SignificanceBar(X, Y, Text, Ax, FontSize=9, TicksDir='down', lw=1, color='k'):
    if TicksDir == 'down':
        from matplotlib.markers import TICKDOWN as Tick
        Yy = max(Y)+(max(Y)*0.02)
    elif TicksDir == 'up':
        from matplotlib.markers import TICKUP as Tick
        Yy = max(Y)-(max(Y)*0.02)
    else: print('TicksDir should be "up" or "down".'); return(None)
    
    Ax.plot(X, Y, color=color, lw=lw, marker=Tick)
    Ax.text(sum(X)/2, Yy, Text, fontsize=FontSize, ha='center', va='center')
    return(None)


def RawCh(Ch, Lines, Cols, XValues=[], Slice=[], Leg=[], FigName='', Colors='', 
          Visible=True, Save=True):
    Params = Set(Params=True)
    from matplotlib import rcParams; rcParams.update(Params)
    from matplotlib import pyplot as plt
    
    PlotNo = len(Ch)
    plt.figure(figsize=(8, PlotNo))
    Axes = [plt.subplot(Lines, Cols, _+1) for _ in range(PlotNo)]
    
    if XValues == []: XValues = range(len(Ch[0]))
    
    for Ind, Ax in enumerate(Axes):
        if Slice: Line = Ax.plot(XValues[Slice[0]:Slice[1]], 
                                 Ch[Ind][Slice[0]:Slice[1]])
        else: Line = Ax.plot(XValues, Ch[Ind])
        
        if Colors: Line[0].set_color(Colors[Ind])
        if Leg: Line[0].set_label(Leg[Ind]); Ax.legend(loc='best')
        Ax.xaxis.set_visible(False)
        Ax.yaxis.set_visible(False)
        Ax.set_title(str(Ind+1))
    
    plt.tight_layout()
    
    if Save: 
        if not FigName: 
            Now = datetime.now().strftime("%Y%m%d%H%M%S")
            FigName = 'RawCh-' + Now + '.svg'
        
        print('Writing to', FigName+'... ', end='')
        plt.savefig(FigName, format='svg')
        print('Done.')
    
    if Visible: plt.show()
    else: plt.close()
    return(None)


def Spectrogram(SxxAx, T, F, Sxx, Colormap='inferno', HighFreqThr=None, 
                Line=None, LineX=None, LineY=None, LineColor='b', LineLim=None, 
                LineYLabel=None):
    if Line:
        if LineY is None: 
            print('If Line=True, LineY should receive plottable data.')
            return(None)
        if LineX is None: LineX = np.arange(len(LineY))
        if LineLim is None: LineLim = [min(LineY), max(LineY)]
    
    if HighFreqThr is None: HighFreqThr = max(F)
    
#        Params = Plot.Set(Params=True)
#        from matplotlib import rcParams; rcParams.update(Params)
#        from matplotlib import pyplot as plt
    
#        Fig, SxxAx = plt.subplots(1,1,figsize=(6, 4))
#        Fig.subplots_adjust(bottom=0.15, right=0.85)
    
    SxxAx.pcolormesh(T, F, Sxx, cmap='inferno')
    SxxAx.set_ylim(0,HighFreqThr)
    SxxAx.yaxis.set_ticks_position('left')
    SxxAx.xaxis.set_ticks_position('bottom')
    SxxAx.set_xlabel('Time [s]'); SxxAx.set_ylabel('Frequency [Hz]')
    
    if Line:
        SAx = SxxAx.twinx()
        SAx.plot(LineX, LineY, LineColor); SAx.set_ylim(LineLim)
        SAx.set_ylabel(LineYLabel, color=LineColor)
        SAx.xaxis.set_ticks_position('bottom')
    
#        if not FigName: FigName = 'Spectrogram.' + Ext
#        if Save: Fig.savefig(FigName, format=Ext)
#        if Visible: plt.show()
#        else: plt.close()


def RTTest(Files, Save=True, Return=False):
    Params = {'backend': 'TkAgg'}
    from matplotlib import rcParams; rcParams.update(Params)
    from matplotlib import pyplot as plt
    
    Jitter = {}
    for File in Files:
        with open(File, 'r') as FHist:
            for LineNo, Line in enumerate(FHist, 1): pass
        
        Jitter[File] = [[] for _ in range(LineNo)]
        with open(File, 'r') as FHist:
            for Ind, Line in enumerate(FHist):
                if Line == 'Avg:\n': Line = '51'
                Jitter[File][Ind] = int(Line)
    
    for F in Jitter.values(): plt.plot(F, 'bo')
    plt.show()
    
    if Return: return(Jitter)
    else: return(None)


def TTLCh(Ch, Std, Threshold, XValues=[]):
    Params = Set(Params=True)
    from matplotlib import rcParams; rcParams.update(Params)
    from matplotlib import pyplot as plt
    
    if not Std: Std = 3
    
    plt.figure(figsize=(8, 1))
    Ax = plt.subplot(111)
    
    if not XValues: Ax.plot(Ch)
    else: Ax.plot(XValues, Ch)
    
    ThrX = [min(Ax.get_xticks()), max(Ax.get_xticks())]
    Ax.plot(ThrX, [Threshold]*2)
    
    Ax.xaxis.set_visible(False); Ax.yaxis.set_visible(False)
    plt.tight_layout(); plt.show()
    return(None)

