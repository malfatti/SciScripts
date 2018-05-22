#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: T. Malfatti <malfatti@disroot.org>
@date: 2017-06-12
@license: GNU GPLv3 <https://raw.githubusercontent.com/malfatti/SciScripts/master/LICENSE>
@homepage: https://github.com/Malfatti/SciScripts
"""
import numpy as np

from datetime import datetime

## Level 0
def GetTicks(Ax, Lim):
    # Override automatic tick formatter
    Step = round((Ax.get_yticks()[1]-Ax.get_yticks()[0])-0.5)
    Ticks = np.arange(min(Lim), max(Lim)+Step, Step)
    
    # Ticks = []
    # while not 0 in Ticks:
    #     Y = None
    #     while not Y:
    #         if not (max(Lim) - min(Lim))%MinTickNo: Y = MinTickNo
    #         else: MinTickNo += 1
    #     
    #     MinTickNo = Y
    #     Step = (max(Lim) - min(Lim))/Y
    #     Ticks = np.arange(min(Lim), max(Lim)+Step, Step)
    
    return(Ticks)


## Level 1
def Set(Backend='Qt5Agg', Ax=(), Fig=(), AxArgs={}, FigTitle='', Params=False, 
        HideControls=False, Tight=True):
    if Params:
        Params = {
          'backend'             : Backend,
          'image.cmap'          : 'cubehelix',
          'savefig.dpi'         : 300,
          'savefig.format'      : 'svg',
          'savefig.transparent' : True,
          
          'xtick.direction'     : 'out',
          'ytick.direction'     : 'out',
          
          'figure.figsize'      : (4, 3),
          'figure.dpi'          : 250,
          
          'svg.fonttype'        : 'none',
          'pdf.fonttype'        : 42,
          'font.family'         : 'sans-serif',
          'font.serif'          : ['DejaVu Serif'],
          'font.sans-serif'     : ['DejaVu Sans'],
          'font.cursive'        : ['Zapf Chancery'],
          'font.monospace'      : ['DejaVu Sans Mono'],
          'font.size'           : 10,
          # 'text.usetex'         : True, 
          # 'text.latex.unicode'  : True,
          # 'text.latex.preamble' : '\\usepackage{siunitx}',
        }
        
        return(Params)
    
    if Ax:
        XLim = Ax.get_xlim(); YLim = Ax.get_ylim()
        
        if 'title' in AxArgs: Ax.set_title(AxArgs['title'])
        
        if 'xlabel' in AxArgs: Ax.set_xlabel(AxArgs['xlabel'])
        if 'xticks' in AxArgs: Ax.set_xticks(AxArgs['xticks'])
        if 'xticklabels' in AxArgs: Ax.set_xticklabels(AxArgs['xticklabels'])
        if 'xlim' in AxArgs: XLim = AxArgs['xlim']
        
        if 'ylabel' in AxArgs: Ax.set_ylabel(AxArgs['ylabel'])
        if 'yticks' in AxArgs: Ax.set_yticks(AxArgs['yticks'])
        if 'yticklabels' in AxArgs: Ax.set_yticklabels(AxArgs['yticklabels'])
        if 'ylim' in AxArgs: YLim = AxArgs['ylim']
        
        Ax.set_xlim(min(XLim), max(XLim))
        Ax.set_ylim(min(YLim), max(YLim))
        
        if (max(XLim) - min(XLim)) % (Ax.get_xticks()[1]-Ax.get_xticks()[0]):
            print('Fixing xticks...')
            Ax.set_xticks(GetTicks(Ax, XLim))
        
        if (max(YLim) - min(YLim)) % (Ax.get_yticks()[1]-Ax.get_yticks()[0]):
            print('Fixing yticks...')
            Ax.set_yticks(GetTicks(Ax, YLim))
        
        # if 'xtickspacing' in AxArgs:
        #     import matplotlib.ticker as ticker
        #     Ax.xaxis.set_major_locator(ticker.MultipleLocator(AxArgs['xtickspacing']))
        
        # if 'ytickspacing' in AxArgs:
        #     import matplotlib.ticker as ticker
        #     Ax.yaxis.set_major_locator(ticker.MultipleLocator(AxArgs['ytickspacing']))
        
        # if 'ylim' in AxArgs:
        #     Ax.spines['left'].set_bounds(Ax.get_yticks()[0], Ax.get_yticks()[-1])
        # else:
        #     Ax.spines['left'].set_bounds(Ax.get_yticks()[1], Ax.get_yticks()[-2])
        
        # if 'xlim' in AxArgs:
        #     Ax.spines['bottom'].set_bounds(Ax.get_xticks()[0], Ax.get_xticks()[-1])
        # else:
        #     Ax.spines['bottom'].set_bounds(Ax.get_xticks()[1], Ax.get_xticks()[-2])
        
        Ax.spines['left'].set_bounds(min(YLim), max(YLim))
        Ax.spines['bottom'].set_bounds(min(XLim), max(XLim))
        Ax.spines['bottom'].set_position(('outward', 5))
        Ax.spines['left'].set_position(('outward', 5))
        
        Ax.tick_params(top='off', right='off')
        Ax.spines['right'].set_visible(False)
        Ax.spines['top'].set_visible(False)
        Ax.patch.set_visible(False)
#        Ax.locator_params(tight=True)
    
    if Fig:
        if FigTitle: Fig.suptitle(FigTitle)
        
        if HideControls:
            if Backend[:2].lower() == 'tk': Fig.canvas.toolbar.pack_forget()
            if Backend[:2].lower() == 'qt': Fig.canvas.toolbar.setVisible(False)
        
        if Tight:
            Fig.tight_layout(pad=0)
            Fig.subplots_adjust(top=0.9)
        
        for Obj in Fig.findobj(): Obj.set_clip_on(False)
        Fig.patch.set_visible(False)
    
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

