#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: T. Malfatti <malfatti@disroot.org>
@year: 2017-06-12
@license: GNU GPLv3 <https://raw.githubusercontent.com/malfatti/SciScripts/master/LICENSE>
@homepage: https://github.com/Malfatti/SciScripts
"""
from DataAnalysis.Plot import Plot


## Level 0
def Traces(GPIASData, XValues, SoundPulseDur, FigName, Ext=['svg'], AxArgs={}, Save=True, Show=True):
    Params = Plot.Set(Params=True)
    from matplotlib import rcParams; rcParams.update(Params)
    from matplotlib import pyplot as plt
    
    print('Plotting...')
    Ind1 = list(XValues).index(0)
    Ind2 = list(XValues).index(int(SoundPulseDur*1000))
    
    YLim = []
    for F, Freq in GPIASData['Trace'].items():
        for G in Freq.values(): YLim.append(max(G)); YLim.append(min(G))
    
    # YLim = [(min(YLim)*1000 + 0.5).round(), (max(YLim)*1000 + 0.5).round()]
    YLim = [(min(YLim)*1000 - 0.5).round(), (max(YLim)*1000 + 0.5).round()]
    XLim = [round(min(XValues)-0.5), round(max(XValues)+0.5)]
    
    PlotNo = len(GPIASData['Trace'].keys())
    Fig = plt.figure(figsize=(6, 1.5*PlotNo))
    Axes = [plt.subplot(PlotNo,1,_+1) for _ in range(PlotNo)]
    
    for F, Freq in enumerate(GPIASData['Trace'].keys()):
        AxArgs['title'] = Freq + ' Hz' + ' Index = ' + str(round(GPIASData['Index'][Freq]['GPIASIndex'], 4))
        LineNoGapLabel = 'No Gap'; LineGapLabel = 'Gap'
        SpanLabel = 'Sound Pulse'
        XLabel = 'time [ms]'; YLabel = 'voltage [µV]'
        
        Axes[F].axvspan(XValues[Ind1], XValues[Ind2], color='k', alpha=0.5, 
                    lw=0, label=SpanLabel)
        Axes[F].plot(XValues, GPIASData['Trace'][Freq]['NoGap']*1000, 
                 color='r', label=LineNoGapLabel, lw=2)
        Axes[F].plot(XValues, GPIASData['Trace'][Freq]['Gap']*1000, 
                 color='b', label=LineGapLabel, lw=2)
        
        AxArgs = {**{'ylabel': YLabel, 'xlabel': XLabel, 
                     'ylim': YLim, 'xlim': XLim}, **AxArgs}
        Plot.Set(Ax=Axes[F], AxArgs=AxArgs)
        
        if F != len(GPIASData['Trace'].keys())-1:
            Axes[F].tick_params(bottom='off')
            Axes[F].spines['bottom'].set_visible(False)
            Axes[F].set_xticklabels([])
            Axes[F].set_xlabel('')
        
        if F != len(GPIASData['Trace'].keys())//2:
            Axes[F].tick_params(left='off')
            Axes[F].spines['left'].set_visible(False)
            Axes[F].set_yticklabels([])
            Axes[F].set_ylabel('')
        else:
            # Axes[F].legend(loc='center left', bbox_to_anchor=(1.05, 0.5),
            #               prop={'size':6})
            Axes[F].legend(loc='upper left', prop={'size':6})
    
    FigTitle = FigName.split('/')[-1]
    Plot.Set(Fig=Fig, FigTitle=FigTitle)
    # Fig.subplots_adjust(right=0.8)
    
    if Save: 
        for E in Ext: Fig.savefig(FigName+'.'+E, format=E)
    if Show: plt.show()
    else: plt.close()
    
    print('Done.')
    return(None)
