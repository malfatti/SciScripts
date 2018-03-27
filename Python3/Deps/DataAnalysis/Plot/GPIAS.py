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
def Traces(GPIASData, XValues, SoundPulseDur, FigName, Ext=['svg'], 
              Save=True, Show=True):
    Params = Plot.Set(Params=True)
    from matplotlib import rcParams; rcParams.update(Params)
    from matplotlib import pyplot as plt
    
    print('Plotting...')
    Ind1 = list(XValues).index(0)
    Ind2 = list(XValues).index(int(SoundPulseDur*1000))
    
    PlotNo = len(GPIASData['Trace'].keys())
    Fig = plt.figure(figsize=(7, 12))
    Axes = [plt.subplot(PlotNo,1,_+1) for _ in range(PlotNo)]
#        
    for FInd, Freq in enumerate(GPIASData['Trace'].keys()):
        SubTitle = Freq + ' Hz' + ' Index = ' + str(round(GPIASData['Index'][Freq]['GPIASIndex'], 4))
        LineNoGapLabel = 'No Gap'; LineGapLabel = 'Gap'
        SpanLabel = 'Sound Pulse'
        XLabel = 'time [ms]'; YLabel = 'voltage [mV]'
        
        Axes[FInd].axvspan(XValues[Ind1], XValues[Ind2], color='k', alpha=0.5, 
                    lw=0, label=SpanLabel)
        Axes[FInd].plot(XValues, GPIASData['Trace'][Freq]['NoGap'], 
                 color='r', label=LineNoGapLabel, lw=2)
        Axes[FInd].plot(XValues, GPIASData['Trace'][Freq]['Gap'], 
                 color='b', label=LineGapLabel, lw=2)
        
        AxArgs = {'title': SubTitle, 'ylabel': YLabel, 'xlabel': XLabel}
        Plot.Set(Ax=Axes[FInd], AxArgs=AxArgs)
        Axes[FInd].legend(loc='center left', bbox_to_anchor=(1.05, 0.5),
                          prop={'size':6})
    
    FigTitle = FigName.split('/')[-1]
    Plot.Set(Fig=Fig, FigTitle=FigTitle)
    Fig.subplots_adjust(right=0.8)
    
    if Save: 
        for E in Ext: Fig.savefig(FigName+'.'+E, format=E)
    if Show: plt.show()
    else: plt.close()
    
    print('Done.')
    return(None)

