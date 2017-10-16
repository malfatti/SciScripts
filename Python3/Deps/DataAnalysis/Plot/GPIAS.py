#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jun 12 15:06:47 2017

@author: malfatti
"""
from DataAnalysis.Plot import Plot


## Level 0
def Traces(GPIASData, XValues, SoundPulseDur, FigName, Ext='svg', 
              Save=True, Visible=False):
    Params = Plot.Set(Params=True)
    from matplotlib import rcParams; rcParams.update(Params)
    from matplotlib import pyplot as plt
    
    print('Plotting...')
    Ind1 = list(XValues).index(0)
    Ind2 = list(XValues).index(int(SoundPulseDur*1000))
    
    PlotNo = len(GPIASData['Trace'].keys())
    Fig, Axes = plt.subplots(PlotNo, 1, figsize=(8, 3*PlotNo), sharex=True)
#        
    for FInd, Freq in enumerate(GPIASData['Trace'].keys()):
        SubTitle = Freq + ' Hz' + ' Index = ' + str(GPIASData['Index'][Freq]['GPIASIndex'])
        LineNoGapLabel = 'No Gap'; LineGapLabel = 'Gap'
        SpanLabel = 'Sound Pulse'
        XLabel = 'time [ms]'; YLabel = 'voltage [mV]'
        
        Axes[FInd].axvspan(XValues[Ind1], XValues[Ind2], color='k', alpha=0.5, 
                    lw=0, label=SpanLabel)
        Axes[FInd].plot(XValues, GPIASData['Trace'][Freq]['NoGap'], 
                 color='r', label=LineNoGapLabel, lw=2)
        Axes[FInd].plot(XValues, GPIASData['Trace'][Freq]['Gap'], 
                 color='b', label=LineGapLabel, lw=2)
        Axes[FInd].legend(loc='best')
        
        AxArgs = {'title': SubTitle, 'ylabel': YLabel, 'xlabel': XLabel}
        Plot.Set(Ax=Axes[FInd], AxArgs=AxArgs)
    
    FigTitle = FigName.split('/')[-1]
    FigName = FigName + '.' + Ext
    Plot.Set(Fig=Fig, FigTitle=FigTitle)
    
    if Save: Fig.savefig(FigName, format=Ext)
    if Visible: plt.show()
    else: plt.close()
    
    print('Done.')
    return(None)

