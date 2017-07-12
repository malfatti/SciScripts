#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jun 12 15:06:47 2017

@author: malfatti
"""
from DataAnalysis.Plot.Plot import Set as PlotSet


## Level 0
def Traces(GPIASData, XValues, SoundPulseDur, FigName, Ext='svg', 
              Save=True, Visible=False):
    Params = PlotSet(Params=True)
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
        Axes[FInd].set_title(SubTitle)
        Axes[FInd].set_ylabel(YLabel); Axes[FInd].set_xlabel(XLabel)

        PlotSet(AxesObj=Axes[FInd], Axes=True)
    
    FigTitle = FigName.split('/')[-1]
    FigName = FigName + '.' + Ext
    PlotSet(FigObj=Fig, FigTitle=FigTitle, Plot=True)
    
    if Save: Fig.savefig(FigName, format=Ext)
    if Visible: plt.show()
    else: plt.close()
    
    print('Done.')
    return(None)

