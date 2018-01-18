#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jul  7 21:13:16 2017

@author: malfatti
"""
import os
import numpy as np

from DataAnalysis.Plot import Plot
from glob import glob
from IO import Hdf5, Txt

def Traces(AnalysisPath, AnalysisFile, InfoFile, FigPath='./Figs', Ext=['svg'], Save=True, Visible=True):
    Data = Hdf5.DataLoad(AnalysisPath, AnalysisFile)[0]
    ABRs, XValues = Data['ABRs'], Data['XValues']
    
    if InfoFile[-4:] == 'hdf5': DataInfo = Hdf5.DictLoad('/DataInfo', InfoFile)
    else: DataInfo = Txt.DictRead(InfoFile)
    
    Params = Plot.Set(Params=True)
    from matplotlib import rcParams; rcParams.update(Params)
    from matplotlib import pyplot as plt
    
    
    print('Plotting...')
    Colormaps = [plt.get_cmap('Reds'), plt.get_cmap('Blues')]
    for S, Stim in ABRs.items():
        for DV, DVCoord in Stim.items():
            for F, Freq in DVCoord.items():
                for T, Trial in Freq.items():
                    YLim = []
                    
                    for ABR in Trial.values():
                        YLim.append(max(ABR)); YLim.append(min(ABR))
                    
                    Intensities = list(Trial.keys())
                    Intensities.sort(reverse=True)
                    TXValues = XValues[S][DV][F][T][:]
                    
                    Fig, Axes = plt.subplots(len(Intensities), sharex=True, 
                                             figsize=(8, 1.5*len(Intensities)))
                    
                    AxArgs = {'ylabel': 'Voltage [mV]', 
                              'xlabel': 'Time [ms]',
                              'ylim': (min(YLim), max(YLim))}
                        
                    for dB, ABR in Trial.items():
                        FigTitle = ' '.join([S, F+'Hz,', DV, 'DV, Trial', T])
                        LineLabel = dB
                        SpanLabel = 'Sound pulse'
                        
                        dBInd = Intensities.index(dB)
                        Colors = [Colormaps[0](255-(dBInd*20)), 
                                  Colormaps[1](255-(dBInd*20))]
                        
                        Ind1 = list(TXValues).index(0)
                        Ind2 = list(TXValues).index(
                                           int(DataInfo['Audio']['SoundPulseDur']*1000))
                        
                        Axes[dBInd].axvspan(TXValues[Ind1], TXValues[Ind2], 
                                           color='k', alpha=0.3, lw=0, 
                                           label=SpanLabel)
                        
                        Axes[dBInd].plot(TXValues, ABR, color=Colors[0], 
                                         label=LineLabel)
                        
                        Plot.Set(Ax=Axes[dBInd], AxArgs=AxArgs)
                        Axes[dBInd].legend(loc='lower right', frameon=False)
                        # Axes[dBInd].spines['bottom'].set_visible(False)
                        # Axes[dBInd].spines['left'].set_bounds(round(0), round(1))
                        # Axes[dBInd].xaxis.set_ticks_position('none')
                        # Axes[dBInd].set_ylabel(YLabel)
                        # Axes[dBInd].set_ylim(min(YLim), max(YLim))
                        
                    Axes[-1].spines['bottom'].set_visible(True)
                    Axes[-1].set_xlabel(AxArgs['xlabel'])
                    Axes[-1].spines['bottom'].set_bounds(round(0), round(1))
                    Plot.Set(Fig=Fig, FigTitle=FigTitle, HideControls=True, Tight=False)
                    
                    if Save:
                        os.makedirs(FigPath, exist_ok=True)    # Figs folder
                        FigName = ''.join([FigPath, '/', 
                                           '-'.join(DataInfo['InfoFile'].split('-')[:-1]), 
                                           '-', S,  '_', DV, 'DV_', F, 'Hz_', T])
                        for E in Ext: Fig.savefig(FigName+'.'+E, format=E)
                        
    
    if Visible: plt.show()
    return(None)


def Triang3D(AnalysisFile, FileName, Azimuth=-110, Elevation=50, Visible=True):
    """ 
    This function will plot the data from ../*Analysis.hdf5. Make sure 
    FileName is a string with the path to only one file.
    
    Also, LaTeX will render all the text in the plots (see SetPlot function). 
    For this, make sure you have a working LaTex installation and dvipng 
    package installed.
    """
    ABRs, XValues = Hdf5.ABRsLoad(AnalysisFile)    
    os.makedirs('Figs', exist_ok=True)    # Figs folder
    
    Params = Plot.Set(Params=True)
    import matplotlib.tri as mtri
    from matplotlib import rcParams; rcParams.update(Params)
    from matplotlib import pyplot as plt
    from mpl_toolkits.mplot3d.axes3d import Axes3D
    
    
    print('Plotting...')
    Colormaps = [plt.get_cmap('Reds'), plt.get_cmap('Blues')]
    for Stim in ABRs.keys():
        for DVCoord in ABRs[Stim].keys():
            for Freq in ABRs[Stim][DVCoord].keys():
                for Trial in ABRs[Stim][DVCoord][Freq].keys():
                    FigTitle = ''.join([Freq, 'Hz, ', DVCoord, 
                                        'DV, trial ', Trial])
                    YLabel = 'Intensity [dBSPL]'; XLabel = 'Time [ms]'
                    ZLabel = 'Voltage [mV]'
                    
                    Fig = plt.figure()
                    Axes = Axes3D(Fig)
                    
                    Intensities = list(ABRs[Stim][DVCoord][Freq][Trial])
                    Intensities.sort(reverse=True)
                    TXValues = XValues[Stim][DVCoord][Freq][Trial][:]
                    
                    for LineIndex in range(len(Intensities)-1):
                        dB0 = Intensities[LineIndex]
                        dB1 = Intensities[LineIndex+1]
                        
                        ABR0 = ABRs[Stim][DVCoord][Freq][Trial][dB0][:]
                        ABR1 = ABRs[Stim][DVCoord][Freq][Trial][dB1][:]
                        
                        X = np.concatenate([TXValues, TXValues])
                        Y = [int(dB0[:-2])]*len(ABR0) + [int(dB1[:-2])]*len(ABR1)
                        Z = np.concatenate([ABR0, ABR1])
                        T = mtri.Triangulation(X, Y)
                        
                        Axes.plot_trisurf(X, Y, Z, triangles=T.triangles, 
                                          cmap=Colormaps[0], edgecolor='none', 
                                          antialiased=False, shade=False)
                    
                    Axes.locator_params(tight=True)
                    Axes.set_xlabel(XLabel); Axes.set_ylabel(YLabel)
                    Axes.set_zlabel(ZLabel)
                    Axes.grid(False)
                    Axes.view_init(Elevation, Azimuth)
                    
                    Fig.suptitle(FigTitle)#; Fig.tight_layout(); 
                    
                    FigName = ''.join(['Figs/', FileName[:-15], '-ABRs3D_', Stim, 
                                       '_', DVCoord, 'DV_', Freq, 'Hz_', Trial, 
                                       '.svg'])
                    Fig.savefig(FigName, format='svg')
    
    if Visible: plt.show()
    
    return(None)


def ABRThresholdsTxtToHdf5(Override={}):
    if 'AnalysisFile' in Override: AnalysisFile =  Override['AnalysisFile']
    else: AnalysisFile = glob('*.hdf5')[0]
    
    if 'FileName' in Override: FileName =  Override['FileName']
    else: FileName = glob('*.txt')[0]
    
    with open(FileName, 'r') as F:
        Lines = [Line for Line in F]
    
    Lines = [Line.split() for Line in Lines]
    Exps = Lines[0][:]; del(Lines[0])
    Freqs = Lines[0][:]; del(Lines[0])
    
    ABRThresholds = {Exp: {} for Exp in Exps}
    for EInd, Exp in enumerate(Exps):
        for FInd, Freq in enumerate(Freqs):
            ABRThresholds[Exp][Freq] = float(Lines[EInd][FInd])
    
    for Key in ABRThresholds:
        Hdf5.DictWrite(ABRThresholds[Key], '/ABRThresholds/'+Key, AnalysisFile, Attrs=False)
    
    return(None)


