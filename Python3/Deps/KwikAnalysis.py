# -*- coding: utf-8 -*-
"""
    Copyright (C) 2015  T. Malfatti
    
    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.
    
    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.
    
    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.

"""

##import h5py
##import Kwik
#from numbers import Number

import Hdf5F
import numpy as np
import os
from datetime import datetime
from glob import glob
from multiprocessing import Process
from scipy import io, signal
from subprocess import call


## Lower-level functions
def CallWaveClus(Rate, Path):
    print('Clustering spikes...')
    
    MLab = '/home/malfatti/Software/Programs/MatLabR2015a/bin/matlab'
    CmdCd = 'cd ' + Path + '; '
    CmdCluster = 'try, SampleRate='+str(int(Rate))+'; Get_spikes;' +\
                 'Do_clustering(SampleRate); end; quit'
    call([MLab, '-r', CmdCd+CmdCluster])
    
    print('Done clustering.')
    return(None)


def GetProc(Raw, Board):
    print('Get proc no. for', Board, 'board... ', end='')
    ProcChs = {Proc: len(Raw[Proc]['data']['0'][1,:]) 
               for Proc in Raw.keys()}
    
    for Proc, Chs in ProcChs.items():
        if Chs == max(ProcChs.values()): OEProc = Proc
        else: RHAProc = Proc
    
    if 'RHAProc' not in locals(): RHAProc = OEProc
    
    if Board == 'OE': Proc = OEProc
    elif Board == 'RHA': Proc = RHAProc
    else: print("Choose Board as 'OE' or 'RHA'."); return(None)
    
    print('Done.')
    return(OEProc, RHAProc, Proc)


def GetProbeChOrder(ProbeTip, ProbeHead, Connector):
    """
    Get probe channels order. It doesn't matter what logic you follow to order 
    your connector channels, but you MUST follow the same logic for your probe 
    head.
    
    I the probe tip channels top-down or bottom-up, the resulting 
    channel map will be ordered accordingly.
    
    Example:
        CustomAdaptor = [5, 6, 7, 8, 9, 10 ,11, 12, 13, 14, 15, 16, 1, 2, 3, 4]
        RHAHeadstage = [16, 15, 14, 13, 12, 11, 10, 9, 1, 2, 3, 4, 5, 6, 7, 8]
        A16 = {'ProbeTip': [9, 8, 10, 7, 13, 4, 12, 5, 15, 2, 16, 1, 14, 3, 11, 6],
               'ProbeHead': [8, 7, 6, 5, 4, 3, 2, 1, 9, 10, 11, 12, 13, 14, 15, 16]}
        
        ChannelMap = GetProbeChOrder(A16['ProbeTip'], A16['ProbeHead'], CustomAdaptor)
    """
    print('Get probe channel order... ', end='')
    ChNo = len(ProbeTip)
    ChMap = [0]*ChNo
    
    for Ch in range(ChNo):
        TipCh = ProbeTip[Ch] # What channel should be the Ch
        HeadCh = ProbeHead.index(TipCh) # Where Ch is in ProbeHead
        ChMap[Ch] = Connector[HeadCh] # Channels in depth order
    
    print('Done.')
    return(ChMap)


def FilterSignal(Signal, Rate, Frequency, FilterOrder=4, Type='bandpass'):
    if Type not in ['bandpass', 'lowpass', 'highpass']:
        print("Choose 'bandpass', 'lowpass' or 'highpass'.")
    
    elif len(Frequency) not in [1, 2]:
        print('Frequency must have 2 elements for bandpass; or 1 element for \
        lowpass or highpass.')
    
    else:
        passband = [_/(Rate/2) for _ in Frequency]
        f2, f1 = signal.butter(FilterOrder, passband, Type)
        Signal = signal.filtfilt(f2, f1, Signal, padtype='odd', padlen=0)
        
        return(Signal)


def GetRecKeys(Raw, Events, AnalogTTLs):
    for Processor in Raw.keys():
        if '0' not in list(Raw[Processor]['data'].keys()):
            print('Rec numbers are wrong. Fixing...')
            for iKey in Raw[Processor].keys():
                Recs = list(Raw[Processor][iKey].keys())
                Recs = [int(_) for _ in Recs]; Min = min(Recs)
                
                for Key in Recs:
                    Raw[Processor][iKey][str(Key-Min)] = \
                        Raw[Processor][iKey].pop(str(Key))
                
                if AnalogTTLs:
                    return(Raw)
                else:
                    EventRec = Events['TTLs']['recording'][:]
                    for _ in range(len(EventRec)): 
                        EventRec[_] = EventRec[_] - Min
                    return(Raw, EventRec)
            
            print('Fixed.')
    
    else:
        if AnalogTTLs: return(Raw)
        else: EventRec = Events['TTLs']['recording']; return(Raw, EventRec)


def GetTTLInfo(Events, EventRec, TTLCh):
    print('Get TTL data...')
    EventID = Events['TTLs']['user_data']['eventID']
    EventCh = Events['TTLs']['user_data']['event_channels']
    EventSample = Events['TTLs']['time_samples']

#    TTLChs = np.nonzero(np.bincount(EventCh))[0]
    TTLRecs = np.nonzero(np.bincount(EventRec))[0]
    TTLRecs = ["{0:02d}".format(_) for _ in TTLRecs]
    TTLsPerRec = {Rec: [EventSample[_] for _ in range(len(EventRec)) 
                         if EventRec[_] == int(Rec)
                         and EventCh[_] == TTLCh-1 
                         and EventID[_] == 1]
                  for Rec in TTLRecs}
#    TTLRising = Kwik.get_rising_edge_times(Files['kwe'], TTLCh-1)
    
    return(TTLsPerRec)


def QuantifyTTLsPerRec(Raw, Rec, AnalogTTLs, ChTTL=-1, Proc='', TTLsPerRec=[], 
                       Rate=[]):
    print('Get TTL timestamps... ', end='')
    if AnalogTTLs:
        TTLCh = Raw[Proc]['data'][Rec][:, ChTTL-1]
        Threshold = max(TTLCh)/2
        TTLs = []
        for _ in range(1, len(TTLCh)):
            if TTLCh[_] > Threshold:
                if TTLCh[_-1] < Threshold: TTLs.append(_)
        
        print('Done.')
        return(TTLs)
    else:
#        TTLNo = [0]
#        for _ in range(1, len(TTLsPerRec)+1):
#            TTLNo = TTLNo + [len(TTLsPerRec[_-1]) + TTLNo[-1]]
#        TTLNo = [0] + [TTLNo[_]-1 for _ in range(1, len(TTLNo))]
#        
#        if Rec == 0:
#            sTTLNo = 0
#        else:
#            sTTLNo = TTLNo[Rec] + 1
#        
#        return(TTLNo, sTTLNo)
        RawTime = [_*Rate for _ in Raw['timestamps'][int(Rec)]]
        TTLs = TTLsPerRec[Rec]
        
        print('Done.')
        return(RawTime, TTLs)


def RemoveDateFromFolderName():
    RenameFolders = input('Rename folders in KwikFiles/* (BE CAREFUL)? [y/N] ')
    if RenameFolders in ['y', 'Y', 'yes', 'Yes', 'YES']:
        DirList = glob('KwikFiles/*'); DirList.sort()
        for FolderName in DirList:
            NewFolderName = ''.join([FolderName[:10], FolderName[21:]])
            NewFolderName = NewFolderName.replace("-", "")
            os.rename(FolderName, NewFolderName)
            print(FolderName, ' moved to ', NewFolderName)
        del(RenameFolders, DirList, FolderName, NewFolderName)    
    
    return(None)


def SepSpksPerCluster(Clusters, Ch):
    Classes = np.unique(Clusters['ClusterClass'])
    Dict = {}                    
#    Dict['Spks'] = [[] for _ in range(len(Classes))]
    Dict['Spks'] = {}
    
    print(Ch+':', str(len(Classes)), 'clusters:')
    for Class in Classes:
        ClassIndex = Clusters['ClusterClass'] == Class
        
        SpkNo = len(Clusters['Spikes'][ClassIndex,:])
        if not SpkNo: continue
        
        Class = "{0:02d}".format(int(Class))
        Dict['Spks'][Class] =  Clusters['Spikes'][ClassIndex,:][:]
    
        print('    Class', Class, '-', str(SpkNo), 'spikes.')
    
    if len(Dict): return(Dict)
    else: return({})


def SetPlot(Backend='TkAgg', AxesObj=(), FigObj=(), FigTitle='', Params=False, 
            Plot=False, Axes=False):
    if Params:
        print('Set plot parameters...')
        Params = {'backend': Backend,
                  'text.usetex': True, 'text.latex.unicode': True,
    #              'text.latex.preamble': '\\usepackage{siunitx}',
                  
                  'font.family': 'serif', 'font.serif': 'Computer Modern Roman',
                  'axes.titlesize': 'medium', 'axes.labelsize': 'medium',
                  'xtick.labelsize': 'small', 'xtick.direction': 'out',
                  'ytick.labelsize': 'small', 'ytick.direction': 'out',
                  'legend.fontsize': 'small', 'legend.labelspacing': 0.4,
                  'figure.titlesize': 'large', 'figure.titleweight': 'normal',
                  
                  'image.cmap': 'cubehelix', 'savefig.transparent': True,
                  'svg.fonttype': 'path'}
        return(Params)
    
    elif Plot:
        print('Set plot figure...')
        FigObj.suptitle(FigTitle); FigObj.tight_layout(); 
        FigObj.subplots_adjust(top=0.925)
    
    elif Axes:
        print('Set plot axes...')
        AxesObj.spines['right'].set_visible(False)
        AxesObj.spines['top'].set_visible(False)
        AxesObj.yaxis.set_ticks_position('left')
        AxesObj.xaxis.set_ticks_position('bottom')
        AxesObj.locator_params(tight=True)
    
    else: print("'Params', 'Plot' or 'Axes' must be True.")
    
    return(None)


def FixTTLs(Array, TTLsToFix):
    for TTL in TTLsToFix:
        nInd = np.random.randint(1, 100)
        while nInd == TTL: nInd = np.random.randint(0, 100)
        
        print('TTL', str(TTL), 'was replaced by', str(nInd))
        Array[TTL] = Array[nInd]
    
    return(Array)


def SliceData(Data, Proc, Rec, TTLs, DataCh, NoOfSamplesBefore, 
              NoOfSamplesAfter, NoOfSamples, AnalogTTLs, RawTime=[]):
    print('Slicing data around TTL...')
    Array = [[0 for _ in range(NoOfSamples)] for _ in range(len(TTLs))]
    TTLsToFix = []
    
    for TTL in range(len(TTLs)):
        if AnalogTTLs: TTLLoc = int(TTLs[TTL])
        else: TTLLoc = int(RawTime.index(TTLs[TTL]))#)/Rate)
        
        Start = TTLLoc-NoOfSamplesBefore
        End = TTLLoc+NoOfSamplesAfter
        
        if Start < 0: Start = 0; End = End+(Start*-1); TTLsToFix.append(TTL)
        
        Array[TTL] = Data[Proc]['data'][Rec][Start:End, DataCh[0]-1] * \
                     Data[Proc]['channel_bit_volts'][Rec][DataCh[0]-1] # in mV
        
        if len(Array[TTL]) != End-Start: TTLsToFix.append(TTL)
            
    Array = FixTTLs(Array, TTLsToFix)
    
    print('Done.')
    return(Array)


def UnitPlotPerCh(ChDict, Ch, XValues, PulseDur, FigName, FigTitle):
    ClusterNo = len(ChDict['Spks'])
    if ClusterNo == 0: print(Ch, 'had no spikes :('); return(None)
    
    Params = SetPlot(Backend='Agg', Params=True)
    from matplotlib import rcParams; rcParams.update(Params)
    from matplotlib import pyplot as plt
    
    Fig, Axes = plt.subplots(ClusterNo,2, figsize=(6, 3*ClusterNo))
    SpksYLabel = 'Voltage [mV]'; SpksXLabel = 'Time [ms]'
    PSTHYLabel = 'Number of spikes in channel'; PSTHXLabel = 'Time [ms]'
    SpanLabel = 'Sound pulse'
    
    for Class in ChDict['Spks'].keys():
        SpkNo = len(ChDict['Spks'][Class])
        print(str(SpkNo), 'Spks in cluster', Class)
        print('Max of', max(ChDict['PSTH'][Class]),  
              'Spks in PSTH')
        
#                        if not SpkNo:
#                            print('No Spk data on cluster', str(Cluster) + '. Skipping...')
#                            Thrash[SKey+FKey+RKey+Ch] = Units[SKey][FKey][RKey][Ch].copy()
#                            continue
        
#                PSTHPeak = max(Units[SKey][FKey][RKey]['PSTH'][Cluster])
#                PSTHMean = np.mean(Units[SKey][FKey][RKey]['PSTH'][Cluster])
#        if max(UnitRec[Key]['PSTH'][Cluster]) < 4: 
#            print('No peaks in PSTH. Skipping cluster', str(Cluster), '...')
#            continue
        
        if SpkNo > 100: SpkNo = np.arange(100); np.random.shuffle(SpkNo)
        else: SpkNo = np.arange(SpkNo)
        
        for Spike in SpkNo:
            if ClusterNo == 1: Axes[0].plot(ChDict['Spks'][Class][Spike], 'r')
            else: Axes[int(Class)-1][0].plot(ChDict['Spks'][Class][Spike], 'r')
        
        if ClusterNo == 1:
#                        Axes[0].set_title('Peak='+str(PSTHPeak)+' Mean='+str(PSTHMean))
            Axes[0].plot(np.mean(ChDict['Spks'][Class], axis=0), 'k')
            Axes[1].bar(XValues, ChDict['PSTH'][Class], 'r')
            
            Ind1 = list(XValues).index(0)
            Ind2 = list(XValues).index(int(PulseDur*1000))
            
            Axes[1].axvspan(XValues[Ind1], XValues[Ind2], color='k', alpha=0.3, 
                            lw=0, label=SpanLabel)
            
            SetPlot(AxesObj=Axes[0], Axes=True)
            SetPlot(AxesObj=Axes[1], Axes=True)
            Axes[0].set_ylabel(SpksYLabel); Axes[0].set_xlabel(SpksXLabel)
            Axes[1].set_ylabel(PSTHYLabel); Axes[1].set_xlabel(PSTHXLabel)
            
        else:
#                    Axes[Cluster][0].set_title('Peak='+str(PSTHPeak)+' Mean='+\
#                                               str(PSTHMean)+' Std='+str(PSTHStd))
            Axes[int(Class)-1][0].plot(np.mean(ChDict['Spks'][Class], axis=0), 'k')
            Axes[int(Class)-1][1].bar(XValues, ChDict['PSTH'][Class], 'r')
            
            Ind1 = list(XValues).index(0)
            Ind2 = list(XValues).index(int(PulseDur*1000))
            
            Axes[int(Class)-1][1].axvspan(XValues[Ind1], XValues[Ind2], 
                                          color='k', alpha=0.3, lw=0, 
                                          label=SpanLabel)
            
            SetPlot(AxesObj=Axes[int(Class)-1][0], Axes=True)
            SetPlot(AxesObj=Axes[int(Class)-1][1], Axes=True)
            Axes[int(Class)-1][0].set_ylabel(SpksYLabel)
            Axes[int(Class)-1][0].set_xlabel(SpksXLabel)
            Axes[int(Class)-1][1].set_ylabel(PSTHYLabel)
            Axes[int(Class)-1][1].set_xlabel(PSTHXLabel)
    
    SetPlot(FigObj=Fig, FigTitle=FigTitle, Plot=True)
    print('Writing to', FigName+'... ', end='')
    Fig.savefig(FigName, format='svg')
    print('Done.')
    return(None)


## Higher-level functions
def ABRAnalysis(FileName, ABRCh=[1], ABRTTLCh=1, ABRTimeBeforeTTL=0, 
                ABRTimeAfterTTL=12, FilterFreq=[300, 3000], FilterOrder=4, 
                StimType='Sound', AnalogTTLs=False, Board='OE', Override={}):
    """
    Analyze ABRs from data in Kwik format. A '*ABRs.hdf5' file will be saved 
    in cwd, containing:
        - ABRs group, where data will be saved as 
          ABRs[Ear][Freq][AmpF][DVCoord][Trial], where:
              Ear = 0 (right) or 1 (left)
              Freq = index of DataInfo['NoiseFrequency']
              AmpF = index of DataInfo['SoundAmpF']['Freq']
              DVCoord = string with DV coordinate at the moment of recording
              Trial = Trial number - 1 (so if it is one trial, Trial=0)
              
        - XValues array, for x axis of the plot;
        
        - DataInfo dict, where all info will be saved.
    
    For this function to work:
        - The Kwik folders must be in 'KwikFiles/';
        - There must be a *Exp.hdf5 file containing all experimental settings 
          (see Python3/SoundBoardControl/SoundAndLaserStimulation.py, 1st cell);
    """
    
    print('Load DataInfo...')
    DirList = glob('KwikFiles/*'); DirList.sort()
    DataInfo = Hdf5F.LoadDict('/DataInfo', FileName)
    
    AnalysisFile = '../' + DataInfo['AnimalName'] + '-Analysis.hdf5'
    Now = datetime.now().strftime("%Y%m%d%H%M%S")
    Here = os.getcwd().split(sep='/')[-1]
    Group = Here + '-ABRs_' + Now
        
    for Stim in StimType:
        if Override != {}: 
            if 'Stim' in Override.keys(): Stim = Override['Stim']
        
        Exps = Hdf5F.LoadExpPerStim(Stim, DirList, FileName)
        
        for RecFolder in Exps:
            ABRs = {}; Info = {}
            ExpInfo = Hdf5F.ExpExpInfo(RecFolder, DirList, FileName)
            
            if AnalogTTLs: 
                Raw, _, Files = Hdf5F.LoadOEKwik(RecFolder, AnalogTTLs)
            else: 
                Raw, Events, _, Files = Hdf5F.LoadOEKwik(RecFolder, AnalogTTLs)
            
            OEProc, RHAProc, ABRProc = GetProc(Raw, Board)
            
            if AnalogTTLs: Raw = GetRecKeys(Raw, [0], AnalogTTLs)
            else:
                Raw, EventRec = GetRecKeys(Raw, Events, AnalogTTLs)
                TTLsPerRec = GetTTLInfo(Events, EventRec, ABRTTLCh)
            
            Rate = Raw[OEProc]['info']['0']['sample_rate']
            NoOfSamplesBefore = ABRTimeBeforeTTL*int(Rate*10**-3)
            NoOfSamplesAfter = ABRTimeAfterTTL*int(Rate*10**-3)
            NoOfSamples = NoOfSamplesBefore + NoOfSamplesAfter
            
            Info['Frequency'] = ''.join([
                            str(DataInfo['NoiseFrequency'][ExpInfo['Hz']][0]),
                            '-',
                            str(DataInfo['NoiseFrequency'][ExpInfo['Hz']][1])])
            
            Info['XValues'] = (range(-NoOfSamplesBefore, 
                                     NoOfSamples-NoOfSamplesBefore)/Rate)*10**3
            
            for Rec in Raw[OEProc]['data'].keys():
                print('Slicing and filtering ABRs Rec ', str(Rec), '...')
                
                if AnalogTTLs:
                    TTLs = QuantifyTTLsPerRec(Raw, Rec, AnalogTTLs, ABRTTLCh, 
                                              OEProc)
                    ABR = SliceData(Raw, ABRProc, Rec, TTLs, ABRCh, 
                                    NoOfSamplesBefore, NoOfSamplesAfter, 
                                    NoOfSamples, AnalogTTLs)
                else:
                    RawTime, TTLs = QuantifyTTLsPerRec(Raw, Rec, AnalogTTLs, 
                                                       TTLsPerRec=TTLsPerRec)
                    ABR = SliceData(Raw, ABRProc, Rec, TTLs, ABRCh, 
                                    NoOfSamplesBefore, NoOfSamplesAfter, 
                                    AnalogTTLs, RawTime)
                
                for TTL in range(len(TTLs)):
                    ABR[TTL] = FilterSignal(ABR[TTL], Rate, [min(FilterFreq)], 
                                            FilterOrder, 'highpass')
                
                ABR = np.mean(ABR, axis=0)
                ABR = FilterSignal(ABR, Rate, [max(FilterFreq)], FilterOrder, 
                                   'lowpass')
                
                dB = str(DataInfo['Intensities'][int(Rec)]) + 'dB'
                ABRs[dB] = ABR[:]; del(ABR)
            
            Path = Stim+'/'+ExpInfo['DVCoord']+'/'+Info['Frequency']
            Hdf5F.WriteABR(ABRs, Info['XValues'], Group, Path, AnalysisFile)
    
    return(None)

    
def ABRPlot(AnalysisFile, FileName, Visible=False):
    """ 
    This function will plot the data from ../*Analysis.hdf5. Make sure 
    FileName is a string with the path to only one file.
    
    Also, LaTeX will render all the text in the plots (see SetPlot function). 
    For this, make sure you have a working LaTex installation and dvipng 
    package installed.
    """
    ABRs, XValues = Hdf5F.LoadABRs(AnalysisFile)
    DataInfo = Hdf5F.LoadDict('/DataInfo', FileName)
    
    os.makedirs('Figs', exist_ok=True)    # Figs folder
    
    Params = SetPlot(Params=True)
    from matplotlib import rcParams; rcParams.update(Params)
    from matplotlib import pyplot as plt
    
    
    print('Plotting...')
    Colormaps = [plt.get_cmap('Reds'), plt.get_cmap('Blues')]
    for Stim in ABRs.keys():
        for DVCoord in ABRs[Stim].keys():
            for Freq in ABRs[Stim][DVCoord].keys():
                for Trial in ABRs[Stim][DVCoord][Freq].keys():
                    YLim = []
                    
                    for ABR in ABRs[Stim][DVCoord][Freq][Trial].values():
                        YLim.append(max(ABR)); YLim.append(min(ABR))
                    
                    Intensities = list(ABRs[Stim][DVCoord][Freq][Trial])
                    Intensities.sort(reverse=True)
                    TXValues = XValues[Stim][DVCoord][Freq][Trial][:]
                    
                    Fig, Axes = plt.subplots(len(Intensities), sharex=True, 
                                             figsize=(8, 1.5*len(Intensities)))
                    
                    for dB, ABR in ABRs[Stim][DVCoord][Freq][Trial].items():
                        FigTitle = ''.join([Freq, 'Hz, ', DVCoord, 
                                            'DV, trial ', Trial])
                        YLabel = 'Voltage [mV]'; XLabel = 'Time [ms]'
                        LineLabel = dB
                        SpanLabel = 'Sound pulse'
                        
                        dBInd = Intensities.index(dB)
                        Colors = [Colormaps[0](255-(dBInd*20)), 
                                  Colormaps[1](255-(dBInd*20))]
                        
                        Ind1 = list(TXValues).index(0)
                        Ind2 = list(TXValues).index(
                                           int(DataInfo['SoundPulseDur']*1000))
                        
                        Axes[dBInd].axvspan(TXValues[Ind1], TXValues[Ind2], 
                                           color='k', alpha=0.3, lw=0, 
                                           label=SpanLabel)
                        
                        Axes[dBInd].plot(TXValues, ABR, color=Colors[0], 
                                         label=LineLabel)
                        
                        SetPlot(AxesObj=Axes[dBInd], Axes=True)
                        Axes[dBInd].legend(loc='lower right', frameon=False)
                        Axes[dBInd].spines['bottom'].set_visible(False)
                        Axes[dBInd].spines['left'].set_bounds(round(0), round(1))
                        Axes[dBInd].xaxis.set_ticks_position('none')
                        Axes[dBInd].set_ylabel(YLabel)
                        Axes[dBInd].set_ylim(min(YLim), max(YLim))
                        
                    Axes[-1].spines['bottom'].set_visible(True)
                    Axes[-1].set_xlabel(XLabel)
                    Axes[-1].spines['bottom'].set_bounds(round(0), round(1))
                    SetPlot(FigObj=Fig, FigTitle=FigTitle, Plot=True)
                    
                    FigName = ''.join(['Figs/', FileName[:-15], '-ABRs_', Stim, 
                                       '_', DVCoord, 'DV_', Freq, 'Hz_', Trial, 
                                       '.svg'])
                    Fig.savefig(FigName, format='svg')
    
    if Visible: plt.show()
    return(None)


def ABRPlot3D(AnalysisFile, FileName, Azimuth=-110, Elevation=50, Visible=True):
    """ 
    This function will plot the data from ../*Analysis.hdf5. Make sure 
    FileName is a string with the path to only one file.
    
    Also, LaTeX will render all the text in the plots (see SetPlot function). 
    For this, make sure you have a working LaTex installation and dvipng 
    package installed.
    """
    ABRs, XValues = Hdf5F.LoadABRs(AnalysisFile)    
    os.makedirs('Figs', exist_ok=True)    # Figs folder
    
    Params = SetPlot(Params=True)
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


def ClusterizeSpks(FileName, StimType=['Sound'], AnalogTTLs=False, Board='OE', 
                   Override={}):
    print('Load DataInfo...')
    DirList = glob('KwikFiles/*'); DirList.sort()
    
    CustomAdaptor = [5, 6, 7, 8, 9, 10 ,11, 12, 13, 14, 15, 16, 1, 2, 3, 4]
    A16 = {'ProbeTip': [9, 8, 10, 7, 13, 4, 12, 5, 15, 2, 16, 1, 14, 3, 11, 6],
           'ProbeHead': [8, 7, 6, 5, 4, 3, 2, 1, 9, 10, 11, 12, 13, 14, 15, 16]}
    ChannelMap = GetProbeChOrder(A16['ProbeTip'], A16['ProbeHead'], CustomAdaptor)
    
    for Stim in StimType:
        if Override != {}: 
            if 'Stim' in Override.keys():
                Stim = Override['Stim']
        
        Exps = Hdf5F.LoadExpPerStim(Stim, DirList, FileName)
        
        for FInd, RecFolder in enumerate(Exps):        
            if AnalogTTLs: Raw, _, Files = Hdf5F.LoadOEKwik(RecFolder, AnalogTTLs)
            else: Raw, Events, _, Files = Hdf5F.LoadOEKwik(RecFolder, AnalogTTLs)
            
            OEProc = GetProc(Raw, Board)[0]
            
            Path = os.getcwd() + '/' + RecFolder +'/SepCh/'
            os.makedirs(Path, exist_ok=True)
            
            Rate = Raw[OEProc]['info']['0']['sample_rate']
            
            Clusters = {}
            for Rec in Raw[OEProc]['data'].keys():
                if Override != {}: 
                    if 'Rec' in Override.keys():
                        Rec = Override['Rec']
                RecS = "{0:02d}".format(int(Rec))
                
                print('Separating channels according to ChannelMap...')
                Data = [Raw[OEProc]['data'][Rec][:, _-1] * 
                        Raw[OEProc]['channel_bit_volts'][Rec][_-1] * 1000
                        for _ in ChannelMap]
                                
                print('Writing files for clustering... ', end='')
                FileList = []
                for Ind, Ch in enumerate(Data):
                    MatName = 'Exp' + Files['100_kwd'][-13:-8] + '_' + \
                              RecS + '-Ch' + "{0:02d}".format(Ind+1) + '.mat'
                    
                    FileList.append(MatName)
                    io.savemat(Path+MatName, {'data': Ch})
                
                TxtFile = open(Path+'Files.txt', 'w')
                for File in FileList: TxtFile.write(File+'\n')
                TxtFile.close()
                print('Done.')
                
                CallWaveClus(Rate, Path)
                
                ClusterList = glob(Path+'times_*'); ClusterList.sort()
                ClusterList = [_ for _ in ClusterList if _[-11:-9] == RecS]
                
                Clusters[RecS] = {}
                for File in ClusterList:
                    Ch = File[-8:-4]
                    
                    ClusterFile = io.loadmat(File)
                    Clusters[RecS][Ch] = {}
                    Clusters[RecS][Ch]['ClusterClass'] = ClusterFile['cluster_class'][:, 0]
                    Clusters[RecS][Ch]['Timestamps'] = ClusterFile['cluster_class'][:, 1]
                    Clusters[RecS][Ch]['Spikes'] = ClusterFile['spikes'][:]
                    
#                    Clusters[RecS][Ch]['Info'] = {}
#                    Clusters[RecS][Ch]['Info']['Parameters'] = ClusterFile['par']
#                    Clusters[RecS][Ch]['Info']['InSpk'] = ClusterFile['inspk'][:]
                
                if Override != {}: 
                    if 'Rec' in Override.keys(): break
            
            Hdf5F.WriteClusters(Clusters, Path[:-6]+'SpkClusters.hdf5')
            ToDelete = glob(Path+'*')
            for File in ToDelete: os.remove(File)
            os.removedirs(Path)


def GPIASAnalysis(RecFolder, FileName, GPIASCh=1, GPIASTTLCh=1, 
                  GPIASTimeBeforeTTL=50, GPIASTimeAfterTTL=150, 
                  FilterFreq=[70, 400], FilterOrder=4, AnalogTTLs=False, 
                  Board='OE'):
    
    print('set paths...')
    DirList = glob('KwikFiles/*'); DirList.sort()
    RecFolder = DirList[RecFolder-1]
    
    DataInfo = Hdf5F.LoadDict('/DataInfo', FileName)    
    AnalysisFile = '../' + DataInfo['AnimalName'] + '-Analysis.hdf5'
    
    for Path in ['Freqs', 'FreqOrder', 'FreqSlot']:
        DataInfo[Path] = Hdf5F.LoadDataset('/DataInfo/'+Path, FileName)
        
    if AnalogTTLs: Raw, _, Files = Hdf5F.LoadOEKwik(RecFolder, AnalogTTLs)
    else: Raw, Events, _, Files = Hdf5F.LoadOEKwik(RecFolder, AnalogTTLs)
    
    if AnalogTTLs: Raw = GetRecKeys(Raw, [0], AnalogTTLs)
    else:
        Raw, EventRec = GetRecKeys(Raw, Events, AnalogTTLs)
        TTLsPerRec = GetTTLInfo(Events, EventRec, GPIASTTLCh)
    
    OEProc = GetProc(Raw, Board)[0]
    
    Rate = Raw[OEProc]['info']['0']['sample_rate']
    NoOfSamplesBefore = int(round((GPIASTimeBeforeTTL*Rate)*10**-3))
    NoOfSamplesAfter = int(round((GPIASTimeAfterTTL*Rate)*10**-3))
    NoOfSamples = NoOfSamplesBefore + NoOfSamplesAfter
    
    XValues = (range(-NoOfSamplesBefore, NoOfSamples-NoOfSamplesBefore)
               /Rate)*10**3
    
    GPIAS = {''.join([str(Freq[0]), '-', str(Freq[1])]): {}
             for Freq in DataInfo['NoiseFrequency']}
    
    for Freq in GPIAS.keys():
        GPIAS[Freq]['NoGap'] = []; GPIAS[Freq]['Gap'] = []
    
    for Rec in Raw[OEProc]['data'].keys():
        print('Slicing and filtering Rec ', Rec, '...')
        Freq = DataInfo['FreqOrder'][int(Rec)][0]; 
        Trial = DataInfo['FreqOrder'][int(Rec)][1];
        
        SFreq = ''.join([str(DataInfo['NoiseFrequency'][Freq][0]), '-', 
                         str(DataInfo['NoiseFrequency'][Freq][1])])
        
        if Trial%2 == 0: STrial = 'NoGap'
        else: STrial = 'Gap'
        
        if AnalogTTLs:
            TTLs = QuantifyTTLsPerRec(Raw, Rec, AnalogTTLs, GPIASTTLCh, OEProc)
            GD = SliceData(Raw, OEProc, Rec, TTLs, GPIASCh, NoOfSamplesBefore, 
                           NoOfSamplesAfter, NoOfSamples, AnalogTTLs)
        else:
            RawTime, TTLs = QuantifyTTLsPerRec(Raw, Rec, AnalogTTLs, 
                                               TTLsPerRec=TTLsPerRec)
            GD = SliceData(Raw, OEProc, Rec, TTLs, GPIASCh, NoOfSamplesBefore, 
                           NoOfSamplesAfter, NoOfSamples, AnalogTTLs, RawTime)
        
        if GPIAS[SFreq][STrial] == []: GPIAS[SFreq][STrial] = GD[:]
        else: 
            # Cumulative moving average
            ElNo = len(GPIAS[SFreq][STrial])
            GPIAS[SFreq][STrial] = ((np.mean(GPIAS[SFreq][STrial], axis=0)
                                     *ElNo) + GD[:])/(ElNo+1)
    
    for Freq in GPIAS.keys():
        # Hilbert
        GPIAS[Freq]['NoGap'] = abs(signal.hilbert(GPIAS[Freq]['NoGap']))
        GPIAS[Freq]['Gap'] = abs(signal.hilbert(GPIAS[Freq]['Gap']))
        
        # Bandpass filter
        GPIAS[Freq]['NoGap'] = FilterSignal(GPIAS[Freq]['NoGap'], Rate, 
                                            FilterFreq, FilterOrder, 
                                            'bandpass')
        GPIAS[Freq]['Gap'] = FilterSignal(GPIAS[Freq]['Gap'], Rate, FilterFreq, 
                                          FilterOrder, 'bandpass')
        
        # SavGol smooth
        GPIAS[Freq]['NoGap'] = signal.savgol_filter(GPIAS[Freq]['NoGap'], 5, 2, 
                                                    mode='nearest')
        GPIAS[Freq]['Gap'] = signal.savgol_filter(GPIAS[Freq]['Gap'], 5, 2, 
                                                  mode='nearest')
        
        # V to mV
        GPIAS[Freq]['NoGap'] = GPIAS[Freq]['NoGap']*1000
        GPIAS[Freq]['Gap'] = GPIAS[Freq]['Gap']*1000
    
    AnalysisFile = '../' + DataInfo['AnimalName'] + '-Analysis.hdf5'
    Hdf5F.WriteGPIAS(GPIAS, RecFolder, XValues, AnalysisFile)
    
    return(None)


def GPIASPlot(AnalysisFile, FileName, Visible=False):
    os.makedirs('Figs', exist_ok=True)    # Figs folder
    
    print('Loading data...')
    ## DataInfo
    DataInfo = Hdf5F.LoadDict('/DataInfo', FileName)
    
    ## GPIAS
    GPIAS, XValues = Hdf5F.LoadGPIAS(AnalysisFile)
    
    Params = SetPlot(Params=True)
    from matplotlib import rcParams; rcParams.update(Params)
    from matplotlib import pyplot as plt
    
    print('Plotting...')
    Ind1 = list(XValues).index(0)
    Ind2 = list(XValues).index(int(DataInfo['SoundLoudPulseDur']*1000))
    
    for Freq in GPIAS.keys():
        FigTitle = Freq + ' Hz'
        LineNoGapLabel = 'No Gap'; LineGapLabel = 'Gap'
        SpanLabel = 'Sound Pulse'
        XLabel = 'time [ms]'; YLabel = 'voltage [mV]'
        
        plt.figure()
        plt.axvspan(XValues[Ind1], XValues[Ind2], color='k', alpha=0.5, 
                    lw=0, label=SpanLabel)
        plt.plot(XValues, GPIAS[Freq]['NoGap'][0,:], 
                 color='r', label=LineNoGapLabel, lw=2)
        plt.plot(XValues, GPIAS[Freq]['Gap'][0,:], 
                 color='b', label=LineGapLabel, lw=2)

        SetPlot(AxesObj=plt.axes(), Axes=True)
        SetPlot(FigObj=plt, FigTitle=FigTitle, Plot=True)
        plt.ylabel(YLabel); plt.xlabel(XLabel)
        plt.legend(loc='lower right')
        
        FigName = 'Figs/' + FileName[:-5] + '-' + Freq + '.svg'
        plt.savefig(FigName, format='svg')
    
    if Visible: plt.show()
    
    print('Done.')
    return(None)


def TTLsLatencyAnalysis(FileName, SoundCh=1, TTLSqCh=2, TTLCh=1, 
                TimeBeforeTTL=5, TimeAfterTTL=8, AnalogTTLs=False):
    print('set paths...')
    RecFolder = glob('KwikFiles/*'); RecFolder = RecFolder[0]
#    SoundCh = 0; TTLSqCh = 1 # Override
    
    TTLsLatency = {}    
    Raw, Events, _, Files = Hdf5F.LoadOEKwik(RecFolder, AnalogTTLs)    
    OEProc = GetProc(Raw, 'OE')[0]
    
    Raw, EventRec = GetRecKeys(Raw, Events, AnalogTTLs)
    TTLsPerRec = GetTTLInfo(Events, EventRec, TTLCh)
    
    Rate = Raw['info']['0']['sample_rate']
    NoOfSamplesBefore = TimeBeforeTTL*int(Rate*10**-3)
    NoOfSamplesAfter = TimeAfterTTL*int(Rate*10**-3)
    NoOfSamples = NoOfSamplesBefore + NoOfSamplesAfter
    
    XValues = list((range((NoOfSamplesBefore)*-1, 0)/Rate)*1000) + \
              list((range(NoOfSamplesAfter)/Rate)*1000)
    
    for Rec in Raw[OEProc]['data'].keys():
        RawTime, TTLs = QuantifyTTLsPerRec(Raw, Rec, AnalogTTLs, 
                                           TTLsPerRec=TTLsPerRec)
        
        SoundPulse = [[0 for _ in range(NoOfSamples)] 
                      for _ in range(len(TTLs))]
        TTLSq = SoundPulse[:]
        
        SoundPulse = SliceData(SoundPulse, Raw, OEProc, Rec, TTLs,SoundCh, 
                               NoOfSamplesBefore, NoOfSamplesAfter, 
                               NoOfSamples, AnalogTTLs, RawTime)
        TTLSq = SliceData(TTLSq, Raw, OEProc, Rec, TTLs,TTLSqCh, 
                          NoOfSamplesBefore, NoOfSamplesAfter, 
                          NoOfSamples, AnalogTTLs, RawTime)
    
        TTLSqDelay = [[0] for _ in range(len(TTLSq))]    
        for TTL in range(len(TTLSq)):
            AnalogTTLs = True
            Peaks = QuantifyTTLsPerRec(Raw, Rec, AnalogTTLs, TTLSqCh, OEProc)
            AnalogTTLs = False
            
            if len(Peaks) != 1: 
                print('Bad TTL detection. Skipping TTL...')
                TTLSqDelay[TTL] = float('NaN')
                continue
            
            TTLSqDelay[TTL] = ((Peaks[0] - NoOfSamplesBefore) / (Rate/1000))*-1
            del(Peaks)
        
        TTLsLatency[Rec] = dict((Name, eval(Name)) 
                                for Name in ['SoundPulse', 'TTLSq', 
                                             'TTLSqDelay'])
    
    Hdf5F.WriteTTLsLatency(TTLsLatency, XValues, FileName)
    
    return(None)


def TTLsLatencyPlot(FileName):
    os.makedirs('Figs', exist_ok=True)    # Figs folder
#    
#    TTLsLatency, XValues = Hdf5F.LoadTTLsLatency(FileName)
#    
#    SetPlot(Params=True)
#    
#    for _ in range(len(TTLsLatency['SoundPulse'])):
#        plt.figure(1); plt.plot(XValues, TTLsLatency['SoundPulse'][_])
#        plt.figure(2); plt.plot(XValues, TTLsLatency['TTLSq'][_])
#    
#    Hist, BinEdges = np.histogram(TTLsLatency['TTLSqDelay'], bins=200)
#    Threshold = (DataInfo['SoundPulseDur']/100)*1000
#    Threshold = 0.08
#    sIndex = min(range(len(BinEdges)), 
#                 key=lambda i: abs(BinEdges[i]-Threshold*-1))
#    eIndex = min(range(len(BinEdges)), 
#                 key=lambda i: abs(BinEdges[i]-Threshold))
#    Sum = sum(Hist); Perc = Sum/len(SoundSq) * 100
#    plt.figure(3); plt.plot(BinEdges[:-1], Hist)
#    plt.axvspan(BinEdges[sIndex], BinEdges[eIndex], color='k', alpha=0.5, lw=0,
#                label=str(Perc) + '\% of pulses with latency $<$ 3µs') 
#    
#    SetPlot(FigObj=plt, FigTitle='TTLs latencies', Plot=True)
#    SetPlot(AxesObj=plt.axes(), Axes=True)
#    plt.legend(loc='upper right')
#    
#    plt.savefig('Figs/SoundTTLLatencies-SoundBoardToOE.svg', format='svg')
    return(None)


def UnitsPSTH(FileName, StimTTLCh=-1, PSTHTimeBeforeTTL=0, 
                 PSTHTimeAfterTTL=300, StimType=['Sound'], AnalogTTLs=False, 
                 Board='OE', Override={}):
    print('Load DataInfo...')
    DirList = glob('KwikFiles/*'); DirList.sort()
    DataInfo = Hdf5F.LoadDict('/DataInfo', FileName)
    
    Units = {}
    for Stim in StimType:
        if Override != {}: 
            if 'Stim' in Override.keys():
                Stim = Override['Stim']
        
        Exps = Hdf5F.LoadExpPerStim(Stim, DirList, FileName)
        Units[Stim] = {}
        
        if isinstance(StimTTLCh, dict): UnitTTLCh = StimTTLCh[Stim]
        else: UnitTTLCh = StimTTLCh
        
        for FInd, RecFolder in enumerate(Exps):
            if AnalogTTLs: Raw, _, Files = Hdf5F.LoadOEKwik(RecFolder, AnalogTTLs)
            else: Raw, Events, _, Files = Hdf5F.LoadOEKwik(RecFolder, AnalogTTLs)
            
            OEProc = GetProc(Raw, Board)[0]
            
            Path = os.getcwd() + '/' + RecFolder
            ClusterFile = Path + '/SpkClusters.hdf5'
            Clusters = Hdf5F.LoadClusters(ClusterFile)
            
            Rate = Raw[OEProc]['info']['0']['sample_rate']
            NoOfSamplesBefore = int(round((PSTHTimeBeforeTTL*Rate)*10**-3))
            NoOfSamplesAfter = int(round((PSTHTimeAfterTTL*Rate)*10**-3))
            NoOfSamples = NoOfSamplesBefore + NoOfSamplesAfter
            XValues = (range(-NoOfSamplesBefore, 
                             NoOfSamples-NoOfSamplesBefore)/Rate)*10**3
            
            FIndS = "{0:02d}".format(FInd)
            Units[Stim][FIndS] = {}
            for Rec in Raw[OEProc]['data'].keys():
                if Override != {}: 
                    if 'Rec' in Override.keys():
                        Rec = Override['Rec']
                RecS = "{0:02d}".format(int(Rec))
                
                TTLs = QuantifyTTLsPerRec(Raw, Rec, AnalogTTLs, UnitTTLCh, 
                                          OEProc)
                
                Units[Stim][FIndS][RecS] = {}
                print('Preparing histograms and spike waveforms...')
                for CKey in Clusters[RecS].keys():
                    Classes = np.unique(Clusters[RecS][CKey]['ClusterClass'])
                    Units[Stim][FIndS][RecS][CKey] = {}
                    Units[Stim][FIndS][RecS][CKey]['PSTH'] = {}
                    
                    for Class in Classes:
                        ClassIndex = Clusters[RecS][CKey]['ClusterClass'] == Class
                        if not len(Clusters[RecS][CKey]['Spikes'][ClassIndex,:]): continue
                        
                        Class = "{0:02d}".format(int(Class))
                        Units[Stim][FIndS][RecS][CKey]['PSTH'][Class] = \
                                                         np.zeros(len(XValues))
                        
                        for TTL in range(len(TTLs)):
                            Firing = Clusters[RecS][CKey]['Timestamps'][ClassIndex] \
                                     - (TTLs[TTL]/(Rate/1000))
                            Firing = Firing[(Firing >= XValues[0]) * 
                                            (Firing < XValues[-1])]
                            SpkCount = np.histogram(Firing, 
                                                    np.hstack((XValues, 300)))[0]
                            
                            Units[Stim][FIndS][RecS][CKey]['PSTH'][Class] = \
                              Units[Stim][FIndS][RecS][CKey]['PSTH'][Class] \
                              + SpkCount
                            
                            del(Firing, SpkCount)
                    
                    print(CKey+':', str(len(Classes)), 'clusters.')
                
                del(TTLs)
                if Override != {}: 
                    if 'Rec' in Override.keys(): break
            
            del(Raw, Clusters)
    
    AnalysisFile = '../' + DataInfo['AnimalName'] + '-Analysis.hdf5'
    Hdf5F.WriteUnits(Units, AnalysisFile, XValues)
    return(None)


def UnitsPSTHInfo(FileName, StimTTLCh=-1, PSTHTimeBeforeTTL=0, 
                  PSTHTimeAfterTTL=300, StimType=['Sound'], AnalogTTLs=False, 
                  Board='OE', Override={}):
#    print('Load DataInfo...')
#    DirList = glob('KwikFiles/*'); DirList.sort()
#    DataInfo = Hdf5F.LoadDict('/DataInfo', FileName)
#    
#    Units = {}
#    for Stim in StimType:
#        if Override != {}: 
#            if 'Stim' in Override.keys():
#                Stim = Override['Stim']
#        
#        Exps = Hdf5F.LoadExpPerStim(Stim, DirList, FileName)
#        Units[Stim] = {}
#        
#        if isinstance(StimTTLCh, dict): UnitTTLCh = StimTTLCh[Stim]
#        else: UnitTTLCh = StimTTLCh
#        
#        for FInd, RecFolder in enumerate(Exps):
#            if AnalogTTLs: Raw, _, Files = Hdf5F.LoadOEKwik(RecFolder, AnalogTTLs)
#            else: Raw, Events, _, Files = Hdf5F.LoadOEKwik(RecFolder, AnalogTTLs)
#            
#            OEProc = GetProc(Raw, Board)[0]
#            
#            Path = os.getcwd() + '/' + RecFolder
#            ClusterFile = Path + '/SpkClusters.hdf5'
#            Clusters = Hdf5F.LoadClusters(ClusterFile)
#            
#            Rate = Raw[OEProc]['info']['0']['sample_rate']
#            NoOfSamplesBefore = int(round((PSTHTimeBeforeTTL*Rate)*10**-3))
#            NoOfSamplesAfter = int(round((PSTHTimeAfterTTL*Rate)*10**-3))
#            NoOfSamples = NoOfSamplesBefore + NoOfSamplesAfter
#            XValues = (range(-NoOfSamplesBefore, 
#                             NoOfSamples-NoOfSamplesBefore)/Rate)*10**3
#            
#            FIndS = "{0:02d}".format(FInd)
#            Units[Stim][FIndS] = {}
#            for Rec in range(len(Raw[OEProc]['data'])):
#                if Override != {}: 
#                    if 'Rec' in Override.keys():
#                        Rec = Override['Rec']
#                RecS = "{0:02d}".format(Rec)
#                
#                TTLs = QuantifyTTLsPerRec(Raw, Rec, AnalogTTLs, UnitTTLCh, 
#                                          OEProc)
#                
#                Units[Stim][FIndS][RecS] = {}
#                print('Preparing histograms and spike waveforms...')
#                for CKey in Clusters[RecS].keys():
#                    Classes = np.unique(Clusters[RecS][CKey]['ClusterClass'])
#                    Units[Stim][FIndS][RecS][CKey] = {}
#                    Units[Stim][FIndS][RecS][CKey]['PSTH_Info'] = [
#                                        np.zeros(len(XValues)) 
#                                        for _ in range(len(Classes))]
#                    
#                    for Cluster in range(len(Classes)):
#                        ClassIndex = Clusters[RecS][CKey]['ClusterClass'] == Cluster
#                        if not len(Clusters[RecS][CKey]['Spikes'][ClassIndex,:]): continue
#                        
#                        for TTL in range(len(TTLs)):
#                            Firing = Clusters[RecS][CKey]['Timestamps'][ClassIndex] \
#                                     - (TTLs[TTL]/(Rate/1000))
#                            Firing = Firing[(Firing >= XValues[0]) * 
#                                            (Firing < XValues[-1])]
#                            SpkCount = np.histogram(Firing, 
#                                                    np.hstack((XValues, 300)))[0]
#                            
#                            Units[Stim][FIndS][RecS][CKey]['PSTH'][Cluster] = \
#                              Units[Stim][FIndS][RecS][CKey]['PSTH'][Cluster] \
#                              + SpkCount
#                            
#                            del(Firing, SpkCount)
#                    
#                    print(CKey+':', str(len(Classes)), 'clusters.')
#                
#                del(TTLs)
#                if Override != {}: 
#                    if 'Rec' in Override.keys(): break
#            
#            del(Raw, Clusters)
#    
#    AnalysisFile = '../' + DataInfo['AnimalName'] + '-Analysis.hdf5'
#    Hdf5F.WriteUnits(Units, AnalysisFile, XValues)
    return(None)


def UnitsSpks(FileName, StimType=['Sound'], Board='OE', Override={}):
    print('Load DataInfo...')
    DirList = glob('KwikFiles/*'); DirList.sort()
    DataInfo = Hdf5F.LoadDict('/DataInfo', FileName)
    
    Units = {}
    for Stim in StimType:
        if Override != {}: 
            if 'Stim' in Override.keys():
                Stim = Override['Stim']
        
        Exps = Hdf5F.LoadExpPerStim(Stim, DirList, FileName)
        Units[Stim] = {}
        
        for FInd, RecFolder in enumerate(Exps):
            Raw = Hdf5F.LoadOEKwik(RecFolder, 'Raw')[0]
            OEProc = GetProc(Raw, Board)[0]
            
            Path = os.getcwd() + '/' + RecFolder 
            ClusterFile = Path + '/SpkClusters.hdf5'
            Clusters = Hdf5F.LoadClusters(ClusterFile)
            
            FIndS = "{0:02d}".format(FInd)
            Units[Stim][FIndS] = {}
            for Rec in Raw[OEProc]['data'].keys():
                if Override != {}: 
                    if 'Rec' in Override.keys():
                        Rec = Override['Rec']
                RecS = "{0:02d}".format(int(Rec))
                
                Units[Stim][FIndS][RecS] = {}
                for Ch in Clusters[RecS].keys():
                    Units[Stim][FIndS][RecS][Ch] = SepSpksPerCluster(Clusters[RecS][Ch], 
                                                                     Ch)
                
                if Override != {}: 
                    if 'Rec' in Override.keys(): break
            
            del(Raw, Clusters)
        
        if Override != {}: 
            if 'Stim' in Override.keys(): break
    
    AnalysisFile = '../' + DataInfo['AnimalName'] + '-Analysis.hdf5'
    Hdf5F.WriteUnits(Units, AnalysisFile)
    return(None)


def UnitsSpksPSTH_ToSVG(AnalysisFile, FileName, Override):    
    Units, XValues = Hdf5F.LoadUnits(AnalysisFile, Override)
    DataInfo = Hdf5F.LoadDict('/DataInfo', FileName)
    
#    Thrash = {}
    for SKey in Units:
        for FKey in Units[SKey]:
            for RKey in Units[SKey][FKey]:
                for Ch in Units[SKey][FKey][RKey]:
                    FigName = 'Figs/' + FileName[:-15] + '-UnitRec_' + SKey + '_Folder' + FKey + '_Rec' + RKey + '_' + Ch + '.svg'
                    FigTitle = SKey.replace('_', '-') + ' ' + Ch
                    UnitPlot = Process(target=UnitPlotPerCh, 
                                       args=(Units[SKey][FKey][RKey][Ch], Ch, XValues, DataInfo['SoundPulseDur'], FigName, FigTitle))
                    UnitPlot.start(); print('PID =', UnitPlot.pid)
                    UnitPlot.join()

