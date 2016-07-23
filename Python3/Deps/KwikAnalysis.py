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

##import datetime
##import h5py
##import Kwik
#from numbers import Number

import Hdf5F
import numpy as np
import os
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
    TTLsPerRec = {_Rec: [EventSample[_] for _ in range(len(EventRec)) 
                         if EventRec[_] == _Rec 
                         and EventCh[_] == TTLCh-1 
                         and EventID[_] == 1]
                  for _Rec in TTLRecs}
#    TTLRising = Kwik.get_rising_edge_times(Files['kwe'], TTLCh-1)
    
    return(TTLsPerRec)


def QuantifyTTLsPerRec(Raw, Rec, AnalogTTLs, ChTTL=-1, Proc='', TTLsPerRec=[], 
                       Rate=[]):
    print('Get TTL timestamps... ', end='')
    if AnalogTTLs:
        TTLCh = Raw[Proc]['data'][str(Rec)][:, ChTTL-1]
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
        RawTime = [_*Rate for _ in Raw['timestamps'][Rec]]
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
        
        Array[TTL] = Data[Proc]['data'][str(Rec)][Start:End, DataCh[0]-1] \
                     * Data[Proc]['channel_bit_volts'][str(Rec)][DataCh[0]-1] \
                     * 1000 # in mV
        
        if len(Array[TTL]) != End-Start: TTLsToFix.append(TTL)
            
    Array = FixTTLs(Array, TTLsToFix)
    
    print('Done.')
    return(Array)


def UnitPlotPerCh(ChDict, Ch, XValues, FigName, FigTitle):
    ClusterNo = len(ChDict['Spks'])
    if ClusterNo == 0: print(Ch, 'had no spikes :('); return(None)
    
    Params = SetPlot(Backend='Agg', Params=True)
    from matplotlib import rcParams; rcParams.update(Params)
    from matplotlib import pyplot as plt
    
    Fig, Axes = plt.subplots(ClusterNo,2, figsize=(6, 3*ClusterNo))

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
            Axes[1].bar(XValues, ChDict['PSTH'][Class])
            SetPlot(AxesObj=Axes[0], Axes=True)
            SetPlot(AxesObj=Axes[1], Axes=True)
        else:
#                    Axes[Cluster][0].set_title('Peak='+str(PSTHPeak)+' Mean='+\
#                                               str(PSTHMean)+' Std='+str(PSTHStd))
            Axes[int(Class)-1][0].plot(np.mean(ChDict['Spks'][Class], axis=0), 'k')
            Axes[int(Class)-1][1].bar(XValues, ChDict['PSTH'][Class])
            SetPlot(AxesObj=Axes[int(Class)-1][0], Axes=True)
            SetPlot(AxesObj=Axes[int(Class)-1][1], Axes=True)
    
    SetPlot(FigObj=Fig, FigTitle=FigTitle, Plot=True)
    print('Writing to', FigName+'... ', end='')
    Fig.savefig(FigName, format='svg')
    print('Done.')
    return(None)


## Higher-level functions
def ABRAnalysis(FileName, ABRCh=[1], ABRTTLCh=1, ABRTimeBeforeTTL=0, 
                ABRTimeAfterTTL=12, FilterFreq=[300, 3000], FilterOrder=4, 
                StimType='Sound', AnalogTTLs=False, Board='OE'):
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
    DataInfo['SoundAmpF'] = Hdf5F.LoadDict('/DataInfo/SoundAmpF', FileName)
    Exps = Hdf5F.LoadExpPerStim('Sound', DirList, FileName)
    
    
    print('Preallocate memory...')
    ABRs = [[0] for _ in range(len(DataInfo['NoiseFrequency']))]
    for Freq in range(len(DataInfo['NoiseFrequency'])):
        Key = str(DataInfo['NoiseFrequency'][Freq][0]) + '-' \
              + str(DataInfo['NoiseFrequency'][Freq][1])
        ABRs[Freq] = [{} for _ in range(len(DataInfo['SoundAmpF'][Key]))]
    del(Freq)
    
    for RecFolder in Exps:
        ExpInfo = Hdf5F.ExpExpInfo(RecFolder, DirList, FileName)
        
        if AnalogTTLs: Raw, _, Files = Hdf5F.LoadOEKwik(RecFolder, AnalogTTLs)
        else: Raw, Events, _, Files = Hdf5F.LoadOEKwik(RecFolder, AnalogTTLs)
        
        OEProc, RHAProc, ABRProc = GetProc(Raw, Board)
        
        if AnalogTTLs: Raw = GetRecKeys(Raw, [0], AnalogTTLs)
        else:
            Raw, EventRec = GetRecKeys(Raw, Events, AnalogTTLs)
            TTLsPerRec = GetTTLInfo(Events, EventRec, ABRTTLCh)
        
        Rate = Raw[OEProc]['info']['0']['sample_rate']
        NoOfSamplesBefore = ABRTimeBeforeTTL*int(Rate*10**-3)
        NoOfSamplesAfter = ABRTimeAfterTTL*int(Rate*10**-3)
        NoOfSamples = NoOfSamplesBefore + NoOfSamplesAfter

        XValues = list((range((NoOfSamplesBefore)*-1, 0)/Rate)*10**3) + \
                  list((range(NoOfSamplesAfter)/Rate)*10**3)
        
        for Rec in range(len(Raw[OEProc]['data'])):
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
            
            if ExpInfo['DVCoord'] not in ABRs[ExpInfo['Hz']][Rec]:
                ABRs[ExpInfo['Hz']][Rec][ExpInfo['DVCoord']] = [ABR]
            else:
                ABRs[ExpInfo['Hz']][Rec][ExpInfo['DVCoord']].append(ABR)
    
    AnalysisFile = '../' + DataInfo['AnimalName'] + '-Analysis.hdf5'
    Hdf5F.WriteABRs(ABRs, XValues, AnalysisFile)
    return(None)

    
def ABRPlot(FileName):
    """ 
    This function will plot the data from *Stim.hdf5. Make sure FileName is a 
    string with the path to only one file.
    
    Also, LaTeX will render all the text in the plots. For this, make sure you 
    have a working LaTex installation and dvipng package installed.
    """
    
    os.makedirs('Figs', exist_ok=True)    # Figs folder
    
    print('Loading data...')
    DataInfo = Hdf5F.LoadDict('/DataInfo', FileName)
    DataInfo['SoundAmpF'] = Hdf5F.LoadDict('/DataInfo/SoundAmpF', FileName)
    ABRs, XValues = Hdf5F.LoadABRs(FileName)
    
    Params = SetPlot(Params=True)
    from matplotlib import rcParams; rcParams.update(Params)
    from matplotlib import pyplot as plt
    
    
    print('Plotting...')
    Colormaps = [plt.get_cmap('Reds'), plt.get_cmap('Blues')]
    Colors = [[Colormaps[0](255-(_*20)), Colormaps[1](255-(_*20))] 
              for _ in range(len(ABRs[0]))]
    
    Keys = list(ABRs[0][0].keys())
    for Key in Keys:
        for Trial in range(len(ABRs[0][0][Key])):
            for Freq in range(len(ABRs)):
                Fig, Axes = plt.subplots(len(ABRs[Freq]), 
                                         sharex=True, 
                                         figsize=(8, 
                                            1.5*len(DataInfo['Intensities'])))
                
                KeyHz = str(DataInfo['NoiseFrequency'][Freq][0]) + '-' \
                        + str(DataInfo['NoiseFrequency'][Freq][1])
                
                if 0.0 in DataInfo['SoundAmpF'][KeyHz]:
                    DataInfo['SoundAmpF'][KeyHz][
                        DataInfo['SoundAmpF'][KeyHz].index(0.0)
                                                ] = 0
#                upYLim = max(ABRs[0][Freq][0][Key][Trial])
#                downYLim = min(ABRs[0][Freq][0][Key][Trial])
                upYLim = [0] * len(ABRs[Freq][0][Key]); downYLim = upYLim[:]
        
                for Trial in range(len(ABRs[Freq][0][Key])):
                    upYLim[Trial] = max(ABRs[Freq][0][Key][Trial])
                    downYLim[Trial] = min(ABRs[Freq][0][Key][Trial])
                
                upYLim = max(upYLim); downYLim = min(downYLim)
                
                for AmpF in range(len(DataInfo['Intensities'])):
                    FigTitle = KeyHz + ' Hz, DV ' + Key + \
                               ', trial ' + str(Trial+1)
                    YLabel = 'Voltage [mV]'; XLabel = 'Time [ms]'
                    LineLabel = str(DataInfo['Intensities'][AmpF]) + ' dB'
                    SpanLabel = 'Sound pulse'
                    
                    Ind1 = list(XValues).index(0)
                    Ind2 = list(XValues).index(3)
                    
                    Axes[AmpF].plot(XValues, ABRs[Freq][AmpF][Key][Trial], 
                                    color=Colors[AmpF][0], label=LineLabel)
                    
                    Axes[AmpF].axvspan(XValues[Ind1], XValues[Ind2], 
                                       color='k', alpha=0.3, lw=0, 
                                       label=SpanLabel)
                    
                    SetPlot(AxesObj=Axes[AmpF], Axes=True)
                    Axes[AmpF].legend(loc='lower right', frameon=False)
                    Axes[AmpF].spines['bottom'].set_visible(False)
                    Axes[AmpF].spines['left'].set_bounds(round(0), round(1))
                    Axes[AmpF].xaxis.set_ticks_position('none')
                    Axes[AmpF].set_ylabel(YLabel)
                    Axes[AmpF].set_ylim(downYLim, upYLim)
                    
                Axes[-1].spines['bottom'].set_visible(True)
                Axes[-1].set_xlabel(XLabel)
                Axes[-1].spines['bottom'].set_bounds(round(0), round(1))
                SetPlot(FigObj=Fig, FigTitle=FigTitle, Plot=True)
                
                FigName = 'Figs/' + FileName[:-15] + '-ABR_DV' +  Key + \
                          '_Trial' + str(Trial) + '_Freq' + KeyHz + '.svg'
                Fig.savefig(FigName, format='svg')


def GPIASAnalysis(RecFolder, FileName, GPIASCh=1, GPIASTTLCh=1, 
                  GPIASTimeBeforeTTL=50, GPIASTimeAfterTTL=150, 
                  FilterFreq=[70, 400], FilterOrder=4, AnalogTTLs=False, 
                  Board='OE'):
    
    print('set paths...')
    DirList = glob('KwikFiles/*'); DirList.sort()
    RecFolder = DirList[RecFolder-1]
    
    DataInfo = Hdf5F.LoadDict('/DataInfo', FileName)
    DataInfo['SoundBackgroundAmpF'] = Hdf5F.LoadDict(
                                         '/DataInfo/SoundBackgroundAmpF', 
                                         FileName, Attrs=False)
    DataInfo['SoundPulseAmpF'] = Hdf5F.LoadDict('/DataInfo/SoundPulseAmpF', 
                                                FileName, Attrs=False)
    
    for Path in ['Freqs', 'FreqOrder', 'FreqSlot']:
        DataInfo[Path] = Hdf5F.LoadDataset('/DataInfo/'+Path, FileName)
    
    print('Preallocate memory...')
    GPIAS = [[0] for _ in range(len(DataInfo['NoiseFrequency']))]
    for Freq in range(len(DataInfo['NoiseFrequency'])):
        GPIAS[Freq] = [[0] for _ in range(round(DataInfo['NoOfTrials']*2))]
    
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
    XValues = (range(-NoOfSamplesBefore, 
                     NoOfSamples-NoOfSamplesBefore)/Rate)*10**3
    
    for Rec in range(len(Raw[OEProc]['data'])):
        print('Slicing and filtering Rec ', str(Rec), '...')
        Freq = DataInfo['FreqOrder'][Rec][0]; 
        Trial = DataInfo['FreqOrder'][Rec][1];
        
        if AnalogTTLs:
            TTLs = QuantifyTTLsPerRec(Raw, Rec, AnalogTTLs, GPIASTTLCh, OEProc)
            GPIAS[Freq][Trial] = SliceData(Raw, OEProc, Rec, TTLs, GPIASCh, 
                                           NoOfSamplesBefore, NoOfSamplesAfter, 
                                           NoOfSamples, AnalogTTLs)
        else:
            RawTime, TTLs = QuantifyTTLsPerRec(Raw, Rec, AnalogTTLs, 
                                               TTLsPerRec=TTLsPerRec)
            GPIAS[Freq][Trial] = SliceData(Raw, OEProc, Rec, TTLs, GPIASCh, 
                                           NoOfSamplesBefore, NoOfSamplesAfter, 
                                           NoOfSamples, AnalogTTLs, RawTime)
        
        GPIAS[Freq][Trial] = GPIAS[Freq][Trial][0][:]
    
    for Freq in range(len(DataInfo['NoiseFrequency'])):
        # Separate Gap/NoGap
        gData = GPIAS[Freq][:]
        NoGapAll = [gData[_] for _ in range(len(gData)) if _%2 == 0]
        GapAll = [gData[_] for _ in range(len(gData)) if _%2 != 0]
        
        # Average
        gData = [0, 0]
        NoGapSum = list(map(sum, zip(*NoGapAll)))
        GapSum = list(map(sum, zip(*GapAll)))
        gData[0] = [_/DataInfo['NoOfTrials'] for _ in NoGapSum]
        gData[1] = [_/DataInfo['NoOfTrials'] for _ in GapSum]
        
        # Hilbert
        gData[0] = abs(signal.hilbert(gData[0]))
        gData[1] = abs(signal.hilbert(gData[1]))
        
        # Bandpass filter
        gData[0] = FilterSignal(gData[0], Rate, FilterFreq, FilterOrder, 
                                                            'bandpass')
        gData[1] = FilterSignal(gData[1], Rate, FilterFreq, FilterOrder, 
                                                            'bandpass')
        
        # SavGol smooth
        gData[0] = signal.savgol_filter(gData[0], 5, 2, mode='nearest')
        gData[1] = signal.savgol_filter(gData[1], 5, 2, mode='nearest')
        
        GPIAS[Freq] = {}
        GPIAS[Freq]['Gap'] = gData[1][:]; GPIAS[Freq]['NoGap'] = gData[0][:]
        
        del(NoGapAll, GapAll, NoGapSum, GapSum, gData)
    
    AnalysisFile = '../' + DataInfo['AnimalName'] + '-Analysis.hdf5'
    Hdf5F.WriteGPIAS(GPIAS, RecFolder, XValues, AnalysisFile)
    
    return(None)


def GPIASPlot(FileName):
    os.makedirs('Figs', exist_ok=True)    # Figs folder
    
    print('Loading data...')
    ## DataInfo
    DataInfo = Hdf5F.LoadDict('/DataInfo', FileName)
    DataInfo['SoundBackgroundAmpF'] = Hdf5F.LoadDict(
                                         '/DataInfo/SoundBackgroundAmpF', 
                                         FileName, Attrs=False)
    DataInfo['SoundPulseAmpF'] = Hdf5F.LoadDict('/DataInfo/SoundPulseAmpF', 
                                                FileName, Attrs=False)
    
    for Path in ['Freqs', 'FreqOrder', 'FreqSlot']:
        DataInfo[Path] = Hdf5F.LoadDataset('/DataInfo/'+Path, FileName)
    
    ## GPIAS
    GPIAS, XValues = Hdf5F.LoadGPIAS(FileName)
    
    Params = SetPlot(Params=True)
    from matplotlib import rcParams; rcParams.update(Params)
    from matplotlib import pyplot as plt
    
    print('Plotting...')
    Ind1 = list(XValues).index(0)
    Ind2 = list(XValues).index(int(DataInfo['SoundLoudPulseDur']*1000))
    
    for Freq in range(len(DataInfo['NoiseFrequency'])):
        FreqKey = str(DataInfo['NoiseFrequency'][Freq][0]) + '-' + \
                  str(DataInfo['NoiseFrequency'][Freq][1])
        
        FigTitle = FreqKey + ' Hz'
        LineNoGapLabel = 'No Gap'; LineGapLabel = 'Gap'
        SpanLabel = 'Sound Pulse'
        XLabel = 'time [ms]'; YLabel = 'voltage [mV]'
        
        plt.figure(Freq)
        plt.plot(XValues, GPIAS[Freq]['NoGap'], 
                 color='r', label=LineNoGapLabel, lw=2)
        plt.plot(XValues, GPIAS[Freq]['Gap'], 
                 color='b', label=LineGapLabel, lw=2)
        plt.axvspan(XValues[Ind1], XValues[Ind2], color='k', alpha=0.5, 
                    lw=0, label=SpanLabel)

        SetPlot(AxesObj=plt.axes(), Axes=True)
        SetPlot(FigObj=plt, FigTitle=FigTitle, Plot=True)
        plt.ylabel(YLabel); plt.xlabel(XLabel)
        plt.legend(loc='lower right')
        
        FigName = 'Figs/' + FileName[:-5] + '-' + FreqKey + '.svg'
        plt.savefig(FigName, format='svg')
        
    print('Done.')


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
    
    for Rec in range(len(Raw[OEProc]['data'])):
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
        
        TTLsLatency[str(Rec)] = dict((Name, eval(Name)) 
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


def ClusterizeAll(FileName, StimType=['Sound'], AnalogTTLs=False, Board='OE', Override={}):
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
            for Rec in range(len(Raw[OEProc]['data'])):
                if Override != {}: 
                    if 'Rec' in Override.keys():
                        Rec = Override['Rec']
                RecS = "{0:02d}".format(Rec)
                
                print('Separating channels according to ChannelMap...')
                Data = [Raw[OEProc]['data'][str(Rec)][:, _-1] * 
                        Raw[OEProc]['channel_bit_volts'][str(Rec)][_-1] * 1000
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
            for Rec in range(len(Raw[OEProc]['data'])):
                if Override != {}: 
                    if 'Rec' in Override.keys():
                        Rec = Override['Rec']
                RecS = "{0:02d}".format(Rec)
                
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
            for Rec in range(len(Raw[OEProc]['data'])):
                if Override != {}: 
                    if 'Rec' in Override.keys():
                        Rec = Override['Rec']
                RecS = "{0:02d}".format(Rec)
                
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
    
#    Thrash = {}
    for SKey in Units:
        for FKey in Units[SKey]:
            for RKey in Units[SKey][FKey]:
                for Ch in Units[SKey][FKey][RKey]:
                    FigName = 'Figs/' + FileName[:-15] + '-UnitRec_' + SKey + '_Folder' + FKey + '_Rec' + RKey + '_' + Ch + '.svg'
                    FigTitle = SKey.replace('_', '-') + ' ' + Ch
                    UnitPlot = Process(target=UnitPlotPerCh, 
                                       args=(Units[SKey][FKey][RKey][Ch], Ch, XValues, FigName, FigTitle))
                    UnitPlot.start(); print('PID =', UnitPlot.pid)
                    UnitPlot.join()

