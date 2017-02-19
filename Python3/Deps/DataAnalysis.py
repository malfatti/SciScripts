#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
    Copyright (C) 2017  T. Malfatti
    
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

import Hdf5F
import h5py
import IntanBin
import numpy as np
import os
from datetime import datetime
from glob import glob
from scipy import io, signal
from subprocess import call

## Level 0
def CallWaveClus(Rate, Path):
    print('Clustering spikes...')
    
#    MLab = '/home/cerebro/Software/Programs/MatLabR2014a/bin/matlab'
    MLab = '/home/malfatti/Software/Programs/MatLabR2015a/bin/matlab'
    CmdCd = 'cd ' + Path + '; '
    CmdCluster = 'try, par.sr=' + str(int(Rate)) + '; ' + \
                 "Get_spikes('Files.txt', 'parallel', true, 'par', par);" + \
                 "Files = dir('./*spikes.mat'); Files = {Files.name}; " + \
                 "Do_clustering(Files, 'make_plots', false); end; quit"
    call([MLab, '-r', CmdCd+CmdCluster])
    
    print('Done clustering.')
    return(None)


def CheckStimulationPattern(TTL):
    if np.mean(TTL) > 1000: 
        print('Sinusoidal stimulation')
        Threshold = (max(TTL) - min(TTL)) / 2
        return(Threshold)
    else:
        print('square pulses stimulation')
        Threshold = ((max(TTL) - min(TTL)) / 2) + 2*(np.std(TTL))
        return(Threshold) 


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


def FixTTLs(Array, TTLsToFix):
    for TTL in TTLsToFix:
        nInd = np.random.randint(1, 100)
        while nInd == TTL: nInd = np.random.randint(0, 100)
        while nInd >= len(Array): nInd = np.random.randint(0, 100)
        
        print('TTL', str(TTL), 'was replaced by', str(nInd))
        Array[TTL] = Array[nInd]
    
    return(Array)


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


def RemapChannels(Tip, Head, Connector):
    """
    Get probe channels order. It doesn't matter what logic you follow to order 
    your connector channels, but you MUST follow the same logic for your probe 
    head.
    
    If the probe tip channels are put top-down or bottom-up, the resulting 
    channel map will be ordered accordingly.
    
    Example:
        CustomAdaptor = [5, 6, 7, 8, 9, 10 ,11, 12, 13, 14, 15, 16, 1, 2, 3, 4]
        RHAHeadstage = [16, 15, 14, 13, 12, 11, 10, 9, 1, 2, 3, 4, 5, 6, 7, 8]
        A16 = {'ProbeTip': [9, 8, 10, 7, 13, 4, 12, 5, 15, 2, 16, 1, 14, 3, 11, 6],
               'ProbeHead': [8, 7, 6, 5, 4, 3, 2, 1, 9, 10, 11, 12, 13, 14, 15, 16]}
        
        ChannelMap = RemapChannels(A16['ProbeTip'], A16['ProbeHead'], CustomAdaptor)
    """
    print('Get channel order... ', end='')
    ChNo = len(Tip)
    ChMap = [0]*ChNo
    
    for Ch in range(ChNo):
        TipCh = Tip[Ch] # What channel should be the Ch
        HeadCh = Head.index(TipCh) # Where Ch is in Head
        ChMap[Ch] = Connector[HeadCh] # Channels in depth order
    
    print('Done.')
    return(ChMap)


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


def SetPlot(Backend='Qt5Agg', AxesObj=(), FigObj=(), FigTitle='', Params=False, 
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
                  'svg.fonttype': 'none'}
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


def WriteUnits(Units, FileName, XValues=[]):
    """ Belongs to Hdf5F """
    print('Writing data to', FileName+'... ', end='')
    Group = 'UnitRec'
    Thrash = []
    with h5py.File(FileName) as F:
        for RKey in Units.keys():
            for Ch in Units[RKey].keys():
                Path = '/'+Group+'/'+RKey+'/'+Ch
                if Path not in F: F.create_group(Path)
                
                for Var in Units[RKey][Ch].keys():
                    if Var == 'Spks':
                        if '/'+Path+'/Spks' in F: del(F[Path]['Spks'])
                        if Path+'/Spks' not in F: 
                            F.create_group(Path+'/Spks')
                        
                        for Class in Units[RKey][Ch]['Spks'].keys():
                            for Key, Spk in enumerate(Units[RKey][Ch]['Spks'][Class]):
                                if not len(Spk):
                                    Thrash.append([Path, Class, Key])
                                    del(Units[RKey][Ch]['Spks'][Class][Key])
                                    continue
                                if Path+'/Spks/'+Class not in F:
                                    F[Path]['Spks'].create_group(Class)
                                F[Path]['Spks'][Class][str(Key)] = Spk
                    
                    elif Var == 'PSTH':
                        if '/'+Path+'/PSTH' in F: del(F[Path]['PSTH'])
                        F[Path].create_group('PSTH')
                        for VarKey in Units[RKey][Ch][Var].keys():
                            F[Path]['PSTH'][VarKey] = Units[RKey][Ch]['PSTH'][VarKey][:]
                    
                    elif Var == 'Spks_Info':
                        for VarKey, VarValue in Units[RKey][Ch][Var].items():
                            F[Path]['Spks'].attrs[VarKey] = VarValue
                    
                    elif Var == 'PSTH_Info':
                        for VarKey, VarValue in Units[RKey][Ch][Var].items():
                            F[Path]['PSTH'].attrs[VarKey] = VarValue
                    
                    else:
                        print(Var, 'not supported')
        
        if len(XValues): F[Group]['XValues'] = XValues
    
    if Thrash: 
        for _ in Thrash: print(_)
    
    return(None)


## Level 1
def ClusterizeSpks(Data, Rate, Path, AnalogTTLs=False, ChannelMap=[], Override={}):
    """ Detect and clusterize spks using WaveClus """
    
    os.makedirs(Path, exist_ok=True)
    
    Clusters = {}
    for Rec in Data.keys():
        if not ChannelMap: ChannelMap = [_ for _ in range(Data[Rec].shape[1])]
        
        if Override != {}: 
            if 'Rec' in Override.keys():
                Rec = Override['Rec']
        
        RecS = "{0:02d}".format(int(Rec))
        
        print('Separating channels according to ChannelMap...')
        Data = [Data[Rec][:, _] for _ in ChannelMap]
                        
        print('Writing files for clustering... ', end='')
        FileList = []
        for Ind, Ch in enumerate(Data):
            MatName = 'Rec' + RecS + '-Ch' + "{0:02d}".format(Ind+1) + '.mat'
            
            FileList.append(MatName)
            io.savemat(Path+MatName, {'data': Ch})
        
        TxtFile = open(Path+'Files.txt', 'w')
        for File in FileList: TxtFile.write(File + '\n')
        TxtFile.close()
        print('Done.')
        
        CallWaveClus(Rate, Path)
        
        ClusterList = glob(Path + 'times_*'); ClusterList.sort()
        ClusterList = [_ for _ in ClusterList if _[-11:-9] == RecS]
        
        Clusters[RecS] = {}
        for File in ClusterList:
            Ch = File[-8:-4]
            
            ClusterFile = io.loadmat(File)
            Clusters[RecS][Ch] = {}
            Clusters[RecS][Ch]['ClusterClass'] = ClusterFile['cluster_class'][:, 0]
            Clusters[RecS][Ch]['Timestamps'] = ClusterFile['cluster_class'][:, 1]
            Clusters[RecS][Ch]['Spikes'] = ClusterFile['spikes'][:]
        
        if Override != {}: 
            if 'Rec' in Override.keys(): break
    
        ToDelete = glob(Path + '*')
        for File in ToDelete: os.remove(File)
    
    os.removedirs(Path)
    Hdf5F.WriteClusters(Clusters, Path[:-6]+'SpkClusters.hdf5')


def QuantifyTTLsPerRec(Data, Rec, AnalogTTLs, TTLCh=0):
    print('Get TTL timestamps... ', end='')
    if AnalogTTLs:
        TTL = Data[Rec][:, TTLCh-1]
        Threshold = CheckStimulationPattern(TTL)
        TTLs = []
        for _ in range(1, len(TTL)):
            if TTL[_] > Threshold:
                if TTL[_-1] < Threshold: TTLs.append(_)
        
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
#        RawTime = [_*Rate for _ in Raw['timestamps'][int(Rec)]]
#        TTLs = TTLsPerRec[Rec]
        
        print('Done.')
#        return(RawTime, TTLs)


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


def UnitsSpks(Data, Path, AnalysisFile, Override={}):
    ClusterFile = Path + '/SpkClusters.hdf5'
    Clusters = Hdf5F.LoadClusters(ClusterFile)
    
    Units = {}
    for Rec in Data.keys():
        if Override != {}: 
            if 'Rec' in Override.keys():
                Rec = Override['Rec']
        
        RecS = "{0:02d}".format(int(Rec))
        
        Units[RecS] = {}
        for Ch in Clusters[RecS].keys():
            Units[RecS][Ch] = SepSpksPerCluster(Clusters[RecS][Ch], 
                                                             Ch)
        
        if Override != {}: 
            if 'Rec' in Override.keys(): break
    
    WriteUnits(Units, AnalysisFile)
    del(Data, Clusters)
    
    return(None)


## Level 2
def UnitsPSTH(Data, Rate, Path, AnalysisFile, TTLCh=0, PSTHTimeBeforeTTL=0, 
                 PSTHTimeAfterTTL=300, AnalogTTLs=False, 
                 Override={}):
    Units = {}
    
    ClusterFile = Path + '/SpkClusters.hdf5'
    Clusters = Hdf5F.LoadClusters(ClusterFile, KeyInd=-1)
    
    for Rec in Data.keys():
        if Override != {}: 
            if 'Rec' in Override.keys():
                Rec = Override['Rec']
        
        TTL = Data[Rec][:, TTLCh-1]
        TTLs = QuantifyTTLsPerRec(Data, Rec, AnalogTTLs, TTLCh)
        StimFreq = 1 / ((TTLs[3] - TTLs[2])/Rate)
        
        if np.mean(TTL) > 1000: 
            URecS = "{0:02d}".format(int(Rec)) + '-Sinusoidal-' + str(int(StimFreq)) + 'Hz'
            NoOfSamplesBefore = int(round((0*Rate)*10**-3))
            NoOfSamplesAfter = int(round((TTLs[3] - TTLs[2])))
        else:
            URecS = "{0:02d}".format(int(Rec)) + '-Squares-' + str(int(StimFreq)) + 'Hz'
            NoOfSamplesBefore = int(round((PSTHTimeBeforeTTL*Rate)*10**-3))
            NoOfSamplesAfter = int(round((PSTHTimeAfterTTL*Rate)*10**-3))
        
        NoOfSamples = NoOfSamplesBefore + NoOfSamplesAfter
        XValues = (range(-NoOfSamplesBefore, 
                         NoOfSamples-NoOfSamplesBefore)/Rate)*10**3
        
        RecS = "{0:02d}".format(int(Rec))
        Units[URecS] = {}
        print('Preparing histograms and spike waveforms...')
        for CKey in Clusters[RecS].keys():
            Classes = np.unique(Clusters[RecS][CKey]['ClusterClass'])
#            Units[URecS][CKey] = {}
            Units[URecS][CKey] = SepSpksPerCluster(Clusters[RecS][CKey], CKey)
            Units[URecS][CKey]['PSTH'] = {}
            
            for Class in Classes:
                ClassIndex = Clusters[RecS][CKey]['ClusterClass'] == Class
                if not len(Clusters[RecS][CKey]['Spikes'][ClassIndex,:]): continue
                
                Class = "{0:02d}".format(int(Class))
                Units[URecS][CKey]['PSTH'][Class] = np.zeros(len(XValues))
                
                for TTL in range(len(TTLs)):
                    Firing = Clusters[RecS][CKey]['Timestamps'][ClassIndex] \
                             - (TTLs[TTL]/(Rate/1000))
                    Firing = Firing[(Firing >= XValues[0]) * 
                                    (Firing < XValues[-1])]
                    SpkCount = np.histogram(Firing, 
                                            np.hstack((XValues, len(XValues)/(Rate/1000))))[0]
                    
                    Units[URecS][CKey]['PSTH'][Class] = Units[URecS][CKey]['PSTH'][Class] + SpkCount
                    
                    del(Firing, SpkCount)
            
            print(CKey+':', str(len(Classes)), 'clusters.')
        
        del(TTLs)
        if Override != {}: 
            if 'Rec' in Override.keys(): break
    
    WriteUnits(Units, AnalysisFile, XValues)
    return(None)


## Level 3
def ABRAnalysis(FileName, ABRCh=[1], ABRTTLCh=1, ABRTimeBeforeTTL=0, 
                ABRTimeAfterTTL=12, FilterFreq=[300, 3000], FilterOrder=4, 
                StimType='Sound', AnalogTTLs=False, Board='OE', Type='kwik',
                Override={}):
    """
    Analyze ABRs from data. A '*ABRs.hdf5' file will be saved 
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
    
    """
    
    print('Load DataInfo...')
    DirList = glob('KwikFiles/*'); DirList.sort()
    DataInfo = Hdf5F.LoadDict('/DataInfo', FileName)
    
    AnalysisFile = './' + DataInfo['AnimalName'] + '-Analysis.hdf5'
    Now = datetime.now().strftime("%Y%m%d%H%M%S")
    Here = os.getcwd().split(sep='/')[-1]
    Group = Here + '-ABRs_' + Now
        
    for Stim in StimType:
        if Override != {}: 
            if 'Stim' in Override.keys(): Stim = Override['Stim']
        
        Exps = Hdf5F.LoadExpPerStim(Stim, DirList, FileName)
        
        for RecFolder in Exps:
            ABRs, Info = ABRCalc()
            Hdf5F.WriteABR(ABRs, Info['XValues'], Group, Info['Path'], AnalysisFile)
        
        
def ABRCalc(FileName, ABRCh=[1], ABRTTLCh=1, ABRTimeBeforeTTL=0, 
                ABRTimeAfterTTL=12, FilterFreq=[300, 3000], FilterOrder=4, 
                StimType='Sound', AnalogTTLs=False, Board='OE', Type='kwik',
                Override={}):
    """
    Analyze ABRs from data. A '*ABRs.hdf5' file will be saved 
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
    
    """
    
    print('Load DataInfo...')
    DirList = glob('KwikFiles/*'); DirList.sort()
    DataInfo = Hdf5F.LoadDict('/DataInfo', FileName)
    
    AnalysisFile = './' + DataInfo['AnimalName'] + '-Analysis.hdf5'
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
            
            Info['Frequency'] = ExpInfo['Hz']
            
            Info['XValues'] = (range(-NoOfSamplesBefore, 
                                     NoOfSamples-NoOfSamplesBefore)/Rate)*10**3
            
            for Rec in Raw[OEProc]['data'].keys():
                print('Slicing and filtering ABRs Rec ', str(Rec), '...')
                
                if len(Raw[OEProc]['data'][Rec]) < 50*Rate:
                    print('Rec', Rec, 'is broken!!!')
                    continue
                
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
            
            Info['Path'] = Stim+'/'+ExpInfo['DVCoord']+'/'+Info['Frequency']
#            AnalysisFile = '../' + DataInfo['AnimalName'] + '-Analysis.hdf5'
            Hdf5F.WriteABR(ABRs, Info['XValues'], Group, Info['Path'], AnalysisFile)
    
    return(None)


#%% Tests
File = 'Intan/Arch45S_153244.int'; ClusterPath = './SepCh/'
Raw = {}
Raw['0'] = IntanBin.IntLoad(File)[0]

RHAHeadstage = [16, 15, 14, 13, 12, 11, 10, 9, 1, 2, 3, 4, 5, 6, 7, 8]
A16 = {'ProbeTip': [9, 8, 10, 7, 13, 4, 12, 5, 15, 2, 16, 1, 14, 3, 11, 6],
       'ProbeHead': [8, 7, 6, 5, 4, 3, 2, 1, 9, 10, 11, 12, 13, 14, 15, 16]}
ChannelMap = RemapChannels(A16['ProbeTip'], A16['ProbeHead'], RHAHeadstage)
