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
#import IntanBin
import numpy as np
import os
from datetime import datetime
from glob import glob
from multiprocessing import Process
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
    call([MLab, '-nosplash', '-r', CmdCd+CmdCluster])
    
    print('Done clustering.')
    return(None)


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


def GenerateTTLVector(TTLs, TTLLen, FullLen):
    TTLVec = np.zeros([FullLen, 1])
    
    for TTL in TTLs:
        TTLVec[TTL:TTL+TTLLen] = np.ones([TTLLen, 1])
    
    return(TTLVec)

def GenerateFakeTTLsRising(Start, Dur, No, PauseBetween):
    Dur = Dur + PauseBetween
    
    FakeTTLs = []
    for Block in Start:
        FakeTTLs = FakeTTLs + list(range(Block, (Dur*No) + Block, Dur))
    
    return(FakeTTLs)


def GetRecXValues(TTLs, Rate, TimeBeforeTTL, TimeAfterTTL):
#    StimFreq = 1 / ((TTLs[3] - TTLs[2])/Rate)
#    
#    if np.mean(TTL) > 1000: 
#        URecS = Rec + '-Sinusoidal-' + str(int(StimFreq)) + 'Hz'
#        NoOfSamplesBefore = int(round((0*Rate)*10**-3))
#        NoOfSamplesAfter = int(round((TTLs[3] - TTLs[2])))
#    else:
#        URecS = Rec + '-Squares-' + str(int(StimFreq)) + 'Hz'
#        NoOfSamplesBefore = int(round((TimeBeforeTTL*Rate)*10**-3))
#        NoOfSamplesAfter = int(round((TimeAfterTTL*Rate)*10**-3))
    NoOfSamplesBefore = int((TimeBeforeTTL*Rate) * 10**-3)
    NoOfSamplesAfter = int((TimeAfterTTL*Rate) * 10**-3)
    NoOfSamples = NoOfSamplesBefore + NoOfSamplesAfter
    XValues = (range(-NoOfSamplesBefore, 
                     NoOfSamples-NoOfSamplesBefore)/Rate) * 10**3
    
    return(XValues)


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


def GetTTLThreshold(TTLCh):
    if np.mean(TTLCh) > 1000: 
        print('Sinusoidal stimulation')
        Threshold = (max(TTLCh) - min(TTLCh)) / 2
        return(Threshold)
    else:
        print('square pulses stimulation')
        Top = max(TTLCh) + abs(min(TTLCh))
        Bot = min(TTLCh) + abs(min(TTLCh))
        Threshold = max(TTLCh) - (Top - Bot)/3 # + 2*(np.std(TTL))
        return(Threshold) 


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


## Level 1
def ClusterizeSpks(Data, Rate, ChannelMap, ClusterPath, AnalysisFile, 
                   AnalysisKey, Rec='0', Override={}, Return=False):
    """ Detect and clusterize spks using WaveClus """
    
    os.makedirs(ClusterPath, exist_ok=True)
    
    Data = [Data[:, _-1] for _ in sorted(ChannelMap)]
    print('Writing files for clustering... ', end='')
    FileList = []
    
    try: Rec = "{0:02d}".format(int(Rec))
    except ValueError: pass
    
    for Ind, Ch in enumerate(Data):
        MatName = 'Rec' + Rec + '-Ch' + "{0:02d}".format(Ind+1) + '.mat'
        
        FileList.append(MatName)
        io.savemat(ClusterPath + '/' + MatName, {'data': Ch})
    
    TxtFile = open(ClusterPath + '/Files.txt', 'w')
    for File in FileList: TxtFile.write(File + '\n')
    TxtFile.close()
    print('Done.')
    
    CallWaveClus(Rate, ClusterPath)
    
    ClusterList = glob(ClusterPath + '/times_*'); ClusterList.sort()
    ClusterList = [_ for _ in ClusterList if _.split('-')[-2].split('_')[-1] == 'Rec'+Rec]
    print(ClusterList)
    Clusters = {Rec: {}}
    for File in ClusterList:
        Ch = File[-8:-4]
        
        ClusterFile = io.loadmat(File)
        Clusters[Rec][Ch] = {}
        Clusters[Rec][Ch]['ClusterClass'] = ClusterFile['cluster_class'][:, 0]
        Clusters[Rec][Ch]['Timestamps'] = ClusterFile['cluster_class'][:, 1]
        Clusters[Rec][Ch]['Spikes'] = ClusterFile['spikes'][:]
    
    Group = AnalysisKey + '/SpkClusters'
    Hdf5F.WriteClusters(Clusters, AnalysisFile, Group)
    
    ToDelete = glob(ClusterPath + '/*')
    for File in ToDelete: 
        if File not in glob(ClusterPath + '/times_*'): os.remove(File)
#    os.removedirs(ClusterPath)
    
    if Return: return(Clusters)
    else: return(None)


def QuantifyTTLsPerRec(AnalogTTLs, Data=[]):
    print('Get TTL timestamps... ', end='')
    if AnalogTTLs:
        Threshold = GetTTLThreshold(Data)
        TTLs = []
        for _ in range(1, len(Data)):
            if Data[_] > Threshold:
                if Data[_-1] < Threshold: TTLs.append(_)
        
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


def SliceData(Data, Rec, TTLs, DataCh, NoOfSamplesBefore, 
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
        
        Array[TTL] = Data[Rec][Start:End, DataCh[0]-1]
        
        if len(Array[TTL]) != End-Start: TTLsToFix.append(TTL)
            
    Array = FixTTLs(Array, TTLsToFix)
    
    print('Done.')
    return(Array)


def UnitsSpks(Clusters, AnalysisFile, AnalysisKey, Rec='0', Override={}):
    try: Rec = "{0:02d}".format(int(Rec))
    except ValueError: pass
    
    Units = {Rec: {}}
    for Ch in Clusters[Rec].keys():
        Units[Rec][Ch] = SepSpksPerCluster(Clusters[Rec][Ch], Ch)
    
    Group = AnalysisKey + '/Units'
    Hdf5F.WriteUnits(Units, AnalysisFile, Group)
    return(None)


## Level 2
def UnitsPSTH(Clusters, TTLCh, Rate, AnalysisFile, AnalysisKey, Rec='0', 
              TimeBeforeTTL=0, TimeAfterTTL=300, AnalogTTLs=True, 
              Override={}):
    try: Rec = "{0:02d}".format(int(Rec))
    except ValueError: pass
        
    if 'TTLs' in Override.keys(): TTLs = Override['TTLs']
    else: TTLs = QuantifyTTLsPerRec(AnalogTTLs, TTLCh)
    
    XValues = GetRecXValues(TTLs, Rate, TimeBeforeTTL, TimeAfterTTL)
    Units = {Rec: {}}
    
    print('Preparing histograms and spike waveforms...')
    for CKey in Clusters[Rec].keys():
        Classes = np.unique(Clusters[Rec][CKey]['ClusterClass'])
#            Units[URecS][CKey] = SepSpksPerCluster(Clusters[Rec][CKey], CKey)
        Units[Rec][CKey] = {'PSTH': {}}
        
        for Class in Classes:
            ClassIndex = Clusters[Rec][CKey]['ClusterClass'] == Class
            if not len(Clusters[Rec][CKey]['Spikes'][ClassIndex,:]): continue
            
            Class = "{0:02d}".format(int(Class))
            
            Hist = np.array([])
            for TTL in range(len(TTLs)):
                Firing = Clusters[Rec][CKey]['Timestamps'][ClassIndex] \
                         - (TTLs[TTL]/(Rate/1000))
                Firing = Firing[(Firing >= XValues[0]) * 
                                (Firing < XValues[-1])]
                
                Hist = np.concatenate((Hist, Firing)); del(Firing)
            
            Units[Rec][CKey]['PSTH'][Class] = Hist[:]
            del(Hist)
        
        print(CKey+':', str(len(Classes)), 'clusters.')
    
    del(TTLs)
    
    Group = AnalysisKey + '/Units'
    Hdf5F.WriteUnits(Units, AnalysisFile, Group, XValues)
    del(Units)
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
            if AnalogTTLs: 
                Raw, _, Files = Hdf5F.LoadOEKwik(RecFolder, AnalogTTLs)
            else: 
                Raw, Events, _, Files = Hdf5F.LoadOEKwik(RecFolder, AnalogTTLs)
            
            ExpInfo = Hdf5F.ExpExpInfo(RecFolder, DirList, FileName)
            OEProc, RHAProc, ABRProc = Hdf5F.GetProc(Raw, Board)
            
            if AnalogTTLs: Raw = Hdf5F.GetRecKeys(Raw, [0], AnalogTTLs)
            else:
                Raw, EventRec = Hdf5F.GetRecKeys(Raw, Events, AnalogTTLs)
#                TTLsPerRec = GetTTLInfo(Events, EventRec, ABRTTLCh)
            
            Rate = Raw['100']['info']['0']['sample_rate']
            ABRs, Info = ABRCalc(Raw, Rate, ExpInfo)
            Hdf5F.WriteABR(ABRs, Info['XValues'], Group, Info['Path'], AnalysisFile)
        

#FileName, ABRCh=[1], ABRTTLCh=1, ABRTimeBeforeTTL=0, 
#                ABRTimeAfterTTL=12, FilterFreq=[300, 3000], FilterOrder=4, 
#                StimType='Sound', AnalogTTLs=False, Board='OE', Type='kwik',
#                Override={}
def ABRCalc(Data, Rate, ExpInfo, DataInfo, Stim, TimeBeforeTTL=3, 
            TimeAfterTTL=12, AnalogTTLs=True, ABRCh=5, TTLCh=17, 
            FilterFreq=[300, 3000], FilterOrder=5, TTLsPerRec=[]):
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
    ABRs = {}; Info = {}
    
    NoOfSamplesBefore = TimeBeforeTTL*int(Rate*10**-3)
    NoOfSamplesAfter = TimeAfterTTL*int(Rate*10**-3)
    NoOfSamples = NoOfSamplesBefore + NoOfSamplesAfter
    
    Info['Frequency'] = ExpInfo['Hz']
    
    Info['XValues'] = (range(-NoOfSamplesBefore, 
                             NoOfSamples-NoOfSamplesBefore)/Rate)*10**3
    
    for Rec in Data.keys():
        print('Slicing and filtering ABRs Rec ', str(Rec), '...')
        
        if len(Data[Rec]) < 50*Rate:
            print('Rec', Rec, 'is broken!!!')
            continue
        
        if AnalogTTLs:
            TTLs = QuantifyTTLsPerRec(Data, Rec, AnalogTTLs, TTLCh)
            ABR = SliceData(Data, Rec, TTLs, ABRCh, NoOfSamplesBefore, 
                            NoOfSamplesAfter, NoOfSamples, AnalogTTLs)
        else:
            RawTime, TTLs = QuantifyTTLsPerRec(Data, Rec, AnalogTTLs, 
                                               TTLsPerRec=TTLsPerRec)
            ABR = SliceData(Data, Rec, TTLs, ABRCh, 
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
    
    return(ABRs, Info)


def GPIASAnalysis(Data, DataInfo, Rate, RecFolderNo, GPIASCh=1, GPIASTTLCh=1, GPIASTimeBeforeTTL=50, 
                  GPIASTimeAfterTTL=150, FilterFreq=[70, 400], FilterOrder=4, 
                  AnalogTTLs=False, Board='OE', Override={}):
    
    print('set paths...')
#    if 'FileName' in Override: FileName =  Override['FileName']
#    else: 
#        FileName = glob('*.hdf5'); FileName.sort()
#        FileName = FileName[RecFolderNo-1]
    
#    DataInfo = Hdf5F.LoadDict('/DataInfo', FileName)
#    
    if 'DirList' in Override: DirList = Override['DirList']
    else: DirList = glob('KwikFiles/*'); DirList.sort()
    
    if 'AnalysisFile' in Override: AnalysisFile =  Override['AnalysisFile']
    else: AnalysisFile = '../' + DataInfo['AnimalName'] + '-Analysis.hdf5'
    
    RecFolder = DirList[RecFolderNo-1]
#    
#    for Path in ['Freqs', 'FreqOrder', 'FreqSlot']:
#        DataInfo[Path] = Hdf5F.LoadDataset('/DataInfo/'+Path, FileName)
        
#    if AnalogTTLs: Raw, _, Files = Hdf5F.LoadOEKwik(RecFolder, AnalogTTLs)
#    else: Raw, Events, _, Files = Hdf5F.LoadOEKwik(RecFolder, AnalogTTLs)
#    
#    if AnalogTTLs: Raw = Hdf5F.GetRecKeys(Raw, [0], AnalogTTLs)
#    else:
#        Raw, EventRec = Hdf5F.GetRecKeys(Raw, Events, AnalogTTLs)
#        TTLsPerRec = GetTTLInfo(Events, EventRec, GPIASTTLCh)
#    
#    OEProc = Hdf5F.GetProc(Raw, Board)[0]
    
#    Rate = Raw[OEProc]['info']['0']['sample_rate']
    NoOfSamplesBefore = int(round((GPIASTimeBeforeTTL*Rate)*10**-3))
    NoOfSamplesAfter = int(round((GPIASTimeAfterTTL*Rate)*10**-3))
    NoOfSamples = NoOfSamplesBefore + NoOfSamplesAfter
    
    XValues = (range(-NoOfSamplesBefore, NoOfSamples-NoOfSamplesBefore)
               /Rate)*10**3
    
    GPIAS = {''.join([str(Freq[0]), '-', str(Freq[1])]): {}
             for Freq in DataInfo['NoiseFrequency']}
    
    for Freq in GPIAS.keys():
        GPIAS[Freq]['NoGap'] = []; GPIAS[Freq]['Gap'] = []
    
    for Rec in Data.keys():
        print('Slicing and filtering Rec ', Rec, '...')
        Freq = DataInfo['FreqOrder'][int(Rec)][0]; 
        Trial = DataInfo['FreqOrder'][int(Rec)][1];
        
        SFreq = ''.join([str(DataInfo['NoiseFrequency'][Freq][0]), '-', 
                         str(DataInfo['NoiseFrequency'][Freq][1])])
        
        if Trial % 2 == 0: STrial = 'NoGap'
        else: STrial = 'Gap'
        
        if AnalogTTLs:
            TTLs = QuantifyTTLsPerRec(Data, Rec, AnalogTTLs, GPIASTTLCh)
            GD = SliceData(Data, Rec, TTLs, GPIASCh, NoOfSamplesBefore, 
                           NoOfSamplesAfter, NoOfSamples, AnalogTTLs)
#        else:
#            RawTime, TTLs = QuantifyTTLsPerRec(Data, Rec, AnalogTTLs, 
#                                               TTLsPerRec=TTLsPerRec)
#            GD = SliceData(Raw, OEProc, Rec, TTLs, GPIASCh, NoOfSamplesBefore, 
#                           NoOfSamplesAfter, NoOfSamples, AnalogTTLs, RawTime)
        
        if GPIAS[SFreq][STrial] == []: GPIAS[SFreq][STrial] = GD[:]
        else: 
            # Cumulative moving average
            ElNo = len(GPIAS[SFreq][STrial])
            GPIAS[SFreq][STrial] = ((np.mean(GPIAS[SFreq][STrial], axis=0)
                                     *ElNo) + GD[:])/(ElNo+1)
    
    for Freq in GPIAS.keys():
        # Fix array location on dict
        GPIAS[Freq]['Gap'] = GPIAS[Freq]['Gap'][0]
        GPIAS[Freq]['NoGap'] = GPIAS[Freq]['NoGap'][0]
        
        # Bandpass filter
        GPIAS[Freq]['Gap'] = FilterSignal(GPIAS[Freq]['Gap'], Rate, FilterFreq, 
                                          FilterOrder, 'bandpass')
        GPIAS[Freq]['NoGap'] = FilterSignal(GPIAS[Freq]['NoGap'], Rate, 
                                            FilterFreq, FilterOrder, 
                                            'bandpass')
        
        # V to mV
        GPIAS[Freq]['Gap'] = GPIAS[Freq]['Gap'] * 1000
        GPIAS[Freq]['NoGap'] = GPIAS[Freq]['NoGap'] * 1000
        
        # Amplitude envelope
        GapAE = abs(signal.hilbert(GPIAS[Freq]['Gap']))
        NoGapAE = abs(signal.hilbert(GPIAS[Freq]['NoGap']))
        
        # RMS
        BGStart = 0; BGEnd = NoOfSamplesBefore - 1
        PulseStart = NoOfSamplesBefore; PulseEnd = len(GapAE) - 1
#        BinSize = XValues[-1] - XValues[-2]
        
#        GapRMSBG = sum(GapAE[BGStart:BGEnd] * BinSize)**0.5
#        GapRMSPulse = sum(GapAE[PulseStart:PulseEnd] * BinSize)**0.5
        GapRMSBG = (np.mean(GapAE[BGStart:BGEnd]**2))**0.5
        GapRMSPulse = (np.mean(GapAE[PulseStart:PulseEnd]**2))**0.5
        GapRMS = GapRMSPulse - GapRMSBG
        
#        NoGapRMSBG = sum(NoGapAE[BGStart:BGEnd] * BinSize)**0.5
#        NoGapRMSPulse = sum(NoGapAE[PulseStart:PulseEnd] * BinSize)**0.5
        NoGapRMSBG = (np.mean(NoGapAE[BGStart:BGEnd]**2))**0.5
        NoGapRMSPulse = (np.mean(NoGapAE[PulseStart:PulseEnd]**2))**0.5
        NoGapRMS = NoGapRMSPulse - NoGapRMSBG
        
        # GPIAS index (How much Gap is different from NoGap)
        GPIAS[Freq]['GPIASIndex'] = (NoGapRMS-GapRMS)/NoGapRMS
    
    if 'DirList' in Override:
        RecExp = RecFolder.split('/')[1]
        Hdf5F.WriteGPIAS(GPIAS, XValues, RecFolder, AnalysisFile, RecExp)
    else:
        Hdf5F.WriteGPIAS(GPIAS, XValues, RecFolder, AnalysisFile)
    
    return(None)


## Classes
class Plot():
    ## Level 0
    def Set(Backend='Qt5Agg', AxesObj=(), FigObj=(), FigTitle='', Params=False, 
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
    
    
    ## Level 1
    def RawCh(Ch, Lines, Cols, XValues=[], Slice=[], Leg=[], Colors='', 
              Visible=True, Save=True, FigName=''):
        Params = Plot.Set(Backend='Qt5Agg', Params=True)
        from matplotlib import rcParams; rcParams.update(Params)
        from matplotlib import pyplot as plt
        
        PlotNo = len(Ch)
        Axes = [plt.subplot(Lines, Cols, _+1) for _ in range(PlotNo)]
        
        if XValues == []: XValues = range(len(Ch[0]))
        
        for Ind, Ax in enumerate(Axes):
            if Slice: Line = Ax.plot(XValues[Slice[0]:Slice[1]], 
                                     Ch[Ind][Slice[0]:Slice[1]])
            else: Line = Ax.plot(XValues, Ch[Ind])
            
            if Colors: Line[0].set_color(Colors[Ind])
            if Leg: Line[0].set_label(Leg[Ind]); Ax.legend(loc='best')
        
        if Save: 
            if not FigName: 
                Now = datetime.now().strftime("%Y%m%d%H%M%S")
                FigName = 'RawCh-' + Now + '.svg'
            
            print('Writing to', FigName+'... ', end='')
            plt.savefig(FigName, format='svg')
            print('Done.')
        
        if Visible: plt.show()
        return(None)
    
    
    def UnitPerCh(ChDict, Ch, XValues, FigName, Ext):
        ClusterNo = len(ChDict['Spks'])
        if ClusterNo == 0: print(Ch, 'had no spikes :('); return(None)
        
        PSTHNo = 0
        for Class in ChDict['PSTH'].keys(): 
            PSTHNo += len(ChDict['PSTH'][Class])
        
        if not PSTHNo:
            print('No Spks in PSTHs of this channel :( Skipping channel...')
            return(None)
        
        Params = Plot.Set(Backend='Agg', Params=True)
        from matplotlib import rcParams; rcParams.update(Params)
        from matplotlib import pyplot as plt
        
        Fig, Axes = plt.subplots(ClusterNo,2, figsize=(8, 3*ClusterNo))
        SpksYLabel = 'Voltage [µV]'; SpksXLabel = 'Time [ms]'
        PSTHYLabel = 'Number of spikes in channel'; PSTHXLabel = 'Time [ms]'
    #    SpanLabel = 'Sound pulse'
        
        for Class in ChDict['Spks'].keys():
            SpkNo = len(ChDict['Spks'][Class])
            print(str(SpkNo), 'Spks in cluster', Class)
            if len(ChDict['PSTH'][Class]) == 0:
                print('No Spks in PSTH. Skipping class...')
                continue
            else:
                print('Max of', len(ChDict['PSTH'][Class]), 'Spks in PSTH.')
            
            if SpkNo > 50: 
                SpkNo = np.arange(SpkNo)
                np.random.shuffle(SpkNo)
                SpkNo = SpkNo[:50]
            else: 
                SpkNo = np.arange(SpkNo)
            
            
            for Spike in SpkNo:
                x = np.arange(len(ChDict['Spks'][Class][Spike])) / 30
                if ClusterNo == 1: Axes[0].plot(x, ChDict['Spks'][Class][Spike], 'r')
                else: Axes[int(Class)-1][0].plot(x, ChDict['Spks'][Class][Spike], 'r')
            
            x = np.arange(len(np.mean(ChDict['Spks'][Class], axis=0))) / 30
            if ClusterNo == 1:
                Axes[0].plot(x, np.mean(ChDict['Spks'][Class], axis=0), 'k')
                Axes[1].hist(ChDict['PSTH'][Class], XValues)
                
    #            Ind1 = list(XValues).index(0)
    #            Ind2 = list(XValues).index(int(PulseDur*1000))
                
    #            Axes[1].axvspan(XValues[Ind1], XValues[Ind2], color='k', alpha=0.3, 
    #                            lw=0, label=SpanLabel)
                
                Plot.Set(AxesObj=Axes[0], Axes=True)
                Plot.Set(AxesObj=Axes[1], Axes=True)
                Axes[0].set_ylabel(SpksYLabel); Axes[0].set_xlabel(SpksXLabel)
                Axes[1].set_ylabel(PSTHYLabel); Axes[1].set_xlabel(PSTHXLabel)
                
            else:
                Axes[int(Class)-1][0].plot(x, np.mean(ChDict['Spks'][Class], axis=0), 'k')
                Axes[int(Class)-1][1].hist(ChDict['PSTH'][Class], XValues)
                
    #            Ind1 = list(XValues).index(0)
    #            Ind2 = list(XValues).index(int(PulseDur*1000))
                
    #            Axes[int(Class)-1][1].axvspan(XValues[Ind1], XValues[Ind2], 
    #                                          color='k', alpha=0.3, lw=0, 
    #                                          label=SpanLabel)
                
                Plot.Set(AxesObj=Axes[int(Class)-1][0], Axes=True)
                Plot.Set(AxesObj=Axes[int(Class)-1][1], Axes=True)
                Axes[int(Class)-1][0].set_ylabel(SpksYLabel)
                Axes[int(Class)-1][0].set_xlabel(SpksXLabel)
                Axes[int(Class)-1][1].set_ylabel(PSTHYLabel)
                Axes[int(Class)-1][1].set_xlabel(PSTHXLabel)
        
        FigTitle = FigName.split('/')[-1][:-4]
        Plot.Set(FigObj=Fig, FigTitle=FigTitle, Plot=True)
        print('Writing to', FigName+'... ', end='')
        Fig.savefig(FigName, format=Ext)
        print('Done.')
        return(None)
    
    
    def UnitTestBinSizeCh(ChDict, Ch, XValuesList, FigName, Ext):
        ClusterNo = len(ChDict['Spks'])
        if ClusterNo == 0: print(Ch, 'had no spikes :('); return(None)
        
        PSTHNo = 0
        for Class in ChDict['PSTH'].keys(): 
            PSTHNo += len(ChDict['PSTH'][Class])
        
        if not PSTHNo:
            print('No Spks in PSTHs of this channel :( Skipping channel...')
            return(None)
        
        Params = Plot.Set(Backend='Agg', Params=True)
        from matplotlib import rcParams; rcParams.update(Params)
        from matplotlib import pyplot as plt
        
        Fig, Axes = plt.subplots(ClusterNo, len(XValuesList), 
                                 figsize=(4*len(XValuesList), 3*ClusterNo))
        
        PSTHYLabel = 'Number of spikes in channel'; PSTHXLabel = 'Time [ms]'
    #    SpanLabel = 'Sound pulse'
        
        for Class in ChDict['Spks'].keys():
            SpkNo = len(ChDict['Spks'][Class])
            print(str(SpkNo), 'Spks in cluster', Class)
            if len(ChDict['PSTH'][Class]) == 0:
                print('No Spks in PSTH. Skipping class...')
                continue
            else:
                print('Max of', len(ChDict['PSTH'][Class]), 'Spks in PSTH.')
            
            if SpkNo > 50: 
                SpkNo = np.arange(SpkNo)
                np.random.shuffle(SpkNo)
                SpkNo = SpkNo[:50]
            else: 
                SpkNo = np.arange(SpkNo)
            
            if ClusterNo == 1:
                for XInd, XValues in enumerate(XValuesList):
                    Axes[XInd].hist(ChDict['PSTH'][Class], XValues)
                    SubTitle = str(XValues[1] - XValues[0]) + ' ms bin size'
                    Plot.Set(AxesObj=Axes[XInd], Axes=True)
                    Axes[XInd].set_ylabel(PSTHYLabel)
                    Axes[XInd].set_xlabel(PSTHXLabel)
                    Axes[XInd].set_title(SubTitle)
                
            else:
                for XInd, XValues in enumerate(XValuesList):
                    Axes[int(Class)-1][XInd].hist(ChDict['PSTH'][Class], XValues)
                    SubTitle = str(XValues[1] - XValues[0]) + ' ms bin size'
                    Plot.Set(AxesObj=Axes[int(Class)-1][XInd], Axes=True)
                    Axes[int(Class)-1][XInd].set_ylabel(PSTHYLabel)
                    Axes[int(Class)-1][XInd].set_xlabel(PSTHXLabel)
                    Axes[int(Class)-1][XInd].set_title(SubTitle)
        
        FigTitle = FigName.split('/')[-1][:-4]
        Plot.Set(FigObj=Fig, FigTitle=FigTitle, Plot=True)
        print('Writing to', FigName+'... ', end='')
        Fig.savefig(FigName, format=Ext)
        print('Done.')
        return(None)
    
    
    ## Level 2
    def UnitsSpksPSTH(AnalysisFile, AnalysisKey, DataFolder='', Ext='svg', 
                      Override={}, Procs=8):
        Units, XValues = Hdf5F.LoadUnits(AnalysisFile, AnalysisKey, Override)
        if 'XValues' in Override: XValues = Override['XValues']
        
#        AnalysisFileName = AnalysisFile.split('/')[-1]
        
        if not DataFolder: DataFolder = '/'.join(AnalysisFile.split('/')[:-1])
        DataFolder += '/Figures'
        
        os.makedirs(DataFolder, exist_ok=True)
        
    #    Thrash = {} [a[s:s+8] for s in range(0,len(a),8)]
        for RKey in Units:
            ChNo = len(Units[RKey])
            ProcLists = [[] for _ in range(0, ChNo, Procs)]
            
            for Ind, Ch in enumerate(Units[RKey]):
                FigName = ''.join([DataFolder, '/', DataFolder.split('/')[0],
                                   '_', DataFolder.split('/')[-1], '_Rec', 
                                   RKey, '_', Ch, '.', Ext])
                ProcLists[Ind//Procs].append(
                        Process(target=Plot.UnitPerCh, 
                                args=(Units[RKey][Ch], Ch, XValues, FigName, Ext))
                        )
            
            for ProcList in ProcLists:
                for Proc in ProcList:
                    Proc.start(); print('PID =', Proc.pid)
                Proc.join()
            
            return(None)
#                UnitPlot.start(); print('PID =', UnitPlot.pid)
#                UnitPlot.join()
            
    
    
    def UnitsPSTHTestBin(XValuesList, AnalysisFile, DataFolder='', Ext='svg', 
                         Override={}):
        Units = Hdf5F.LoadUnits(AnalysisFile, Override)[0]
        
        AnalysisFileName = AnalysisFile.split('/')[-1]
        
        if not DataFolder: DataFolder = '/'.join(AnalysisFile.split('/')[:-1])
        DataFolder += '/Figures/'
        
        os.makedirs(DataFolder, exist_ok=True)
        
    #    Thrash = {}
        for RKey in Units:
            for Ch in Units[RKey]:
                FigName = ''.join([DataFolder, AnalysisFileName[:-5], '_Rec', 
                                   RKey, '_', Ch, '-BinSizeTest.', Ext])
                UnitPlot = Process(target=Plot.UnitTestBinSizeCh, 
                                   args=(Units[RKey][Ch], Ch, XValuesList, 
                                         FigName, Ext))
                UnitPlot.start(); print('PID =', UnitPlot.pid)
                UnitPlot.join()


#%% Tests
#Here = os.getcwd(); Path = 'SepCh/'; ClusterPath = Here + '/' + Path
#
#RHAHeadstage = [16, 15, 14, 13, 12, 11, 10, 9, 1, 2, 3, 4, 5, 6, 7, 8]
#A16 = {'Tip': [9, 8, 10, 7, 13, 4, 12, 5, 15, 2, 16, 1, 14, 3, 11, 6],
#       'Head': [8, 7, 6, 5, 4, 3, 2, 1, 9, 10, 11, 12, 13, 14, 15, 16]}
#
#ChannelMap = RemapChannels(A16['Tip'], A16['Head'], RHAHeadstage)
#
#Rate = np.array([25000]); AnalogTTLs = True; Override={}
#
#Files = [['IntanSound/Arch3n2-43S_184905.int', [154934, 4121425]],
#         ['IntanSound/Arch3n2-35S_175902.int', [1702747]]]
#
#AnalysisFile = './CaMKIIaArch3_02.hdf5'
#
#for File in Files:
#    Data = {}; Rec = File[0].split('-')[-1].split('_')[0]
#    Data[Rec] = IntanBin.IntLoad(File)[0]
#    Override['RecS'] = Rec
#    
#    Clusters = ClusterizeSpks(Data, Rate, ClusterPath, AnalysisFile, 
#                              ChannelMap, Override, Return=True)
#    
##    Clusters = Hdf5F.LoadClusters(AnalysisFile)
##    UnitsSpks(Clusters, AnalysisFile, Override)
#    
#    
#    Starts = QuantifyTTLsPerRec(AnalogTTLs, Data[Rec][:, 16])
#    Override['TTLs'] = GenerateFakeTTLsRising(Starts, int(0.003*Rate), 200, int(0.090*Rate))
#    TimeBeforeTTL = 0; TimeAfterTTL = 300
#    
#    UnitsPSTH(Clusters, [], Rate, AnalysisFile, TimeBeforeTTL, TimeAfterTTL, 
#              AnalogTTLs, Override)
#
#    Plot.UnitsSpksPSTH(AnalysisFile, 'svg', Override)
#
#
#
#
#
