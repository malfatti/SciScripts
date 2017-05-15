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
import os
import numpy as np
from datetime import datetime
from glob import glob
from multiprocessing import Process
from rpy2.robjects import packages as RPackages
from rpy2 import robjects as RObj
from scipy import io, signal
from subprocess import call

## Level 0
def CumulativeMA(Base, Add, ElNo):
    if len(Base) == 0: 
        Base = Add
    else:
        Base = ((Base * ElNo) + Add)/(ElNo+1)
    
    return(Base)


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


def FindPeaks(Data, Dist, LowThr=None, HighThr=None):
    if not LowThr: LowThr = np.mean(Data) + (np.std(Data)/(len(Data)**0.5))
    if not HighThr: HighThr = np.mean(Data) + (2*np.std(Data))
    
    Range = np.arange(0, len(Data), Dist, dtype=int)
    
    PVal = [max(Data[Range[i-1]:Range[i]])
             for i in range(1, len(Range)) 
             if max(Data[Range[i-1]:Range[i]]) > LowThr 
             and max(Data[Range[i-1]:Range[i]]) < HighThr]
#    PInd = [np.where(Data == max(Data[Range[i-1]:Range[i]]))[0][0]
#             for i in range(1, len(Range)) if max(Data[Range[i-1]:Range[i]]) > Thr]
    PInd = [np.where(Data == Val)[0][0] for Val in PVal]
    
    return(PInd, PVal)


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
#        Top = max(TTLCh) + abs(min(TTLCh))
#        Bot = min(TTLCh) + abs(min(TTLCh))
#        Threshold = max(TTLCh) - (Top - Bot)/3 + 2*(np.std(TTLCh))
        Amp = max(TTLCh) - min(TTLCh)
        Threshold = max(TTLCh) - (Amp/4)
        return(Threshold) 


def PSD(Data, Rate, Scaling='density'):
    Window = signal.hanning(len(Data)//(Rate/1000))
    F, PxxSp = signal.welch(Data, Rate, Window, nperseg=len(Window), 
                            noverlap=0, scaling=Scaling)
    
    return(F, PxxSp)


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


def StrRange(Start='a', End='e', Step=1):
    if max(len(Start), len(End)) > 1: 
        print('Only 1-char length strings are accepted.')
        return(None)
    else:
        Range = map(chr, range(ord(Start), ord(End), Step))
        return(Range)


## Level 1
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


def SignalIntensity(Data, Rate, FreqBand, Ref):
    Intensity = {}
    
    if Data.shape[-1] == 2: F, PxxSp = PSD(Data[:,0], Rate)
    else: F, PxxSp = PSD(Data, Rate)
    
    Start = np.where(F > FreqBand[0])[0][0]-1
    End = np.where(F > FreqBand[1])[0][0]-1
#    BinSize = F[1] - F[0]
#    RMS = (sum(PxxSp) * BinSize)**0.5
    RMS = (sum(PxxSp[Start:End]))**0.5
    dB = 20*(np.log10(RMS/Ref)) + 94
    
    Intensity['PSD'] = [F, PxxSp]
    Intensity['RMS'] = RMS
    Intensity['dB'] = dB
    
    return(Intensity)


def SliceData(Data, TTLs, NoOfSamplesBefore, 
              NoOfSamplesAfter, NoOfSamples, AnalogTTLs, RawTime=[]):
    print('Slicing data around TTL...')
#    Array = [[0 for _ in range(NoOfSamples)] for _ in range(len(TTLs))]
    Array = np.zeros((len(TTLs), NoOfSamples))
    TTLsToFix = []
    
    for TTL in range(len(TTLs)):
        if AnalogTTLs: TTLLoc = int(TTLs[TTL])
        else: TTLLoc = int(RawTime.index(TTLs[TTL]))#)/Rate)
        
        Start = TTLLoc-NoOfSamplesBefore
        End = TTLLoc+NoOfSamplesAfter
        
        if Start < 0: Start = 0; End = End+(Start*-1); TTLsToFix.append(TTL)
        
        Array[TTL] = Data[Start:End]
        
        if len(Array[TTL]) != End-Start: TTLsToFix.append(TTL)
            
    if TTLsToFix: Array = FixTTLs(Array, TTLsToFix)
    
    print('Done.')
    return(Array)


## Level 3
#def ABRAnalysis(FileName, ABRCh=[1], ABRTTLCh=1, ABRTimeBeforeTTL=0, 
#                ABRTimeAfterTTL=12, FilterFreq=[300, 3000], FilterOrder=4, 
#                StimType='Sound', AnalogTTLs=False, Board='OE', Type='kwik',
#                Override={}):
#    print('Load DataInfo...')
#    DirList = glob('KwikFiles/*'); DirList.sort()
#    DataInfo = Hdf5F.LoadDict('/DataInfo', FileName)
#    
#    AnalysisFile = './' + DataInfo['AnimalName'] + '-Analysis.hdf5'
#    Now = datetime.now().strftime("%Y%m%d%H%M%S")
#    Here = os.getcwd().split(sep='/')[-1]
#    Group = Here + '-ABRs_' + Now
#        
#    for Stim in StimType:
#        if Override != {}: 
#            if 'Stim' in Override.keys(): Stim = Override['Stim']
#        
#        Exps = Hdf5F.LoadExpPerStim(Stim, DirList, FileName)
#        
#        for RecFolder in Exps:
#            if AnalogTTLs: 
#                Raw, _, Files = Hdf5F.LoadOEKwik(RecFolder, AnalogTTLs)
#            else: 
#                Raw, Events, _, Files = Hdf5F.LoadOEKwik(RecFolder, AnalogTTLs)
#            
#            ExpInfo = Hdf5F.ExpExpInfo(RecFolder, DirList, FileName)
#            OEProc, RHAProc, ABRProc = Hdf5F.GetProc(Raw, Board)
#            
#            if AnalogTTLs: Raw = Hdf5F.GetRecKeys(Raw, [0], AnalogTTLs)
#            else:
#                Raw, EventRec = Hdf5F.GetRecKeys(Raw, Events, AnalogTTLs)
##                TTLsPerRec = GetTTLInfo(Events, EventRec, ABRTTLCh)
#            
#            Rate = Raw['100']['info']['0']['sample_rate']
#            ABRs, Info = ABRCalc(Raw, Rate, ExpInfo)
#            Hdf5F.WriteABR(ABRs, Info['XValues'], Group, Info['Path'], AnalysisFile)
#        

#FileName, ABRCh=[1], ABRTTLCh=1, ABRTimeBeforeTTL=0, 
#                ABRTimeAfterTTL=12, FilterFreq=[300, 3000], FilterOrder=4, 
#                StimType='Sound', AnalogTTLs=False, Board='OE', Type='kwik',
#                Override={}
def ABRCalc(Data, Rate, ExpInfo, DataInfo, Stim, TimeBeforeTTL=3, 
            TimeAfterTTL=12, AnalogTTLs=True, ABRCh=5, TTLCh=17, 
            FilterFreq=[300, 3000], FilterOrder=5, TTLsPerRec=[]):
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
            TTLs = QuantifyTTLsPerRec(AnalogTTLs, Data[Rec][:,TTLCh-1])
            ABR = SliceData(Data[Rec][:, ABRCh-1], TTLs, NoOfSamplesBefore, 
                            NoOfSamplesAfter, NoOfSamples, AnalogTTLs)
#        else:
#            RawTime, TTLs = QuantifyTTLsPerRec(Data, Rec, AnalogTTLs, 
#                                               TTLsPerRec=TTLsPerRec)
#            ABR = SliceData(Data[Rec][:,ABRCh-1], TTLs,NoOfSamplesBefore, 
#                            NoOfSamplesAfter, NoOfSamples, AnalogTTLs, RawTime)
        
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



## Classes
class GPIAS():
    ## Level 0
    def CheckGPIASRecs(Data, SizeLimits, Plot=False):
        ToCheck = [Rec for Rec in Data.keys() 
                       if len(Data[Rec])<min(SizeLimits)
                       or len(Data[Rec])>max(SizeLimits)]
        
        if ToCheck:
            if Plot:
                Params = {'backend': 'TkAgg'}
                from matplotlib import rcParams; rcParams.update(Params)
                import matplotlib.pyplot as plt
                
                for Rec in ToCheck:
                    print('Showing Rec', Rec+', size', Data[Rec].shape[0])
                    plt.plot(Data[Rec])
                    plt.show()
                
            return(ToCheck)
        else:
            print('All recs within expected size.')
            return(None)

    def IndexCalc(Data, Keys, PulseSampleStart, SliceSize):
        Index = {}
        for Key in Keys:
            BGStart = 0; BGEnd = SliceSize
            PulseStart = PulseSampleStart; PulseEnd = PulseSampleStart + SliceSize
            
            ResRMSBG = (np.mean(Data[Key[0]][BGStart:BGEnd]**2))**0.5
            ResRMSPulse = (np.mean(Data[Key[0]][PulseStart:PulseEnd]**2))**0.5
#            ResRMS = ResRMSPulse
            if ResRMSPulse < ResRMSBG: ResRMS = ResRMSPulse
            else: ResRMS = ResRMSPulse - ResRMSBG
            
            RefRMSBG = (np.mean(Data[Key[1]][BGStart:BGEnd]**2))**0.5
            RefRMSPulse = (np.mean(Data[Key[1]][PulseStart:PulseEnd]**2))**0.5
#            RefRMS = RefRMSPulse
            if RefRMSPulse < RefRMSBG: RefRMS = RefRMSPulse
            else: RefRMS = RefRMSPulse - RefRMSBG
            
            # GPIAS index (How much Res is different from Ref)
            Index[Key[2]] = (RefRMS-ResRMS)/RefRMS
        
        return(Index)
    
    
    def PreallocateDict(DataInfo, PrePostFreq):
        Dict = {}
        for Key in ['Trace', 'Index']:
            Dict[Key] = {''.join([str(Freq[0]), '-', str(Freq[1])]): {}
                         for Freq in DataInfo['NoiseFrequency']}
        
        Dict['Trace'][PrePostFreq]['Pre'] = []
        Dict['Trace'][PrePostFreq]['Post'] = []
        Dict['Index'][PrePostFreq]['Pre'] = []
        Dict['Index'][PrePostFreq]['Post'] = []
        
        for Freq in Dict['Trace'].keys():
            Dict['Trace'][Freq]['NoGap'] = []; Dict['Trace'][Freq]['Gap'] = []
            Dict['Index'][Freq]['NoGap'] = []; Dict['Index'][Freq]['Gap'] = []
        
        return(Dict)
    
    
    def OrganizeRecs(Dict, Data, DataInfo, AnalogTTLs, NoOfSamplesBefore, 
                     NoOfSamplesAfter, NoOfSamples):
        for Rec in Data.keys():
            print('Slicing and filtering Rec ', Rec, '...')
            Freq = DataInfo['FreqOrder'][int(Rec)][0]; 
            Trial = DataInfo['FreqOrder'][int(Rec)][1];
            
            SFreq = ''.join([str(DataInfo['NoiseFrequency'][Freq][0]), '-', 
                             str(DataInfo['NoiseFrequency'][Freq][1])])
            
            if Trial == -1: STrial = 'Pre'
            elif Trial == -2: STrial = 'Post'
            elif Trial % 2 == 0: STrial = 'NoGap'
            else: STrial = 'Gap'
            
            if AnalogTTLs:
                TTLs = QuantifyTTLsPerRec(AnalogTTLs, Data[Rec][:,DataInfo['TTLCh']-1])
    #            if len(TTLs) > 1: TTLs = [TTLs[0]]
                GD = SliceData(Data[Rec][:,DataInfo['PiezoCh'][0]-1], TTLs, 
                               NoOfSamplesBefore, NoOfSamplesAfter, NoOfSamples, 
                               AnalogTTLs)
    #        else:
    #            RawTime, TTLs = QuantifyTTLsPerRec(Data, Rec, AnalogTTLs, 
    #                                               TTLsPerRec=TTLsPerRec)
    #            GD = SliceData(Raw, OEProc, Rec, TTLs, GPIASCh, NoOfSamplesBefore, 
    #                           NoOfSamplesAfter, NoOfSamples, AnalogTTLs, RawTime)
            
            Dict['Index'][SFreq][STrial].append(GD[0])
            Dict['Trace'][SFreq][STrial].append(GD[0])
        
        return(Dict)
        
    
    ## Level 1    
    def Analysis(Data, DataInfo, Rate, AnalysisFile, AnalysisKey, 
                 GPIASTimeBeforeTTL=50, GPIASTimeAfterTTL=150, 
                 FilterFreq=[70, 400], FilterOrder=4, SliceSize=100,
                 AnalogTTLs=True, Return=False):
        
        NoOfSamplesBefore = int(round((GPIASTimeBeforeTTL*Rate)*10**-3))
        NoOfSamplesAfter = int(round((GPIASTimeAfterTTL*Rate)*10**-3))
        NoOfSamples = NoOfSamplesBefore + NoOfSamplesAfter
        
        XValues = (range(-NoOfSamplesBefore, NoOfSamples-NoOfSamplesBefore)
                   /Rate)*10**3
        
        PrePostFreq = DataInfo['FreqOrder'][0][0]
        PrePostFreq = '-'.join([str(DataInfo['NoiseFrequency'][PrePostFreq][0]),
                                str(DataInfo['NoiseFrequency'][PrePostFreq][1])])
        
        GPIASData = GPIAS.PreallocateDict(DataInfo, PrePostFreq)
        GPIASData = GPIAS.OrganizeRecs(GPIASData, Data, DataInfo, AnalogTTLs,
                                       NoOfSamplesBefore, NoOfSamplesAfter, 
                                       NoOfSamples)
        
        SliceSize = int(SliceSize * (Rate/1000))
        
        for Freq in GPIASData['Index'].keys():
            for Key in GPIASData['Index'][Freq].keys():
                # Average trials for traces
                GPIASData['Trace'][Freq][Key] = np.mean(GPIASData['Trace'][Freq][Key], axis=0)
                
                # Bandpass filter
                GPIASData['Trace'][Freq][Key] = FilterSignal(GPIASData['Trace'][Freq][Key], 
                                                           Rate, FilterFreq, FilterOrder, 
                                                           'bandpass')
                
                for Tr in range(len(GPIASData['Index'][Freq][Key])):
                    # Bandpass filter
                    AE = FilterSignal(GPIASData['Index'][Freq][Key][Tr], Rate, 
                                         FilterFreq, FilterOrder, 'bandpass')
                    
                    # Amplitude envelope
                    AE = abs(signal.hilbert(AE))
                    
                    GPIASData['Index'][Freq][Key][Tr] = AE
                
                GPIASData['Index'][Freq][Key] = np.mean(GPIASData['Index'][Freq][Key], axis=0)
                
            # RMS
            if Freq == PrePostFreq: Keys = [['Gap', 'NoGap', 'GPIASIndex'], 
                                            ['Post', 'Pre', 'PrePost']]
            else: Keys = [['Gap', 'NoGap', 'GPIASIndex']]
            
            GPIASData['Index'][Freq] = GPIAS.IndexCalc(
                                           GPIASData['Index'][Freq], Keys, 
                                           NoOfSamplesBefore, SliceSize)
        
        Hdf5F.WriteGPIAS(GPIASData, AnalysisKey, AnalysisFile, XValues)
        
        if Return: return(GPIASData, XValues)
        else: return(None)


class Plot():
    ## Level 0
    def Set(Backend='TkAgg', AxesObj=(), FigObj=(), FigTitle='', Params=False, 
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
    
    
    def SignificanceBar(X, Y, Text, Ax, TicksDir='down', lw=1, color='k'):
        if TicksDir == 'down':
            from matplotlib.markers import TICKDOWN as Tick
            Yy = max(Y)+(max(Y)*0.05)
        elif TicksDir == 'up':
            from matplotlib.markers import TICKUP as Tick
            Yy = max(Y)-(max(Y)*0.01)
        else: print('TicksDir should be "up" or "down".'); return(None)
        
        Ax.plot(X, Y, color=color, lw=lw, marker=Tick)
        Ax.text(sum(X)/2, Yy, Text, fontsize=6, ha='center', va='center')
        return(None)
    
    
    ## Level 1
    def GPIAS(GPIASData, XValues, SoundPulseDur, FigName, Ext='svg', 
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
            Axes[FInd].set_title(SubTitle)
            Axes[FInd].set_ylabel(YLabel); Axes[FInd].set_xlabel(XLabel)
    
            Plot.Set(AxesObj=Axes[FInd], Axes=True)
        
        FigTitle = FigName.split('/')[-1]
        FigName = FigName + '.' + Ext
        Plot.Set(FigObj=Fig, FigTitle=FigTitle, Plot=True)
        
        if Save: Fig.savefig(FigName, format=Ext)
        if Visible: plt.show()
        else: plt.close()
        
        print('Done.')
        return(None)
    
    
    def RawCh(Ch, Lines, Cols, XValues=[], Slice=[], Leg=[], FigName='', Colors='', 
              Visible=True, Save=True):
        Params = Plot.Set(Params=True)
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
        else: plt.close()
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
    def UnitsSpksPSTH(Units, XValues, FigName, Mode='SpksPSTH', Ext='svg', Procs=8):
        for RKey in Units:
            ChNo = len(Units[RKey])
            ProcLists = [[] for _ in range(0, ChNo, Procs)]
            
            for Ind, Ch in enumerate(Units[RKey]):
                if Mode == 'BinSizeTest':
                    FigName = FigName + '_' + Ch + '-BinSizeTest.', Ext
                    ProcLists[Ind//Procs].append(
                        Process(target=Plot.UnitTestBinSizeCh, 
                                args=(Units[RKey][Ch], Ch, XValues, FigName, Ext))
                    )
                elif Mode == 'SpksPSTH':
                    FigName = FigName + '_' + Ch + '-UnitsPSTH.' + Ext
                    ProcLists[Ind//Procs].append(
                        Process(target=Plot.UnitPerCh, 
                                args=(Units[RKey][Ch], Ch, XValues, FigName, Ext))
                    )
            
            for ProcList in ProcLists:
                for Proc in ProcList:
                    Proc.start(); print('PID =', Proc.pid)
                Proc.join()
            
            return(None)


class Stats():
    ## Level 0
    def RCheckPackage(Packages):
        RPacksToInstall = [Pack for Pack in Packages 
                           if not RPackages.isinstalled(Pack)]
        if len(RPacksToInstall) > 0:
            print(str(RPacksToInstall), 'not installed. Install now?')
            Ans = input('[y/N]: ')
            
            if Ans.lower() in ['y', 'yes']:
                from rpy2.robjects.vectors import StrVector as RStrVector
                
                RUtils = RPackages.importr('utils')
                RUtils.chooseCRANmirror(ind=1)
                
                RUtils.install_packages(RStrVector(RPacksToInstall))
            
            else: print('Aborted.')
        
        else: print('Packages', str(Packages), 'installed.')
        
        return(None)
    
    
    def AdjustNaNs(Array):
        NaN = RObj.NA_Real
        
        for I, A in enumerate(Array):
            if A != A: Array[I] = NaN
            
        return(Array)
    
    
    ## Level 1
    def RAnOVa(DataA, DataB):
        Raov = RObj.r['aov']
        
        Data = DataA + DataB
        Trtmnt = ['A']*len(DataA) + ['B']*len(DataB)
        Formula = RObj.Formula('Data ~ Trtmnt'); FEnv = Formula.environment
        FEnv['Data'] = RObj.FloatVector(Data)
        FEnv['Trtmnt'] = RObj.StrVector(Trtmnt)
        
        Results = Raov(Formula)
        print(RObj.r['summary'](Results))
        
        return(Results)
    
    
    def RPwrAnOVa(GroupNo=RObj.NULL, SampleSize=RObj.NULL, Power=RObj.NULL, 
               SigLevel=RObj.NULL, EffectSize=RObj.NULL):
        Stats.RCheckPackage(['pwr']); Rpwr = RPackages.importr('pwr')
        
        Results = Rpwr.pwr_anova_test(k=GroupNo, power=Power, sig_level=SigLevel, 
                                      f=EffectSize, n=SampleSize)
        
        print('Calculating', Results.rx('method')[0][0] + '... ', end='')
        AnOVaResults = {}
        for Key, Value in {'k': 'GroupNo', 'n': 'SampleSize', 'f': 'EffectSize', 
                           'power':'Power', 'sig.level': 'SigLevel'}.items():
            AnOVaResults[Value] = Results.rx(Key)[0][0]
        
        print('Done.')
        return(AnOVaResults)
    
    
    def RTTest(DataA, DataB, Paired=True, Alt='less', Confidence=0.95):
        Rttest = RObj.r['t.test']
        
        DataA = Stats.AdjustNaNs(DataA); DataB = Stats.AdjustNaNs(DataB)
        
        Results = Rttest(RObj.FloatVector(DataA), RObj.FloatVector(DataB), 
                         paired=Paired, var_equal=False, alternative=Alt, 
                         conf_level=RObj.FloatVector([Confidence]), 
                         na_action=RObj.r['na.omit'])
        
        print('Calculating', Results.rx('method')[0][0] + '... ', end='')
        TTestResults = {}; Names = list(Results.names)
        for Name in Names:
            TTestResults[Name] = Results.rx(Name)[0][0]
        
        print('Done.')
        return(TTestResults)



class Units():
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
    
    
    def PSTH(Clusters, TTLCh, Rate, AnalysisFile, AnalysisKey, Rec='0', 
                  TimeBeforeTTL=0, TimeAfterTTL=300, AnalogTTLs=True, 
                  Override={}):
        try: Rec = "{0:02d}".format(int(Rec))
        except ValueError: pass
            
        if 'TTLs' in Override.keys(): TTLs = Override['TTLs']
        else: TTLs = QuantifyTTLsPerRec(AnalogTTLs, TTLCh)
        
        XValues = GetRecXValues(TTLs, Rate, TimeBeforeTTL, TimeAfterTTL)
        UnitsData = {Rec: {}}
        
        print('Preparing histograms and spike waveforms...')
        for CKey in Clusters[Rec].keys():
            Classes = np.unique(Clusters[Rec][CKey]['ClusterClass'])
    #            UnitsData[URecS][CKey] = SepSpksPerCluster(Clusters[Rec][CKey], CKey)
            UnitsData[Rec][CKey] = {'PSTH': {}}
            
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
                
                UnitsData[Rec][CKey]['PSTH'][Class] = Hist[:]
                del(Hist)
            
            print(CKey+':', str(len(Classes)), 'clusters.')
        
        del(TTLs)
        
        Group = AnalysisKey + '/Units'
        Hdf5F.WriteUnits(UnitsData, AnalysisFile, Group, XValues)
        del(UnitsData)
        return(None)
    
    
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
        
        Units.CallWaveClus(Rate, ClusterPath)
        
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
    
    
    def Spks(Clusters, AnalysisFile, AnalysisKey, Rec='0', Override={}):
        try: Rec = "{0:02d}".format(int(Rec))
        except ValueError: pass
        
        UnitsData = {Rec: {}}
        for Ch in Clusters[Rec].keys():
            UnitsData[Rec][Ch] = Units.SepSpksPerCluster(Clusters[Rec][Ch], Ch)
        
        Group = AnalysisKey + '/Units'
        Hdf5F.WriteUnits(UnitsData, AnalysisFile, Group)
        return(None)

