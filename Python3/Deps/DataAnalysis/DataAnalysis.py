#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: T. Malfatti <malfatti@disroot.org>
@year: 2017
@license: GNU GPLv3 <https://raw.githubusercontent.com/malfatti/SciScripts/master/LICENSE>
@homepage: https://github.com/Malfatti/SciScripts
"""
import numpy as np

from itertools import tee
from multiprocessing import Process
from scipy import signal


## Level 0
def BitsToVolts(Data, ChInfo):
    for R in Data.keys():
        for C, Ch in enumerate(sorted(ChInfo.keys())):
            Data[R][:,C] = Data[R][:,C] * float(ChInfo[Ch]['gain'])
    
    return(Data)


def CumulativeMA(Base, Add, ElNo):
    if len(Base) == 0: 
        Base = Add
    else:
        Base = ((Base * ElNo) + Add)/(ElNo+1)
    
    return(Base)


def FilterSignal(Signal, Rate, Frequency, FilterOrder=4, Coeff='butter', Type='bandpass'):
    if Coeff == 'butter':
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
    
    elif Coeff == 'fir':
        Freqs = np.arange(1,(Rate/2)+1)
        DesiredFreqs = np.zeros(int(Rate/2))
        DesiredFreqs[min(Frequency):max(Frequency)] = 1
        
        o = FilterOrder + ((FilterOrder%2)*-1) +1
        a = signal.firls(o, Freqs, DesiredFreqs, nyq=Rate/2)
        Signal = signal.filtfilt(a, 1.0, Signal, padtype='odd', padlen=0)
        
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


def GenTTLVector(TTLs, TTLLen, FullLen):
    TTLVec = np.zeros([FullLen, 1])
    
    for TTL in TTLs:
        TTLVec[TTL:TTL+TTLLen] = np.ones([TTLLen, 1])
    
    return(TTLVec)


def GenFakeTTLsRising(Rate, PulseDur, PauseBefore, PauseAfter, SampleStart, PulseNo):
#    Dur = Dur + PauseBetween
#    
#    FakeTTLs = []
#    for Block in Start:
#        FakeTTLs = FakeTTLs + list(range(Block, (Dur*No) + Block, Dur))
    BlockDur = PauseBefore+PulseDur+PauseAfter
    FakeTTLs = [int(((PauseBefore)+(_*BlockDur))*Rate)+SampleStart for _ in range(PulseNo)]
    return(FakeTTLs)


def GetEventEdges(EventsDict, Ch, Proc, Rec, Edge='rise'):
    if Edge.lower() == 'rise': Id = 1
    elif Edge.lower() == 'fall': Id = 0
    else: print('Edge should be "rise" or "fall"')
    
    Events = EventsDict['sampleNum'][(EventsDict['channel'] == Ch) * 
                                     (EventsDict['recordingNumber'] == Rec) * 
                                     (EventsDict['nodeId'] == Proc) *
                                     (EventsDict['eventId'] == Id) * 
                                     (EventsDict['eventType'] == 3)]
    
    return(Events)


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


def GetTTLThreshold(TTLCh, StdNo=3):
    if not StdNo: StdNo = 3
    
    if np.mean(TTLCh) > 1000: 
        print('Sinusoidal stimulation')
        Threshold = (max(TTLCh) - min(TTLCh)) / 2
        return(Threshold)
    else:
        print('square pulses stimulation')
#        Top = max(TTLCh) + abs(min(TTLCh))
#        Bot = min(TTLCh) + abs(min(TTLCh))
#        Threshold = max(TTLCh) - (Top - Bot)/3 + 2*(np.std(TTLCh))
#         Amp = max(TTLCh) - min(TTLCh)
#         Threshold = max(TTLCh) - (Amp/4)
        Threshold = np.mean(TTLCh) + StdNo*(np.std(TTLCh))
        return(Threshold) 


def MultiProcess(Function, Args, Procs=8):
    TotalNo = len(Args)
    ProcLists = [[] for _ in range(0, TotalNo, Procs)]
    
    for A, Arg in enumerate(Args):
        ProcLists[A//Procs].append(Process(target=Function, args=Arg))
    
    for ProcList in ProcLists:
        for Proc in ProcList:
            Proc.start(); print('PID =', Proc.pid)
        Proc.join()
    
    return(None)


def NestedClean(Nest):
    if 'numpy' in str(type(Nest)):
        if not Nest.size: return(None)
        else: return(Nest)
    else:
        if len(Nest) == 0: return(None)
    
    ToDel = []
    if type(Nest) == dict:
        for K, Key in Nest.items(): 
            Nest[K] = NestedClean(Key)
            
            if type(Nest[K]) == np.ndarray: 
                if not Nest[K].size: ToDel.append(K)
            else:
                if not Nest[K]: ToDel.append(K)
        
        for K in ToDel: del(Nest[K])
        return(Nest)
    
    elif type(Nest) in [list, tuple]:
        for E, El in enumerate(Nest): 
            Nest[E] = NestedClean(El)
            
            if type(Nest[E]) == np.ndarray: 
                if not Nest[E].size: ToDel.append(E)
            else:
                if not Nest[E]: ToDel.append(E)
        
        for E in ToDel: del(Nest[E])
        return(Nest)
    
    else:
        return(Nest)


def Pairwise(iterable):
    """ from https://docs.python.org/3.6/library/itertools.html#itertools-recipes
    s -> (s0,s1), (s1,s2), (s2, s3), ..."""
    
    a, b = tee(iterable)
    next(b, None)
    return zip(a, b)


def PSD(Data, Rate, Scaling='density', WindowSize=4096):
#    Window = signal.hanning(round(len(Data)/(Rate/10000)), sym=False)
    Window = signal.hanning(WindowSize, sym=False)
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
        A16OM16 = [13, 12, 14, 11, 15, 10, 16, 9, 5, 4, 6, 3, 7, 2, 8, 1]
        A16 = {'Tip': [9, 8, 10, 7, 13, 4, 12, 5, 15, 2, 16, 1, 14, 3, 11, 6],
               'Head': [8, 7, 6, 5, 4, 3, 2, 1, 9, 10, 11, 12, 13, 14, 15, 16]}
        
        ChannelMap = RemapChannels(A16['Tip'], A16['Head'], CustomAdaptor)
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


def Spectrogram(Data, Rate, HigherFreq):
    print('Analyzing data spectrum... ', end='')
    NFFT = int((1/HigherFreq)*10*Rate)
    
    Window = signal.hanning(NFFT)
    F, T, Sxx = signal.spectrogram(Data, Rate, window=Window, 
                                   nperseg=int(NFFT), noverlap=int(NFFT*0.5), 
                                   nfft=int(NFFT))
    
    print('Done.')
    return(F, T, Sxx)


def StrRange(Start='a', End='e', Step=1):
    if max(len(Start), len(End)) > 1: 
        print('Only 1-char length strings are accepted.')
        return(None)
    else:
        Range = map(chr, range(ord(Start), ord(End), Step))
        return(Range)


def UniqueStr(List, KeepOrder=False):
    if KeepOrder:
        used = set()
        UniqueList = [x for x in List if x not in used and (used.add(x) or True)]
    else:
        UniqueList = set(List); UniqueList = sorted(list(UniqueList))
    
    return(UniqueList)


## Level 1
def QuantifyTTLsPerRec(AnalogTTLs, Data=[], StdNo=3, EventsDict={}, TTLCh=None, 
                       Proc=None, Rec=None, Edge='rise'):
    print('Get TTL timestamps... ', end='')
    if AnalogTTLs:
        Threshold = GetTTLThreshold(Data, StdNo); print('TTL threshold:', Threshold)
        
        if Edge == 'rise':
            TTLs = np.where((Data[:-1] < Threshold)*(Data[1:] > Threshold))[0]
        elif Edge == 'fall':
            TTLs = np.where((Data[:-1] > Threshold)*(Data[1:] < Threshold))[0]
        else:
            print('"Edge" should be "rise" of "fall".')
            return(None)
        
        print('Done.')
        return(TTLs)
    else:
        TTLs = GetEventEdges(EventsDict, TTLCh, Proc, Rec, Edge='rise')
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
        return(TTLs)


def RemapCh(Probe, Adaptor):
    Probes = {
        'A16': {'Tip': [9, 8, 10, 7, 13, 4, 12, 5, 15, 2, 16, 1, 14, 3, 11, 6],
                'Head': [8, 7, 6, 5, 4, 3, 2, 1, 9, 10, 11, 12, 13, 14, 15, 16]}
    }
    
    Adaptors = {
        'CustomAdaptor': [5, 6, 7, 8, 9, 10 ,11, 12, 13, 14, 15, 16, 1, 2, 3, 4],
        'RHAHeadstage': [16, 15, 14, 13, 12, 11, 10, 9, 1, 2, 3, 4, 5, 6, 7, 8],
        'A16OM16': [13, 12, 14, 11, 15, 10, 16, 9, 5, 4, 6, 3, 7, 2, 8, 1]
    }
    
    if Probe not in Probes or Adaptor not in Adaptors:
        print('Unknown probe and/or adaptor.')
        print('Known probes:')
        for P in Probes.keys(): print('    ' + P)
        print('Known adaptors:')
        for A in Adaptors.keys(): print('    ' + A)
        return(None)
    
    Map = RemapChannels(Probes[Probe]['Tip'], Probes[Probe]['Head'], Adaptors[Adaptor])
    return(Map)


def SignalIntensity(Data, Rate, FreqBand, Ref, NoiseRMS=None, WindowSize=None):
    Intensity = {}
    
    if not WindowSize: WindowSize = len(Data)
    
    F, PxxSp = PSD(Data, Rate, WindowSize=WindowSize)
    
    Range = (F > FreqBand[0])*(F < FreqBand[1])
    BinSize = F[1] - F[0]
    
    RMS = (sum(PxxSp[Range]) * BinSize)**0.5
    if NoiseRMS: RMS = RMS - NoiseRMS
    
    dB = 20*(np.log10((RMS/Ref)/0.00002))
    
    IntensityPSD = [F, PxxSp]
    Intensity['RMS'] = RMS
    Intensity['dB'] = dB
    
    return(Intensity, IntensityPSD)


def SliceData(Data, TTLs, NoOfSamplesBefore, 
              NoOfSamplesAfter, NoOfSamples, AnalogTTLs, RawTime=[]):
    print('Slicing data around TTL...')
#    Array = [[0 for _ in range(NoOfSamples)] for _ in range(len(TTLs))]
    Array = np.zeros((len(TTLs), NoOfSamples))
#    TTLsToFix = []
    
    for TTL in range(len(TTLs)):
        if AnalogTTLs: TTLLoc = int(TTLs[TTL])
        else: TTLLoc = int(RawTime.index(TTLs[TTL]))#)/Rate)
        
        Start = TTLLoc-NoOfSamplesBefore
        End = TTLLoc+NoOfSamplesAfter
        
#        if Start < 0: Start = 0; End = End+(Start*-1); TTLsToFix.append(TTL)
        if Start < 0 or End > len(Data):
            print('TTL too close to the edge. Skipping...')
            continue
        
        Array[TTL,:] = Data[Start:End]
        
#        if len(Array[TTL]) != End-Start: TTLsToFix.append(TTL)
            
#    if TTLsToFix: Array = FixTTLs(Array, TTLsToFix)
    
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
#def ABRCalc(Data, Rate, ExpInfo, DataInfo, Stim, TimeBeforeTTL=3, 
#            TimeAfterTTL=12, AnalogTTLs=True, ABRCh=5, TTLCh=17, 
#            FilterFreq=[300, 3000], FilterOrder=5, TTLsPerRec=[]):
#    ABRs = {}; Info = {}
#    
#    NoOfSamplesBefore = TimeBeforeTTL*int(Rate*10**-3)
#    NoOfSamplesAfter = TimeAfterTTL*int(Rate*10**-3)
#    NoOfSamples = NoOfSamplesBefore + NoOfSamplesAfter
#    
#    Info['Frequency'] = ExpInfo['Hz']
#    
#    Info['XValues'] = (range(-NoOfSamplesBefore, 
#                             NoOfSamples-NoOfSamplesBefore)/Rate)*10**3
#    
#    for Rec in Data.keys():
#        print('Slicing and filtering ABRs Rec ', str(Rec), '...')
#        
#        if len(Data[Rec]) < 50*Rate:
#            print('Rec', Rec, 'is broken!!!')
#            continue
#        
#        if AnalogTTLs:
#            TTLs = QuantifyTTLsPerRec(AnalogTTLs, Data[Rec][:,TTLCh-1])
#            ABR = SliceData(Data[Rec][:, ABRCh-1], TTLs, NoOfSamplesBefore, 
#                            NoOfSamplesAfter, NoOfSamples, AnalogTTLs)
##        else:
##            RawTime, TTLs = QuantifyTTLsPerRec(Data, Rec, AnalogTTLs, 
##                                               TTLsPerRec=TTLsPerRec)
##            ABR = SliceData(Data[Rec][:,ABRCh-1], TTLs,NoOfSamplesBefore, 
##                            NoOfSamplesAfter, NoOfSamples, AnalogTTLs, RawTime)
#        
#        for TTL in range(len(TTLs)):
#            ABR[TTL] = FilterSignal(ABR[TTL], Rate, [min(FilterFreq)], 
#                                    FilterOrder, 'butter', 'highpass')
#        
#        ABR = np.mean(ABR, axis=0)
#        ABR = FilterSignal(ABR, Rate, [max(FilterFreq)], FilterOrder, 
#                           'butter', 'lowpass')
#        
#        dB = str(DataInfo['Intensities'][int(Rec)]) + 'dB'
#        ABRs[dB] = ABR[:]; del(ABR)
#    
#    Info['Path'] = Stim+'/'+ExpInfo['DV']+'/'+Info['Frequency']
#    
#    return(ABRs, Info)

