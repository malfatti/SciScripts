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

import datetime
import glob
import h5py
import Kwik
import Hdf5F
from matplotlib import rcParams
from matplotlib import pyplot as plt
from numbers import Number
import numpy as np
import os
from scipy import signal


## Lower-level functions
def GetProc(Raw, Board):
    ProcChs = {Proc: len(Raw[Proc]['data']['0'][1,:]) 
               for Proc in Raw.keys()}
    
    for Proc, Chs in ProcChs.items():
        if Chs == max(ProcChs.values()): OEProc = Proc
        else: RHAProc = Proc
    
    if 'RHAProc' not in locals(): RHAProc = OEProc
    
    if Board == 'OE': Proc = OEProc
    elif Board == 'RHA': Proc = RHAProc
    else: print("Choose Board as 'OE' or 'RHA'."); return(None)
    
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
        
        ChannelMap = GetChOrder(A16['ProbeTip'], A16['ProbeHead'], CustomAdaptor)
    """
    ChNo = len(ProbeTip)
    ChMap = [0]*ChNo
    
    for Ch in range(ChNo):
        TipCh = ProbeTip[Ch] # What channel should be the Ch
        HeadCh = ProbeHead.index(TipCh) # Where Ch is in ProbeHead
        ChMap[Ch] = Connector[HeadCh] # Channels in depth order
    
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


def QuantifyTTLsPerRec(Raw, Rec, AnalogTTLs, ChTTL=-1, Proc='', TTLsPerRec=[]):
    if AnalogTTLs:
        TTLCh = Raw[Proc]['data'][str(Rec)][:, ChTTL-1]
        Threshold = max(TTLCh)/2
        TTLs = []
        for _ in range(1, len(TTLCh)):
            if TTLCh[_] > Threshold:
                if TTLCh[_-1] < Threshold: TTLs.append(_)
        
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
        return(RawTime, TTLs)


def RemoveDateFromFolderName():
    RenameFolders = input('Rename folders in KwikFiles/* (BE CAREFUL)? [y/N] ')
    if RenameFolders in ['y', 'Y', 'yes', 'Yes', 'YES']:
        DirList = glob.glob('KwikFiles/*'); DirList.sort()
        for FolderName in DirList:
            NewFolderName = ''.join([FolderName[:10], FolderName[21:]])
            NewFolderName = NewFolderName.replace("-", "")
            os.rename(FolderName, NewFolderName)
            print(FolderName, ' moved to ', NewFolderName)
        del(RenameFolders, DirList, FolderName, NewFolderName)    
    
    return(None)


def SetPlot(AxesObj=(), FigObj=(), FigTitle='', Params=False, Plot=False, 
            Axes=False):
    if Params:
        print('Set plot parameters...')
        Params = {'backend': 'TkAgg',
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
        rcParams.update(Params)
    
    elif Plot:
        FigObj.suptitle(FigTitle); FigObj.tight_layout(); 
        FigObj.subplots_adjust(top=0.925)
    
    elif Axes:
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
              NoOfSamplesAfter, AnalogTTLs, RawTime=[]):
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
        
    return(Array)


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
    DirList = glob.glob('KwikFiles/*'); DirList.sort()
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
        
        if AnalogTTLs: Raw, _, Files = Hdf5F.LoadOEFiles(RecFolder, AnalogTTLs)
        else: Raw, Events, _, Files = Hdf5F.LoadOEFiles(RecFolder, AnalogTTLs)
        
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
                                AnalogTTLs)
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
    
    Hdf5F.WriteABRs(ABRs, XValues, FileName)
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
    
    SetPlot(Params=True)
    
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
    DirList = glob.glob('KwikFiles/*'); DirList.sort()
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
    
    if AnalogTTLs: Raw, _, Files = Hdf5F.LoadOEFiles(RecFolder, AnalogTTLs)
    else: Raw, Events, _, Files = Hdf5F.LoadOEFiles(RecFolder, AnalogTTLs)
    
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
                                           AnalogTTLs)
        else:
            RawTime, TTLs = QuantifyTTLsPerRec(Raw, Rec, AnalogTTLs, 
                                               TTLsPerRec=TTLsPerRec)
            GPIAS[Freq][Trial] = SliceData(Raw, OEProc, Rec, TTLs, GPIASCh, 
                                           NoOfSamplesBefore, NoOfSamplesAfter, 
                                           AnalogTTLs, RawTime)
        
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
    
    Hdf5F.WriteGPIAS(GPIAS, RecFolder, XValues, FileName)
    return(None)


def GPIASPlot(FileName):
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
    
    SetPlot(Params=True)
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


def TTLsLatency(FileName, SoundCh=1, TTLSqCh=2, TTLCh=1, 
                TimeBeforeTTL=5, TimeAfterTTL=8, AnalogTTLs=False):
    print('set paths...')
    os.makedirs('Figs', exist_ok=True)    # Figs folder
    RecFolder = glob.glob('KwikFiles/*'); RecFolder = RecFolder[0]
#    SoundCh = 0; TTLSqCh = 1 # Override
    
    TTLsLatency = {}    
    Raw, Events, _, Files = Hdf5F.LoadOEFiles(RecFolder, AnalogTTLs)    
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
                               AnalogTTLs, RawTime)
        TTLSq = SliceData(TTLSq, Raw, OEProc, Rec, TTLs,TTLSqCh, 
                          NoOfSamplesBefore, NoOfSamplesAfter, 
                          AnalogTTLs, RawTime)
    
        TTLSqDelay = [[0] for _ in range(len(TTLSq))]    
        for TTL in range(len(TTLSq)):
            Peaks = QuantifyTTLsPerRec(Raw, Rec, AnalogTTLs=True, TTLSqCh, 
                                       OEProc)
            
            if len(Peaks) != 1: 
                print('Bad TTL detection. Skipping TTL...')
                TTLSqDelay[TTL] = float('NaN')
                continue
            
            TTLSqDelay[TTL] = ((Peaks[0] - NoOfSamplesBefore) / (Rate/1000))*-1
            del(Peaks)
        
        TTLsLatency[str(Rec)] = dict((Name, eval(Name)) 
                                     for Name in ['SoundPulse', 'TTLSq', 
                                                  'TTLSqSoundDelay'])
    
    Hdf5F.WriteTTLsLatency(TTLsLatency, XValues, FileName)
    
    return(None)


def PlotTTLsLatency(FileName):
    DataInfo = {}
    with h5py.File(FileName) as F:
        SoundPulse = F['Analysis']['SoundPulse'][:]
        SoundSq = F['Analysis']['SoundSq'][:]
        SoundSqDelay = F['Analysis']['SoundSqDelay'][:]
        XValues = list(F['Analysis']['XValues'])
    
        for Key, Value in F['DataInfo'].items():
            DataInfo['SoundAmpF'] = {}
            for aKey, aValue in F['DataInfo']['SoundAmpF'].items():
                DataInfo['SoundAmpF'][aKey] = aValue[:]
        
        for bKey, bValue in F['DataInfo'].attrs.items():
            if isinstance(bValue, Number):
                DataInfo[bKey] = float(bValue)
            else:
                DataInfo[bKey] = bValue
    
    SetPlot(Params=True)
    
    for _ in range(len(SoundPulse)):
        plt.figure(1); plt.plot(XValues, SoundPulse[_])
        plt.figure(2); plt.plot(XValues, SoundSq[_])
    
    Hist, BinEdges = np.histogram(SoundSqDelay, bins=200)
    Threshold = (DataInfo['SoundPulseDur']/100)*1000
    Threshold = 0.08
    sIndex = min(range(len(BinEdges)), 
                 key=lambda i: abs(BinEdges[i]-Threshold*-1))
    eIndex = min(range(len(BinEdges)), 
                 key=lambda i: abs(BinEdges[i]-Threshold))
    Sum = sum(Hist); Perc = Sum/len(SoundSq) * 100
    plt.figure(3); plt.plot(BinEdges[:-1], Hist)
    plt.axvspan(BinEdges[sIndex], BinEdges[eIndex], color='k', alpha=0.5, lw=0,
                label=str(Perc) + '\% of pulses with latency $<$ 3µs') 
    plt.legend(loc='upper right')
    plt.locator_params(tight=True)
    plt.axes().spines['right'].set_visible(False)
    plt.axes().spines['top'].set_visible(False)
    plt.axes().yaxis.set_ticks_position('left')
    plt.axes().xaxis.set_ticks_position('bottom')
    plt.savefig('Figs/SoundTTLLatencies-SoundBoardToOE.svg', format='svg')
