# -*- coding: utf-8 -*-
"""
Just drafts
"""
#%% GPIASMatToHdf5


GPIASTimeBeforeTTL = 20    # in ms
GPIASTimeAfterTTL = 100    # in ms
FilterFreq = [70, 400]     # frequency for filter
FilterOrder = 3       # butter order
RecFolder = 1

#%% ReAnalysis
import Hdf5F
import KwikAnalysis
import MatF
from glob import glob

Animals = ['CaMKIIahM4Dn06', 'CaMKIIahM4Dn07', 'CaMKIIahM4Dn08', 'CaMKIIahM4Dn09']
#AlreadRun = []
AlreadyRun = ['CaMKIIahM4Dn06/CaMKIIahM4Dn06-20160217-GPIASn01Screening',
              'CaMKIIahM4Dn06/CaMKIIahM4Dn06-20160623-GPIAS',
              'CaMKIIahM4Dn06/CaMKIIahM4Dn06-20160629-GPIAS',
              'CaMKIIahM4Dn06/CaMKIIahM4Dn06-20160630-GPIAS',
              'CaMKIIahM4Dn07/CaMKIIahM4Dn07-20160217-GPIASn01Screening',
              'CaMKIIahM4Dn07/CaMKIIahM4Dn07-20160218-GPIASn01Retest'
              'CaMKIIahM4Dn08/CaMKIIahM4Dn08-20160218-GPIASn02Screening',
              'CaMKIIahM4Dn08/CaMKIIahM4Dn08-20160624-GPIAS',
              'CaMKIIahM4Dn08/CaMKIIahM4Dn08-20160629-GPIAS',
              'CaMKIIahM4Dn08/CaMKIIahM4Dn08-20160630-GPIAS',
              'CaMKIIahM4Dn08/CaMKIIahM4Dn08-20160701-GPIAS',
              'CaMKIIahM4Dn08/CaMKIIahM4Dn08-20160702-GPIAS'
              'CaMKIIahM4Dn09/CaMKIIahM4Dn09-20160218-GPIASn02Screening',
              'CaMKIIahM4Dn09/CaMKIIahM4Dn09-20160624-GPIAS']
Exp = 'GPIAS'

GPIASCh = [1]
GPIASTTLCh = 2
GPIASTimeBeforeTTL = 20    # in ms
GPIASTimeAfterTTL = 100    # in ms
FilterFreq = [70, 400]     # frequency for filter
FilterOrder = 3       # butter order
AnalogTTLs = True
RecFolderNo = 1
Override = {}

for Animal in Animals:
    Override['AnalysisFile'] = glob(Animal+'/*.hdf5')[0]
    Exps = [F for F in glob(Animal+'/*') if Exp in F]; Exps.sort()
    
    for RecExp in Exps:
        if RecExp in AlreadyRun:
            continue
        
        print('Running', RecExp + '...')
        DataFolders = glob(RecExp + '/*Files')
        
        for DataFolder in DataFolders:
            DataType = DataFolder.split('/')[-1]
            
            if DataType == 'KwikFiles':
                Override['FileName'] = glob(RecExp + '/*.hdf5')
                Override['FileName'].sort(); 
                Override['FileName'] = Override['FileName'][RecFolderNo-1]
                
                Override['DirList'] = glob(RecExp + '/KwikFiles/*')
                Override['DirList'].sort()
                
                KwikAnalysis.GPIASAnalysis(RecFolderNo, GPIASCh, GPIASTTLCh, 
                                           GPIASTimeBeforeTTL, 
                                           GPIASTimeAfterTTL, FilterFreq, 
                                           FilterOrder, AnalogTTLs, 
                                           Override=Override)
                
            elif DataType == 'MatFiles':
                Override['FileName'] = glob(RecExp + '/*.mat')[0]
                
                Override['DirList'] = glob(RecExp + '/MatFiles/*')
                Override['DirList'].sort()
                
                MatF.GPIASAnalysisMat(RecFolderNo, GPIASTimeBeforeTTL, 
                                      GPIASTimeAfterTTL, FilterFreq, 
                                      FilterOrder, Override)
                    
            else: 
                print(DataType, 'not supported (yet).')
                continue
        
        AlreadyRun.append(RecExp)
                
        
#%% GPIASIndex
Params = KwikAnalysis.SetPlot(Backend='TkAgg', Params=True)
from matplotlib import rcParams; rcParams.update(Params)
from matplotlib import pyplot as plt

Animals = ['CaMKIIahM4Dn06', 'CaMKIIahM4Dn07', 'CaMKIIahM4Dn08', 'CaMKIIahM4Dn09']
GPIASIndex = {}

for Animal in Animals:
    AnalysisFile = glob(Animal+'/*.hdf5')[0]
    GPIASIndex[Animal] = {}
    
    with h5py,File(AnalysisFile, 'r') as F:
         Keys = [Key for Key in F.keys() if 'GPIAS_' in Key]; Keys.sort()
         
         for Key in Keys:
             
             GPIASIndex[Animal]
             
    GPIAS, XValues = Hdf5F.LoadGPIAS(AnalysisFile)
    
    for Freq in GPIAS:
        GPIASIndex[Animal]
        


Group = {}

YLabel = 'GPIAS index'
ExpLabels = ['Before ANT', '2d after ANT', '1w after ANT (NaCl)', 
             '1w after ANT (CNO)']

Fig, Axes = plt.subplots(1,4, sharey=True)

for AInd, Animal in enumerate(Animals):
    AnalysisFile = glob(Animal+'/*.hdf5')[0]
    GPIAS, XValues = Hdf5F.LoadGPIAS(AnalysisFile)
    
    Freqs = list(GPÌAS.keys())
    for FInd, Freq in enumerate(Freqs):
        
        
    
#%% ABRThresholds
import Hdf5F
import numpy as np
from glob import glob

FileName = glob('*.txt')[0]

with open(FileName, 'r') as F:
    Lines = [Line for Line in F]

Lines = [Line.split() for Line in Lines]
Exps = Lines[0][:]; del(Lines[0])
Freqs = Lines[0][:]; del(Lines[0])

ABRThresholds = {Exps[_]: np.array(Lines[_], dtype='f') 
                 for _ in range(len(Lines))}

AnalysisFile = glob('*.hdf5')[0]
Hdf5F.WriteDict(ABRThresholds, '/ABRThresholds', AnalysisFile, Attrs=False)

del(ABRThresholds, AnalysisFile, Exps, FileName, Freqs, F, Lines)


#%% Plot ABRThresholds
import KwikAnalysis
import Hdf5F
import numpy as np
from glob import glob

Animals = ['CaMKIIahM4Dn06', 'CaMKIIahM4Dn07', 'CaMKIIahM4Dn08', 'CaMKIIahM4Dn09']

Params = KwikAnalysis.SetPlot(Backend='TkAgg', Params=True)
from matplotlib import rcParams; rcParams.update(Params)
from matplotlib import pyplot as plt

YLabel = 'Intensity [dBSPL]'
Freqs = ['8$\sim$10', '9$\sim$11', '10$\sim$12', '12$\sim$14', '14$\sim$16 [KHz]']
Colors = ['--ro', '--bx', '--g^', '--ms']
When = ['Before ANT', '48h after ANT', '2w after ANT (NaCl)', '2w after ANT (CNO)']

Fig, Axes = plt.subplots(1,4, sharey=True)
for AInd, Animal in enumerate(Animals):
    AnalysisFile = glob(Animal + '/*.hdf5')[0]
    ABRThresholds = Hdf5F.LoadDict('/ABRThresholds', AnalysisFile, Attrs=False)
    Exps = [Exp for Exp in ABRThresholds.keys()]; Exps.sort()
    
    for EInd, Exp in enumerate(Exps):
        Axes[EInd].plot(ABRThresholds[Exp], Colors[AInd], label=Animal)
        KwikAnalysis.SetPlot(AxesObj=Axes[EInd], Axes=True)
        Axes[EInd].xaxis.set_ticks(range(len(ABRThresholds[Exp])))
        Axes[EInd].set_xticklabels(Freqs)
        Axes[EInd].set_xlabel(When[EInd])
        Axes[EInd].legend(loc='upper left')
        Axes[EInd].spines['bottom'].set_position(('outward', 5))
        if EInd == 0: Axes[EInd].spines['left'].set_position(('outward', 5))
        else: 
            Axes[EInd].spines['left'].set_visible(False)
            Axes[EInd].yaxis.set_ticks_position('none')
#        
#    adjust_spines(Axes[0], ['left, bottom'])
    Axes[0].set_ylabel(YLabel)
plt.show()

CharNo = [len(Animal) for Animal in Animals]

PrefixList = []; Prefix = ''; SuffixList = []; SuffixA = ''; SuffixB = ''
if len(np.unique(CharNo)) == 1:
    for Animal in range(1, len(Animals)):
        for Char in range(CharNo[0]):
            if Animals[Animal-1][Char] == Animals[Animal][Char]:
                Prefix = ''.join([Prefix, Animals[Animal][Char]])
            else: 
                SuffixA = ''.join([SuffixA, Animals[Animal-1][Char:]])
                SuffixB = ''.join([SuffixB, Animals[Animal][Char:]])
                break
        PrefixList.append(Prefix[:]); Prefix = ''
        SuffixList.append([SuffixA[:], SuffixB[:]]); SuffixA = ''; SuffixB = ''

Prefix = np.unique(PrefixList)[0]
FigName = Prefix + Suffix[:-1]

#%% Tar+Wave revolution

ABRCh = [10]         # [RightChannel, LeftChannel], if order matters
ABRTTLCh = 17            # TTL ch for ABR
ABRTimeBeforeTTL = 3    # in ms
ABRTimeAfterTTL = 9    # in ms
FilterFreq = [300, 3000]         # frequeutter order
AnalogTTLs = True
FilterOrder = 5         # b
StimType = ['Sound_NaCl', 'Sound_CNO']
#StimType = ['Sound']
Board = 'OE'
Override = {}

import Hdf5F
import numpy as np
import os
import struct
import tarfile
import wave
from datetime import datetime
from glob import glob
from scipy import signal

FileName = glob('*.hdf5')[0]

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


def WriteABRTar(ABRs, XValues, Path, FileName):
    print('Writing data to', FileName+'... ', end='')    
    FileList = glob(Path + '/**/*.*', recursive=True); FileList.sort()
    
    with tarfile.open(FileName, 'a') as F:
        for File in FileList:
            F.add(File)
    
    print('Done.')
    return(None)


def WriteWav(Data, ChNo, DataSize, Rate, FileName):
    print('Writing data to', FileName+'... ', end='')
    with wave.open(Path+'/'+FileName, 'w') as F:
        F.setparams((ChNo, 4, Rate, DataSize*ChNo, 'NONE', 'uncompressed'))
        F.writeframes(Data)
    
    print('Done.')
    return(None)

def WriteABRWave(ABRs, XValues, Rate, Group, Path):
    Path = Group + '/' + Path
    Trial = len(glob(Path+'/*')); Trial = "{0:02d}".format(Trial)
    Path = Path + '/' + Trial
    os.makedirs(Path, exist_ok=True)
    
    Keys = list(ABRs.keys()); Keys.sort()
    ChNo = len(ABRs); DataSize = len(ABRs[Keys[0]])
    
    Data = [ABRs[dB][Sample] for Sample in range(DataSize) for dB in Keys]
    Data = np.array('f', Data); Data = bytes(Data)
    
    FileName = Path.split(sep='/'); del(FileName[0])
    FileName = 'ABRs-' + '_'.join(FileName)
    WriteWav(Data, ChNo, DataSize, Rate, FileName)
    
    XValues = np.array('f', XValues); XValues = bytes(XValues)
    
    FileName = Path.split(sep='/'); del(FileName[0])
    FileName = 'XValues-' + '_'.join(FileName)
    WriteWav(XValues, 1, len(XValues), Rate, FileName)


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



#%%
Backend = 'Qt5Agg'
from matplotlib import rcParams; rcParams.update({'backend': Backend})
from matplotlib import pyplot as plt
from random import random
from time import time

a = list(range(7500)); b = [random() for _ in range(7500)]

Time = []
for aa in range(5):
    Start = time(); plt.bar(a,b); Time.append(time() - Start)
    print(str(aa), '=', str(Time[aa]))

Time = sum(Time)/5
print(Backend, '=', str(Time))

Tk = 8.832437372207641
Qt4 =  27.928475093841552


#%%
import h5py
with h5py.File(FileName) as F:
    for key in F['ExpInfo'].keys():
        N = "{0:02d}".format(int(key))
        F['ExpInfo'][N] = F['ExpInfo'][key]
        del(F['ExpInfo'][key])

for ind, key in enumerate(F['ExpInfo'].keys()):
    if ind <5: F['ExpInfo'][key].attrs['StimType'] = np.string_(['Sound_NaCl'])
    else: F['ExpInfo'][key].attrs['StimType'] = np.string_(['Sound_CNO'])



#%%
#import h5py
#from glob import glob

Files = glob('*/*.hdf5'); Files.sort()
AnalysisFile = Files[0][:14] + '-Analysis.hdf5'

Analysis = h5py.File(AnalysisFile)

for Ind, File in enumerate(Files):
    F = h5py.File(File)
    for Key in F.keys():
        for Exp in ['ABRs-', 'GPIAS']:
            if Exp in Key:
                F.copy('/'+Key, Analysis)
                del(F[Key])
    F.close()

print(list(Analysis.keys()))
Analysis.close()


#%%
def GPIASAAA(FileName, GPIASTimeBeforeTTL=50, GPIASTimeAfterTTL=150, FilterLow=3, 
          FilterHigh=300, FilterOrder=4, GPIASTTLCh=2, PiezoCh=1):
    
    print('set paths...')
    os.makedirs('Figs', exist_ok=True)    # Figs folder
    DirList = glob.glob('KwikFiles/*'); DirList.sort()
    
    for RecFolder in DirList:
        print('Load files...')
        FilesList = glob.glob(''.join([RecFolder, '/*']))
        Files = {}
        for File in FilesList:
            if '.kwd' in File:
                try:
                    Raw = Kwik.load(File, 'all')
                    Files['kwd'] = File
                except OSError:
                        print('File ', File, " is corrupted :'(")
                
            elif '.kwe' in File:
                try:
                    Events = Kwik.load(File)
                    Files['kwe'] = File
                except OSError:
                    print('File ', File, " is corrupted :'(")
                
            elif '.kwik' in File:
                try:
                    Events = Kwik.load(File)
                    Files['kwik'] = File
                except OSError:
                    print('File ', File, " is corrupted :'(")
            
#            elif '.db' in File:
#                with shelve.open(File[:-3]) as Shelve: 
#                    DataInfo = Shelve['DataInfo']
        
        DataInfo = Hdf5F.LoadDict('/DataInfo', FileName)
        DataInfo['SoundBackgroundAmpF'] = Hdf5F.LoadDict(
                                             '/DataInfo/SoundBackgroundAmpF', 
                                             FileName, Attrs=False)
        DataInfo['SoundPulseAmpF'] = Hdf5F.LoadDict('/DataInfo/SoundPulseAmpF', 
                                                    FileName, Attrs=False)
        
        print('Check if files are ok...')
        if 'Raw' not in locals():
            print('.kwd file is corrupted. Skipping dataset...')
            continue
        
        if 'Events' not in locals():
            print('.kwe/.kwik file is corrupted. Skipping dataset...')
            continue
        
        if 'DataInfo' not in locals():
            print('No data info. Skipping dataset...')
            continue
        
        print('Data from ', RecFolder, ' loaded.')
        
        print('Preallocate memory...')
        GPIAS = [[0] for _ in range(len(DataInfo['NoiseFrequency']))]
        for Freq in range(len(DataInfo['NoiseFrequency'])):
            GPIAS[Freq] = [[0] for _ in range(round(DataInfo['NoOfTrials']*2))]
        
#        AllTTLs = GPIAS.copy()
        
        print('Get TTL data...')
        EventID = Events['TTLs']['user_data']['eventID']
        EventCh = Events['TTLs']['user_data']['event_channels']
        EventRec = Events['TTLs']['recording']
        EventSample = Events['TTLs']['time_samples']
        
        TTLChs = np.nonzero(np.bincount(EventCh))[0]
        TTLRecs = np.nonzero(np.bincount(EventRec))[0]
        TTLsPerRec = {_Rec: [EventSample[_] for _ in range(len(EventRec)) 
                             if EventRec[_] == _Rec 
                             and EventCh[_] ==  GPIASTTLCh-1 
                             and EventID[_] == 1]
                      for _Rec in TTLRecs}
        TTLsDPerRec = {_Rec: [EventSample[_] for _ in range(len(EventRec)) 
                             if EventRec[_] == _Rec 
                             and EventCh[_] ==  GPIASTTLCh-1 
                             and EventID[_] == 0]
                      for _Rec in TTLRecs}
        TTLRising = Kwik.get_rising_edge_times(Files['kwe'], GPIASTTLCh-1)
        
        Rate = Raw['info']['0']['sample_rate']
        NoOfSamplesBefore = int(round((GPIASTimeBeforeTTL*Rate)*10**-3))
        NoOfSamplesAfter = int(round((GPIASTimeAfterTTL*Rate)*10**-3))
        NoOfSamples = NoOfSamplesBefore + NoOfSamplesAfter
        XValues = (range(-NoOfSamplesBefore, 
                         NoOfSamples-NoOfSamplesBefore)/Rate)*10**3
        
        print('Set filter...')
        passband = [FilterLow/(Rate/2), FilterHigh/(Rate/2)]
        f2, f1 = signal.butter(FilterOrder, passband, 'bandpass')
        
        for Rec in range(len(Raw['data'])):
            TTLNo = [0]
            for _ in range(1, len(TTLsPerRec)+1):
                TTLNo = TTLNo + [len(TTLsPerRec[_-1]) + TTLNo[-1]]

            TTLNo = [0] + [TTLNo[_]-1 for _ in range(1, len(TTLNo))]
            
            RawTime = list(Raw['timestamps'][str(Rec)])
            
            if Rec == 0:
                sTTLNo = 0
            else:
                sTTLNo = TTLNo[Rec] + 1
            
            print('Slicing and filtering Rec ', str(Rec), '...')
            for TTL in range(sTTLNo, TTLNo[Rec+1]):
                TTLLoc = RawTime.index(TTLRising[TTL])
                    
                Start = TTLLoc-NoOfSamplesBefore
                End = TTLLoc+NoOfSamplesAfter
                
                Index = TTL-sTTLNo
                
                gData = Raw['data'][str(Rec)][Start:End, PiezoCh-1]
#                gData = [float(gData[_]) for _ in range(NoOfSamples)]
                gData = abs(signal.hilbert(gData))
                gData = signal.filtfilt(f2, f1, gData, padtype='odd', padlen=0)
                
                Freq = DataInfo['FreqOrder'][Rec][0]; 
                Trial = DataInfo['FreqOrder'][Rec][1];
                GPIAS[Freq][Trial] = gData[:]
#                AllTTLs[Freq][Trial] = [TTLStart, TTLEnd]
                
                del(TTLLoc, Start, End, Freq, Trial, gData)
        
        for Freq in range(len(DataInfo['NoiseFrequency'])):
            gData = GPIAS[Freq][:]
            NoGapAll = [gData[_] for _ in range(len(gData)) if _%2 == 0]
            GapAll = [gData[_] for _ in range(len(gData)) if _%2 != 0]
            NoGapSum = list(map(sum, zip(*NoGapAll)))
            GapSum = list(map(sum, zip(*GapAll)))
            
            gData = [0, 0]
            gData[0] = [_/DataInfo['NoOfTrials'] for _ in NoGapSum]
            gData[1] = [_/DataInfo['NoOfTrials'] for _ in GapSum]
            gData[0] = signal.savgol_filter(gData[0], 5, 2, mode='nearest')
            gData[1] = signal.savgol_filter(gData[1], 5, 2, mode='nearest')
            GPIAS[Freq] = gData[:]
            
#            tData = AllTTLs[Freq][:]
#            TTLNoGapAll = [tData[_] for _ in range(len(tData)) if _%2 == 0]
#            TTLGapAll = [tData[_] for _ in range(len(tData)) if _%2 != 0]
#            TTLNoGapSum = list(map(sum, zip(*TTLNoGapAll)))
#            TTLGapSum = list(map(sum, zip(*TTLGapAll)))
#            
#            tData = [0, 0]
#            tData[0] = [int(round(_/DataInfo['NoOfTrials'])) for _ in TTLNoGapSum]
#            tData[1] = [int(round(_/DataInfo['NoOfTrials'])) for _ in TTLGapSum]
#            AllTTLs[Freq] = tData[:]
            
            del(NoGapAll, GapAll, NoGapSum, GapSum, gData)
#            del(TTLNoGapAll, TTLGapAll, TTLNoGapSum, TTLGapSum, tData)
            
        FileName = DataInfo['AnimalName'] + '-GPIAS-' + RecFolder[10:]
        print('Saving data to ' + FileName)
#        with shelve.open(FileName) as Shelve:
#            Shelve['GPIAS'] = GPIAS
#            Shelve['AllTTLs'] = AllTTLs
#            Shelve['XValues'] = XValues
#            Shelve['RecFolder'] = RecFolder
#            Shelve['DataInfo'] = DataInfo
        print('Done.')
    
    print('Finished.')
    return(None)

#%% SavGol
        # SavGol smooth
#        GPIAS[Freq]['NoGap'] = signal.savgol_filter(GPIAS[Freq]['NoGap'], 5, 2, 
#                                                    mode='nearest')
#        GPIAS[Freq]['Gap'] = signal.savgol_filter(GPIAS[Freq]['Gap'], 5, 2, 
#                                                  mode='nearest')
        

#%% RenameFiles
import os
from glob import glob

ModifyFreq = ['20160217-GPIASn01Screening']

for Folder in ModifyFreq:
    Files = glob(Folder + '/*.mat')
    
    for File in Files:
        NewFile = File.replace('  ', '_')
        os.rename(File, NewFile)


#%%
import os
from glob import glob

ModifyName = glob('*GPIAS*'); ModifyName.sort()

for Folder in ModifyName:
    Files = glob(Folder + '/*.mat')
    Start = len(Folder) + 1
        
    for File in Files:
        NewFile = File[Start:Start+14] + '_' + File[Start+15:]
        NewFile = NewFile.split('_')
        Date = NewFile[0]; Animal = NewFile[1]
        DateIndex = NewFile.index(Date); AnimalIndex = NewFile.index(Animal)
        NewFile[0], NewFile[1] = NewFile[AnimalIndex], NewFile[DateIndex]
        NewFile = File[:Start] + '_'.join(NewFile)
        os.rename(File, NewFile)
#        print(File + '\n', NewFile)

#%%
import os
from glob import glob

ModifyFreq = ['20160217-GPIASn01Screening']

for Folder in ModifyFreq:
    Files = glob(Folder + '/*.mat')
    Start = len(Folder) + 1
    
    for File in Files:
        if File[Start:Start+4] == '1OD-': NewFile = File[:Start] + '1R_' + File[Start+4:]
        elif File[Start:Start+4] == '1OE-': NewFile = File[:Start] + '1L_' + File[Start+4:]
        elif File[Start:Start+4] == '2OD-': NewFile = File[:Start] + '2R_' + File[Start+4:]
        elif File[Start:Start+4] == '2OE-': NewFile = File[:Start] + '2L_' + File[Start+4:]
        elif File[Start:Start+5] == '1ODE-': NewFile = File[:Start] + '1LR_' + File[Start+5:]
        
        os.rename(File, NewFile)
