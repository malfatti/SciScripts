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

import Hdf5F
import MatF
import os
from glob import glob


def RemoveDateFromFolderName(Type):
    RenameFolders = input('Rename folders in' + Type + '/* (BE CAREFUL)? [y/N] ')
    if RenameFolders in ['y', 'Y', 'yes', 'Yes', 'YES']:
        DirList = glob(Type+'/*'); DirList.sort()
        for FolderName in DirList:
            NewFolderName = ''.join([FolderName[:10], FolderName[21:]])
            NewFolderName = NewFolderName.replace("-", "")
            os.rename(FolderName, NewFolderName)
            print(FolderName, ' moved to ', NewFolderName)
        del(RenameFolders, DirList, FolderName, NewFolderName)    
    
    return(None)


def GPIASAnalysisGroup(RecFolderNo, GPIASCh, GPIASTTLCh, GPIASTimeBeforeTTL, 
                       GPIASTimeAfterTTL, FilterFreq, FilterOrder, AnalogTTLs, 
                       Animals, Exp, AlreadyRun=[], Override={}, Visible=False):
    for Animal in Animals:
        Override['AnalysisFile'] = glob(Animal+'/*.hdf5')[0]
        Exps = [F for F in glob(Animal+'/*') if Exp in F]; Exps.sort()
        
        for RecExp in Exps:
            if RecExp in AlreadyRun:
                continue
            
            Override['ExpPath'] = RecExp
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
                    
                    GPIASAnalysis(RecFolderNo, GPIASCh, GPIASTTLCh, 
                                  GPIASTimeBeforeTTL, GPIASTimeAfterTTL, 
                                  FilterFreq, FilterOrder, AnalogTTLs, 
                                  Override=Override)
                    
                    GPIASPlot(RecFolderNo, Visible=Visible, Override=Override)
                    
                elif DataType == 'MatFiles':
                    Override['FileName'] = glob(RecExp + '/*.mat')[0]
                    
                    Override['DirList'] = glob(RecExp + '/MatFiles/*')
                    Override['DirList'].sort()
                    
                    MatF.GPIASAnalysisMat(RecFolderNo, GPIASTimeBeforeTTL, 
                                          GPIASTimeAfterTTL, FilterFreq, 
                                          FilterOrder, Override)
                    
                    GPIASPlot(RecFolderNo, Visible=Visible, Override=Override)
                    
                else: 
                    print(DataType, 'not supported (yet).')
                    continue
            
            AlreadyRun.append(RecExp)
    
    return(AlreadyRun)


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

