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
#%% Set experiment details

PiezoCh = 1
GPIASTTLCh = 1
GPIASTimeBeforeTTL = 50    # in ms
GPIASTimeAfterTTL = 150    # in ms
FilterLow = 3       # High-pass frequency for bandpass filter
FilterHigh = 300     # Low-pass frequency
FilterOrder = 3       # butter order

import KwikAnalysis
import glob

FileList = glob.glob('*.db'); FileList.sort()

KwikAnalysis.GPIAS(GPIASTimeBeforeTTL, GPIASTimeAfterTTL, FilterLow, 
                   FilterHigh, FilterOrder, GPIASTTLCh, PiezoCh)
KwikAnalysis.PlotGPIAS(FileList)


######################
import glob
import Kwik
import matplotlib.pyplot as plt
import numpy as np
import os
import shelve
from scipy import signal


with shelve.open(FileName) as Shelve:
    DataInfo = Shelve['DataInfo']
    Freqs = Shelve['Freqs']
    FreqOrder = Shelve['FreqOrder']
    FreqSlot = Shelve['FreqSlot']

for Key, Value in DataInfo.items():
    exec(str(Key) + '=' + 'Value')
del(Key, Value)

print('Preallocate memory...')
GPIAS = [[0] for _ in range(len(NoiseFrequency))]
for Freq in range(len(NoiseFrequency)):
    GPIAS[Freq] = [[0] for _ in range(NoOfTrials*2)]

AllTTLs = [[0] for _ in range(len(NoiseFrequency))]
for Freq in range(len(NoiseFrequency)):
    AllTTLs[Freq] = [[0] for _ in range(NoOfTrials*2)]

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
    
    print('Check if files are ok...')
    if 'Raw' not in globals():
        print('.kwd file is corrupted. Skipping dataset...')
        continue
    
    if 'Events' not in globals():
        print('.kwe/.kwik file is corrupted. Skipping dataset...')
        continue
    
    print('Data from ', RecFolder, ' loaded.')
    
    print('Get TTL data...')
    EventID = Events['TTLs']['user_data']['eventID']
    EventCh = Events['TTLs']['user_data']['event_channels']
    EventRec = Events['TTLs']['recording']
    EventSample = Events['TTLs']['time_samples']

    TTLChs = np.nonzero(np.bincount(EventCh))[0]
    TTLNo = {_: (sum(np.squeeze(EventCh) == _))//2 for _ in TTLChs}
#    TTLTimesRise = {_: Kwik.get_rising_edge_times(Files['kwe'], _)
#                    for _ in TTLChs}
#    TTLTimesFall = {_: Kwik.get_falling_edge_times(Files['kwe'], _) 
#                    for _ in TTLChs}
    
    for Rec in range(len(Raw['data'])):
        Rate = Raw['info'][str(Rec)]['sample_rate']
        RawTime = [int(round(Raw['timestamps'][str(Rec)][_]*Rate)) 
                   for _ in range(len(Raw['timestamps'][str(Rec)]))]
        
        NoOfSamplesBefore = int(round((TimeBeforeTTL*Rate)*10**-3))
        NoOfSamplesAfter = int(round((TimeAfterTTL*Rate)*10**-3))
        NoOfSamples = NoOfSamplesBefore + NoOfSamplesAfter
        XValues = (range(NoOfSamples)/Rate)*10**3
        
        print('Set filter...')
        passband = [FilterLow/(Rate/2), FilterHigh/(Rate/2)]
        f2, f1 = signal.butter(FilterOrder, passband, 'bandpass')
        
        print('Find TTL...')
        TTLLoc = [EventSample[_] for _ in range(len(EventID)) 
                  if EventCh[_] == GPIASTTLCh-1 and EventRec[_] == Rec]
#        TTLLoc = [int(round(TTLTimesRise[GPIASTTLCh-1][_]*Rate))
#                  for _ in range(len(TTLTimesRise[GPIASTTLCh-1])) 
#                  if EventRec[_] == Rec]
#        TTLLoc = TTLLoc[0] + (NoOfSamplesBefore+NoOfSamplesAfter)
        Start = int(round(TTLLoc[0]-NoOfSamplesBefore))
        End = int(round(TTLLoc[1]+NoOfSamplesAfter))
        
        Start = RawTime.index(Start)
        End = RawTime.index(End)
        
        print('Slicing and filtering...')
        Freq = FreqOrder[Rec][0]; Trial = FreqOrder[Rec][1]
        
        GPIAS[Freq][Trial] = Raw['data'][str(Rec)][Start:End,PiezoCh-1]        
        GPIAS[Freq][Trial] = [float(GPIAS[Freq][Trial][_]) 
                            for _ in range(NoOfSamples)]
        GPIAS[Freq][Trial] = abs(signal.hilbert(GPIAS[Freq][Trial]))
        GPIAS[Freq][Trial] = signal.filtfilt(f2, f1, GPIAS[Freq][Trial], 
                                           padtype='odd', padlen=0)
        
        AllTTLs[Freq][Trial] = [RawTime.index(TTLLoc[0])-Start, 
                                RawTime.index(TTLLoc[1])-RawTime.index(TTLLoc[0])]
        
        del(TTLLoc, Start, End, Freq, Trial)
    
    for Freq in range(len(NoiseFrequency)):
        NoGapAll = [GPIAS[Freq][_] for _ in range(len(GPIAS[Freq])) if _%2 == 0]
        GapAll = [GPIAS[Freq][_] for _ in range(len(GPIAS[Freq])) if _%2 != 0]
        NoGapSum = list(map(sum, zip(*NoGapAll)))
        GapSum = list(map(sum, zip(*GapAll)))
        
        GPIAS[Freq] = [0, 0]
        GPIAS[Freq][0] = [_/NoOfTrials for _ in NoGapSum]
        GPIAS[Freq][1] = [_/NoOfTrials for _ in GapSum]
        GPIAS[Freq][0] = signal.savgol_filter(GPIAS[Freq][0], 5, 2, mode='nearest')
        GPIAS[Freq][1] = signal.savgol_filter(GPIAS[Freq][1], 5, 2, mode='nearest')
        
        TTLNoGapAll = [AllTTLs[Freq][_] for _ in range(len(AllTTLs[Freq])) if _%2 == 0]
        TTLGapAll = [AllTTLs[Freq][_] for _ in range(len(AllTTLs[Freq])) if _%2 != 0]
        TTLNoGapSum = list(map(sum, zip(*TTLNoGapAll)))
        TTLGapSum = list(map(sum, zip(*TTLGapAll)))
        
        AllTTLs[Freq] = [0, 0]
        AllTTLs[Freq][0] = [int(round(_/NoOfTrials)) for _ in TTLNoGapSum]
        AllTTLs[Freq][1] = [int(round(_/NoOfTrials)) for _ in TTLGapSum]
        
        del(NoGapAll, GapAll, NoGapSum, GapSum)
        del(TTLNoGapAll, TTLGapAll, TTLNoGapSum, TTLGapSum)
                
    print('Plotting...')
    for Freq in range(len(NoiseFrequency)):
        plt.figure(Freq)
        plt.plot(XValues, GPIAS[Freq][0], color='r', label='No Gap')
        plt.plot(XValues, GPIAS[Freq][1], color='b', label='Gap')
        plt.axvspan(XValues[AllTTLs[Freq][0][0]], XValues[AllTTLs[Freq][0][1]], 
                    color='r', alpha=0.5, lw=0, label='Sound pulse')
#        plt.axvspan(XValues[AllTTLs[Freq][1][0]], XValues[AllTTLs[Freq][1][1]], 
#                    color='b', alpha=0.5, lw=0, label='Sound pulse (Gap)')
        plt.ylabel('voltage [mV]'); plt.xlabel('time [ms]')
        plt.legend(loc='best', frameon=False)
        plt.locator_params(tight=True)
        plt.tick_params(direction='out')
        plt.axes().spines['right'].set_visible(False)
        plt.axes().spines['top'].set_visible(False)
        plt.axes().yaxis.set_ticks_position('left')
        plt.axes().xaxis.set_ticks_position('bottom')
        plt.savefig('Figs/' + AnimalName + '-' + RecFolder[10:] + '-' + 
                    str(NoiseFrequency[Freq][0]) + '_' + 
                    str(NoiseFrequency[Freq][1]), transparent=True)
    print('Done.')
