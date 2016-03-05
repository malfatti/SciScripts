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

FileName = '20160305105342-GPIAS-TestSetup02'

PiezoCh = 1
GPIASTTLCh = 1
TimeBeforeTTL = 50    # in ms
TimeAfterTTL = 190    # in ms
FilterLow = 300       # High-pass frequency for bandpass filter
FilterHigh = 3000     # Low-pass frequency
FilterOrder = 4       # butter order


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
GPIASs = [[0] for _ in range(len(NoiseFrequency))]
for Freq in range(len(NoiseFrequency)):
    GPIASs[Freq] = [[0], [0]]

Gap = [[0] for _ in range(len(NoiseFrequency))]
for Freq in range(len(NoiseFrequency)):
    Gap[Freq] = [[0] for _ in range(NoOfTrials*2)]
NoGap = Gap[:]

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
    TTLTimesRise = {_: Kwik.get_rising_edge_times(Files['kwe'], _)
                    for _ in TTLChs}
    TTLTimesFall = {_: Kwik.get_falling_edge_times(Files['kwe'], _) 
                    for _ in TTLChs}
    
    for Rec in range(len(Raw['data'])):
        Rate = Raw['info'][str(Rec)]['sample_rate']
#        RawTime = [int(round(Raw['timestamps'][str(Rec)][_]*Rate)) 
#                   for _ in range(len(Raw['timestamps'][str(Rec)]))]
#        
        NoOfSamplesBefore = int(round((TimeBeforeTTL*Rate)*10**-3))
        NoOfSamplesAfter = int(round((TimeAfterTTL*Rate)*10**-3))
        NoOfSamples = NoOfSamplesBefore + NoOfSamplesAfter
        XValues = (range(NoOfSamples)/Rate)*10**3
        
        print('Set filter...')
        passband = [FilterLow/(Rate/2), FilterHigh/(Rate/2)]
        f2, f1 = signal.butter(FilterOrder, passband, 'bandpass')
        
        print('Find TTL...')
        TTLLoc = [int(round(TTLTimesRise[GPIASTTLCh-1][_]*Rate))
                  for _ in range(len(TTLTimesRise[GPIASTTLCh-1])) 
                  if EventRec[_] == Rec]
#        TTLLoc = TTLLoc[0] + (NoOfSamplesBefore+NoOfSamplesAfter)
        Start = TTLLoc-NoOfSamplesBefore
        End = TTLLoc+NoOfSamplesAfter
        
        print('Slicing and filtering...')
        Freq = FreqOrder[Rec][0]; Trial = FreqOrder[Rec][1]
        Gap[Freq][Trial] = Raw['data'][str(Rec)][Start:End,PiezoCh-1]        
        Gap[Freq][Trial] = [float(Gap[Freq][Trial][_]) 
                            for _ in range(NoOfSamples)]
        Gap[Freq][Trial] = abs(signal.hilbert(Gap[Freq][Trial]))
        Gap[Freq][Trial] = signal.filtfilt(f2, f1, Gap[Freq][Trial], 
                                           padtype='odd', padlen=0)
    
    for Freq in range(len(NoiseFrequency)):
        NoGapAll = [Gap[Freq][_] for _ in range(len(Gap[Freq])) if _%2 == 0]
        GapAll = [Gap[Freq][_] for _ in range(len(Gap[Freq])) if _%2 != 0]
    
    NoGapSum = list(map(sum, zip(*NoGapAll)))
    GapSum = list(map(sum, zip(*GapAll)))
        
        
        
        del(TTLLoc, Start, End, Freq, Trial)