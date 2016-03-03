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

Experiment: stimulation of the brainstem using light and sound, recording with
a silicon probe (16 channels) + 2 tungsten wires + reference screw.
"""
#%% Set experiment details

FileName = '20160303101319-TestSetup01-SoundStim'

# ABR
ABRCh = [1, 16]         # [RightChannel, LeftChannel], if order matters
ABRTimeBeforeTTL = 0    # in ms
ABRTimeAfterTTL = 12    # in ms
ABRTTLCh = 1            # TTL ch for ABR

#==========#==========#==========#==========#

import glob
import Kwik
import matplotlib.pyplot as plt
import numpy as np
import os
import shelve
from scipy import signal


with shelve.open(FileName) as Shelve:
    DataInfo = Shelve['DataInfo']

for Key, Value in DataInfo.items():
    exec(str(Key) + '=' + 'Value')
del(Key, Value)


## Remove date from folder name
RenameFolders = input('Rename folders (BE CAREFUL)? [y/N] ')
if RenameFolders in ['y', 'Y', 'yes', 'Yes', 'YES']:
    DirList = glob.glob('KwikFiles/*'); DirList.sort()
    for FolderName in DirList:
        NewFolderName = ''.join([FolderName[:10], FolderName[21:]])
        NewFolderName = NewFolderName.replace("-", "")
        os.rename(FolderName, NewFolderName)
        print(FolderName, ' moved to ', NewFolderName)
    del(RenameFolders, DirList, FolderName, NewFolderName)


## ABR traces

print('Preallocate memory...')
ABRs = [[], []]
ABRs[0] = [[0] for _ in range(len(NoiseFrequency))]
ABRs[1] = [[0] for _ in range(len(NoiseFrequency))]
for Freq in range(len(NoiseFrequency)):
    ABRs[0][Freq] = [[[0]] for _ in range(len(SoundAmpF))]
    ABRs[1][Freq] = [[[0]] for _ in range(len(SoundAmpF))]
del(Freq)

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
                Raw = Kwik.load(File)
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
            
        elif '.kwx' in File:
            try:
                Spks = Kwik.load(File)
                Files['kwx'] = File
            except OSError:
                print('File ', File, " is corrupted :'(")
        
        elif 'dB.db' in File:
            with shelve.open(File[:-3]) as Shelve:
                ExpInfo = Shelve['ExpInfo']
            
            Files['db'] = File
            for Key, Value in ExpInfo.items():
                exec(str(Key) + '=' + 'Value')
            del(Key, Value)
    
    print('Check if files are ok...')
    if 'Raw' not in globals():
        print('.kwd file is corrupted. Skipping dataset...')
        continue
    
    if 'Events' not in globals():
        print('.kwe/.kwik file is corrupted. Skipping dataset...')
        continue
    
    print('Data from ', RecFolder, ' loaded.')
    
    print('Get rate and calculate No of samples...')
    Rate = Raw['info']['sample_rate']
    NoOfSamplesBefore = int((ABRTimeBeforeTTL*Rate)*10**-3)
    NoOfSamplesAfter = int((ABRTimeAfterTTL*Rate)*10**-3)
    NoOfSamples = NoOfSamplesBefore + NoOfSamplesAfter
    XValues = (range(NoOfSamples)/Rate)*10**3
    
    print('Set filter...')
    passband = [300/(Rate/2), 3000/(Rate/2)]
    f2, f1 = signal.butter(4, passband, 'bandpass')
    
    print('Get TTL data...')
    EventID = Events['TTLs']['user_data']['eventID']
    EventCh = Events['TTLs']['user_data']['event_channels']
    EventSample = Events['TTLs']['time_samples']

    TTLChs = np.nonzero(np.bincount(EventCh))[0]
    TTLNo = {_: sum(np.squeeze(EventCh) == _) for _ in TTLChs}
    TTLTimes = {_: Kwik.get_rising_edge_times(Files['kwe'], _) for _ in TTLChs}
    
    rABR = [[0]*NoOfSamples]*TTLNo[ABRTTLCh-1]; lABR = rABR[:]
    
    print('Slicing and filtering ABRs...')
    for TTL in range(TTLNo[ABRTTLCh-1]):
        TTLLoc = int(TTLTimes[ABRTTLCh-1][TTL])
        Start = TTLLoc-NoOfSamplesBefore
        End = TTLLoc+NoOfSamplesAfter
        
        rData = Raw['data'][Start:End, ABRCh[0]-1]
        lData = Raw['data'][Start:End, ABRCh[1]-1]
        
        rABR[TTL] = [float(rData[_]) for _ in range(NoOfSamples)]
        lABR[TTL] = [float(lData[_]) for _ in range(NoOfSamples)]
        
        rABR[TTL] = signal.filtfilt(f2, f1, rABR[TTL], padtype='odd', padlen=0)
        lABR[TTL] = signal.filtfilt(f2, f1, lABR[TTL], padtype='odd', padlen=0)
        
        del(TTLLoc, Start, End, rData, lData)
    
    print('Saving ABRs...')
    ABRs[0][Freq][AmpF][len(ABRs[0][Freq][AmpF])-1] = np.mean(rABR, axis=0)
    ABRs[1][Freq][AmpF][len(ABRs[0][Freq][AmpF])-1] = np.mean(lABR, axis=0)
    ABRs[0][Freq][AmpF].append(0); ABRs[1][Freq][AmpF].append(0)

print('Done. Plotting...')

for Trial in range(len(ABRs[0][0][0])-1):
    Fig, Axes = plt.subplots(len(NoiseFrequency), 2, sharex=True, figsize=(12, 12))
    BigAx = Fig.add_subplot(111)
    BigAx.set_axis_bgcolor('none')
    BigAx.spines['top'].set_color('none'); BigAx.spines['bottom'].set_color('none')
    BigAx.spines['left'].set_color('none'); BigAx.spines['right'].set_color('none')
    BigAx.tick_params(labelcolor='none', top='off', 
                      bottom='off', left='off', right='off')
    BigAx.set_xlabel('Time [ms]')
    BigAx.set_ylabel('Voltage [µv]')
    
    for Ear in range(2):
        for Freq in range(len(NoiseFrequency)):
            for AmpF in range(len(SoundAmpF)):
                if Ear == 0:
                    Axes[Freq][Ear].plot(XValues, 
                                         ABRs[Ear][Freq][AmpF][Trial], 
                                         color='r')
                else:
                    Axes[Freq][Ear].plot(XValues, 
                                         ABRs[Ear][Freq][AmpF][Trial], 
                                         color='b')

