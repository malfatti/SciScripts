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

FileName = '20160308170928-TestABR01-SoundStim'

# ABR
ABRCh = [1, 16]         # [RightChannel, LeftChannel], if order matters
ABRTimeBeforeTTL = 0    # in ms
ABRTimeAfterTTL = 12    # in ms
ABRTTLCh = 1            # TTL ch for ABR
FilterLow = 300         # High-pass frequency for bandpass filter
FilterHigh = 3000       # Low-pass frequency
FilterOrder = 4         # butter order

#==========#==========#==========#==========#

import KwikAnalysis
import shelve

with shelve.open(FileName) as Shelve: DataInfo = Shelve['DataInfo']
for Key, Value in DataInfo.items():
    exec(str(Key) + '=' + 'Value')
del(Shelve, Key, Value)

ABRs, XValues = KwikAnalysis.ABR(FileName, NoiseFrequency, SoundAmpF, ABRCh, 
                                 ABRTimeBeforeTTL, ABRTimeAfterTTL, ABRTTLCh, 
                                 FilterLow, FilterHigh, FilterOrder)
KwikAnalysis.PlotABR(ABRs, XValues, FileName, NoiseFrequency, SoundAmpF)

#######
import glob
import Kwik
import matplotlib.pyplot as plt
import numpy as np
import os
import shelve
from scipy import signal


with shelve.open(FileName) as Shelve: DataInfo = Shelve['DataInfo']
for Key, Value in DataInfo.items():
    exec(str(Key) + '=' + 'Value')
del(Shelve, Key, Value)


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
            
        elif '.db' in File:
            with shelve.open(File[:-3]) as Shelve: ExpInfo = Shelve['ExpInfo']
            for Key, Value in ExpInfo.items():
                exec(str(Key) + '=' + 'Value')
            del(Shelve, Key, Value)
    
    print('Check if files are ok...')
    if 'Raw' not in locals():
        print('.kwd file is corrupted. Skipping dataset...')
        continue
    
    if 'Events' not in locals():
        print('.kwe/.kwik file is corrupted. Skipping dataset...')
        continue
    
    print('Data from ', RecFolder, ' loaded.')
    
    print('Get TTL data...')
    EventID = Events['TTLs']['user_data']['eventID']
    EventCh = Events['TTLs']['user_data']['event_channels']
    EventRec = Events['TTLs']['recording']
    EventSample = Events['TTLs']['time_samples']

    TTLChs = np.nonzero(np.bincount(EventCh))[0]
    TTLRecs = np.nonzero(np.bincount(EventRec))[0]
    TTLsPerRec = {_Rec: [EventSample[_] for _ in range(len(EventRec)) 
                         if EventRec[_] == _Rec 
                         and EventCh[_] ==  ABRTTLCh-1 
                         and EventID[_] == 1]
                  for _Rec in TTLRecs}
    
    
    Rate = Raw['info']['0']['sample_rate']
    NoOfSamplesBefore = int(round((ABRTimeBeforeTTL*Rate)*10**-3))
    NoOfSamplesAfter = int(round((ABRTimeAfterTTL*Rate)*10**-3))
    NoOfSamples = NoOfSamplesBefore + NoOfSamplesAfter
    XValues = (range(NoOfSamples)/Rate)*10**3
    
    print('Set filter...')
    passband = [FilterLow/(Rate/2), FilterHigh/(Rate/2)]
    f2, f1 = signal.butter(FilterOrder, passband, 'bandpass')
    
    for Rec in range(len(Raw['data'])):
        RawTime = [int(round(Raw['timestamps'][str(Rec)][_]*Rate)) 
                   for _ in range(len(Raw['timestamps'][str(Rec)]))]
        
        rABR = [[0 for _ in range(NoOfSamples)] 
                for _ in range(len(TTLsPerRec[Rec]))]
        lABR = [[0 for _ in range(NoOfSamples)] 
                for _ in range(len(TTLsPerRec[Rec]))]
        
        print('Slicing and filtering ABRs...')
        for TTL in range(len(TTLsPerRec[Rec])):
            TTLLoc = TTLsPerRec[Rec][TTL]
            Start = int(round(TTLLoc-NoOfSamplesBefore))
            End = int(round(TTLLoc+NoOfSamplesAfter))
            Start = RawTime.index(Start)
            End = RawTime.index(End)
            
            rData = Raw['data'][str(Rec)][Start:End, ABRCh[0]-1]
            lData = Raw['data'][str(Rec)][Start:End, ABRCh[1]-1]
            
            rABR[TTL] = [float(rData[_]) for _ in range(NoOfSamples)]
            lABR[TTL] = [float(lData[_]) for _ in range(NoOfSamples)]
            
            rABR[TTL] = signal.filtfilt(f2, f1, rABR[TTL], padtype='odd', padlen=0)
            lABR[TTL] = signal.filtfilt(f2, f1, lABR[TTL], padtype='odd', padlen=0)
            
            del(TTLLoc, Start, End, rData, lData)
    
        print('Saving ABRs...')
        ABRs[0][Freq][Rec][len(ABRs[0][Freq][Rec])-1] = np.mean(rABR, axis=0)
        ABRs[1][Freq][Rec][len(ABRs[1][Freq][Rec])-1] = np.mean(lABR, axis=0)
        ABRs[0][Freq][Rec].append(0); ABRs[1][Freq][Rec].append(0)


# Remove extra zero
for Freq in range(len(ABRs[0])):
    for Rec in range(len(ABRs[0][Freq])):    
        del(ABRs[0][Freq][Rec][-1])
        del(ABRs[1][Freq][Rec][-1])

print('Saving data to ABRs.db...')
with shelve.open('ABRs') as Shelve:
    Shelve['ABRs'] = ABRs
    Shelve['DataInfo'] = DataInfo
    Shelve['XValues'] = XValues
print('Done.')


print('Plotting...')
with shelve.open(FileName) as Shelve: SoundIntensity = Shelve['SoundIntensity']

Colormaps = [plt.get_cmap('Reds'), plt.get_cmap('Blues')]
Colors = [[Colormaps[0](255-(_*20)), Colormaps[1](255-(_*20))] 
          for _ in range(len(SoundAmpF))]

for Trial in range(len(ABRs[0][0][0])):
    Fig, Axes = plt.subplots(len(NoiseFrequency), 2, sharex=True, 
                             figsize=(12, 12))
#    BigAx = Fig.add_subplot(111)
#    BigAx.set_axis_bgcolor('none')
#    BigAx.spines['top'].set_color('none')
#    BigAx.spines['bottom'].set_color('none')
#    BigAx.spines['left'].set_color('none')
#    BigAx.spines['right'].set_color('none')
#    BigAx.tick_params(labelcolor='none', top='off', 
#                      bottom='off', left='off', right='off')
#    BigAx.set_xlabel('Time [ms]')
#    BigAx.set_ylabel('Voltage [µV]')
#    BigAx.set_title('Auditory brainstem responses')
    
    for Ear in range(2):
        for Freq in range(len(NoiseFrequency)):
            for AmpF in range(len(SoundAmpF)):
                if Ear == 0:
                    Axes[Freq][Ear].plot(XValues, 
                                         ABRs[Ear][Freq][AmpF][Trial], 
                                         color=Colors[AmpF][Ear], 
                                         label=str(round(SoundIntensity[Freq]
                                                     [str(SoundAmpF[AmpF])]))
                                               + ' dB')
                else:
                    Axes[Freq][Ear].plot(XValues, 
                                         ABRs[Ear][Freq][AmpF][Trial], 
                                         color=Colors[AmpF][Ear],
                                         label=str(round(SoundIntensity[Freq]
                                                     [str(SoundAmpF[AmpF])]))
                                               + ' dB')
                
                Axes[Freq][Ear].legend(loc='best', frameon=False)
                Axes[Freq][Ear].spines['right'].set_visible(False)
                Axes[Freq][Ear].spines['top'].set_visible(False)
                Axes[Freq][Ear].yaxis.set_ticks_position('left')
                Axes[Freq][Ear].xaxis.set_ticks_position('bottom')
                Axes[Freq][Ear].set_title(str(NoiseFrequency[Freq]))
                Axes[Freq][Ear].set_ylabel('voltage [µV]')
#                Axes[Freq][Ear].tick_params(direction='out')
    Axes[-1][0].set_xlabel('time [ms]'); Axes[-1][1].set_xlabel('Time [ms]')
    Fig.tight_layout()
