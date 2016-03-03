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

FileName = '20160302190255-TestSetup01-SoundStim'

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


## Create folders
os.makedirs('Figs', exist_ok=True)    # Figs folder


## ABR traces
ABRs = [[], []]
ABRs[0] = [[0] for _ in range(len(NoiseFrequency))]
ABRs[1] = [[0] for _ in range(len(NoiseFrequency))]
for Freq in range(len(NoiseFrequency)):
    ABRs[0][Freq] = [[[0]] for _ in range(len(SoundAmpF))]
    ABRs[1][Freq] = [[[0]] for _ in range(len(SoundAmpF))]
del(Freq)


DirList = glob.glob('KwikFiles/*'); DirList.sort()
for RecFolder in DirList:
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
        
        elif '.db' in File:
            with shelve.open(File) as Shelve:
                ExpInfo = Shelve['ExpInfo']
            
            Files['db'] = File
            for Key, Value in ExpInfo.items():
                exec(str(Key) + '=' + 'Value')
            del(Key, Value)
    
    if 'Raw' not in globals():
        print('.kwd file is corrupted. Skipping dataset...')
        continue
    
    if 'Events' not in globals():
        print('.kwe/.kwik file is corrupted. Skipping dataset...')
        continue
    
    print('Data from ', RecFolder, ' loaded.')
    
    # Getting data...
    Rate = Raw['info']['sample_rate']
    NoOfSamplesBefore = int((ABRTimeBeforeTTL*Rate)*10**-3)
    NoOfSamplesAfter = int((ABRTimeAfterTTL*Rate)*10**-3)
    NoOfSamples = NoOfSamplesBefore + NoOfSamplesAfter
    XValues = (range(NoOfSamples)/Rate)*10**3
    
    EventID = Events['TTLs']['user_data']['eventID']
    EventCh = Events['TTLs']['user_data']['event_channels']
    EventSample = Events['TTLs']['time_samples']

    TTLChs = np.nonzero(np.bincount(EventCh))[0]
    TTLNo = {_: sum(np.squeeze(EventCh) == _) for _ in TTLChs}
    TTLTimes = {_: Kwik.get_rising_edge_times(Files['kwe'], _) for _ in TTLChs}
    
    rABR = [[0]*NoOfSamples]*TTLNo[ABRTTLCh-1]; lABR = rABR[:]
    
    for TTL in range(TTLNo[ABRTTLCh-1]):
        TTLLoc = int(TTLTimes[ABRTTLCh-1][TTL])
        rABR[TTL] = [float(Raw['data'][TTLLoc-NoOfSamplesBefore:
                                       TTLLoc+NoOfSamplesAfter, ABRCh[0]-1][_]) 
                     for _ in range(NoOfSamples)]
        lABR[TTL] = [float(Raw['data'][TTLLoc-NoOfSamplesBefore:
                                       TTLLoc+NoOfSamplesAfter, ABRCh[1]-1][_]) 
                     for _ in range(NoOfSamples)]
        
        passband = [300/(Rate/2), 3000/(Rate/2)]
        f2, f1 = signal.butter(4, passband, 'bandpass')
        rABR[TTL] = signal.filtfilt(f2, f1, rABR[TTL], padtype='odd', padlen=0)
        lABR[TTL] = signal.filtfilt(f2, f1, lABR[TTL], padtype='odd', padlen=0)
    
    ABRs[0][Freq][AmpF][len(ABRs[0][Freq][AmpF])-1] = np.mean(rABR, axis=0)
    ABRs[1][Freq][AmpF][len(ABRs[0][Freq][AmpF])-1] = np.mean(lABR, axis=0)
    ABRs[0][Freq][AmpF].append(0); ABRs[1][Freq][AmpF].append(0)


plt.figure()
plt.plot(XValues, np.mean(rABR[:62], axis=0), color='r', label='Right')
plt.plot(XValues, np.mean(lABR[:62], axis=0), color='b', label='Left')
plt.ylabel('Voltage [µV]'); plt.xlabel('Time [ms]')
plt.legend(loc='best', frameon=False)
plt.locator_params(tight=True)
plt.tick_params(direction='out')
plt.axes().spines['right'].set_visible(False)
plt.axes().spines['top'].set_visible(False)
plt.axes().yaxis.set_ticks_position('left')
plt.axes().xaxis.set_ticks_position('bottom')
plt.show()
#input('Press enter to save figure 1.')
plt.savefig('Figs/', RecFolder[10:], '-ABR.pdf', transparent=True)
