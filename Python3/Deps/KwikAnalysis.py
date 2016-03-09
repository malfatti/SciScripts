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

import glob
import Kwik
import matplotlib.pyplot as plt
import numpy as np
import os
import shelve
from scipy import signal

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


def ABR(FileName, NoiseFrequency, SoundAmpF, ABRCh=[1, 16], ABRTimeBeforeTTL=0, 
        ABRTimeAfterTTL=12, ABRTTLCh=1, FilterLow=300, FilterHigh=3000, 
        FilterOrder=4):
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
                with shelve.open(File[:-3]) as Shelve: 
                    ExpInfo = Shelve['ExpInfo']
                locals().update(ExpInfo)
        
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
        Shelve['XValues'] = XValues
    print('Done.')
    
    return(ABRs, XValues)


def PlotABR(AnimalName, ABRs, XValues, FileName, NoiseFrequency, SoundAmpF):
    print('Plotting...')
    with shelve.open(FileName) as Shelve: SoundIntensity = Shelve['SoundIntensity']
    
    Colormaps = [plt.get_cmap('Reds'), plt.get_cmap('Blues')]
    Colors = [[Colormaps[0](255-(_*20)), Colormaps[1](255-(_*20))] 
              for _ in range(len(SoundAmpF))]
    
    for Trial in range(len(ABRs[0][0][0])):
        Fig, Axes = plt.subplots(len(NoiseFrequency), 2, sharex=True, 
                                 figsize=(12, 12))

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
        Axes[-1][0].set_xlabel('time [ms]'); Axes[-1][1].set_xlabel('Time [ms]')
        Fig.tight_layout()


def GPIAS(GPIASTimeBeforeTTL=50, GPIASTimeAfterTTL=150, FilterLow=3, 
          FilterHigh=300, FilterOrder=4, GPIASTTLCh=1, PiezoCh=1):
    
    print('set paths...')
    os.makedirs('Figs', exist_ok=True)    # Figs folder
    DirList = glob.glob('KwikFiles/*'); DirList.sort()
    
    GPIAS = [0]
    AllTTLs = [0]
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
                with shelve.open(File[:-3]) as Shelve: 
                    DataInfo = Shelve['DataInfo']
                    Freqs = Shelve['Freqs']
                    FreqOrder = Shelve['FreqOrder']
                    FreqSlot = Shelve['FreqSlot']
        
        print('Check if files are ok...')
        if 'Raw' not in locals():
            print('.kwd file is corrupted. Skipping dataset...')
            continue
        
        if 'Events' not in locals():
            print('.kwe/.kwik file is corrupted. Skipping dataset...')
            continue
        
        print('Data from ', RecFolder, ' loaded.')
        
        print('Preallocate memory...')
        GPIASTemp = [[0] for _ in range(len(DataInfo['NoiseFrequency']))]
        for Freq in range(len(DataInfo['NoiseFrequency'])):
            GPIASTemp[Freq] = [[0] for _ in range(DataInfo['NoOfTrials']*2)]
        GPIAS[len(GPIAS)-1] = [GPIASTemp]
        
        AllTTLsTemp = [[0] for _ in range(len(DataInfo['NoiseFrequency']))]
        for Freq in range(len(DataInfo['NoiseFrequency'])):
            AllTTLsTemp[Freq] = [[0] for _ in range(DataInfo['NoOfTrials']*2)]
        AllTTLs[len(AllTTLs)-1] = [AllTTLsTemp]
        
        print('Get TTL data...')
        EventID = Events['TTLs']['user_data']['eventID']
        EventCh = Events['TTLs']['user_data']['event_channels']
        EventRec = Events['TTLs']['recording']
        EventSample = Events['TTLs']['time_samples']
        
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
            RawTime = [int(round(Raw['timestamps'][str(Rec)][_]*Rate)) 
                       for _ in range(len(Raw['timestamps'][str(Rec)]))]
            
            print('Find TTL...')
            TTLLoc = [EventSample[_] for _ in range(len(EventID)) 
                      if EventCh[_] == GPIASTTLCh-1 and EventRec[_] == Rec]
            
            Start = int(round(TTLLoc[0]-NoOfSamplesBefore))
            End = int(round(TTLLoc[1]+NoOfSamplesAfter))
            Start = RawTime.index(Start)
            End = RawTime.index(End)
            
            TTLStart = RawTime.index(TTLLoc[0])-Start
            TTLEnd = RawTime.index(TTLLoc[1])-RawTime.index(TTLLoc[0])
            
            print('Slicing and filtering...')            
            gData = Raw['data'][str(Rec)][Start:End,PiezoCh-1]        
            gData = [float(gData[_]) for _ in range(NoOfSamples)]
            gData = abs(signal.hilbert(gData))
            gData = signal.filtfilt(f2, f1, gData, padtype='odd', padlen=0)
            
            Freq = FreqOrder[Rec][0]; Trial = FreqOrder[Rec][1]
            AllTTLs[len(AllTTLs)-1][Freq][Trial] = [TTLStart, TTLEnd]
            GPIAS[len(GPIAS)-1][Freq][Trial] = gData[:]
            
            del(TTLLoc, Start, End, TTLStart, TTLEnd, Freq, Trial, gData)
        
        for Freq in range(len(DataInfo['NoiseFrequency'])):
            gData = GPIAS[len(GPIAS)-1][Freq][:]
            NoGapAll = [gData[_] for _ in range(len(gData)) if _%2 == 0]
            GapAll = [gData[_] for _ in range(len(gData)) if _%2 != 0]
            NoGapSum = list(map(sum, zip(*NoGapAll)))
            GapSum = list(map(sum, zip(*GapAll)))
            
            gData = [0, 0]
            gData[0] = [_/DataInfo['NoOfTrials'] for _ in NoGapSum]
            gData[1] = [_/DataInfo['NoOfTrials'] for _ in GapSum]
            gData[0] = signal.savgol_filter(gData[0], 5, 2, mode='nearest')
            gData[1] = signal.savgol_filter(gData[1], 5, 2, mode='nearest')
            GPIAS[len(GPIAS)-1][Freq] = gData[:]
            
            tData = AllTTLs[len(AllTTLs)-1][Freq][:]
            TTLNoGapAll = [tData[_] for _ in range(len(tData)) if _%2 == 0]
            TTLGapAll = [tData[_] for _ in range(len(tData)) if _%2 != 0]
            TTLNoGapSum = list(map(sum, zip(*TTLNoGapAll)))
            TTLGapSum = list(map(sum, zip(*TTLGapAll)))
            
            tData = [0, 0]
            tData[0] = [int(round(_/DataInfo['NoOfTrials'])) for _ in TTLNoGapSum]
            tData[1] = [int(round(_/DataInfo['NoOfTrials'])) for _ in TTLGapSum]
            AllTTLs[len(AllTTLs)-1][Freq] = tData[:]
            
            del(NoGapAll, GapAll, NoGapSum, GapSum, gData)
            del(TTLNoGapAll, TTLGapAll, TTLNoGapSum, TTLGapSum)
            
        
        print('Saving data to ' + DataInfo['AnimalName'] + '-GPIAS-' + 
              RecFolder[10:])
        with shelve.open('GPIAS') as Shelve:
            Shelve['GPIAS'] = GPIAS
            Shelve['AllTTLs'] = AllTTLs
            Shelve['XValues'] = XValues
            Shelve['RecFolder'] = RecFolder
        print('Done.')
        GPIAS.append(0)
        AllTTLs.append(0)
    
    del(GPIAS[-1], AllTTLs[-1])
    print('Finished.')

    return(GPIAS, AllTTLs, XValues)

def PlotGPIAS(AnimalName, NoiseFrequency):
    if 'GPIAS' not in globals():
        if os.path.isfile('GPIAS.db'):
            print('Loading data from GPIAS.db ...')
            with shelve.open('GPIAS') as Shelve:
                GPIAS = Shelve['GPIAS']
                AllTTLs = Shelve['AllTTLs']
                XValues = Shelve['XValues']
        else:
            print('Data not found')
            return None
    else:
        global GPIAS, XValues
    
    print('Plotting...')
    for Freq in range(len(NoiseFrequency)):
        plt.figure(Freq)
        plt.plot(XValues, GPIAS[Freq][0], color='r', label='No Gap')
        plt.plot(XValues, GPIAS[Freq][1], color='b', label='Gap')
        plt.axvspan(XValues[AllTTLs[Freq][0][0]], XValues[AllTTLs[Freq][0][1]], 
                    color='k', alpha=0.5, lw=0, label='Sound pulse')
    #        plt.axvspan(XValues[AllTTLs[Freq][1][0]], XValues[AllTTLs[Freq][1][1]], 
    #                    color='b', alpha=0.5, lw=0, label='Sound pulse (Gap)')
        plt.ylabel('Voltage [mV]'); plt.xlabel('Time [ms]')
        plt.legend(loc='best', frameon=False)
        plt.locator_params(tight=True)
        plt.tick_params(direction='out')
        plt.axes().spines['right'].set_visible(False)
        plt.axes().spines['top'].set_visible(False)
        plt.axes().yaxis.set_ticks_position('left')
        plt.axes().xaxis.set_ticks_position('bottom')
        plt.savefig('Figs/' + AnimalName + '-GPIAS-' + RecFolder[10:] + 
                    '-' + str(NoiseFrequency[Freq][0]) + '_' + 
                    str(NoiseFrequency[Freq][1]), transparent=True)
    print('Done.')
