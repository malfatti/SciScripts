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
import LoadHdf5Files
from matplotlib import rcParams
from matplotlib import pyplot as plt
from numbers import Number
import numpy as np
import os
from scipy import signal


## Lower-level functions

def SetLaTexPlot():
    print('Set plot...')
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


## Higher-level functions

def ABR(FileName, ABRCh=[1, 16], ABRTimeBeforeTTL=0, ABRTimeAfterTTL=12, 
        ABRTTLCh=1, FilterLow=300, FilterHigh=3000, FilterOrder=4, 
        StimType='Sound'):
    """
    Analyze ABRs from data recorded with OpenEphys. A '*ABRs.hdf5' file will be 
    saved in cwd, containing:
        - ABRs dict, where data will be saved as 
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
    print('set paths...')
    os.makedirs('Figs', exist_ok=True)    # Figs folder
    DirList = glob.glob('KwikFiles/*'); DirList.sort()
    
    print('Load DataInfo...')
    DataInfo, Exps = LoadHdf5Files.ExpDataInfo(FileName, DirList, StimType, 
                                               Var='Exps')
    
    print('Preallocate memory...')
    ABRs = [[], []]
    ABRs[0] = [[0] for _ in range(len(DataInfo['NoiseFrequency']))]
    ABRs[1] = [[0] for _ in range(len(DataInfo['NoiseFrequency']))]
    for Freq in range(len(DataInfo['NoiseFrequency'])):
        Key = str(DataInfo['NoiseFrequency'][Freq][0]) + '-' \
              + str(DataInfo['NoiseFrequency'][Freq][1])
        ABRs[0][Freq] = [{} for _ in range(len(DataInfo['SoundAmpF'][Key]))]
        ABRs[1][Freq] = [{} for _ in range(len(DataInfo['SoundAmpF'][Key]))]
    del(Freq)
    
    for RecFolder in Exps:
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
        
        ExpInfo = LoadHdf5Files.ExpExpInfo(FileName, DirList, RecFolder)
        
        print('Check if files are ok...')
        if 'Raw' not in locals():
            print('.kwd file is corrupted. Skipping dataset...')
            continue
        
        if 'Events' not in locals():
            print('.kwe/.kwik file is corrupted. Skipping dataset...')
            continue
        
        print('Data from ', RecFolder, ' loaded.')
        
        if '0' not in list(Raw['data'].keys()):
            print('Rec numbers are wrong. Fixing...')
            for iKey in Raw.keys():
                Recs = list(Raw[iKey].keys())
                Recs = [int(_) for _ in Recs]; Min = min(Recs)
                
                EventRec = Events['TTLs']['recording'][:]
                for _ in range(len(EventRec)): EventRec[_] = EventRec[_] - Min
                
                for Key in Recs:
                    Raw[iKey][str(Key-Min)] = Raw[iKey].pop(str(Key))
                
            print('Fixed.')
        else:
            EventRec = Events['TTLs']['recording']
        
        print('Get TTL data...')
        EventID = Events['TTLs']['user_data']['eventID']
        EventCh = Events['TTLs']['user_data']['event_channels']
        EventSample = Events['TTLs']['time_samples']
    
        TTLChs = np.nonzero(np.bincount(EventCh))[0]
        TTLRecs = np.nonzero(np.bincount(EventRec))[0]
        TTLsPerRec = {_Rec: [EventSample[_] for _ in range(len(EventRec)) 
                             if EventRec[_] == _Rec 
                             and EventCh[_] ==  ABRTTLCh-1 
                             and EventID[_] == 1]
                      for _Rec in TTLRecs}
        TTLRising = Kwik.get_rising_edge_times(Files['kwe'], ABRTTLCh-1)
        
        Rate = Raw['info']['0']['sample_rate']
#        NoOfSamplesBefore = int(round((ABRTimeBeforeTTL*Rate)*10**-3))
#        NoOfSamplesAfter = int(round((ABRTimeAfterTTL*Rate)*10**-3))
#        NoOfSamples = NoOfSamplesBefore + NoOfSamplesAfter
#        XValues = (range(NoOfSamples)/Rate)*10**3
        NoOfSamplesBefore = ABRTimeBeforeTTL*int(Rate*10**-3)
        NoOfSamplesAfter = ABRTimeAfterTTL*int(Rate*10**-3)
        NoOfSamples = NoOfSamplesBefore + NoOfSamplesAfter
#        XValues = (range(NoOfSamples)/Rate)*10**3
        XValues = list((range((NoOfSamplesBefore)*-1, 0)/Rate)*10**3) + \
                  list((range(NoOfSamplesAfter)/Rate)*10**3)
        
#        print('Set filter...')
#        passband = [FilterLow/(Rate/2), FilterHigh/(Rate/2)]
#        f2, f1 = signal.butter(FilterOrder, passband, 'bandpass')
        Lf2, Lf1 = signal.butter(FilterOrder, FilterLow/(Rate/2), 'highpass')
        Hf2, Hf1 = signal.butter(FilterOrder, FilterHigh/(Rate/2), 'lowpass')
        
        for Rec in range(len(Raw['data'])):
#            RawTime = [int(round(Raw['timestamps'][str(Rec)][_]*Rate)) 
#                       for _ in range(len(Raw['timestamps'][str(Rec)]))]
            TTLNo = [0]
            for _ in range(1, len(TTLsPerRec)+1):
                TTLNo = TTLNo + [len(TTLsPerRec[_-1]) + TTLNo[-1]]

            TTLNo = [0] + [TTLNo[_]-1 for _ in range(1, len(TTLNo))]
            
            RawTime = list(Raw['timestamps'][str(Rec)])
            
            rABR = [[0 for _ in range(NoOfSamples)] 
                    for _ in range(len(TTLsPerRec[Rec]))]
            lABR = [[0 for _ in range(NoOfSamples)] 
                    for _ in range(len(TTLsPerRec[Rec]))]
            
            if Rec == 0:
                sTTLNo = 0
            else:
                sTTLNo = TTLNo[Rec] + 1
            
            print('Slicing and filtering ABRs Rec ', str(Rec), '...')
            for TTL in range(sTTLNo, TTLNo[Rec+1]):
#                TTLLoc = int(TTLsPerRec[Rec][TTL])
#                Start = TTLLoc-NoOfSamplesBefore
#                End = TTLLoc+NoOfSamplesAfter
#                try:
#                    Start = RawTime.index(Start)
#                    End = RawTime.index(End)
#                except ValueError:
#                    print('ValueError: Timestamp is messed up in Freq', 
#                          str(ExpInfo['Hz']), ' AmpF ', str(Rec), ' :(')
#                    break
                TTLLoc = RawTime.index(TTLRising[TTL])
                    
                Start = TTLLoc-NoOfSamplesBefore
                End = TTLLoc+NoOfSamplesAfter
                
                Index = TTL-sTTLNo
                
                rABR[Index] = Raw['data'][str(Rec)][Start:End, ABRCh[0]-1]
                lABR[Index] = Raw['data'][str(Rec)][Start:End, ABRCh[1]-1]
                
                rABR[Index] = signal.filtfilt(Lf2, Lf1, rABR[Index], padlen=0)
#                lABR[Index] = signal.filtfilt(Lf2, Lf1, lABR[Index], padlen=0)

                del(TTLLoc, Start, End)
            
#            # CrazyMA
#            rData = rABR[0]; lData = lABR[0]
#            for _ in range(1, len(rABR)):
#                rData = (rData + rABR[_])/2
#                lData = (lData + lABR[_])/2
#            rABR = rData; lABR = lData
            
            # Mean
            rABR = np.mean(rABR, axis=0)
#            lABR = np.mean(lABR, axis=0)
#            rData = np.mean(rABR, axis=0)
#            lData = np.mean(lABR, axis=0)
            
            rABR = signal.filtfilt(Hf2, Hf1, rABR, padlen=0)
#            lABR = signal.filtfilt(Hf2, Hf1, lABR, padlen=0)
#            rData = signal.filtfilt(Hf2, Hf1, rData, padlen=0)
#            lData = signal.filtfilt(Hf2, Hf1, lData, padlen=0)

#            rABR = signal.savgol_filter(rABR, 9, 2, mode='nearest')
#            lABR = signal.savgol_filter(lABR, 9, 2, mode='nearest')
#            rABR = (rABR/10) + DataInfo['Intensities'][Rec]
#            lABR = (lABR/10) + DataInfo['Intensities'][Rec]
            
            if ExpInfo['DVCoord'] not in ABRs[0][ExpInfo['Hz']][Rec]:
                ABRs[0][ExpInfo['Hz']][Rec][ExpInfo['DVCoord']] = [rABR]
                ABRs[1][ExpInfo['Hz']][Rec][ExpInfo['DVCoord']] = [lABR]
            else:
                ABRs[0][ExpInfo['Hz']][Rec][ExpInfo['DVCoord']].append(rABR)
                ABRs[1][ExpInfo['Hz']][Rec][ExpInfo['DVCoord']].append(lABR)
    
    print('Saving data to ' + FileName)
    GroupName = 'ABRs-' + datetime.datetime.now().strftime("%Y%m%d%H%M%S")
    with h5py.File(FileName) as F:
        if 'ABRs' in F.keys():
            F[GroupName] = F['ABRs']; del(F['ABRs'])
        
        F.create_group('ABRs')
        F['ABRs'].attrs['XValues'] = XValues
        F['ABRs'].create_group('Right'); F['ABRs'].create_group('Left')
            
        for Freq in range(len(ABRs[0])):
            for AmpF in range(len(ABRs[0][Freq])):
                for DV in ABRs[0][Freq][AmpF].keys():
                    Path = str(Freq) + '/' + \
                           str(AmpF) + '/' + \
                           DV
                    
                    F['ABRs']['Right'].create_group(Path)
                    F['ABRs']['Left'].create_group(Path)
                    del(Path)
                    
                    for Trial in range(len(ABRs[0][Freq][AmpF][DV])):
                        F['ABRs']['Right'][str(Freq)][str(AmpF)][DV][str(Trial)] = \
                            ABRs[0][Freq][AmpF][DV][Trial][:]
                        F['ABRs']['Left'][str(Freq)][str(AmpF)][DV][str(Trial)] = \
                            ABRs[1][Freq][AmpF][DV][Trial][:]
    
    print('Done.')
    return(None)

    
def ABRAnalogTTLs(FileName, ABRCh=[1, 16], ABRTimeBeforeTTL=0, ABRTimeAfterTTL=12, 
        ABRTTLCh=1, FilterLow=300, FilterHigh=3000, FilterOrder=4, 
        StimType='Sound'):
    """
    Analyze ABRs from data recorded with OpenEphys. A '*ABRs.hdf5' file will be 
    saved in cwd, containing:
        - ABRs dict, where data will be saved as 
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
    print('set paths...')
    os.makedirs('Figs', exist_ok=True)    # Figs folder
    DirList = glob.glob('KwikFiles/*'); DirList.sort()
    
    print('Load DataInfo...')
    DataInfo, Exps = LoadHdf5Files.ExpDataInfo(FileName, DirList, StimType, 
                                               Var='Exps')
    
    print('Preallocate memory...')
    ABRs = [[], []]
    ABRs[0] = [[0] for _ in range(len(DataInfo['NoiseFrequency']))]
    ABRs[1] = [[0] for _ in range(len(DataInfo['NoiseFrequency']))]
    for Freq in range(len(DataInfo['NoiseFrequency'])):
        Key = str(DataInfo['NoiseFrequency'][Freq][0]) + '-' \
              + str(DataInfo['NoiseFrequency'][Freq][1])
        ABRs[0][Freq] = [{} for _ in range(len(DataInfo['SoundAmpF'][Key]))]
        ABRs[1][Freq] = [{} for _ in range(len(DataInfo['SoundAmpF'][Key]))]
    del(Freq)
    
    for RecFolder in Exps:
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
        
        ExpInfo = LoadHdf5Files.ExpExpInfo(FileName, DirList, RecFolder)
        
        print('Check if files are ok...')
        if 'Raw' not in locals():
            print('.kwd file is corrupted. Skipping dataset...')
            continue
        
        print('Data from ', RecFolder, ' loaded.')
        
        if '0' not in list(Raw['data'].keys()):
            print('Rec numbers are wrong. Fixing...')
            for iKey in Raw.keys():
                Recs = list(Raw[iKey].keys())
                Recs = [int(_) for _ in Recs]; Min = min(Recs)
                
                for Key in Recs:
                    Raw[iKey][str(Key-Min)] = Raw[iKey].pop(str(Key))
                
            print('Fixed.')
        
        Rate = Raw['info']['0']['sample_rate']

        NoOfSamplesBefore = ABRTimeBeforeTTL*int(Rate*10**-3)
        NoOfSamplesAfter = ABRTimeAfterTTL*int(Rate*10**-3)
        NoOfSamples = NoOfSamplesBefore + NoOfSamplesAfter
        XValues = list((range((NoOfSamplesBefore)*-1, 0)/Rate)*10**3) + \
                  list((range(NoOfSamplesAfter)/Rate)*10**3)
        
        Lf2, Lf1 = signal.butter(FilterOrder, FilterLow/(Rate/2), 'highpass')
        Hf2, Hf1 = signal.butter(FilterOrder, FilterHigh/(Rate/2), 'lowpass')
        
        for Rec in range(len(Raw['data'])):
            RecChNo = Raw['data'][str(Rec)].shape[1]
#            TTLCh = Raw['data'][str(Rec)][:, ABRTTLCh + (RecChNo-9)]
            TTLCh = Raw['data'][str(Rec)][:, -1]
            Threshold = max(TTLCh)/5
            TTLs = []
            for _ in range(1, len(TTLCh)):
                if TTLCh[_] > Threshold:
                    if TTLCh[_-1] < Threshold:
                        TTLs.append(_)
            
            rABR = [[0 for _ in range(NoOfSamples)] 
                    for _ in range(len(TTLs))]
            lABR = [[0 for _ in range(NoOfSamples)] 
                    for _ in range(len(TTLs))]
            sTTL = [[0 for _ in range(NoOfSamples)] 
                    for _ in range(len(TTLs))]
            
            print('Slicing and filtering ABRs Rec ', str(Rec), '...')
            for TTL in range(len(TTLs)):
                TTLLoc = int(TTLs[TTL])
                Start = TTLLoc-NoOfSamplesBefore
                End = TTLLoc+NoOfSamplesAfter
                
                rABR[TTL] = Raw['data'][str(Rec)][Start:End, ABRCh[0]-1] * \
                            Raw['channel_bit_volts'][str(Rec)][ABRCh[0]-1]
                
                lABR[TTL] = Raw['data'][str(Rec)][Start:End, ABRCh[1]-1] * \
                            Raw['channel_bit_volts'][str(Rec)][ABRCh[1]-1]
                
                rABR[TTL] = signal.filtfilt(Lf2, Lf1, rABR[TTL], padlen=0)
                lABR[TTL] = signal.filtfilt(Lf2, Lf1, lABR[TTL], padlen=0)
                
                del(TTLLoc, Start, End)
            
            # Mean
            rABR = np.mean(rABR, axis=0)
            lABR = np.mean(lABR, axis=0)
#            rData = np.mean(rABR, axis=0)
#            lData = np.mean(lABR, axis=0)
            
            rABR = signal.filtfilt(Hf2, Hf1, rABR, padlen=0)
            lABR = signal.filtfilt(Hf2, Hf1, lABR, padlen=0)
#            rData = signal.filtfilt(Hf2, Hf1, rData, padlen=0)
#            lData = signal.filtfilt(Hf2, Hf1, lData, padlen=0)
            
            if ExpInfo['DVCoord'] not in ABRs[0][ExpInfo['Hz']][Rec]:
                ABRs[0][ExpInfo['Hz']][Rec][ExpInfo['DVCoord']] = [rABR]
                ABRs[1][ExpInfo['Hz']][Rec][ExpInfo['DVCoord']] = [lABR]
            else:
                ABRs[0][ExpInfo['Hz']][Rec][ExpInfo['DVCoord']].append(rABR)
                ABRs[1][ExpInfo['Hz']][Rec][ExpInfo['DVCoord']].append(lABR)
    
    if 'XValues' in locals():
        print('Saving data to ' + FileName)
        GroupName = 'ABRs-' + datetime.datetime.now().strftime("%Y%m%d%H%M%S")
        with h5py.File(FileName) as F:
            if 'ABRs' in F.keys():
                F[GroupName] = F['ABRs']; del(F['ABRs'])
            
            F.create_group('ABRs')
            F['ABRs'].attrs['XValues'] = XValues
            F['ABRs'].create_group('Right'); F['ABRs'].create_group('Left')
                
            for Freq in range(len(ABRs[0])):
                for AmpF in range(len(ABRs[0][Freq])):
                    for DV in ABRs[0][Freq][AmpF].keys():
                        Path = str(Freq) + '/' + \
                               str(AmpF) + '/' + \
                               DV
                        
                        F['ABRs']['Right'].create_group(Path)
                        F['ABRs']['Left'].create_group(Path)
                        del(Path)
                        
                        for Trial in range(len(ABRs[0][Freq][AmpF][DV])):
                            F['ABRs']['Right'][str(Freq)][str(AmpF)][DV][str(Trial)] = \
                                ABRs[0][Freq][AmpF][DV][Trial][:]
                            F['ABRs']['Left'][str(Freq)][str(AmpF)][DV][str(Trial)] = \
                                ABRs[1][Freq][AmpF][DV][Trial][:]
    
    print('Done.')
    return(None)


def PlotABR(FileName):
    """ 
    This function will plot the data from *Stim.hdf5. Make sure FileName is a 
    string with the path to only one file.
    
    Also, LaTeX will render all the text in the plots. For this, make sure you 
    have a working LaTex installation and dvipng package installed.
    """
    
    print('Loading data...')
    DataInfo = {}
    with h5py.File(FileName) as F:
        for Key, Value in F['DataInfo'].items():
            DataInfo['SoundAmpF'] = {}
            for aKey, aValue in F['DataInfo']['SoundAmpF'].items():
                DataInfo['SoundAmpF'][aKey] = aValue[:]
        
        for bKey, bValue in F['DataInfo'].attrs.items():
            if isinstance(bValue, Number):
                DataInfo[bKey] = float(bValue)
            else:
                DataInfo[bKey] = bValue
        
        ABRs = [[], []]
        ABRs[0] = [0]*len(F['ABRs']['Right'])
        ABRs[1] = [0]*len(F['ABRs']['Left'])
        for Freq in range(len(F['ABRs']['Right'])):
            ABRs[0][Freq] = [0]*len(F['ABRs']['Right'][str(Freq)])
            ABRs[1][Freq] = [0]*len(F['ABRs']['Left'][str(Freq)])
            
            for AmpF in range(len(F['ABRs']['Right'][str(Freq)])):
                ABRs[0][Freq][AmpF] = {}; ABRs[1][Freq][AmpF] = {}
                
                for DV in F['ABRs']['Right'][str(Freq)][str(AmpF)].keys():
                    ABRs[0][Freq][AmpF][DV] = [0]*len(F['ABRs']['Right'][str(Freq)][str(AmpF)][DV])
                    ABRs[1][Freq][AmpF][DV] = [0]*len(F['ABRs']['Left'][str(Freq)][str(AmpF)][DV])
                    
                    for Trial in range(len(F['ABRs']['Right'][str(Freq)][str(AmpF)][DV])):
                        ABRs[0][Freq][AmpF][DV][Trial] = F['ABRs']['Right'][str(Freq)][str(AmpF)][DV][str(Trial)][:]
                        ABRs[1][Freq][AmpF][DV][Trial] = F['ABRs']['Left'][str(Freq)][str(AmpF)][DV][str(Trial)][:]
        
        XValues = F['ABRs'].attrs['XValues'][:]
    
    SoundIntensity = LoadHdf5Files.SoundMeasurement(DataInfo['CalibrationFile'], 'SoundIntensity')
    
    SetLaTexPlot()
    print('Plotting...')
    Colormaps = [plt.get_cmap('Reds'), plt.get_cmap('Blues')]
    Colors = [[Colormaps[0](255-(_*20)), Colormaps[1](255-(_*20))] 
              for _ in range(len(ABRs[0][0]))]
    
    Keys = list(ABRs[0][0][0].keys())
    for Key in Keys:
        for Trial in range(len(ABRs[0][0][0][Key])):
            Fig, Axes = plt.subplots(len(DataInfo['NoiseFrequency']), 2, 
                                         sharex=True, figsize=(12, 12))
    
            for Ear in range(2):
                for Freq in range(len(DataInfo['NoiseFrequency'])):
                    KeyHz = str(DataInfo['NoiseFrequency'][Freq][0]) + '-' \
                            + str(DataInfo['NoiseFrequency'][Freq][1])
                    
                    for AmpF in range(len(DataInfo['SoundAmpF'][KeyHz])):
                        FigTitle = Key + ' DV, trial ' + str(Trial+1)
                        AxTitle = KeyHz
                        YLabel = 'Voltage [\si{\micro}V]'
                        XLabel = 'Time [ms]'
                        if 0.0 in DataInfo['SoundAmpF'][KeyHz]:
                            DataInfo['SoundAmpF'][KeyHz][
                                DataInfo['SoundAmpF'][KeyHz].index(0.0)
                                                        ] = 0
                        
                        AmpStr = str(DataInfo['SoundAmpF'][KeyHz][AmpF])
                        LineLabel = str(round(SoundIntensity[KeyHz][AmpStr]
                                              )) + ' dB'
                        
                        if Ear == 0:
                            Axes[Freq][Ear].plot(
                                XValues, ABRs[Ear][Freq][AmpF][Key][Trial], 
                                color=Colors[AmpF][Ear], 
                                label=LineLabel)
                        else:
                            Axes[Freq][Ear].plot(
                                XValues, ABRs[Ear][Freq][AmpF][Key][Trial], 
                                color=Colors[AmpF][Ear], 
                                label=LineLabel)
                        
#                        Axes[Freq][Ear].legend(loc='lower right')#, frameon=False)
                        Axes[Freq][Ear].spines['right'].set_visible(False)
                        Axes[Freq][Ear].spines['top'].set_visible(False)
                        Axes[Freq][Ear].spines['left'].set_visible(False)
                        Axes[Freq][Ear].yaxis.set_ticks_position('none')
                        Axes[Freq][Ear].xaxis.set_ticks_position('bottom')
                        Axes[Freq][Ear].set_title(AxTitle)
                        Axes[Freq][Ear].set_ylabel(YLabel)
                        Axes[Freq][Ear].locator_params(tight=True)
            Axes[-1][0].set_xlabel(XLabel)
            Axes[-1][1].set_xlabel(XLabel)
            Fig.suptitle(FigTitle)
            Fig.tight_layout()
            Fig.subplots_adjust(top=0.95)
            Fig.savefig('Figs/' + FileName + '-DV' +  Key + '-Trial' + 
                        str(Trial) + '.svg', format='svg')


def PlotABR2(FileName):
    """ 
    This function will plot the data from *Stim.hdf5. Make sure FileName is a 
    string with the path to only one file.
    
    Also, LaTeX will render all the text in the plots. For this, make sure you 
    have a working LaTex installation and dvipng package installed.
    """
    
    print('Loading data...')
    DataInfo = {}
    with h5py.File(FileName) as F:
        for Key, Value in F['DataInfo'].items():
            DataInfo['SoundAmpF'] = {}
            for aKey, aValue in F['DataInfo']['SoundAmpF'].items():
                DataInfo['SoundAmpF'][aKey] = aValue[:]
        
        for bKey, bValue in F['DataInfo'].attrs.items():
            if isinstance(bValue, Number):
                DataInfo[bKey] = float(bValue)
            else:
                DataInfo[bKey] = bValue
        
        ABRs = [[], []]
        ABRs[0] = [0]*len(F['ABRs']['Right'])
        ABRs[1] = [0]*len(F['ABRs']['Left'])
        for Freq in range(len(F['ABRs']['Right'])):
            ABRs[0][Freq] = [0]*len(F['ABRs']['Right'][str(Freq)])
            ABRs[1][Freq] = [0]*len(F['ABRs']['Left'][str(Freq)])
            
            for AmpF in range(len(F['ABRs']['Right'][str(Freq)])):
                ABRs[0][Freq][AmpF] = {}; ABRs[1][Freq][AmpF] = {}
                
                for DV in F['ABRs']['Right'][str(Freq)][str(AmpF)].keys():
                    ABRs[0][Freq][AmpF][DV] = [0]*len(F['ABRs']['Right'][str(Freq)][str(AmpF)][DV])
                    ABRs[1][Freq][AmpF][DV] = [0]*len(F['ABRs']['Left'][str(Freq)][str(AmpF)][DV])
                    
                    for Trial in range(len(F['ABRs']['Right'][str(Freq)][str(AmpF)][DV])):
                        ABRs[0][Freq][AmpF][DV][Trial] = F['ABRs']['Right'][str(Freq)][str(AmpF)][DV][str(Trial)][:]
                        ABRs[1][Freq][AmpF][DV][Trial] = F['ABRs']['Left'][str(Freq)][str(AmpF)][DV][str(Trial)][:]
        
        XValues = F['ABRs'].attrs['XValues'][:]
    
    SetLaTexPlot()
    print('Plotting...')
    Colormaps = [plt.get_cmap('Reds'), plt.get_cmap('Blues')]
    Colors = [[Colormaps[0](255-(_*20)), Colormaps[1](255-(_*20))] 
              for _ in range(len(ABRs[0][0]))]
    
    Keys = list(ABRs[0][0][0].keys())
    for Key in Keys:
        for Trial in range(len(ABRs[0][0][0][Key])):
            for Freq in range(len(ABRs[0])):
                Fig, Axes = plt.subplots(len(ABRs[0][0]), 
                                     sharex=True, figsize=(8, 12))
                
                KeyHz = str(DataInfo['NoiseFrequency'][Freq][0]) + '-' \
                        + str(DataInfo['NoiseFrequency'][Freq][1])
                
                if 0.0 in DataInfo['SoundAmpF'][KeyHz]:
                    DataInfo['SoundAmpF'][KeyHz][
                        DataInfo['SoundAmpF'][KeyHz].index(0.0)
                                                ] = 0
                upYLim = max(ABRs[0][Freq][0][Key][Trial])
                downYLim = min(ABRs[0][Freq][0][Key][Trial])
                for AmpF in range(len(DataInfo['SoundAmpF'][KeyHz])):
                    FigTitle = KeyHz + ' Hz, trial ' + str(Trial+1)
                    YLabel = 'Voltage [mV]'
                    XLabel = 'Time [ms]'
                    LineLabel = str(DataInfo['Intensities'][AmpF]) + ' dB'
                    SpanLabel = 'Sound pulse'
                    
                    Ind1 = list(XValues).index(0)
                    Ind2 = list(XValues).index(3)
                    
                    Axes[AmpF].plot(XValues, ABRs[0][Freq][AmpF][Key][Trial], 
                                    color=Colors[AmpF][0], label=LineLabel)
                    
                    Axes[AmpF].axvspan(XValues[Ind1], XValues[Ind2], 
                                       color='k', alpha=0.3, lw=0, 
                                       label=SpanLabel)
                    
                    Axes[AmpF].legend(loc='lower right', frameon=False)
                    Axes[AmpF].spines['right'].set_visible(False)
                    Axes[AmpF].spines['top'].set_visible(False)
                    Axes[AmpF].spines['bottom'].set_visible(False)
                    Axes[AmpF].spines['left'].set_bounds(round(0), round(1))
                    Axes[AmpF].yaxis.set_ticks_position('left')
                    Axes[AmpF].xaxis.set_ticks_position('none')
                    Axes[AmpF].set_ylabel(YLabel)
                    Axes[AmpF].set_ylim(downYLim, upYLim)
                    Axes[AmpF].locator_params(tight=True)
                    
                Axes[-1].spines['bottom'].set_visible(True)
                Axes[-1].set_xlabel(XLabel)
                Axes[-1].spines['bottom'].set_bounds(round(0), round(1))
                Fig.suptitle(FigTitle)
                Fig.tight_layout()
                Fig.subplots_adjust(top=0.95)
                
                FigName = 'Figs/' + FileName[:-15] + '-ABR_DV' +  Key + \
                          '_Trial' + str(Trial) + '_Freq' + KeyHz + '.svg'
                Fig.savefig(FigName, format='svg')


def GPIAS(FileName, GPIASTimeBeforeTTL=50, GPIASTimeAfterTTL=150, FilterLow=3, 
          FilterHigh=300, FilterOrder=4, GPIASTTLCh=1, PiezoCh=1):
    
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
        
        DataInfo = LoadHdf5Files.GPIASDataInfo(FileName, DirList)
        
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


def GPIASAnalogTTLs(RecFolder, FileName, GPIASTimeBeforeTTL=50, 
                    GPIASTimeAfterTTL=150, FilterLow=3, FilterHigh=300, 
                    FilterOrder=4, GPIASTTLCh=1, PiezoCh=1):
    
    print('set paths...')
    os.makedirs('Figs', exist_ok=True)    # Figs folder
    DirList = glob.glob('KwikFiles/*'); DirList.sort()
    RecFolder = DirList[RecFolder-1]
    
    print('Load files...')
    FilesList = glob.glob(''.join([RecFolder, '/*']))
    Files = {}
    for File in FilesList:
        if '.kwd' in File:
            try:
                Raw = Kwik.load(File, 'all')
                Files['kwd'] = File
            except OSError:
                    print('File', File, "is corrupted :'(")
            except KeyError: 
                # Old OE versions do not have channel_bit_volts
                print('No channel_bit_volts in the file!')
                with h5py.File(File, 'r') as F:
                    Raw = {}
                    Raw['info'] = {Rec: F['recordings'][Rec].attrs 
                                   for Rec in F['recordings'].keys()}
                    Raw['data'] = {Rec: F['recordings'][Rec]['data']
                                   for Rec in F['recordings'].keys()}
                    Raw['timestamps'] = {Rec: ((
                                    np.arange(0,Raw['data'][Rec].shape[0])
                                    + Raw['info'][Rec]['start_time'])
                                   / Raw['info'][Rec]['sample_rate'])
                                         for Rec in F['recordings']}
            
        elif '.kwe' in File:
            try:
                Events = Kwik.load(File)
                Files['kwe'] = File
            except OSError:
                print('File', File, "is corrupted :'(")
            
        elif '.kwik' in File:
            try:
                Events = Kwik.load(File)
                Files['kwik'] = File
            except OSError:
                print('File', File, "is corrupted :'(")
        
    
    DataInfo = LoadHdf5Files.GPIASDataInfo(FileName)
    
    print('Check if files are ok...')
    if 'Raw' not in locals():
        print('.kwd file is corrupted. Skipping dataset...')
        raise SystemExit
    
    if 'Events' not in locals():
        print('.kwe/.kwik file is corrupted. Skipping dataset...')
        raise SystemExit
    
    if 'DataInfo' not in locals():
        print('No data info. Skipping dataset...')
        raise SystemExit
    
    print('Data from ', RecFolder, ' loaded.')
    
    print('Preallocate memory...')
    GPIAS = [[0] for _ in range(len(DataInfo['NoiseFrequency']))]
    for Freq in range(len(DataInfo['NoiseFrequency'])):
        GPIAS[Freq] = [[0] for _ in range(round(DataInfo['NoOfTrials']*2))]
    
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
        TTLCh = Raw['data'][str(Rec)][:, -1]
        if len(TTLCh) == 0:
            print('Bad recording. Skipping...')
            continue
        
        Threshold = max(TTLCh)/5
        TTLs = []
        for _ in range(1, len(TTLCh)):
            if TTLCh[_] > Threshold:
                if TTLCh[_-1] < Threshold:
                    TTLs.append(_)
        
        print('Slicing and filtering Rec ', str(Rec), '...')
        for TTL in range(len(TTLs)):
            TTLLoc = int(TTLs[TTL])
            Start = TTLLoc-NoOfSamplesBefore
            End = TTLLoc+NoOfSamplesAfter
            
            gData = Raw['data'][str(Rec)][Start:End, PiezoCh-1]
#                gData = [float(gData[_]) for _ in range(NoOfSamples)]
            gData = abs(signal.hilbert(gData))
            gData = signal.filtfilt(f2, f1, gData, padtype='odd', padlen=0)
            
            Freq = DataInfo['FreqOrder'][Rec][0]; 
            Trial = DataInfo['FreqOrder'][Rec][1];
            GPIAS[Freq][Trial] = gData[:]
            
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
        
        del(NoGapAll, GapAll, NoGapSum, GapSum, gData)
    
    if 'XValues' in locals():
        print('Saving data to ' + FileName)
        GroupName = 'GPIAS-' + RecFolder[10:] + '-' + \
                    datetime.datetime.now().strftime("%Y%m%d%H%M%S")
        with h5py.File(FileName) as F:
            F.create_group(GroupName)
            F[GroupName].attrs['XValues'] = XValues

            for Freq in range(len(GPIAS)):
                F[GroupName].create_group(str(Freq))
                
                F[GroupName][str(Freq)]['NoGap'] = GPIAS[Freq][0]
                F[GroupName][str(Freq)]['Gap'] = GPIAS[Freq][1]
        
        print('Done.')
    
    print('All data saved.')
    return(None)


def PlotGPIAS(FileList):
    for File in FileList:
        print('Loading data from GPIAS.db ...')
#        with shelve.open(File[:-3]) as Shelve:
#            GPIAS = Shelve['GPIAS']
#            AllTTLs = Shelve['AllTTLs']
#            XValues = Shelve['XValues']
#            DataInfo = Shelve['DataInfo']
        
        print('Plotting...')
        for Freq in range(len(DataInfo['NoiseFrequency'])):
            FigTitle = str(DataInfo['NoiseFrequency'][Freq]) + '\ Hz'
            Line0Label = 'No\ Gap'; Line1Label = 'Gap'
            SpanLabel = 'Sound\ Pulse'
            XLabel = 'time\ [ms]'; YLabel = 'voltage\ [mV]'
            
            plt.figure(Freq)
            plt.plot(XValues, GPIAS[Freq][0], 
                     color='r', label='$'+Line0Label+'$')
            plt.plot(XValues, GPIAS[Freq][1], 
                     color='b', label='$'+Line1Label+'$')
            plt.axvspan(XValues[AllTTLs[Freq][0][0]], XValues[AllTTLs[Freq][0][1]], 
                        color='k', alpha=0.5, lw=0, label='$'+SpanLabel+'$')
#            plt.axvspan(XValues[AllTTLs[Freq][1][0]], XValues[AllTTLs[Freq][1][1]], 
#                        color='b', alpha=0.5, lw=0, label='Sound pulse (Gap)')
            plt.suptitle('$'+FigTitle+'$')
            plt.ylabel('$'+YLabel+'$'); plt.xlabel('$'+XLabel+'$')
            plt.legend(loc='lower right')
            plt.locator_params(tight=True)
            plt.axes().spines['right'].set_visible(False)
            plt.axes().spines['top'].set_visible(False)
            plt.axes().yaxis.set_ticks_position('left')
            plt.axes().xaxis.set_ticks_position('bottom')
            plt.savefig('Figs/' + File[:-3] + '-' + 
                        str(DataInfo['NoiseFrequency'][Freq][0]) + '_' + 
                        str(DataInfo['NoiseFrequency'][Freq][1]) + '.svg', 
                        format='svg')
        print('Done.')


def PlotGPIAS2(FileName):
    print('Loading data...')
    DataInfo = LoadHdf5Files.GPIASDataInfo(FileName)
    with h5py.File(FileName) as F:
        Keys = list(F.keys())
    
    Keys.remove('DataInfo')
    
    SetLaTexPlot()
    
#    for Key in Keys:
    Key = Keys[-1]
    with h5py.File(FileName) as F:
        XValues = F[Key].attrs['XValues']
        
        GPIAS = [[0] for _ in range(len(DataInfo['NoiseFrequency']))]
        for Freq in range(len(GPIAS)):
            GPIAS[Freq] = [[], []]
            GPIAS[Freq][0] = F[Key][str(Freq)]['NoGap'][:]
            GPIAS[Freq][1] = F[Key][str(Freq)]['Gap'][:]
    
    print('Plotting...')
    Ind1 = list(XValues).index(0)
    Ind2 = list(XValues).index(int(DataInfo['SoundLoudPulseDur']*1000))
    
    for Freq in range(len(DataInfo['NoiseFrequency'])):
        FigTitle = str(DataInfo['NoiseFrequency'][Freq]) + ' Hz'
        Line0Label = 'No Gap'; Line1Label = 'Gap'
        SpanLabel = 'Sound Pulse'
        XLabel = 'time [ms]'; YLabel = 'voltage [mV]'
        
        plt.figure(Freq)
        plt.plot(XValues, GPIAS[Freq][0], 
                 color='r', label=Line0Label)
        plt.plot(XValues, GPIAS[Freq][1], 
                 color='b', label=Line1Label)
        plt.axvspan(XValues[Ind1], XValues[Ind2], color='k', alpha=0.5, 
                    lw=0, label=SpanLabel)

        plt.suptitle(FigTitle)
        plt.ylabel(YLabel); plt.xlabel(XLabel)
        plt.legend(loc='lower right')
        plt.locator_params(tight=True)
        plt.axes().spines['right'].set_visible(False)
        plt.axes().spines['top'].set_visible(False)
        plt.axes().yaxis.set_ticks_position('left')
        plt.axes().xaxis.set_ticks_position('bottom')
        plt.savefig('Figs/' + FileName[:-9] + '-' + 
                    str(DataInfo['NoiseFrequency'][Freq][0]) + '_' + 
                    str(DataInfo['NoiseFrequency'][Freq][1]) + '.svg', 
                    format='svg')
    print('Done.')


def TTLsLatency(FileName, SoundCh=1, SoundSqCh=2, SoundTTLCh=1, 
                TimeBeforeTTL=5, TimeAfterTTL=8):
    print('set paths...')
    os.makedirs('Figs', exist_ok=True)    # Figs folder
    RecFolder = glob.glob('KwikFiles/*'); RecFolder = RecFolder[0]
    SoundCh = 0; SoundSqCh = 1
    
    print('Load DataInfo...')
    DataInfo= LoadHdf5Files.ExpDataInfo(FileName, list(RecFolder), 'Sound')
    
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
    if 'Raw' not in locals():
        print('.kwd file is corrupted  >:(')
        raise SystemExit
    
    if 'Events' not in locals():
        print('.kwd file is corrupted  >:(')
        raise SystemExit
    
    print('Data from ', RecFolder, ' loaded.')
    
    if '0' not in list(Raw['data'].keys()):
        print('Rec numbers are wrong. Fixing...')
        for iKey in Raw.keys():
            Recs = list(Raw[iKey].keys())
            Recs = [int(_) for _ in Recs]; Min = min(Recs)
            
            EventRec = Events['TTLs']['recording'][:]
            for _ in range(len(EventRec)): EventRec[_] = EventRec[_] - Min
            
            for Key in Recs:
                Raw[iKey][str(Key-Min)] = Raw[iKey].pop(str(Key))
            
        print('Fixed.')
    else:
        EventRec = Events['TTLs']['recording']
    
    print('Get TTL data...')
    EventID = Events['TTLs']['user_data']['eventID']
    EventCh = Events['TTLs']['user_data']['event_channels']
    EventSample = Events['TTLs']['time_samples']

    TTLChs = np.nonzero(np.bincount(EventCh))[0]
    TTLRecs = np.nonzero(np.bincount(EventRec))[0]
    TTLsPerRec = {_Rec: [EventSample[_] for _ in range(len(EventRec)) 
                         if EventRec[_] == _Rec 
                         and EventCh[_] == SoundTTLCh-1 
                         and EventID[_] == 1]
                  for _Rec in TTLRecs}
    TTLRising = Kwik.get_rising_edge_times(Files['kwe'], SoundTTLCh-1)
    
    Rate = Raw['info']['0']['sample_rate']
    NoOfSamplesBefore = TimeBeforeTTL*int(Rate*10**-3)
    NoOfSamplesAfter = TimeAfterTTL*int(Rate*10**-3)
    NoOfSamples = NoOfSamplesBefore + NoOfSamplesAfter
    
    XValues = list((range((NoOfSamplesBefore)*-1, 0)/Rate)*10**3) + \
              list((range(NoOfSamplesAfter)/Rate)*10**3)
    
    TTLNo = [0]
    for _ in range(1, len(TTLsPerRec)+1):
        TTLNo = TTLNo + [len(TTLsPerRec[_-1]) + TTLNo[-1]]
    TTLNo = [0] + [TTLNo[_]-1 for _ in range(1, len(TTLNo))]
    
    RawTime = list(Raw['timestamps']['0'])
    
    SoundPulse = [[0 for _ in range(NoOfSamples)] 
                  for _ in range(len(TTLsPerRec[0]))]
    SoundSq = SoundPulse[:]
    
    print('Slicing data...')
    for TTL in range(len(TTLsPerRec[0])):
        TTLLoc = RawTime.index(TTLRising[TTL])
        Start = TTLLoc-NoOfSamplesBefore
        End = TTLLoc+NoOfSamplesAfter
        try:
            SoundPulse[TTL] = Raw['data']['0'][Start:End, SoundCh]
            SoundSq[TTL] = Raw['data']['0'][Start:End, SoundSqCh]
        except ValueError:
            print('ValueError: Slice outside data limits. Skipping', 
                  str(TTL), '...')
            continue
        
        del(TTLLoc, Start, End)
    
    SoundSqDelay = [[0] for _ in range(len(SoundSq))]    
    for TTL in range(len(SoundSq)):
        Threshold = max(SoundSq[TTL])/5
        Peaks = []
        for _ in range(len(SoundSq[TTL])):
            if SoundSq[TTL][_] > Threshold:
                if SoundSq[TTL][_-1] < Threshold:
                    Peaks.append(_)
        
        SoundSqDelay[TTL] = ((Peaks[0] - NoOfSamplesBefore) / (Rate/1000))*-1
        del(Peaks)
    
    print('Saving data to ' + FileName + ' ...')
    GroupName = 'Analysis-' + datetime.datetime.now().strftime("%Y%m%d%H%M%S")
    with h5py.File(FileName) as F:
        if 'Analysis' in F.keys():
            F[GroupName] = F['Analysis']; del(F['Analysis'])
        
        F.create_group('Analysis')
        F['Analysis']['SoundPulse'] = SoundPulse
        F['Analysis']['SoundSq'] = SoundSq
        F['Analysis']['SoundSqDelay'] = SoundSqDelay
        F['Analysis']['XValues'] = XValues
    
    print('Done.')
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
    
    SetLaTexPlot()
    
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
