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

Functions for manipulating specific .mat files.
"""
import numpy as np
import os

from DataAnalysis.DataAnalysis import FilterSignal
from IO import Hdf5F
from glob import glob
from scipy import io, signal


def GPIASAnalysis(RecFolderNo, GPIASTimeBeforeTTL=20, GPIASTimeAfterTTL=100, 
                     FilterFreq=[70, 400], FilterOrder=3, Override={}):
    if 'DirList' in Override: DirList = Override['DirList']
    else: DirList = glob('MatFiles/*'); DirList.sort()
    
    if 'FileName' in Override: FileName =  Override['FileName']
    else: 
        FileName = glob('*.mat')[0]
    
    DataInfo = io.loadmat(FileName)
    DataInfo['AnimalName'] = os.getcwd().split('/')[-2]
    
    if 'AnalysisFile' in Override: AnalysisFile =  Override['AnalysisFile']
    else: AnalysisFile = '../' + DataInfo['AnimalName'] + '-Analysis.hdf5'
    
    RecFolder = DirList[RecFolderNo-1]
    
    Rate = np.array(int(DataInfo['Rate'])); PulseStart = DataInfo['PulseStart']
    Freqs = glob(RecFolder + '/*.mat')
    GPIAS = {}
    
    for Freq in Freqs:
        SFreq = Freq[-9:-4]; SFreq = SFreq.split('_')
        try:
            SFreq = str(int(SFreq[0])*1000) + '-' + str(int(SFreq[1])*1000)
        except ValueError:
            SFreq = str(int(SFreq[0][1:])*1000) + '-' + str(int(SFreq[1])*1000)
        
        if SFreq not in GPIAS.keys(): GPIAS[SFreq] = {}
        
        NoOfSamplesBefore = int(round((GPIASTimeBeforeTTL*Rate)*10**-3))
        NoOfSamplesAfter = int(round((GPIASTimeAfterTTL*Rate)*10**-3))
        NoOfSamples = NoOfSamplesBefore + NoOfSamplesAfter
        
        XValues = (range(-NoOfSamplesBefore, NoOfSamples-NoOfSamplesBefore)
                   /Rate)*10**3
        
        Start = int(PulseStart-NoOfSamplesBefore)
        End = int(PulseStart+NoOfSamplesAfter)
        
        Data = io.loadmat(Freq)                
        Data['Gap'] = Data['Gap'][0,Start:End] * 1000 # in mV
        Data['NoGap'] = Data['NoGap'][0,Start:End] * 1000 # in mV
        
        Data['Gap'] = FilterSignal(Data['Gap'], Rate, 
                                                FilterFreq, 
                                                FilterOrder, 
                                                'bandpass')
        Data['NoGap'] = FilterSignal(Data['NoGap'], Rate, 
                                                  FilterFreq, 
                                                  FilterOrder, 
                                                  'bandpass')
        
         # Amplitude envelope
        GapAE = abs(signal.hilbert(Data['Gap']))
        NoGapAE = abs(signal.hilbert(Data['NoGap']))
        
        # RMS
        Half = len(GapAE)//2
        BGStart = 0; BGEnd = Half - 1
        PulseStart = Half; PulseEnd = len(GapAE) - 1
        BinSize = XValues[-1] - XValues[-2]
        
        GapRMSBG = sum(GapAE[BGStart:BGEnd] * BinSize)**0.5
        GapRMSPulse = sum(GapAE[PulseStart:PulseEnd] * BinSize)**0.5
        GapRMS = GapRMSPulse - GapRMSBG
        
        NoGapRMSBG = sum(NoGapAE[BGStart:BGEnd] * BinSize)**0.5
        NoGapRMSPulse = sum(NoGapAE[PulseStart:PulseEnd] * BinSize)**0.5
        NoGapRMS = NoGapRMSPulse - NoGapRMSBG
        
        # GPIAS index (How much Gap is different from NoGap)
        Data['GPIASIndex'] = (NoGapRMS-GapRMS)/NoGapRMS
        
        GPIAS[SFreq]['Gap'] = Data['Gap'][:]
        GPIAS[SFreq]['NoGap'] = Data['NoGap'][:]
        GPIAS[SFreq]['GPIASIndex']= Data['GPIASIndex']
        del(Data)
    
    if 'DirList' in Override:
        RecExp = RecFolder.split('/')[1]
        Hdf5F.GPIASWrite(GPIAS, XValues, RecFolder, AnalysisFile, RecExp)
    else:
        Hdf5F.GPIASWrite(GPIAS, XValues, RecFolder, AnalysisFile)


def WriteDataToMMSS(FileName, StimType=['Sound'], Override={}):
    DirList = glob('KwikFiles/*'); DirList.sort()
    for Stim in StimType:
        if Override != {}: 
            if 'Stim' in Override.keys():
                Stim = Override['Stim']
        
        Exps = Hdf5F.ExpPerStimLoad(Stim, DirList, FileName)
        
        for FInd, RecFolder in enumerate(Exps): 
            Path = os.getcwd() + '/' + RecFolder
            ClusterFile = Path + '/SpkClusters.hdf5'
            Clusters = Hdf5F.ClustersLoad(ClusterFile)
            
            Path = os.getcwd() + '/' + RecFolder +'/SepRec/'
            os.makedirs(Path, exist_ok=True)
            
            for Rec in Clusters.keys():
                RecS = "{0:02d}".format(int(Rec))
                        
                print('Writing files for clustering... ', end='')
                WF = []; ST = []; ChList = []
                for Ch in Clusters[RecS].keys():
                    WF.append(Clusters[RecS][Ch]['Spikes'][:])
                    ST.append(Clusters[RecS][Ch]['Timestamps'][:])
                    ChList.append(Ch)
                
                data = {'waveforms': np.array(WF, dtype='object'),
                        'spiketimes': np.array(ST, dtype='object'),
                        'ChList': np.string_(ChList)}
                
                MatName = ''.join(['Exp_', RecS, '.mat'])
                io.savemat(Path+MatName, {'data': data})
                print('Done.')
    return(None)

