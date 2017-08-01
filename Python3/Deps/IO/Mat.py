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
from IO import Hdf5
from glob import glob
from scipy import io, signal

def GPIASAnalysis(Folders, InfoFiles, AnalysisFile, GPIASTimeBeforeTTL=100, GPIASTimeAfterTTL=100, 
                     FilterFreq=[70, 400], FilterOrder=3, Filter = 'butter', SliceSize=100):
    
    SliceSizeMS = [SliceSize][0]
    for F, Folder in enumerate(Folders):
        DataInfo = io.loadmat(InfoFiles[F])
        DataInfo['AnimalName'] = Folder.split('_')[-1]
        
        Rate = np.array(int(DataInfo['Rate'])); TTL = DataInfo['PulseStart']
        Freqs = glob(Folder + '/*.mat')
        GPIAS = {'Trace': {}, 'Index':{}}
        
        SliceSize = int(SliceSizeMS * (Rate/1000))
        
        for Freq in Freqs:
            SFreq = Freq[-9:-4]; SFreq = SFreq.split('_')
            try:
                SFreq = str(int(SFreq[0])*1000) + '-' + str(int(SFreq[1])*1000)
            except ValueError:
                SFreq = str(int(SFreq[0][1:])*1000) + '-' + str(int(SFreq[1])*1000)
            
            for Key in GPIAS.keys():
                if SFreq not in GPIAS[Key].keys(): GPIAS[Key][SFreq] = {}
            
            NoOfSamplesBefore = int(round((GPIASTimeBeforeTTL*Rate)*10**-3))
            NoOfSamplesAfter = int(round((GPIASTimeAfterTTL*Rate)*10**-3))
            NoOfSamples = NoOfSamplesBefore + NoOfSamplesAfter
            
            XValues = (range(-NoOfSamplesBefore, NoOfSamples-NoOfSamplesBefore)
                       /Rate)*10**3
            
            Start = int(TTL-NoOfSamplesBefore)
            End = int(TTL+NoOfSamplesAfter)
            
            Data = io.loadmat(Freq)                
            Data['Gap'] = Data['Gap'][0,Start:End] * 1000 # in mV
            Data['NoGap'] = Data['NoGap'][0,Start:End] * 1000 # in mV
            
            Data['Gap'] = FilterSignal(Data['Gap'], Rate, 
                                                    FilterFreq, 
                                                    FilterOrder, Filter,
                                                    'bandpass')
            Data['NoGap'] = FilterSignal(Data['NoGap'], Rate, 
                                                      FilterFreq, 
                                                      FilterOrder, Filter,
                                                      'bandpass')
            
             # Amplitude envelope
            GapAE = abs(signal.hilbert(Data['Gap']))
            NoGapAE = abs(signal.hilbert(Data['NoGap']))
            
            # RMS
            BGStart = 0; BGEnd = SliceSize
            PulseStart = int(GPIASTimeBeforeTTL*(Rate/1000)); PulseEnd = PulseStart + SliceSize
            
            GapRMSBG = (np.mean(GapAE[BGStart:BGEnd]**2))**0.5
            GapRMSPulse = (np.mean(GapAE[PulseStart:PulseEnd]**2))**0.5
            if GapRMSPulse < GapRMSBG: GapRMS = GapRMSPulse
            else: GapRMS = GapRMSPulse - GapRMSBG
            
            NoGapRMSBG = (np.mean(NoGapAE[BGStart:BGEnd]**2))**0.5
            NoGapRMSPulse = (np.mean(NoGapAE[PulseStart:PulseEnd]**2))**0.5
            if NoGapRMSPulse < NoGapRMSBG: NoGapRMS = NoGapRMSPulse
            else: NoGapRMS = NoGapRMSPulse - NoGapRMSBG
            
            # GPIAS index (How much Gap is different from NoGap)
            Data['GPIASIndex'] = (NoGapRMS-GapRMS)/NoGapRMS
            
            GPIAS['Trace'][SFreq]['Gap'] = Data['Gap'][:]
            GPIAS['Trace'][SFreq]['NoGap'] = Data['NoGap'][:]
            GPIAS['Index'][SFreq]['GPIASIndex']= Data['GPIASIndex']
            del(Data)
        
        AnalysisKey = InfoFiles[F].split('/')[1].split('-')
        AnalysisKey[0] = AnalysisKey[0]+'000000'
        AnalysisKey[1] = DataInfo['AnimalName']
        AnalysisKey = '-'.join(AnalysisKey) + '/Sound'
        Hdf5.DataWrite({'GPIAS': GPIAS, 'XValues': XValues}, AnalysisKey, AnalysisFile)


def WriteDataToMMSS(FileName, StimType=['Sound'], Override={}):
    DirList = glob('KwikFiles/*'); DirList.sort()
    for Stim in StimType:
        if Override != {}: 
            if 'Stim' in Override.keys():
                Stim = Override['Stim']
        
        Exps = Hdf5.ExpPerStimLoad(Stim, DirList, FileName)
        
        for FInd, RecFolder in enumerate(Exps): 
            Path = os.getcwd() + '/' + RecFolder
            ClusterFile = Path + '/SpkClusters.hdf5'
            Clusters = Hdf5.ClustersLoad(ClusterFile)
            
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

