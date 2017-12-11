# -*- coding: utf-8 -*-
"""
@author: T. Malfatti <malfatti@disroot.org>
@year: 2015
@license: GNU GPLv3 <https://raw.githubusercontent.com/malfatti/SciScripts/master/LICENSE>
@homepage: https://github.com/Malfatti/SciScripts

Functions for manipulating specific .mat files.
"""
import numpy as np
import os

from DataAnalysis.DataAnalysis import FilterSignal
from DataAnalysis import GPIAS
from IO import Asdf, Hdf5
from glob import glob
from scipy import io#, signal

def GPIASAnalysis(Folders, InfoFiles, AnalysisFolder, GPIASTimeBeforeTTL=100, GPIASTimeAfterTTL=100, 
                     FilterFreq=[70, 400], FilterOrder=3, Filter = 'butter', SliceSize=100):
    
    SliceSizeMS = [SliceSize][0]
    for F, Folder in enumerate(Folders):
        DataInfo = io.loadmat(InfoFiles[F])
        DataInfo['AnimalName'] = Folder.split('_')[-1]
        
        Rate = np.array(int(DataInfo['Rate'])); TTL = DataInfo['PulseStart']
        Freqs = glob(Folder + '/*.mat')
        GPIASRec = {'Trace': {}, 'Index':{}}
        
        SliceSize = int(SliceSizeMS * (Rate/1000))
        
        for Freq in Freqs:
            SFreq = Freq[-9:-4]; SFreq = SFreq.split('_')
            try:
                SFreq = str(int(SFreq[0])*1000) + '-' + str(int(SFreq[1])*1000)
            except ValueError:
                SFreq = str(int(SFreq[0][1:])*1000) + '-' + str(int(SFreq[1])*1000)
            
            for Key in GPIASRec.keys():
                if SFreq not in GPIASRec[Key].keys(): GPIASRec[Key][SFreq] = {}
            
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
            
            GPIASRec['Trace'][SFreq]['Gap'] = Data['Gap'][:]
            GPIASRec['Trace'][SFreq]['NoGap'] = Data['NoGap'][:]
            GPIASRec['Index'][SFreq]['Gap'] = Data['Gap'][:]
            GPIASRec['Index'][SFreq]['NoGap'] = Data['NoGap'][:]
            
            Keys = [['Gap', 'NoGap', 'GPIASIndex']]
            GPIASRec['Index'][SFreq] = GPIAS.IndexCalc(
                                       GPIASRec['Index'][SFreq], Keys, 
                                       NoOfSamplesBefore, SliceSize)
            
            del(Data)
        
        AnalysisKey = InfoFiles[F].split('/')[1].split('-')
        AnalysisKey[0] = AnalysisKey[0]+'000000'
        AnalysisKey[1] = DataInfo['AnimalName']
        AnalysisKey = '-'.join(AnalysisKey) + '-Sound-Recovery_' + 'GPIAS'
#        Hdf5.DataWrite({'GPIAS': GPIASRec, 'XValues': XValues}, AnalysisKey, AnalysisFile)
        Asdf.Write({'GPIAS': GPIASRec, 'XValues': XValues}, '/', AnalysisFolder + '/' + AnalysisKey+'.asdf')


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

