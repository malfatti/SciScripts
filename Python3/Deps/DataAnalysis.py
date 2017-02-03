#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
    Copyright (C) 2017  T. Malfatti
    
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

import Hdf5F
import h5py
import IntanBin
import numpy as np
import os
from glob import glob
from scipy import io
from subprocess import call

## Level 0
def CallWaveClus(Rate, Path):
    print('Clustering spikes...')
    
#    MLab = '/home/cerebro/Software/Programs/MatLabR2014a/bin/matlab'
    MLab = '/home/malfatti/Software/Programs/MatLabR2015a/bin/matlab'
    CmdCd = 'cd ' + Path + '; '
    CmdCluster = 'try, par.sr=' + str(int(Rate)) + '; ' + \
                 "Get_spikes('Files.txt', 'parallel', true, 'par', par);" + \
                 "Files = dir('./*spikes.mat'); Files = {Files.name}; " + \
                 "Do_clustering(Files, 'make_plots', false); end; quit"
    call([MLab, '-r', CmdCd+CmdCluster])
    
    print('Done clustering.')
    return(None)


def CheckStimulationPattern(TTL):
    if np.mean(TTL) > 1000: 
        print('Sinusoidal stimulation')
        Threshold = (max(TTL) - min(TTL)) / 2
        return(Threshold)
    else:
        print('square pulses stimulation')
        Threshold = ((max(TTL) - min(TTL)) / 2) + 2*(np.std(TTL))
        return(Threshold) 


def RemapChannels(Tip, Head, Connector):
    """
    Get probe channels order. It doesn't matter what logic you follow to order 
    your connector channels, but you MUST follow the same logic for your probe 
    head.
    
    If the probe tip channels are put top-down or bottom-up, the resulting 
    channel map will be ordered accordingly.
    
    Example:
        CustomAdaptor = [5, 6, 7, 8, 9, 10 ,11, 12, 13, 14, 15, 16, 1, 2, 3, 4]
        RHAHeadstage = [16, 15, 14, 13, 12, 11, 10, 9, 1, 2, 3, 4, 5, 6, 7, 8]
        A16 = {'ProbeTip': [9, 8, 10, 7, 13, 4, 12, 5, 15, 2, 16, 1, 14, 3, 11, 6],
               'ProbeHead': [8, 7, 6, 5, 4, 3, 2, 1, 9, 10, 11, 12, 13, 14, 15, 16]}
        
        ChannelMap = RemapChannels(A16['ProbeTip'], A16['ProbeHead'], CustomAdaptor)
    """
    print('Get channel order... ', end='')
    ChNo = len(Tip)
    ChMap = [0]*ChNo
    
    for Ch in range(ChNo):
        TipCh = Tip[Ch] # What channel should be the Ch
        HeadCh = Head.index(TipCh) # Where Ch is in Head
        ChMap[Ch] = Connector[HeadCh] # Channels in depth order
    
    print('Done.')
    return(ChMap)


def SepSpksPerCluster(Clusters, Ch):
    Classes = np.unique(Clusters['ClusterClass'])
    Dict = {}                    
#    Dict['Spks'] = [[] for _ in range(len(Classes))]
    Dict['Spks'] = {}
    
    print(Ch+':', str(len(Classes)), 'clusters:')
    for Class in Classes:
        ClassIndex = Clusters['ClusterClass'] == Class
        
        SpkNo = len(Clusters['Spikes'][ClassIndex,:])
        if not SpkNo: continue
        
        Class = "{0:02d}".format(int(Class))
        Dict['Spks'][Class] =  Clusters['Spikes'][ClassIndex,:][:]
    
        print('    Class', Class, '-', str(SpkNo), 'spikes.')
    
    if len(Dict): return(Dict)
    else: return({})


def WriteUnits(Units, FileName, XValues=[]):
    """ Belongs to Hdf5F """
    print('Writing data to', FileName+'... ', end='')
    Group = 'UnitRec'
    Thrash = []
    with h5py.File(FileName) as F:
        for RKey in Units.keys():
            for Ch in Units[RKey].keys():
                Path = '/'+Group+'/'+RKey+'/'+Ch
                if Path not in F: F.create_group(Path)
                
                for Var in Units[RKey][Ch].keys():
                    if Var == 'Spks':
                        if '/'+Path+'/Spks' in F: del(F[Path]['Spks'])
                        if Path+'/Spks' not in F: 
                            F.create_group(Path+'/Spks')
                        
                        for Class in Units[RKey][Ch]['Spks'].keys():
                            for Key, Spk in enumerate(Units[RKey][Ch]['Spks'][Class]):
                                if not len(Spk):
                                    Thrash.append([Path, Class, Key])
                                    del(Units[RKey][Ch]['Spks'][Class][Key])
                                    continue
                                if Path+'/Spks/'+Class not in F:
                                    F[Path]['Spks'].create_group(Class)
                                F[Path]['Spks'][Class][str(Key)] = Spk
                    
                    elif Var == 'PSTH':
                        if '/'+Path+'/PSTH' in F: del(F[Path]['PSTH'])
                        F[Path].create_group('PSTH')
                        for VarKey in Units[RKey][Ch][Var].keys():
                            F[Path]['PSTH'][VarKey] = Units[RKey][Ch]['PSTH'][VarKey][:]
                    
                    elif Var == 'Spks_Info':
                        for VarKey, VarValue in Units[RKey][Ch][Var].items():
                            F[Path]['Spks'].attrs[VarKey] = VarValue
                    
                    elif Var == 'PSTH_Info':
                        for VarKey, VarValue in Units[RKey][Ch][Var].items():
                            F[Path]['PSTH'].attrs[VarKey] = VarValue
                    
                    else:
                        print(Var, 'not supported')
        
        if len(XValues): F[Group]['XValues'] = XValues
    
    if Thrash: 
        for _ in Thrash: print(_)
    
    return(None)


## Level 1
def ClusterizeSpks(Data, Rate, Path, AnalogTTLs=False, ChannelMap=[], Override={}):
    """ Detect and clusterize spks using WaveClus """
    
    os.makedirs(Path, exist_ok=True)
    
    Clusters = {}
    for Rec in Data.keys():
        if not ChannelMap: ChannelMap = [_ for _ in range(Data[Rec].shape[1])]
        
        if Override != {}: 
            if 'Rec' in Override.keys():
                Rec = Override['Rec']
        
        RecS = "{0:02d}".format(int(Rec))
        
        print('Separating channels according to ChannelMap...')
        Data = [Data[Rec][:, _] for _ in ChannelMap]
                        
        print('Writing files for clustering... ', end='')
        FileList = []
        for Ind, Ch in enumerate(Data):
            MatName = 'Rec' + RecS + '-Ch' + "{0:02d}".format(Ind+1) + '.mat'
            
            FileList.append(MatName)
            io.savemat(Path+MatName, {'data': Ch})
        
        TxtFile = open(Path+'Files.txt', 'w')
        for File in FileList: TxtFile.write(File + '\n')
        TxtFile.close()
        print('Done.')
        
        CallWaveClus(Rate, Path)
        
        ClusterList = glob(Path + 'times_*'); ClusterList.sort()
        ClusterList = [_ for _ in ClusterList if _[-11:-9] == RecS]
        
        Clusters[RecS] = {}
        for File in ClusterList:
            Ch = File[-8:-4]
            
            ClusterFile = io.loadmat(File)
            Clusters[RecS][Ch] = {}
            Clusters[RecS][Ch]['ClusterClass'] = ClusterFile['cluster_class'][:, 0]
            Clusters[RecS][Ch]['Timestamps'] = ClusterFile['cluster_class'][:, 1]
            Clusters[RecS][Ch]['Spikes'] = ClusterFile['spikes'][:]
        
        if Override != {}: 
            if 'Rec' in Override.keys(): break
    
        ToDelete = glob(Path + '*')
        for File in ToDelete: os.remove(File)
    
    os.removedirs(Path)
    Hdf5F.WriteClusters(Clusters, Path[:-6]+'SpkClusters.hdf5')


def QuantifyTTLsPerRec(Data, Rec, AnalogTTLs, TTLCh=0):
    print('Get TTL timestamps... ', end='')
    if AnalogTTLs:
        TTL = Data[Rec][:, TTLCh-1]
        Threshold = CheckStimulationPattern(TTL)
        TTLs = []
        for _ in range(1, len(TTL)):
            if TTL[_] > Threshold:
                if TTL[_-1] < Threshold: TTLs.append(_)
        
        print('Done.')
        return(TTLs)          


## Level 2
def UnitsPSTH(Data, Rate, Path, AnalysisFile, TTLCh=0, PSTHTimeBeforeTTL=0, 
                 PSTHTimeAfterTTL=300, AnalogTTLs=False, 
                 Override={}):
    Units = {}
    
    ClusterFile = Path + '/SpkClusters.hdf5'
    Clusters = Hdf5F.LoadClusters(ClusterFile, KeyInd=-1)
    
    for Rec in Data.keys():
        if Override != {}: 
            if 'Rec' in Override.keys():
                Rec = Override['Rec']
        
        TTL = Data[Rec][:, TTLCh-1]
        TTLs = QuantifyTTLsPerRec(Data, Rec, AnalogTTLs, TTLCh)
        StimFreq = 1 / ((TTLs[3] - TTLs[2])/Rate)
        
        if np.mean(TTL) > 1000: 
            URecS = "{0:02d}".format(int(Rec)) + '-Sinusoidal-' + str(int(StimFreq)) + 'Hz'
            NoOfSamplesBefore = int(round((0*Rate)*10**-3))
            NoOfSamplesAfter = int(round((TTLs[3] - TTLs[2])))
        else:
            URecS = "{0:02d}".format(int(Rec)) + '-Squares-' + str(int(StimFreq)) + 'Hz'
            NoOfSamplesBefore = int(round((PSTHTimeBeforeTTL*Rate)*10**-3))
            NoOfSamplesAfter = int(round((PSTHTimeAfterTTL*Rate)*10**-3))
        
        NoOfSamples = NoOfSamplesBefore + NoOfSamplesAfter
        XValues = (range(-NoOfSamplesBefore, 
                         NoOfSamples-NoOfSamplesBefore)/Rate)*10**3
        
        RecS = "{0:02d}".format(int(Rec))
        Units[URecS] = {}
        print('Preparing histograms and spike waveforms...')
        for CKey in Clusters[RecS].keys():
            Classes = np.unique(Clusters[RecS][CKey]['ClusterClass'])
#            Units[URecS][CKey] = {}
            Units[URecS][CKey] = SepSpksPerCluster(Clusters[RecS][CKey], CKey)
            Units[URecS][CKey]['PSTH'] = {}
            
            for Class in Classes:
                ClassIndex = Clusters[RecS][CKey]['ClusterClass'] == Class
                if not len(Clusters[RecS][CKey]['Spikes'][ClassIndex,:]): continue
                
                Class = "{0:02d}".format(int(Class))
                Units[URecS][CKey]['PSTH'][Class] = np.zeros(len(XValues))
                
                for TTL in range(len(TTLs)):
                    Firing = Clusters[RecS][CKey]['Timestamps'][ClassIndex] \
                             - (TTLs[TTL]/(Rate/1000))
                    Firing = Firing[(Firing >= XValues[0]) * 
                                    (Firing < XValues[-1])]
                    SpkCount = np.histogram(Firing, 
                                            np.hstack((XValues, len(XValues)/(Rate/1000))))[0]
                    
                    Units[URecS][CKey]['PSTH'][Class] = Units[URecS][CKey]['PSTH'][Class] + SpkCount
                    
                    del(Firing, SpkCount)
            
            print(CKey+':', str(len(Classes)), 'clusters.')
        
        del(TTLs)
        if Override != {}: 
            if 'Rec' in Override.keys(): break
    
    WriteUnits(Units, AnalysisFile, XValues)
    return(None)


def UnitsSpks(Data, Path, AnalysisFile, Override={}):
    ClusterFile = Path + '/SpkClusters.hdf5'
    Clusters = Hdf5F.LoadClusters(ClusterFile)
    
    Units = {}
    for Rec in Data.keys():
        if Override != {}: 
            if 'Rec' in Override.keys():
                Rec = Override['Rec']
        
        RecS = "{0:02d}".format(int(Rec))
        
        Units[RecS] = {}
        for Ch in Clusters[RecS].keys():
            Units[RecS][Ch] = SepSpksPerCluster(Clusters[RecS][Ch], 
                                                             Ch)
        
        if Override != {}: 
            if 'Rec' in Override.keys(): break
    
    WriteUnits(Units, AnalysisFile)
    del(Data, Clusters)
    
    return(None)



#%% Tests
File = 'Intan/Arch45S_153244.int'; ClusterPath = './SepCh/'
Raw = {}
Raw['0'] = IntanBin.IntLoad(File)[0]

RHAHeadstage = [16, 15, 14, 13, 12, 11, 10, 9, 1, 2, 3, 4, 5, 6, 7, 8]
A16 = {'ProbeTip': [9, 8, 10, 7, 13, 4, 12, 5, 15, 2, 16, 1, 14, 3, 11, 6],
       'ProbeHead': [8, 7, 6, 5, 4, 3, 2, 1, 9, 10, 11, 12, 13, 14, 15, 16]}
ChannelMap = RemapChannels(A16['ProbeTip'], A16['ProbeHead'], RHAHeadstage)
