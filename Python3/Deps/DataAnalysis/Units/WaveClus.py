#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jun 12 14:12:37 2017

@author: malfatti
"""
import numpy as np
import os

from DataAnalysis.DataAnalysis import GetRecXValues, QuantifyTTLsPerRec
from glob import glob
from IO import Hdf5, WaveClus
from scipy import io


## Level 0
def PSTH(Clusters, TTLCh, Rate, AnalysisFile, AnalysisKey, Rec='0', 
              TimeBeforeTTL=0, TimeAfterTTL=300, AnalogTTLs=True, 
              Override={}):
    try: Rec = "{0:02d}".format(int(Rec))
    except ValueError: pass
        
    if 'TTLs' in Override.keys(): TTLs = Override['TTLs']
    else: TTLs = QuantifyTTLsPerRec(AnalogTTLs, TTLCh)
    
    XValues = GetRecXValues(TTLs, Rate, TimeBeforeTTL, TimeAfterTTL)
    UnitsData = {Rec: {}}
    
    print('Preparing histograms and spike waveforms...')
    for CKey in Clusters[Rec].keys():
        Ids = np.unique(Clusters[Rec][CKey]['Id'])
#            UnitsData[URecS][CKey] = SepSpksPerCluster(Clusters[Rec][CKey], CKey)
        UnitsData[Rec][CKey] = {'PSTH': {}}
        
        for Id in Ids:
            IdIndex = Clusters[Rec][CKey]['Id'] == Id
            if not len(Clusters[Rec][CKey]['Spikes'][IdIndex,:]): continue
            
            Id = "{0:02d}".format(int(Id))
            
            Hist = np.array([])
            for TTL in range(len(TTLs)):
                Firing = Clusters[Rec][CKey]['Timestamps'][IdIndex] \
                         - (TTLs[TTL]/(Rate/1000))
                Firing = Firing[(Firing >= XValues[0]) * 
                                (Firing < XValues[-1])]
                
                Hist = np.concatenate((Hist, Firing)); del(Firing)
            
            UnitsData[Rec][CKey]['PSTH'][Id] = Hist[:]
            del(Hist)
        
        print(CKey+':', str(len(Ids)), 'clusters.')
    
    del(TTLs)
    
    Group = AnalysisKey + '/Units'
    Hdf5.WriteUnits(UnitsData, AnalysisFile, Group, XValues)
    del(UnitsData)
    return(None)


def SepSpksPerCluster(Clusters, Ch):
    Ids = np.unique(Clusters['Id'])
    Dict = {}                    
#    Dict['Spks'] = [[] for _ in range(len(Ids))]
    Dict['Spks'] = {}
    
    print(Ch+':', str(len(Ids)), 'clusters:')
    for Id in Ids:
        IdIndex = Clusters['Id'] == Id
        
        SpkNo = len(Clusters['Spikes'][IdIndex,:])
        if not SpkNo: continue
        
        Id = "{0:02d}".format(int(Id))
        Dict['Spks'][Id] =  Clusters['Spikes'][IdIndex,:][:]
    
        print('    Id', Id, '-', str(SpkNo), 'spikes.')
    
    if len(Dict): return(Dict)
    else: return({})


## Level 1
def ClusterizeSpks(Data, Rate, ChannelMap, ClusterPath, AnalysisFile, 
                   AnalysisKey, Rec='0', Override={}, Return=False):
    """ Detect and clusterize spks using WaveClus """
    
    os.makedirs(ClusterPath, exist_ok=True)
    
    Data = [Data[:, _-1] for _ in sorted(ChannelMap)]
    print('Writing files for clustering... ', end='')
    FileList = []
    
    try: Rec = "{0:02d}".format(int(Rec))
    except ValueError: pass
    
    for Ind, Ch in enumerate(Data):
        MatName = 'Rec' + Rec + '-Ch' + "{0:02d}".format(Ind+1) + '.mat'
        
        FileList.append(MatName)
        io.savemat(ClusterPath + '/' + MatName, {'data': Ch})
    
    TxtFile = open(ClusterPath + '/Files.txt', 'w')
    for File in FileList: TxtFile.write(File + '\n')
    TxtFile.close()
    print('Done.')
    
    WaveClus.Run(Rate, ClusterPath)
    
    ClusterList = glob(ClusterPath + '/times_*'); ClusterList.sort()
    ClusterList = [_ for _ in ClusterList if _.split('-')[-2].split('_')[-1] == 'Rec'+Rec]
#        print(ClusterList)
    Clusters = {}
    Clusters[Rec] = WaveClus.LoadClusters(ClusterList, 'WaveClus')
    
    Group = AnalysisKey + '/SpkClusters'
    Hdf5.WriteClusters(Clusters, AnalysisFile, Group)
    
    ToDelete = glob(ClusterPath + '/*')
    for File in ToDelete: 
        if File not in glob(ClusterPath + '/times_*'): os.remove(File)
#    os.removedirs(ClusterPath)
    
    if Return: return(Clusters)
    else: return(None)


def Spks(Clusters, AnalysisFile, AnalysisKey, Rec='0', Override={}):
    try: Rec = "{0:02d}".format(int(Rec))
    except ValueError: pass
    
    UnitsData = {Rec: {}}
    for Ch in Clusters[Rec].keys():
        UnitsData[Rec][Ch] = SepSpksPerCluster(Clusters[Rec][Ch], Ch)
    
    Group = AnalysisKey + '/Units'
    Hdf5.WriteUnits(UnitsData, AnalysisFile, Group)
    return(None)

