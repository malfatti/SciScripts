#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jun 12 15:06:47 2017

@author: malfatti
"""
import numpy as np
from multiprocessing import Process
from Plot.Plot import Set


## Level 0
def UnitPerCh(ChDict, Ch, XValues, FigName, Ext):
    ClusterNo = len(ChDict['Spks'])
    if ClusterNo == 0: print(Ch, 'had no spikes :('); return(None)
    
    PSTHNo = 0
    for Id in ChDict['PSTH'].keys(): 
        PSTHNo += len(ChDict['PSTH'][Id])
    
    if not PSTHNo:
        print('No Spks in PSTHs of this channel :( Skipping channel...')
        return(None)
    
    Params = Set(Backend='Agg', Params=True)
    from matplotlib import rcParams; rcParams.update(Params)
    from matplotlib import pyplot as plt
    
    Fig, Axes = plt.subplots(ClusterNo,2, figsize=(8, 3*ClusterNo))
    SpksYLabel = 'Voltage [µV]'; SpksXLabel = 'Time [ms]'
    PSTHYLabel = 'Number of spikes in channel'; PSTHXLabel = 'Time [ms]'
#    SpanLabel = 'Sound pulse'
    
    for Id in ChDict['Spks'].keys():
        SpkNo = len(ChDict['Spks'][Id])
        print(str(SpkNo), 'Spks in cluster', Id)
        if len(ChDict['PSTH'][Id]) == 0:
            print('No Spks in PSTH. Skipping Id...')
            continue
        else:
            print('Max of', len(ChDict['PSTH'][Id]), 'Spks in PSTH.')
        
        if SpkNo > 50: 
            SpkNo = np.arange(SpkNo)
            np.random.shuffle(SpkNo)
            SpkNo = SpkNo[:50]
        else: 
            SpkNo = np.arange(SpkNo)
        
        
        for Spike in SpkNo:
            x = np.arange(len(ChDict['Spks'][Id][Spike])) / 30
            if ClusterNo == 1: Axes[0].plot(x, ChDict['Spks'][Id][Spike], 'r')
            else: Axes[int(Id)-1][0].plot(x, ChDict['Spks'][Id][Spike], 'r')
        
        x = np.arange(len(np.mean(ChDict['Spks'][Id], axis=0))) / 30
        if ClusterNo == 1:
            Axes[0].plot(x, np.mean(ChDict['Spks'][Id], axis=0), 'k')
            Axes[1].hist(ChDict['PSTH'][Id], XValues)
            
#            Ind1 = list(XValues).index(0)
#            Ind2 = list(XValues).index(int(PulseDur*1000))
            
#            Axes[1].axvspan(XValues[Ind1], XValues[Ind2], color='k', alpha=0.3, 
#                            lw=0, label=SpanLabel)
            
            Set(AxesObj=Axes[0], Axes=True)
            Set(AxesObj=Axes[1], Axes=True)
            Axes[0].set_ylabel(SpksYLabel); Axes[0].set_xlabel(SpksXLabel)
            Axes[1].set_ylabel(PSTHYLabel); Axes[1].set_xlabel(PSTHXLabel)
            
        else:
            Axes[int(Id)-1][0].plot(x, np.mean(ChDict['Spks'][Id], axis=0), 'k')
            Axes[int(Id)-1][1].hist(ChDict['PSTH'][Id], XValues)
            
#            Ind1 = list(XValues).index(0)
#            Ind2 = list(XValues).index(int(PulseDur*1000))
            
#            Axes[int(Id)-1][1].axvspan(XValues[Ind1], XValues[Ind2], 
#                                          color='k', alpha=0.3, lw=0, 
#                                          label=SpanLabel)
            
            Set(AxesObj=Axes[int(Id)-1][0], Axes=True)
            Set(AxesObj=Axes[int(Id)-1][1], Axes=True)
            Axes[int(Id)-1][0].set_ylabel(SpksYLabel)
            Axes[int(Id)-1][0].set_xlabel(SpksXLabel)
            Axes[int(Id)-1][1].set_ylabel(PSTHYLabel)
            Axes[int(Id)-1][1].set_xlabel(PSTHXLabel)
    
    FigTitle = FigName.split('/')[-1][:-4]
    Set(FigObj=Fig, FigTitle=FigTitle, Plot=True)
    print('Writing to', FigName+'... ', end='')
    Fig.savefig(FigName, format=Ext)
    print('Done.')
    return(None)


def TestBinSizeCh(ChDict, Ch, XValuesList, FigName, Ext):
    ClusterNo = len(ChDict['Spks'])
    if ClusterNo == 0: print(Ch, 'had no spikes :('); return(None)
    
    PSTHNo = 0
    for Id in ChDict['PSTH'].keys(): 
        PSTHNo += len(ChDict['PSTH'][Id])
    
    if not PSTHNo:
        print('No Spks in PSTHs of this channel :( Skipping channel...')
        return(None)
    
    Params = Set(Backend='Agg', Params=True)
    from matplotlib import rcParams; rcParams.update(Params)
    from matplotlib import pyplot as plt
    
    Fig, Axes = plt.subplots(ClusterNo, len(XValuesList), 
                             figsize=(4*len(XValuesList), 3*ClusterNo))
    
    PSTHYLabel = 'Number of spikes in channel'; PSTHXLabel = 'Time [ms]'
#    SpanLabel = 'Sound pulse'
    
    for Id in ChDict['Spks'].keys():
        SpkNo = len(ChDict['Spks'][Id])
        print(str(SpkNo), 'Spks in cluster', Id)
        if len(ChDict['PSTH'][Id]) == 0:
            print('No Spks in PSTH. Skipping Id...')
            continue
        else:
            print('Max of', len(ChDict['PSTH'][Id]), 'Spks in PSTH.')
        
        if SpkNo > 50: 
            SpkNo = np.arange(SpkNo)
            np.random.shuffle(SpkNo)
            SpkNo = SpkNo[:50]
        else: 
            SpkNo = np.arange(SpkNo)
        
        if ClusterNo == 1:
            for XInd, XValues in enumerate(XValuesList):
                Axes[XInd].hist(ChDict['PSTH'][Id], XValues)
                SubTitle = str(XValues[1] - XValues[0]) + ' ms bin size'
                Set(AxesObj=Axes[XInd], Axes=True)
                Axes[XInd].set_ylabel(PSTHYLabel)
                Axes[XInd].set_xlabel(PSTHXLabel)
                Axes[XInd].set_title(SubTitle)
            
        else:
            for XInd, XValues in enumerate(XValuesList):
                Axes[int(Id)-1][XInd].hist(ChDict['PSTH'][Id], XValues)
                SubTitle = str(XValues[1] - XValues[0]) + ' ms bin size'
                Set(AxesObj=Axes[int(Id)-1][XInd], Axes=True)
                Axes[int(Id)-1][XInd].set_ylabel(PSTHYLabel)
                Axes[int(Id)-1][XInd].set_xlabel(PSTHXLabel)
                Axes[int(Id)-1][XInd].set_title(SubTitle)
    
    FigTitle = FigName.split('/')[-1][:-4]
    Set(FigObj=Fig, FigTitle=FigTitle, Plot=True)
    print('Writing to', FigName+'... ', end='')
    Fig.savefig(FigName, format=Ext)
    print('Done.')
    return(None)


## Level 1
def SpksPSTH(Units, XValues, FigName, Mode='SpksPSTH', Ext='svg', Procs=8):
    for RKey in Units:
        ChNo = len(Units[RKey])
        ProcLists = [[] for _ in range(0, ChNo, Procs)]
        
        for Ind, Ch in enumerate(Units[RKey]):
            if Mode == 'BinSizeTest':
                FigName = FigName + '_' + Ch + '-BinSizeTest.', Ext
                ProcLists[Ind//Procs].append(
                    Process(target=TestBinSizeCh, 
                            args=(Units[RKey][Ch], Ch, XValues, FigName, Ext))
                )
            elif Mode == 'SpksPSTH':
                FigName = FigName + '_' + Ch + '-UnitsPSTH.' + Ext
                ProcLists[Ind//Procs].append(
                    Process(target=UnitPerCh, 
                            args=(Units[RKey][Ch], Ch, XValues, FigName, Ext))
                )
        
        for ProcList in ProcLists:
            for Proc in ProcList:
                Proc.start(); print('PID =', Proc.pid)
            Proc.join()
        
        return(None)

