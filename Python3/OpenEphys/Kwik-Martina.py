#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Nov 25 13:29:22 2016

@author: cerebro
"""
import Hdf5F
import h5py
import numpy as np
import os
#from datetime import datetime
from glob import glob
from multiprocessing import Process
from scipy import io
from subprocess import call

def CallWaveClus(Rate, Path):
    print('Clustering spikes...')
    
    MLab = '/home/cerebro/Software/MatLabR2014a/bin/matlab'
    CmdCd = 'cd ' + Path + '; '
    CmdCluster = 'try, par.sr=' + str(int(Rate)) + '; ' + \
                 "Get_spikes('Files.txt', 'parallel', true, 'par', par);" + \
                 "Files = dir('./*spikes.mat'); Files = {Files.name}; " + \
                 "Do_clustering(Files, 'make_plots', false); end; quit"
    call([MLab, '-r', CmdCd+CmdCluster])
    
    print('Done clustering.')
    return(None)


def CheckStimulationPattern(TTLCh):
    if np.mean(TTLCh) > 1000: 
        print('Sinusoidal stimulation')
        Threshold = (max(TTLCh) - min(TTLCh)) / 2
        return(Threshold)
    else:
        print('square pulses stimulation')
        Threshold = ((max(TTLCh) - min(TTLCh)) / 2) + 2*(np.std(TTLCh))
        return(Threshold) 


def GetProc(Raw, Board):
    print('Get proc no. for', Board, 'board... ', end='')
    ProcChs = {Proc: len(Raw[Proc]['data']['0'][1,:]) 
               for Proc in Raw.keys()}
    
    for Proc, Chs in ProcChs.items():
        if Chs == max(ProcChs.values()): OEProc = Proc
        else: RHAProc = Proc
    
    if 'RHAProc' not in locals(): RHAProc = OEProc
    
    if Board == 'OE': Proc = OEProc
    elif Board == 'RHA': Proc = RHAProc
    else: print("Choose Board as 'OE' or 'RHA'."); return(None)
    
    print('Done.')
    return(OEProc, RHAProc, Proc)


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


## Medium level
def QuantifyTTLsPerRec(Raw, Rec, AnalogTTLs, ChTTL=0, Proc=''):
    print('Get TTL timestamps... ', end='')
    if AnalogTTLs:
        TTLCh = Raw[Proc]['data'][Rec][:, ChTTL-1]
        Threshold = CheckStimulationPattern(TTLCh)
        TTLs = []
        for _ in range(1, len(TTLCh)):
            if TTLCh[_] > Threshold:
                if TTLCh[_-1] < Threshold: TTLs.append(_)
        
        print('Done.')
        return(TTLs)          


## Higher level functions
def ClusterizeSpks(Folder, AnalogTTLs=False, Board='OE', Override={}):
    print('Load DataInfo...')
    if AnalogTTLs: Raw, _, Files = Hdf5F.LoadOEKwik(Folder, AnalogTTLs)
    else: Raw, Events, _, Files = Hdf5F.LoadOEKwik(Folder, AnalogTTLs)
    
    OEProc = GetProc(Raw, Board)[0]
    
    Path = os.getcwd() + '/' + Folder +'/SepCh/'
    os.makedirs(Path, exist_ok=True)
    
    Rate = Raw[OEProc]['info']['0']['sample_rate']
    
    Clusters = {}
    for Rec in Raw[OEProc]['data'].keys():
        if Override != {}: 
            if 'Rec' in Override.keys():
                Rec = Override['Rec']
        
        RecS = "{0:02d}".format(int(Rec))
        
        print('Separating channels according to ChannelMap...')
        Data = [Raw[OEProc]['data'][Rec][:, _] * 
                Raw[OEProc]['channel_bit_volts'][Rec][_]
                for _ in range(len(Raw[OEProc]['data'][Rec][0,:-1]))]
                        
        print('Writing files for clustering... ', end='')
        FileList = []
        for Ind, Ch in enumerate(Data):
            MatName = 'Exp' + Files[OEProc+'_kwd'][-13:-8] + '_' + \
                      RecS + '-Ch' + "{0:02d}".format(Ind+1) + '.mat'
            
            FileList.append(MatName)
            io.savemat(Path+MatName, {'data': Ch})
        
        TxtFile = open(Path+'Files.txt', 'w')
        for File in FileList: TxtFile.write(File+'\n')
        TxtFile.close()
        print('Done.')
        
        CallWaveClus(Rate, Path)
        
        ClusterList = glob(Path+'times_*'); ClusterList.sort()
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
    
        ToDelete = glob(Path+'*')
        for File in ToDelete: os.remove(File)
    
    os.removedirs(Path)
    Hdf5F.WriteClusters(Clusters, Path[:-6]+'SpkClusters.hdf5')


def UnitsSpks(Folder, AnalysisFile, Board='OE', Override={}):
    print('Load DataInfo...')
    Raw = Hdf5F.LoadOEKwik(Folder, 'Raw')[0]
    OEProc = GetProc(Raw, Board)[0]
    
    Path = os.getcwd() + '/' + Folder 
    ClusterFile = Path + '/SpkClusters.hdf5'
    Clusters = Hdf5F.LoadClusters(ClusterFile)
    
    Units = {}
    for Rec in Raw[OEProc]['data'].keys():
        if Override != {}: 
            if 'Rec' in Override.keys():
                Rec = Override['Rec']
        
        RecS = "{0:02d}".format(int(Rec))
        
        TTLCh = Raw[OEProc]['data'][Rec][:,16]
        if np.mean(TTLCh) > 1000: 
            print('Sinusoidal stimulation')
        
        Units[RecS] = {}
        for Ch in Clusters[RecS].keys():
            Units[RecS][Ch] = SepSpksPerCluster(Clusters[RecS][Ch], 
                                                             Ch)
        
        if Override != {}: 
            if 'Rec' in Override.keys(): break
        
    del(Raw, Clusters)
        
    WriteUnits(Units, AnalysisFile)
    return(None)


def UnitsPSTH(Folder, AnalysisFile, TTLChNo=0, PSTHTimeBeforeTTL=0, 
                 PSTHTimeAfterTTL=300, AnalogTTLs=False, 
                 Board='OE', Override={}):
    print('Load DataInfo...')
    if AnalogTTLs: Raw, _, Files = Hdf5F.LoadOEKwik(Folder, AnalogTTLs)
    else: Raw, Events, _, Files = Hdf5F.LoadOEKwik(Folder, AnalogTTLs)
    
    Units = {}
    
    Path = os.getcwd() + '/' + Folder
    ClusterFile = Path + '/SpkClusters.hdf5'
    Clusters = Hdf5F.LoadClusters(ClusterFile)
    
    OEProc = GetProc(Raw, Board)[0]
    
    for Rec in Raw[OEProc]['data'].keys():
        if Override != {}: 
            if 'Rec' in Override.keys():
                Rec = Override['Rec']
        
        Rate = Raw[OEProc]['info'][Rec]['sample_rate']
        TTLCh = Raw[OEProc]['data'][Rec][:, TTLChNo-1]
        TTLs = QuantifyTTLsPerRec(Raw, Rec, AnalogTTLs, TTLChNo, OEProc)
        StimFreq = 1 / ((TTLs[3] - TTLs[2])/Rate)
        
        if np.mean(TTLCh) > 1000: 
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


def UnitPlotPerCh(ChDict, Ch, XValues, PulseDur, FigName, FigTitle):
    ClusterNo = len(ChDict['Spks'])
    if ClusterNo == 0: print(Ch, 'had no spikes :('); return(None)
    
    Params = {'backend': 'Agg'}
    from matplotlib import rcParams; rcParams.update(Params)
    from matplotlib import pyplot as plt
    
    Fig, Axes = plt.subplots(ClusterNo,2, figsize=(6, 3*ClusterNo))
    SpksYLabel = 'Voltage [mV]'; SpksXLabel = 'Time [ms]'
    PSTHYLabel = 'Number of spikes in channel'; PSTHXLabel = 'Time [ms]'
    SpanLabel = 'Sound pulse'
    
    for Class in ChDict['Spks'].keys():
        SpkNo = len(ChDict['Spks'][Class])
        print(str(SpkNo), 'Spks in cluster', Class)
        print('Max of', max(ChDict['PSTH'][Class]),  
              'Spks in PSTH')
        
#                        if not SpkNo:
#                            print('No Spk data on cluster', str(Cluster) + '. Skipping...')
#                            Thrash[SKey+FKey+RKey+Ch] = Units[SKey][FKey][RKey][Ch].copy()
#                            continue
        
#                PSTHPeak = max(Units[SKey][FKey][RKey]['PSTH'][Cluster])
#                PSTHMean = np.mean(Units[SKey][FKey][RKey]['PSTH'][Cluster])
#        if max(UnitRec[Key]['PSTH'][Cluster]) < 4: 
#            print('No peaks in PSTH. Skipping cluster', str(Cluster), '...')
#            continue
        
        if SpkNo > 100: SpkNo = np.arange(100); np.random.shuffle(SpkNo)
        else: SpkNo = np.arange(SpkNo)
        
        for Spike in SpkNo:
            if ClusterNo == 1: Axes[0].plot(ChDict['Spks'][Class][Spike], 'r')
            else: Axes[int(Class)-1][0].plot(ChDict['Spks'][Class][Spike], 'r')
        
        if ClusterNo == 1:
#                        Axes[0].set_title('Peak='+str(PSTHPeak)+' Mean='+str(PSTHMean))
            Axes[0].plot(np.mean(ChDict['Spks'][Class], axis=0), 'k')
            Axes[1].bar(XValues, ChDict['PSTH'][Class], 'r')
            
            Ind1 = list(XValues).index(0)
            Ind2 = list(XValues).index(int(PulseDur*1000))
            
            Axes[1].axvspan(XValues[Ind1], XValues[Ind2], color='k', alpha=0.3, 
                            lw=0, label=SpanLabel)
            
            SetPlot(AxesObj=Axes[0], Axes=True)
            SetPlot(AxesObj=Axes[1], Axes=True)
            Axes[0].set_ylabel(SpksYLabel); Axes[0].set_xlabel(SpksXLabel)
            Axes[1].set_ylabel(PSTHYLabel); Axes[1].set_xlabel(PSTHXLabel)
            
        else:
#                    Axes[Cluster][0].set_title('Peak='+str(PSTHPeak)+' Mean='+\
#                                               str(PSTHMean)+' Std='+str(PSTHStd))
            Axes[int(Class)-1][0].plot(np.mean(ChDict['Spks'][Class], axis=0), 'k')
            Axes[int(Class)-1][1].bar(XValues, ChDict['PSTH'][Class], 'r')
            
            Ind1 = list(XValues).index(0)
            Ind2 = list(XValues).index(int(PulseDur*1000))
            
            Axes[int(Class)-1][1].axvspan(XValues[Ind1], XValues[Ind2], 
                                          color='k', alpha=0.3, lw=0, 
                                          label=SpanLabel)
            
            SetPlot(AxesObj=Axes[int(Class)-1][0], Axes=True)
            SetPlot(AxesObj=Axes[int(Class)-1][1], Axes=True)
            Axes[int(Class)-1][0].set_ylabel(SpksYLabel)
            Axes[int(Class)-1][0].set_xlabel(SpksXLabel)
            Axes[int(Class)-1][1].set_ylabel(PSTHYLabel)
            Axes[int(Class)-1][1].set_xlabel(PSTHXLabel)
    
    SetPlot(FigObj=Fig, FigTitle=FigTitle, Plot=True)
    print('Writing to', FigName+'... ', end='')
    Fig.savefig(FigName, format='svg')
    print('Done.')
    return(None)


#%% Run
TTLChNo = 17
PSTHTimeBeforeTTL = 150
PSTHTimeAfterTTL = 150
AnalogTTLs = True
Override = {}

Folders = glob('Animal*')
for Folder in Folders:
    Subfolders = glob(Folder+'/2016*')
    
    for Subfolder in Subfolders:
        ClusterizeSpks(Subfolder, AnalogTTLs)
        AnalysisFile = '-'.join(Subfolder.split('/'))
        UnitsPSTH(Subfolder, AnalysisFile, TTLChNo, PSTHTimeBeforeTTL, 
                  PSTHTimeAfterTTL, AnalogTTLs)
        
#AnimalName = os.getcwd().split('/')[-1]
#AnalysisFile = Folder + '/' + AnimalName + '-Analysis.hdf5'


#%%
import Hdf5F
import numpy as np

Params = {'backend': 'TkAgg'}
from matplotlib import rcParams; rcParams.update(Params)
from matplotlib import pyplot as plt

for s in range(len(Units['00-Sinusoidal-2Hz']['Ch09']['Spks']['01'][:,0])): 
    plt.plot(Units['00-Sinusoidal-2Hz']['Ch09']['Spks']['01'][s,:])
plt.show()

Folder = '2016-06-08_21-43-09'

Raw, Spks, Files = Hdf5F.LoadOEKwik(Folder, True)
Means = [0] * len(Raw['102']['data'])
for _ in range(len(Raw['102']['data'])):
    print(_)
    TTLCh = Raw['102']['data'][str(_)][15*30000:-15*30000,16]
    Means[_] = np.mean(TTLCh)

print(Means)
T = np.arange(0,len(TTLCh)); T = T/30000

plt.plot(T, TTLCh); plt.show()
