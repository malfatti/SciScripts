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


## Higher level functions
def ClusterizeSpks(Folder, AnalogTTLs=False, Board='OE', Override={}):
    print('Load DataInfo...')
    DirList = glob(Folder+'/*'); DirList.sort()
    
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
                for _ in range(len(Raw[OEProc]['data'][Rec][0,:]))]
                        
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
    
    Raw, Spks, Files = Hdf5F.LoadOEKwik(Folder, True)
    OEProc = GetProc(Raw, Board)[0]
    TTLCh = Raw[OEProc]['data']['0'][:,16]
    del(Raw)
    
    if np.mean(TTLCh) > 1000: 
        print('Sinusoidal stimulation')
        
        Units = {}
        Raw = Hdf5F.LoadOEKwik(Folder, 'Raw')[0]
        OEProc = GetProc(Raw, Board)[0]
        
        Path = os.getcwd() + '/' + Folder 
        ClusterFile = Path + '/SpkClusters.hdf5'
        Clusters = Hdf5F.LoadClusters(ClusterFile)
        
        for Rec in Raw[OEProc]['data'].keys():
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
        
        del(Raw, Clusters)
        
    WriteUnits(Units, AnalysisFile)
    return(None)


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
        
        if len(XValues): F[Group].attrs['XValues'] = XValues
    
    if Thrash: 
        for _ in Thrash: print(_)
    
    return(None)


#%% Clustering
AnalogTTLs = True
Override = {}

Folders = glob('./Animal*')
for Folder in Folders:
    Subfolders = glob(Folder+'/2016*')
    
    for Subfolder in Subfolders:
        ClusterizeSpks(Subfolder, AnalogTTLs)

AnimalName = os.getcwd().split('/')[-1]
AnalysisFile = Folder + '/' + AnimalName + '-Analysis.hdf5'
UnitsSpks(Folder, AnalysisFile)


#%%
import Hdf5F
import numpy as np

Params = {'backend': 'TkAgg'}
from matplotlib import rcParams; rcParams.update(Params)
from matplotlib import pyplot as plt

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
