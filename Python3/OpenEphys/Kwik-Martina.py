#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Nov 25 13:29:22 2016

@author: cerebro
"""
import Hdf5F
import h5py
import Kwik
import numpy as np
import os
#from datetime import datetime
from glob import glob
from multiprocessing import Process
from numbers import Number
from scipy import io
from subprocess import call

# Level 0
def CallWaveClus(Rate, Path):
    print('Clustering spikes...')
    
    MLab = '/home/cerebro/Software/Programs/MatLabR2014a/bin/matlab'
    CmdCd = 'cd ' + Path + '; '
    CmdCluster = 'try, par.sr=' + str(int(Rate)) + '; ' + \
                 "Get_spikes('Files.txt', 'parallel', true, 'par', par);" + \
                 "Files = dir('./*spikes.mat'); Files = {Files.name}; " + \
                 "Do_clustering(Files, 'make_plots', false); end; quit"
    call([MLab, '-r', CmdCd+CmdCluster])
    
    print('Done clustering.')
    return(None)


def CheckStimulationPattern(TTLCh, StimType=''):
    if StimType == 'Sin':
        print('Sinusoidal stimulation')
        Threshold = (max(TTLCh) - min(TTLCh)) / 2
        return(Threshold)
    elif StimType == 'Sq':
        print('square pulses stimulation')
        Threshold = ((max(TTLCh) - min(TTLCh)) / 2) + 2*(np.std(TTLCh))
        return(Threshold) 
    else:
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


def GetExpKeys(ExpStr, OpenedFile, KeyInd=None, KeyStr=None):
    Keys = [K for K in OpenedFile.keys() if ExpStr in K]; Keys.sort()
    
    if KeyStr != None:
        Key = [K for K in Keys if KeyStr in K]; Key = Key[0]
        
    elif KeyInd != None:
        Key = Keys[KeyInd]
        
    else: 
        if len(Keys) > 1:
            print('Choose dataset:')
            for Ind, K in enumerate(Keys):
                print(str(Ind), '=' , K)
            Key = input(': ')
            Key = Keys[int(Key)]
        else:
            Key = Keys[0]
    
    return(Key)


def LoadOEKwik(RecFolder, AnalogTTLs, Unit='uV'):
    if AnalogTTLs: RawRO, _, Spks, Files = Kwik.load_all_files(RecFolder)
    else: RawRO, Events, Spks, Files = Kwik.load_all_files(RecFolder)
    
    print('Check if files are ok...')
    if 'RawRO' not in locals():
        print('.kwd file is corrupted. Skipping dataset...')
        return(None)
    
    if not AnalogTTLs:
        if 'Events' not in locals():
            print('.kwe/.kwik file is corrupted. Skipping dataset...')
            return(None)
    
    Raw = {}
    for Proc in RawRO.keys():
        Raw[Proc] = {}
        Raw[Proc]['data'] = {Rec: RawRO[Proc]['data'][Rec][:, :]
                             for Rec in RawRO[Proc]['data'].keys()}
        
        Raw[Proc]['channel_bit_volts'] = RawRO[Proc]['channel_bit_volts']
        Raw[Proc]['info'] = RawRO[Proc]['info'].copy()
        Raw[Proc]['timestamps'] = RawRO[Proc]['timestamps'].copy()
    
    del(RawRO)
    
    if Unit == 'uV': Raw = Hdf5F.BitsToVolts(Raw)
            
    print('Data from ', RecFolder, ' loaded.')
    
    if AnalogTTLs: return(Raw, Spks, Files)
    else: return(Raw, Events, Spks, Files)


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


def SetPlot(Backend='Qt5Agg', AxesObj=(), FigObj=(), FigTitle='', Params=False, 
            Plot=False, Axes=False):
    if Params:
        print('Set plot parameters...')
        Params = {'backend': Backend,
                  'text.usetex': True, 'text.latex.unicode': True,
    #              'text.latex.preamble': '\\usepackage{siunitx}',
                  
                  'font.family': 'serif', 'font.serif': 'Computer Modern Roman',
                  'axes.titlesize': 'medium', 'axes.labelsize': 'medium',
                  'xtick.labelsize': 'small', 'xtick.direction': 'out',
                  'ytick.labelsize': 'small', 'ytick.direction': 'out',
                  'legend.fontsize': 'small', 'legend.labelspacing': 0.4,
                  'figure.titlesize': 'large', 'figure.titleweight': 'normal',
                  
                  'image.cmap': 'cubehelix', 'savefig.transparent': True,
                  'svg.fonttype': 'none'}
        return(Params)
    
    elif Plot:
        print('Set plot figure...')
        FigObj.suptitle(FigTitle); FigObj.tight_layout(); 
        FigObj.subplots_adjust(top=0.925)
    
    elif Axes:
        print('Set plot axes...')
        AxesObj.spines['right'].set_visible(False)
        AxesObj.spines['top'].set_visible(False)
        AxesObj.yaxis.set_ticks_position('left')
        AxesObj.xaxis.set_ticks_position('bottom')
        AxesObj.locator_params(tight=True)
    
    else: print("'Params', 'Plot' or 'Axes' must be True.")
    
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
        
            F[Group][RKey].attrs['XValues'] = XValues
    
    if Thrash: 
        for _ in Thrash: print(_)
    
    return(None)


## Level 1
def ResampleData(NewRate, Folder, AnalogTTLs=False, Board='OE', Override={}):
    print('Load DataInfo...')
    if AnalogTTLs: Raw, _, Files = LoadOEKwik(Folder, AnalogTTLs)
    else: Raw, Events, _, Files = LoadOEKwik(Folder, AnalogTTLs)
    
    OEProc = GetProc(Raw, Board)[0]
    
    Path = os.getcwd() + '/' + Folder
    os.makedirs(Path, exist_ok=True)
    
    for Rec in Raw[OEProc]['data'].keys():
        if Override != {}: 
            if 'Rec' in Override.keys():
                Rec = Override['Rec']
        
        Rate = Raw[OEProc]['info'][Rec]['sample_rate']
        RecS = "{0:02d}".format(int(Rec))
        
        print('Resampling data...')
        Data = Raw[OEProc]['data'][Rec][::Rate/NewRate, :-1]
#        Data = [Raw[OEProc]['data'][Rec][:, _] * 
#                Raw[OEProc]['channel_bit_volts'][Rec][_]
#                for _ in range(len(Raw[OEProc]['data'][Rec][0,:-1]))]
        MatName = 'Exp' + Files[OEProc+'_kwd'][-13:-8] + '_' + RecS  + '.mat'
        
        io.savemat(Path+'/'+MatName, {'Data1000Hz': Data})
    
    return(None)


def ClusterizeSpks(Folder, AnalogTTLs=False, Board='OE', Override={}):
    print('Load DataInfo...')
    if AnalogTTLs: Raw, _, Files = LoadOEKwik(Folder, AnalogTTLs)
    else: Raw, Events, _, Files = LoadOEKwik(Folder, AnalogTTLs)
    
    if Raw == {}: print('No data in', Folder, '!'); return(None)
    
    OEProc = GetProc(Raw, Board)[0]
#    OEProc = list(Raw.keys())[0]
    
    Path = os.getcwd() + '/' + Folder +'/SepCh/'
    os.makedirs(Path, exist_ok=True)
    
    Clusters = {}
    for Rec in Raw[OEProc]['data'].keys():
        if Override != {}: 
            if 'Rec' in Override.keys():
                Rec = Override['Rec']
        
        RecS = "{0:02d}".format(int(Rec))
        
        Rate = Raw[OEProc]['info'][Rec]['sample_rate']
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
        ClusterList = [_ for _ in ClusterList if _.split('-')[-2].split('_')[-1] == RecS]
        
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


def LoadUnits(FileName, Override={}):
    with h5py.File(FileName, 'r') as F:
        Key = GetExpKeys('UnitRec', F)
        
        Units, XValues = {}, {}
        for RKey in F[Key].keys():
            if 'XValues' in F[Key][RKey].attrs: 
                XValues[RKey] = F[Key][RKey].attrs['XValues'][:]
                
            elif RKey == 'XValues': XValues[RKey] = F[Key][RKey][:]; continue
            else: XValues[RKey] = []
            
            if Override != {}: 
                if 'Rec' in Override.keys():
                    RKey = "{0:02d}".format(Override['Rec'])
            
            Units[RKey] = {}
            
            for Ch in F[Key][RKey].keys():
                Units[RKey][Ch] = {}
                Path = '/'+Key+'/'+RKey+'/'+Ch
                print('Loading', Path+'...')
                
                if 'Spks' in F[Path].keys():
                    Units[RKey][Ch]['Spks'] = {}
                    for Cluster in F[Path]['Spks'].keys():
                        SpkNo = len(list(F[Path]['Spks'][Cluster].keys()))
                        Units[RKey][Ch]['Spks'][Cluster] = [[] for _ in range(SpkNo)]
                        
                        for Spk in F[Path]['Spks'][Cluster].keys():
                            Units[RKey][Ch]['Spks'][Cluster][int(Spk)] = F[Path]['Spks'][Cluster][Spk][:]
                    
                    if F[Path]['Spks'].attrs.keys():
                        Units[RKey][Ch]['Spks_Info'] = {}
                        for VarKey, VarValue in F[Path]['Spks'].attrs.items():
                            if isinstance(VarValue, Number): 
                                Units[RKey][Ch]['Spks_Info'][VarKey] = float(VarValue)
                            else: 
                                Units[RKey][Ch]['Spks_Info'][VarKey] = VarValue
                
                if 'PSTH' in F[Path].keys():
                    Units[RKey][Ch]['PSTH'] = {}
                    
                    for VarKey in F[Path]['PSTH'].keys():
                        Units[RKey][Ch]['PSTH'][VarKey] = F[Path]['PSTH'][VarKey][:]
                    
                    if F[Path]['PSTH'].attrs.keys():
                        Units[RKey][Ch]['PSTH_Info'] = {}
                        for VarKey, VarValue in F[Path]['PSTH'].attrs.items():
                            if isinstance(VarValue, Number): 
                                Units[RKey][Ch]['PSTH_Info'][VarKey] = float(VarValue)
                            else: 
                                Units[RKey][Ch]['PSTH_Info'][VarKey] = VarValue
            
            if Override != {}: 
                if 'Rec' in Override.keys(): break
    
    return(Units, XValues)


def QuantifyTTLsPerRec(Raw, Rec, AnalogTTLs, ChTTL=0, Proc='', StimType=''):
    print('Get TTL timestamps... ', end='')
    if AnalogTTLs:
        TTLCh = Raw[Proc]['data'][Rec][:, ChTTL-1]
        Threshold = CheckStimulationPattern(TTLCh, StimType)
        TTLs = []
        for _ in range(1, len(TTLCh)):
            if TTLCh[_] > Threshold:
                if TTLCh[_-1] < Threshold: TTLs.append(_)
        
        print('Done.')
        return(TTLs)          


def UnitsSpks(Folder, AnalysisFile, Board='OE', Override={}):
    print('Load DataInfo...')
    Raw = LoadOEKwik(Folder, 'Raw')[0]
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
        
#        TTLCh = Raw[OEProc]['data'][Rec][:,16]
#        
#        if np.mean(TTLCh) > 1000: 
#            print('Sinusoidal stimulation')
#        
        Units[RecS] = {}
        for Ch in Clusters[RecS].keys():
            Units[RecS][Ch] = SepSpksPerCluster(Clusters[RecS][Ch], Ch)
        
        WriteUnits(Units, AnalysisFile)
        if Override != {}: 
            if 'Rec' in Override.keys(): break
        
    del(Raw, Clusters)
    
    return(None)


## Level 2
def UnitsPSTH(Folder, AnalysisFile, TTLChNo=0, PSTHTimeBeforeTTL=0, 
                 PSTHTimeAfterTTL=300, BinSize=1, AnalogTTLs=False, 
                 Board='OE', Override={}):
    print('Load DataInfo...')
    if AnalogTTLs: Raw, _, Files = LoadOEKwik(Folder, AnalogTTLs)
    else: Raw, Events, _, Files = LoadOEKwik(Folder, AnalogTTLs)
    
#    Units = {}
    
    Path = os.getcwd() + '/' + Folder
    ClusterFile = Path + '/SpkClusters.hdf5'
    Clusters = Hdf5F.LoadClusters(ClusterFile, KeyInd=-1)
    
    OEProc = GetProc(Raw, Board)[0]
    Freqs = [2, 4, 8, 16, 32, 40]
    RecNo = len(list(Raw[OEProc]['data'].keys()))
    StimType = Folder.split('-')[-1]
    print('[ StimType', StimType, ']')
    
    for Rec in Raw[OEProc]['data'].keys():
        print('[ Rec', Rec, ']')
        if Override != {}: 
            if 'Rec' in Override.keys():
                Rec = Override['Rec']
        
        Rate = Raw[OEProc]['info'][Rec]['sample_rate']
        TTLCh = Raw[OEProc]['data'][Rec][:, TTLChNo-1]
        TTLs = QuantifyTTLsPerRec(Raw, Rec, AnalogTTLs, TTLChNo, OEProc, StimType)
        if RecNo == 6: StimFreq = Freqs[int(Rec)]
        else: print('Guessing frequency...'); StimFreq = 1 / ((TTLs[3] - TTLs[2])/Rate)
        
        if StimType in ['Sin', 'Sq']:
            if StimType == 'Sin':
                URecS = "{0:02d}".format(int(Rec)) + '-Sin-' + str(int(StimFreq)) + 'Hz'
                NoOfSamplesBefore = 0
#                NoOfSamplesAfter = int(round((TTLs[3] - TTLs[2])))
                NoOfSamplesAfter = int(round(1000/StimFreq))
            else:
                URecS = "{0:02d}".format(int(Rec)) + '-Sq-' + str(int(StimFreq)) + 'Hz'
#                NoOfSamplesBefore = int(round((PSTHTimeBeforeTTL*Rate)*10**-3))
#                NoOfSamplesAfter = int(round((PSTHTimeAfterTTL*Rate)*10**-3))
                NoOfSamplesBefore = PSTHTimeBeforeTTL
                NoOfSamplesAfter = PSTHTimeAfterTTL
        else:
            if np.mean(TTLCh) > 1000: 
                URecS = "{0:02d}".format(int(Rec)) + '-Sin-' + str(int(StimFreq)) + 'Hz'
                NoOfSamplesBefore = 0
#                NoOfSamplesAfter = int(round((TTLs[3] - TTLs[2])))
                NoOfSamplesAfter = int(round(1000/StimFreq))
            else:
                URecS = "{0:02d}".format(int(Rec)) + '-Sq-' + str(int(StimFreq)) + 'Hz'
#                NoOfSamplesBefore = int(round((PSTHTimeBeforeTTL*Rate)*10**-3))
#                NoOfSamplesAfter = int(round((PSTHTimeAfterTTL*Rate)*10**-3))
                NoOfSamplesBefore = PSTHTimeBeforeTTL
                NoOfSamplesAfter = PSTHTimeAfterTTL
        
        XValues = np.arange(-NoOfSamplesBefore, NoOfSamplesAfter, BinSize)
        print('[ Expected XValues len:', str(NoOfSamplesBefore+NoOfSamplesAfter))
        print('  Actual value:', str(len(XValues)), ']')
#        NoOfSamples = NoOfSamplesBefore + NoOfSamplesAfter
#        XValues = (range(-NoOfSamplesBefore, 
#                         NoOfSamples-NoOfSamplesBefore)/Rate)*10**3
        
        RecS = "{0:02d}".format(int(Rec))
#        Units[URecS] = {}
        Units = {}; Units[URecS] = {}
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
#                Units[URecS][CKey]['PSTH'][Class] = np.zeros(len(XValues))
                Hist = np.array([])
                for TTL in range(len(TTLs)):
                    Firing = Clusters[RecS][CKey]['Timestamps'][ClassIndex] \
                             - (TTLs[TTL]/(Rate/1000))
                    Firing = Firing[(Firing >= XValues[0]) * 
                                    (Firing < XValues[-1])]
#                    SpkCount = np.histogram(Firing, 
#                                            np.hstack((XValues, len(XValues)/(Rate/1000))))[0]
#                    
#                    Units[URecS][CKey]['PSTH'][Class] = Units[URecS][CKey]['PSTH'][Class] + SpkCount
#                    
#                    del(Firing, SpkCount)
                    Hist = np.concatenate((Hist, Firing)); del(Firing)
                
                Units[URecS][CKey]['PSTH'][Class] = Hist[:]
                del(Hist)
            
            print(CKey+':', str(len(Classes)), 'clusters.')
        
        del(TTLs)
        WriteUnits(Units, AnalysisFile, XValues)
        
        if Override != {}: 
            if 'Rec' in Override.keys(): break
    
    return(None)


def UnitPlotPerCh(ChDict, Ch, XValues, FigName, Ext):
    ClusterNo = len(ChDict['Spks'])
    if ClusterNo == 0: print(Ch, 'had no spikes :('); return(None)
    
    PSTHNo = 0
    for Class in ChDict['PSTH'].keys(): 
        PSTHNo += len(ChDict['PSTH'][Class])
    
    if not PSTHNo:
        print('No Spks in PSTHs of this channel :( Skipping channel...')
        return(None)
    
    Params = SetPlot(Backend='Agg', Params=True)
    from matplotlib import rcParams; rcParams.update(Params)
    from matplotlib import pyplot as plt
    
    Fig, Axes = plt.subplots(ClusterNo,2, figsize=(8, 3*ClusterNo))
    SpksYLabel = 'Voltage [µV]'; SpksXLabel = 'Time [ms]'
    PSTHYLabel = 'Number of spikes in channel'; PSTHXLabel = 'Time [ms]'
#    SpanLabel = 'Sound pulse'
    
    for Class in ChDict['Spks'].keys():
        SpkNo = len(ChDict['Spks'][Class])
        print(str(SpkNo), 'Spks in cluster', Class)
        if len(ChDict['PSTH'][Class]) == 0:
            print('No Spks in PSTH. Skipping class...')
            continue
        else:
            print('Max of', len(ChDict['PSTH'][Class]), 'Spks in PSTH.')
        
        
#                        if not SpkNo:
#                            print('No Spk data on cluster', str(Cluster) + '. Skipping...')
#                            Thrash[SKey+FKey+RKey+Ch] = Units[SKey][FKey][RKey][Ch].copy()
#                            continue
        
#                PSTHPeak = max(Units[SKey][FKey][RKey]['PSTH'][Cluster])
#                PSTHMean = np.mean(Units[SKey][FKey][RKey]['PSTH'][Cluster])
#        if max(ChDict['PSTH'][Class]) < 4: 
#            print('No peaks in PSTH. Skipping cluster', Class, '...')
#            continue
        
        if SpkNo > 50: 
            SpkNo = np.arange(SpkNo)
            np.random.shuffle(SpkNo)
            SpkNo = SpkNo[:50]
        else: 
            SpkNo = np.arange(SpkNo)
        
        
        for Spike in SpkNo:
            x = np.arange(len(ChDict['Spks'][Class][Spike])) / 30
            if ClusterNo == 1: Axes[0].plot(x, ChDict['Spks'][Class][Spike], 'r')
            else: Axes[int(Class)-1][0].plot(x, ChDict['Spks'][Class][Spike], 'r')
        
        x = np.arange(len(np.mean(ChDict['Spks'][Class], axis=0))) / 30
        if ClusterNo == 1:
            Axes[0].plot(x, np.mean(ChDict['Spks'][Class], axis=0), 'k')
            Axes[1].hist(ChDict['PSTH'][Class], XValues)
#            
#            if len(XValues) != len(ChDict['PSTH'][Class]):
#                x = np.arange(len(ChDict['PSTH'][Class])) / 30
#                Axes[1].bar(x, ChDict['PSTH'][Class])
#            else:
#                Axes[1].bar(XValues, ChDict['PSTH'][Class])
            
#            Ind1 = list(XValues).index(0)
#            Ind2 = list(XValues).index(int(PulseDur*1000))
            
#            Axes[1].axvspan(XValues[Ind1], XValues[Ind2], color='k', alpha=0.3, 
#                            lw=0, label=SpanLabel)
            
            SetPlot(AxesObj=Axes[0], Axes=True)
            SetPlot(AxesObj=Axes[1], Axes=True)
            Axes[0].set_ylabel(SpksYLabel); Axes[0].set_xlabel(SpksXLabel)
            Axes[1].set_ylabel(PSTHYLabel); Axes[1].set_xlabel(PSTHXLabel)
            
        else:
#                    Axes[Cluster][0].set_title('Peak='+str(PSTHPeak)+' Mean='+\
#                                               str(PSTHMean)+' Std='+str(PSTHStd))
            Axes[int(Class)-1][0].plot(x, np.mean(ChDict['Spks'][Class], axis=0), 'k')
            Axes[int(Class)-1][1].hist(ChDict['PSTH'][Class], XValues)
            
#            if len(XValues) != len(ChDict['PSTH'][Class]):
#                x = np.arange(len(ChDict['PSTH'][Class])) / 30
#                Axes[int(Class)-1][1].bar(x, ChDict['PSTH'][Class])
#            else:
#                Axes[int(Class)-1][1].bar(XValues, ChDict['PSTH'][Class])
            
#            Ind1 = list(XValues).index(0)
#            Ind2 = list(XValues).index(int(PulseDur*1000))
            
#            Axes[int(Class)-1][1].axvspan(XValues[Ind1], XValues[Ind2], 
#                                          color='k', alpha=0.3, lw=0, 
#                                          label=SpanLabel)
            
            SetPlot(AxesObj=Axes[int(Class)-1][0], Axes=True)
            SetPlot(AxesObj=Axes[int(Class)-1][1], Axes=True)
            Axes[int(Class)-1][0].set_ylabel(SpksYLabel)
            Axes[int(Class)-1][0].set_xlabel(SpksXLabel)
            Axes[int(Class)-1][1].set_ylabel(PSTHYLabel)
            Axes[int(Class)-1][1].set_xlabel(PSTHXLabel)
    
    FigTitle = FigName.split('/')[-1][:-4]
    FigTitle = r'''\verb+''' + FigTitle + "+"
    SetPlot(FigObj=Fig, FigTitle=FigTitle, Plot=True)
    print('Writing to', FigName+'... ', end='')
    Fig.savefig(FigName, format=Ext)
    print('Done.')
    return(None)


def UnitsSpksPSTH_ToVector(AnalysisFile, Ext='svg', Override={}):    
    Units, XValues = LoadUnits(AnalysisFile, Override)
    if 'XValues' in Override: XValues = Override['XValues']
    
    AnalysisFileName = AnalysisFile.split('/')[-1]
    Folder = '/'.join(AnalysisFile.split('/')[:-1]) + '/Figures/'
    os.makedirs(Folder, exist_ok=True)
    
#    Thrash = {}
    for RKey in Units:
        for Ch in Units[RKey]:
            FigName = Folder + AnalysisFileName[:-5] + '_Rec' + RKey + '_' + Ch + '.' + Ext
            UnitPlot = Process(target=UnitPlotPerCh, 
                               args=(Units[RKey][Ch], Ch, XValues[RKey], FigName, Ext))
            UnitPlot.start(); print('PID =', UnitPlot.pid)
            UnitPlot.join()


def VerifyStimType():
    pass


#%% Run
TTLChNo = 17
PSTHTimeBeforeTTL = 150
PSTHTimeAfterTTL = 150
BinSize = 1
AnalogTTLs = True
Override = {}
Ext='pdf'

Folders = glob('./DIO*'); Done = []; Errors = []; ErrorsLog = []; Skipped = []
for Folder in Folders:
    Subfolders = glob(Folder+'/2017*')
    
    for Subfolder in Subfolders:
        if Subfolder in Done or Subfolder in Errors: 
            print('!!!==========')
            print('Already done. Skipping...')
            print('!!!==========')
            continue
        
        try:
            ClusterizeSpks(Subfolder, AnalogTTLs)
#            ResampleData(1000, Subfolder, AnalogTTLs)
#            Done.append(Subfolder)
            
            if Subfolder+'/SpkClusters.hdf5' not in glob(Subfolder + '/*.hdf5'):
                print('!!!==========')
                print('No Clusters!')
                print('!!!==========')
#                Skipped.append(Subfolder)
                continue
            
            AnalysisFile = '-'.join(Subfolder.split('/')[-2:]) + '.hdf5'
            AnalysisFile = Subfolder + '/' + AnalysisFile
            UnitsPSTH(Subfolder, AnalysisFile, TTLChNo, PSTHTimeBeforeTTL, 
                      PSTHTimeAfterTTL, BinSize, AnalogTTLs)
            
            UnitsSpksPSTH_ToVector(AnalysisFile, Ext, Override)
            Done.append(Subfolder)
            
        except Exception as e:
            Errors.append(Subfolder)
            ErrorsLog.append(e)
            print('!!!==========')
            print(e)
            print('!!!==========')


#Ext='pdf'
#Override = {}
#RKey = '01-Sinusoidal-4Hz'
#Ch = 'Ch09'
#Subfolder = 'Animal_899/2016-06-08_21-43-09'
#AnalysisFile = '-'.join(Subfolder.split('/')) + '.hdf5'
#AnalysisFile = Subfolder + '/' + AnalysisFile

#%% Folders with no data
Errors = ['/home/cerebro/Martina/Data/201702/DIOChR2_93/2017-02-22_22-42-59_1700-Sin',
 '/home/cerebro/Martina/Data/201702/DIOChR2_93/2017-02-22_22-41-41_1700-Sin',
 '/home/cerebro/Martina/Data/201702/DIOChR2_93/2017-02-22_22-42-11_1700-Sin',
 '/home/cerebro/Martina/Data/201702/DIOChR2_93/2017-02-22_21-20-21',
 '/home/cerebro/Martina/Data/201702/DIOChR2_93/2017-02-22_22-40-44_1700-Sin',
 '/home/cerebro/Martina/Data/201702/DIOChR2_19/2017-02-23_21-33-15',
 '/home/cerebro/Martina/Data/201702/DIOChR2_91/2017-02-24_23-38-33',
 '/home/cerebro/Martina/Data/201702/DIOChR2_92/2017-02-22_18-15-58_Test',
 '/home/cerebro/Martina/Data/201702/DIOChR2_92/2017-02-22_17-29-16_Test',
 '/home/cerebro/Martina/Data/201702/DIOChR2_92/2017-02-22_18-15-11_Test',
 '/home/cerebro/Martina/Data/201702/DIOChR2_92/2017-02-22_18-14-17_Test',
 '/home/cerebro/Martina/Data/201702/DIOChR2_9/2017-02-21_16-55-17_1500-Sin',
 '/home/cerebro/Martina/Data/201702/DIOChR2_97/2017-02-24_19-14-39',
 '/home/cerebro/Martina/Data/201702/DIOChR2_97/2017-02-24_19-19-37',
 '/home/cerebro/Martina/Data/201702/DIOChR2_97/2017-02-24_19-15-01',
 '/home/cerebro/Martina/Data/201702/DIOChR2_18/2017-02-21_21-56-57_1500-Sin',
 '/home/cerebro/Martina/Data/201702/DIOChR2_18/2017-02-21_21-58-07_1500-Sin',
 '/home/cerebro/Martina/Data/201702/DIOChR2_18/2017-02-22_02-29-37_1554-Sin',
 '/home/cerebro/Martina/Data/201702/DIOChR2_95/2017-02-23_14-49-46_1700-Sq',
 '/home/cerebro/Martina/Data/201702/DIOChR2_95/2017-02-23_14-49-20_1700-Sq',
 '/home/cerebro/Martina/Data/201702/DIOChR2_90/2017-02-24_14-04-14_900-Sq',
 './DIOChR2_91/2017-02-25_01-13-24_2300-Sq_3freq',
 './DIOChR2_95/2017-02-23_18-35-23_1800-Sq_5freq',
 './DIOChR2_90/2017-02-24_15-48-02_1700-Sin_5freq']
