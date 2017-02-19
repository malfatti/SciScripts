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
from numbers import Number
from scipy import io
from subprocess import call

# Level 0
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
        
        if len(XValues): F[Group]['XValues'] = XValues
    
    if Thrash: 
        for _ in Thrash: print(_)
    
    return(None)


## Level 1
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


def LoadUnits(FileName, Override={}):
    XValues = []
    with h5py.File(FileName, 'r') as F:
        Key = GetExpKeys('UnitRec', F)
        
        Units = {}
        for RKey in F[Key].keys():
            if RKey == 'XValues':
                XValues = F[Key][RKey][:]
                continue
            
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


## Level 2
def UnitsPSTH(Folder, AnalysisFile, TTLChNo=0, PSTHTimeBeforeTTL=0, 
                 PSTHTimeAfterTTL=300, AnalogTTLs=False, 
                 Board='OE', Override={}):
    print('Load DataInfo...')
    if AnalogTTLs: Raw, _, Files = Hdf5F.LoadOEKwik(Folder, AnalogTTLs)
    else: Raw, Events, _, Files = Hdf5F.LoadOEKwik(Folder, AnalogTTLs)
    
    Units = {}
    
    Path = os.getcwd() + '/' + Folder
    ClusterFile = Path + '/SpkClusters.hdf5'
    Clusters = Hdf5F.LoadClusters(ClusterFile, KeyInd=-1)
    
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


def UnitPlotPerCh(ChDict, Ch, XValues, FigName, Ext):
    ClusterNo = len(ChDict['Spks'])
    if ClusterNo == 0: print(Ch, 'had no spikes :('); return(None)
    
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
        print('Max of', max(ChDict['PSTH'][Class]),  
              'Spks in PSTH')
        
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
            
            if len(XValues) != len(ChDict['PSTH'][Class]):
                x = np.arange(len(ChDict['PSTH'][Class])) / 30
                Axes[1].bar(x, ChDict['PSTH'][Class])
            else:
                Axes[1].bar(XValues, ChDict['PSTH'][Class])
            
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
            
            if len(XValues) != len(ChDict['PSTH'][Class]):
                x = np.arange(len(ChDict['PSTH'][Class])) / 30
                Axes[int(Class)-1][1].bar(x, ChDict['PSTH'][Class])
            else:
                Axes[int(Class)-1][1].bar(XValues, ChDict['PSTH'][Class])
            
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
    AnalysisFileName = AnalysisFile.split('/')[-1]
    Folder = '/'.join(AnalysisFile.split('/')[:-1]) + '/Figures/'
    os.makedirs(Folder, exist_ok=True)
    
#    Thrash = {}
    for RKey in Units:
        for Ch in Units[RKey]:
            FigName = Folder + AnalysisFileName[:-5] + '_Rec' + RKey + '_' + Ch + '.' + Ext
            UnitPlot = Process(target=UnitPlotPerCh, 
                               args=(Units[RKey][Ch], Ch, XValues, FigName, Ext))
            UnitPlot.start(); print('PID =', UnitPlot.pid)
            UnitPlot.join()


def VerifyStimType():
    pass


#%% Run
TTLChNo = 17
PSTHTimeBeforeTTL = 150
PSTHTimeAfterTTL = 150
AnalogTTLs = True
Override = {}
Ext='pdf'

Folders = glob('Animal*'); Done = []; Errors = []; ErrorsLog = []
for Folder in Folders:
    Subfolders = glob(Folder+'/2016*')
    
    for Subfolder in Subfolders:
        if Subfolder in Done or Subfolder in Errors: 
            print('!!!==========')
            print('Already done. Skipping...')
            print('!!!==========')
            continue
        
        if glob(Subfolder + '/SpkClusters.hdf5') == []:
            print('!!!==========')
            print('No SpkClusters.hdf5 file. Skipping...')
            print('!!!==========')
            continue
        
        try:
#            ClusterizeSpks(Subfolder, AnalogTTLs)
#            Done.append(Subfolder)
            
            AnalysisFile = '-'.join(Subfolder.split('/')) + '.hdf5'
            AnalysisFile = Subfolder + '/' + AnalysisFile
#            UnitsPSTH(Subfolder, AnalysisFile, TTLChNo, PSTHTimeBeforeTTL, 
#                      PSTHTimeAfterTTL, AnalogTTLs)
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

###
plt.figure()
for s in range(len(Units['00-Sinusoidal-2Hz']['Ch09']['Spks']['01'])):
    plt.plot(Units['00-Sinusoidal-2Hz']['Ch09']['Spks']['01'][s])

plt.figure(); plt.hist(XValues, Units['00-Sinusoidal-2Hz']['Ch09']['PSTH']['01'])
plt.show()

#%% Keeping track

ClusterDone = ['Animal_366/2016-06-12_17-07-03',
 'Animal_366/2016-06-12_16-49-22',
 'Animal_366/2016-06-12_17-32-18',
 'Animal_366/2016-06-12_17-21-19',
 'Animal_366/2016-06-12_16-30-17',
 'Animal_366/2016-06-12_15-49-34',
 'Animal_366/2016-06-12_15-37-47',
 'Animal_366/2016-06-12_16-16-35',
 'Animal_366/2016-06-12_15-24-57',
 'Animal_366/2016-06-12_16-42-17',
 'Animal_366/2016-06-12_16-02-49',
 'Animal_369/2016-06-02_18-34-28SquareDV191',
 'Animal_369/2016-06-02_19-30-13SquareDV217',
 'Animal_369/2016-06-02_18-52-18SinDV191',
 'Animal_369/2016-06-02_19-44-12SinDV217',
 'Animal_689/2016-06-04_21-34-31',
 'Animal_689/2016-06-04_16-26-21',
 'Animal_689/2016-06-04_19-38-58',
 'Animal_689/2016-06-04_21-21-12',
 'Animal_689/2016-06-04_17-59-09',
 'Animal_689/2016-06-04_19-52-00',
 'Animal_689/2016-06-04_20-31-51',
 'Animal_689/2016-06-04_20-19-10',
 'Animal_897/2016-06-09_19-35-36',
 'Animal_897/2016-06-09_19-06-50',
 'Animal_897/2016-06-09_20-11-32',
 'Animal_897/2016-06-09_16-40-33',
 'Animal_897/2016-06-09_18-51-34',
 'Animal_897/2016-06-09_16-51-46',
 'Animal_897/2016-06-09_18-32-44',
 'Animal_897/2016-06-09_16-46-55',
 'Animal_897/2016-06-09_16-17-58',
 'Animal_897/2016-06-09_18-27-40',
 'Animal_897/2016-06-09_17-50-02',
 'Animal_897/2016-06-09_20-00-43',
 'Animal_897/2016-06-09_19-25-09',
 'Animal_897/2016-06-09_18-56-43',
 'Animal_897/2016-06-09_17-04-22',
 'Animal_900/2016-06-05_22-56-25',
 'Animal_900/2016-06-06_00-06-35',
 'Animal_900/2016-06-05_23-07-50',
 'Animal_900/2016-06-05_21-36-23',
 'Animal_900/2016-06-05_23-33-56',
 'Animal_900/2016-06-05_23-23-02',
 'Animal_875/2016-06-05_00-59-08',
 'Animal_875/2016-06-05_00-02-10',
 'Animal_875/2016-06-05_00-19-04',
 'Animal_875/2016-06-05_00-30-24',
 'Animal_875/2016-06-04_22-22-34',
 'Animal_875/2016-06-04_23-11-41',
 'Animal_875/2016-06-04_23-46-13',
 'Animal_875/2016-06-04_23-57-12',
 'Animal_875/2016-06-04_23-36-57',
 'Animal_875/2016-06-04_23-23-09',
 'Animal_875/2016-06-04_22-02-32',
 'Animal_875/2016-06-04_22-15-21',
 'Animal_870/2016-06-11_20-16-15',
 'Animal_870/2016-06-11_19-53-24',
 'Animal_870/2016-06-11_22-22-23',
 'Animal_870/2016-06-11_19-54-45',
 'Animal_870/2016-06-11_22-04-57',
 'Animal_870/2016-06-11_21-55-06',
 'Animal_870/2016-06-11_19-54-20',
 'Animal_870/2016-06-11_23-07-51',
 'Animal_870/2016-06-11_20-18-48',
 'Animal_870/2016-06-11_20-30-03',
 'Animal_870/2016-06-11_19-09-43',
 'Animal_870/2016-06-11_19-34-34',
 'Animal_899/2016-06-09_00-04-31',
 'Animal_899/2016-06-09_01-21-05',
 'Animal_899/2016-06-09_01-15-06',
 'Animal_899/2016-06-08_20-46-01',
 'Animal_899/2016-06-08_22-15-44',
 'Animal_899/2016-06-09_01-31-04',
 'Animal_899/2016-06-08_21-05-28',
 'Animal_899/2016-06-08_19-45-14',
 'Animal_899/2016-06-09_00-55-04',
 'Animal_899/2016-06-08_23-10-50',
 'Animal_899/2016-06-09_00-18-19',
 'Animal_899/2016-06-09_00-44-21',
 'Animal_899/2016-06-08_22-31-43',
 'Animal_899/2016-06-08_21-43-09',
 'Animal_899/2016-06-08_22-08-35',
 'Animal_830/2016-06-09_22-40-30',
 'Animal_830/2016-06-09_22-30-11',
 'Animal_830/2016-06-09_22-19-25',
 'Animal_363/2016-06-03_22-30-28',
 'Animal_363/2016-06-03_19-26-04',
 'Animal_363/2016-06-03_21-48-29',
 'Animal_363/2016-06-03_22-18-55',
 'Animal_363/2016-06-03_21-35-56',
 'Animal_363/2016-06-03_17-07-19',
 'Animal_363/2016-06-03_19-49-07',
 'Animal_868/2016-06-12_03-03-56',
 'Animal_868/2016-06-12_04-09-28',
 'Animal_868/2016-06-12_03-36-44',
 'Animal_868/2016-06-12_03-45-13',
 'Animal_868/2016-06-12_04-11-09',
 'Animal_868/2016-06-12_05-49-45',
 'Animal_868/2016-06-12_03-15-29',
 'Animal_868/2016-06-12_03-05-32',
 'Animal_868/2016-06-12_03-35-06',
 'Animal_868/2016-06-12_04-25-47',
 'Animal_868/2016-06-12_06-22-18',
 'Animal_873/2016-06-12_00-11-46',
 'Animal_873/2016-06-12_00-43-36',
 'Animal_873/2016-06-12_00-13-51',
 'Animal_873/2016-06-12_00-31-20',
 'Animal_873/2016-06-11_23-31-43',
 'Animal_873/2016-06-11_23-11-47',
 'Animal_873/2016-06-12_00-01-26',
 'Animal_873/2016-06-12_00-41-58',
 'Animal_873/2016-06-11_23-57-54',
 'Animal_872/2016-06-12_08-03-29',
 'Animal_872/2016-06-12_07-30-39',
 'Animal_364/2016-06-10_16-36-57',
 'Animal_364/2016-06-10_16-35-33',
 'Animal_364/2016-06-10_16-40-53',
 'Animal_364/2016-06-10_17-49-40',
 'Animal_364/2016-06-10_18-02-32',
 'Animal_364/2016-06-10_18-25-41',
 'Animal_364/2016-06-10_18-35-11',
 'Animal_364/2016-06-10_18-28-11',
 'Animal_365/2016-06-11_01-38-37',
 'Animal_365/2016-06-11_00-39-22',
 'Animal_365/2016-06-11_01-11-59',
 'Animal_871/2016-06-12_13-34-25',
 'Animal_871/2016-06-12_11-57-12',
 'Animal_871/2016-06-12_12-24-56',
 'Animal_871/2016-06-12_13-41-43',
 'Animal_871/2016-06-12_13-46-24',
 'Animal_871/2016-06-12_10-06-36',
 'Animal_871/2016-06-12_12-32-38',
 'Animal_871/2016-06-12_10-17-31',
 'Animal_871/2016-06-12_12-42-52',
 'Animal_829/2016-06-10_01-26-43',
 'Animal_829/2016-06-10_03-36-04',
 'Animal_829/2016-06-10_02-23-47',
 'Animal_829/2016-06-10_01-55-28',
 'Animal_829/2016-06-10_03-24-53',
 'Animal_829/2016-06-10_02-06-06',
 'Animal_829/2016-06-10_00-26-42',
 'Animal_829/2016-06-10_02-34-23',
 'Animal_898/2016-06-05_15-51-20',
 'Animal_898/2016-06-05_17-39-08',
 'Animal_898/2016-06-05_15-40-07',
 'Animal_898/2016-06-05_16-38-09',
 'Animal_898/2016-06-05_18-39-44',
 'Animal_898/2016-06-05_15-06-47',
 'Animal_898/2016-06-05_17-36-11',
 'Animal_898/2016-06-05_19-13-55',
 'Animal_898/2016-06-05_18-11-33',
 'Animal_898/2016-06-05_14-15-44',
 'Animal_898/2016-06-05_18-18-39',
 'Animal_898/2016-06-05_19-51-34',
 'Animal_898/2016-06-05_16-19-01',
 'Animal_898/2016-06-05_17-01-53',
 'Animal_898/2016-06-05_13-43-18',
 'Animal_898/2016-06-05_16-07-10']
 
ClusterErrors = ['Animal_366/2016-06-12_14-46-12',
 'Animal_689/2016-06-04_18-40-41',
 'Animal_689/2016-06-04_15-58-22',
 'Animal_897/2016-06-09_14-48-53',
 'Animal_897/2016-06-09_13-38-32',
 'Animal_897/2016-06-09_14-36-49',
 'Animal_897/2016-06-09_14-04-10',
 'Animal_897/2016-06-09_13-34-22',
 'Animal_897/2016-06-09_18-27-07',
 'Animal_897/2016-06-09_18-01-40',
 'Animal_897/2016-06-09_13-39-27',
 'Animal_897/2016-06-09_18-40-24',
 'Animal_897/2016-06-09_14-41-25',
 'Animal_897/2016-06-09_13-36-33',
 'Animal_900/2016-06-05_21-35-07',
 'Animal_900/2016-06-05_21-21-05',
 'Animal_900/2016-06-06_00-58-46',
 'Animal_900/2016-06-05_21-09-15',
 'Animal_900/2016-06-06_00-04-20',
 'Animal_875/2016-06-04_23-00-39',
 'Animal_875/2016-06-04_23-02-57',
 'Animal_875/2016-06-04_22-55-46',
 'Animal_870/2016-06-11_19-06-56',
 'Animal_870/2016-06-11_22-53-26',
 'Animal_870/2016-06-11_19-05-19',
 'Animal_870/2016-06-11_22-33-10',
 'Animal_899/2016-06-08_19-18-18',
 'Animal_899/2016-06-08_19-40-38',
 'Animal_899/2016-06-08_19-20-30',
 'Animal_899/2016-06-08_19-30-11',
 'Animal_899/2016-06-08_19-12-59',
 'Animal_899/2016-06-08_19-17-05',
 'Animal_899/2016-06-08_19-10-35',
 'Animal_899/2016-06-08_19-13-20',
 'Animal_830/2016-06-09_21-38-34',
 'Animal_830/2016-06-09_21-52-37',
 'Animal_830/2016-06-09_23-04-14',
 'Animal_830/2016-06-09_21-41-06',
 'Animal_830/2016-06-09_21-36-00',
 'Animal_830/2016-06-09_21-48-24',
 'Animal_830/2016-06-09_21-58-31',
 'Animal_830/2016-06-09_21-46-50',
 'Animal_363/2016-06-03_17-06-07',
 'Animal_868/2016-06-12_02-35-48',
 'Animal_868/2016-06-12_04-28-26',
 'Animal_868/2016-06-12_05-32-25',
 'Animal_873/2016-06-11_22-33-59',
 'Animal_873/2016-06-12_01-08-44',
 'Animal_873/2016-06-11_22-04-13',
 'Animal_873/2016-06-12_00-08-33',
 'Animal_873/2016-06-12_00-07-54',
 'Animal_873/2016-06-11_23-41-53',
 'Animal_873/2016-06-11_22-44-47',
 'Animal_872/2016-06-12_07-13-19',
 'Animal_364/2016-06-10_16-00-49',
 'Animal_364/2016-06-10_15-39-44',
 'Animal_364/2016-06-10_17-05-57',
 'Animal_364/2016-06-10_16-51-37',
 'Animal_364/2016-06-10_15-36-24',
 'Animal_365/2016-06-11_02-10-09',
 'Animal_365/2016-06-11_00-25-57',
 'Animal_365/2016-06-11_02-10-46',
 'Animal_365/2016-06-11_00-22-11',
 'Animal_365/2016-06-11_01-24-10',
 'Animal_871/2016-06-12_13-20-08',
 'Animal_871/2016-06-12_10-45-40',
 'Animal_871/2016-06-12_13-08-23',
 'Animal_871/2016-06-12_08-57-10',
 'Animal_871/2016-06-12_09-35-09',
 'Animal_871/2016-06-12_13-06-49',
 'Animal_871/2016-06-12_09-17-33',
 'Animal_871/2016-06-12_10-59-49',
 'Animal_871/2016-06-12_11-47-08',
 'Animal_871/2016-06-12_09-47-55',
 'Animal_829/2016-06-10_00-14-44',
 'Animal_829/2016-06-10_00-22-23',
 'Animal_898/2016-06-05_13-42-18',
 'Animal_898/2016-06-05_13-39-05',
 'Animal_898/2016-06-05_17-08-23']
 
ClusterErrorsLog = [OSError("Can't read data (Wrong b-tree signature)"),
 UnboundLocalError("local variable 'OEProc' referenced before assignment"),
 UnboundLocalError("local variable 'OEProc' referenced before assignment"),
 UnboundLocalError("local variable 'OEProc' referenced before assignment"),
 UnboundLocalError("local variable 'OEProc' referenced before assignment"),
 UnboundLocalError("local variable 'OEProc' referenced before assignment"),
 UnboundLocalError("local variable 'OEProc' referenced before assignment"),
 UnboundLocalError("local variable 'OEProc' referenced before assignment"),
 UnboundLocalError("local variable 'OEProc' referenced before assignment"),
 UnboundLocalError("local variable 'OEProc' referenced before assignment"),
 UnboundLocalError("local variable 'OEProc' referenced before assignment"),
 UnboundLocalError("local variable 'OEProc' referenced before assignment"),
 UnboundLocalError("local variable 'OEProc' referenced before assignment"),
 UnboundLocalError("local variable 'OEProc' referenced before assignment"),
 UnboundLocalError("local variable 'OEProc' referenced before assignment"),
 UnboundLocalError("local variable 'OEProc' referenced before assignment"),
 ValueError('Did not fully consume compressed contents of an miCOMPRESSED element. This can indicate that the .mat file is corrupted.'),
 UnboundLocalError("local variable 'OEProc' referenced before assignment"),
 UnboundLocalError("local variable 'OEProc' referenced before assignment"),
 UnboundLocalError("local variable 'OEProc' referenced before assignment"),
 UnboundLocalError("local variable 'OEProc' referenced before assignment"),
 UnboundLocalError("local variable 'OEProc' referenced before assignment"),
 UnboundLocalError("local variable 'OEProc' referenced before assignment"),
 UnboundLocalError("local variable 'OEProc' referenced before assignment"),
 UnboundLocalError("local variable 'OEProc' referenced before assignment"),
 UnboundLocalError("local variable 'OEProc' referenced before assignment"),
 UnboundLocalError("local variable 'OEProc' referenced before assignment"),
 UnboundLocalError("local variable 'OEProc' referenced before assignment"),
 UnboundLocalError("local variable 'OEProc' referenced before assignment"),
 UnboundLocalError("local variable 'OEProc' referenced before assignment"),
 UnboundLocalError("local variable 'OEProc' referenced before assignment"),
 UnboundLocalError("local variable 'OEProc' referenced before assignment"),
 UnboundLocalError("local variable 'OEProc' referenced before assignment"),
 UnboundLocalError("local variable 'OEProc' referenced before assignment"),
 UnboundLocalError("local variable 'OEProc' referenced before assignment"),
 UnboundLocalError("local variable 'OEProc' referenced before assignment"),
 UnboundLocalError("local variable 'OEProc' referenced before assignment"),
 UnboundLocalError("local variable 'OEProc' referenced before assignment"),
 UnboundLocalError("local variable 'OEProc' referenced before assignment"),
 UnboundLocalError("local variable 'OEProc' referenced before assignment"),
 UnboundLocalError("local variable 'OEProc' referenced before assignment"),
 UnboundLocalError("local variable 'OEProc' referenced before assignment"),
 UnboundLocalError("local variable 'OEProc' referenced before assignment"),
 UnboundLocalError("local variable 'OEProc' referenced before assignment"),
 KeyError('Unable to open object (Unknown flag for message)'),
 OSError("Can't read data (Wrong b-tree signature)"),
 OSError("Can't read data (Wrong b-tree signature)"),
 UnboundLocalError("local variable 'OEProc' referenced before assignment"),
 OSError("Can't read data (Wrong b-tree signature)"),
 UnboundLocalError("local variable 'OEProc' referenced before assignment"),
 UnboundLocalError("local variable 'OEProc' referenced before assignment"),
 OSError("Can't read data (Wrong b-tree signature)"),
 KeyError('Unable to open object (Unknown flag for message)'),
 UnboundLocalError("local variable 'OEProc' referenced before assignment"),
 UnboundLocalError("local variable 'OEProc' referenced before assignment"),
 UnboundLocalError("local variable 'OEProc' referenced before assignment"),
 ValueError('Index (1) out of range (0--1)'),
 OSError("Can't read data (Wrong b-tree signature)"),
 UnboundLocalError("local variable 'OEProc' referenced before assignment"),
 UnboundLocalError("local variable 'OEProc' referenced before assignment"),
 UnboundLocalError("local variable 'OEProc' referenced before assignment"),
 UnboundLocalError("local variable 'OEProc' referenced before assignment"),
 UnboundLocalError("local variable 'OEProc' referenced before assignment"),
 OSError("Can't read data (Wrong b-tree signature)"),
 UnboundLocalError("local variable 'OEProc' referenced before assignment"),
 OSError("Can't read data (Wrong b-tree signature)"),
 OSError("Can't read data (Wrong b-tree signature)"),
 UnboundLocalError("local variable 'OEProc' referenced before assignment"),
 OSError("Can't read data (Wrong b-tree signature)"),
 OSError("Can't read data (Wrong b-tree signature)"),
 UnboundLocalError("local variable 'OEProc' referenced before assignment"),
 OSError("Can't read data (Wrong b-tree signature)"),
 OSError("Can't read data (Wrong b-tree signature)"),
 OSError("Can't read data (Wrong b-tree signature)"),
 UnboundLocalError("local variable 'OEProc' referenced before assignment"),
 UnboundLocalError("local variable 'OEProc' referenced before assignment"),
 UnboundLocalError("local variable 'OEProc' referenced before assignment"),
 UnboundLocalError("local variable 'OEProc' referenced before assignment"),
 UnboundLocalError("local variable 'OEProc' referenced before assignment")]

PSTHDone = ['Animal_366/2016-06-12_17-07-03',
 'Animal_366/2016-06-12_17-32-18',
 'Animal_366/2016-06-12_17-21-19',
 'Animal_366/2016-06-12_16-30-17',
 'Animal_366/2016-06-12_15-49-34',
 'Animal_366/2016-06-12_15-37-47',
 'Animal_366/2016-06-12_16-16-35',
 'Animal_366/2016-06-12_16-02-49',
 'Animal_369/2016-06-02_18-34-28SquareDV191',
 'Animal_369/2016-06-02_19-30-13SquareDV217',
 'Animal_369/2016-06-02_18-52-18SinDV191',
 'Animal_369/2016-06-02_19-44-12SinDV217',
 'Animal_689/2016-06-04_21-34-31',
 'Animal_689/2016-06-04_16-26-21',
 'Animal_689/2016-06-04_19-38-58',
 'Animal_689/2016-06-04_21-21-12',
 'Animal_689/2016-06-04_17-59-09',
 'Animal_689/2016-06-04_19-52-00',
 'Animal_689/2016-06-04_20-31-51',
 'Animal_689/2016-06-04_20-19-10',
 'Animal_897/2016-06-09_19-35-36',
 'Animal_897/2016-06-09_19-06-50',
 'Animal_897/2016-06-09_20-11-32',
 'Animal_897/2016-06-09_16-40-33',
 'Animal_897/2016-06-09_18-51-34',
 'Animal_897/2016-06-09_16-51-46',
 'Animal_897/2016-06-09_18-32-44',
 'Animal_897/2016-06-09_16-46-55',
 'Animal_897/2016-06-09_18-27-40',
 'Animal_897/2016-06-09_17-50-02',
 'Animal_897/2016-06-09_20-00-43',
 'Animal_897/2016-06-09_19-25-09',
 'Animal_897/2016-06-09_18-56-43',
 'Animal_897/2016-06-09_17-04-22',
 'Animal_900/2016-06-05_22-56-25',
 'Animal_900/2016-06-06_00-06-35',
 'Animal_900/2016-06-05_23-07-50',
 'Animal_900/2016-06-05_21-36-23',
 'Animal_900/2016-06-05_23-33-56',
 'Animal_900/2016-06-05_23-23-02',
 'Animal_875/2016-06-05_00-59-08',
 'Animal_875/2016-06-05_00-02-10',
 'Animal_875/2016-06-05_00-19-04',
 'Animal_875/2016-06-05_00-30-24',
 'Animal_875/2016-06-04_22-22-34',
 'Animal_875/2016-06-04_23-11-41',
 'Animal_875/2016-06-04_23-46-13',
 'Animal_875/2016-06-04_23-57-12',
 'Animal_875/2016-06-04_23-36-57',
 'Animal_875/2016-06-04_23-23-09',
 'Animal_875/2016-06-04_22-02-32',
 'Animal_870/2016-06-11_20-16-15',
 'Animal_870/2016-06-11_19-53-24',
 'Animal_870/2016-06-11_22-22-23',
 'Animal_870/2016-06-11_22-04-57',
 'Animal_870/2016-06-11_21-55-06',
 'Animal_870/2016-06-11_19-54-20',
 'Animal_870/2016-06-11_23-07-51',
 'Animal_870/2016-06-11_20-18-48',
 'Animal_870/2016-06-11_20-30-03',
 'Animal_870/2016-06-11_19-09-43',
 'Animal_870/2016-06-11_19-34-34',
 'Animal_899/2016-06-09_00-04-31',
 'Animal_899/2016-06-09_01-21-05',
 'Animal_899/2016-06-09_01-15-06',
 'Animal_899/2016-06-08_20-46-01',
 'Animal_899/2016-06-08_22-15-44',
 'Animal_899/2016-06-09_01-31-04',
 'Animal_899/2016-06-08_21-05-28',
 'Animal_899/2016-06-09_00-55-04',
 'Animal_899/2016-06-09_00-18-19',
 'Animal_899/2016-06-09_00-44-21',
 'Animal_899/2016-06-08_22-31-43',
 'Animal_899/2016-06-08_21-43-09',
 'Animal_899/2016-06-08_22-08-35',
 'Animal_830/2016-06-09_22-40-30',
 'Animal_830/2016-06-09_22-30-11',
 'Animal_830/2016-06-09_22-19-25',
 'Animal_363/2016-06-03_22-30-28',
 'Animal_363/2016-06-03_19-26-04',
 'Animal_363/2016-06-03_21-48-29',
 'Animal_363/2016-06-03_21-35-56',
 'Animal_363/2016-06-03_17-07-19',
 'Animal_363/2016-06-03_19-49-07',
 'Animal_868/2016-06-12_03-03-56',
 'Animal_868/2016-06-12_04-09-28',
 'Animal_868/2016-06-12_03-36-44',
 'Animal_868/2016-06-12_04-11-09',
 'Animal_868/2016-06-12_05-49-45',
 'Animal_868/2016-06-12_03-15-29',
 'Animal_868/2016-06-12_03-05-32',
 'Animal_868/2016-06-12_03-35-06',
 'Animal_868/2016-06-12_04-25-47',
 'Animal_868/2016-06-12_06-22-18',
 'Animal_873/2016-06-12_00-11-46',
 'Animal_873/2016-06-12_00-43-36',
 'Animal_873/2016-06-12_00-13-51',
 'Animal_873/2016-06-12_00-31-20',
 'Animal_873/2016-06-11_23-31-43',
 'Animal_873/2016-06-12_00-01-26',
 'Animal_873/2016-06-12_00-41-58',
 'Animal_873/2016-06-11_23-57-54',
 'Animal_872/2016-06-12_08-03-29',
 'Animal_872/2016-06-12_07-30-39',
 'Animal_364/2016-06-10_16-36-57',
 'Animal_364/2016-06-10_16-35-33',
 'Animal_364/2016-06-10_16-40-53',
 'Animal_364/2016-06-10_17-49-40',
 'Animal_364/2016-06-10_18-02-32',
 'Animal_364/2016-06-10_18-25-41',
 'Animal_364/2016-06-10_18-35-11',
 'Animal_364/2016-06-10_18-28-11',
 'Animal_365/2016-06-11_01-38-37',
 'Animal_365/2016-06-11_00-39-22',
 'Animal_365/2016-06-11_01-11-59',
 'Animal_871/2016-06-12_13-34-25',
 'Animal_871/2016-06-12_11-57-12',
 'Animal_871/2016-06-12_12-24-56',
 'Animal_871/2016-06-12_13-41-43',
 'Animal_871/2016-06-12_13-46-24',
 'Animal_871/2016-06-12_10-06-36',
 'Animal_871/2016-06-12_12-32-38',
 'Animal_871/2016-06-12_12-42-52',
 'Animal_829/2016-06-10_01-26-43',
 'Animal_829/2016-06-10_03-36-04',
 'Animal_829/2016-06-10_02-23-47',
 'Animal_829/2016-06-10_01-55-28',
 'Animal_829/2016-06-10_03-24-53',
 'Animal_829/2016-06-10_02-06-06',
 'Animal_829/2016-06-10_00-26-42',
 'Animal_829/2016-06-10_02-34-23',
 'Animal_898/2016-06-05_15-51-20',
 'Animal_898/2016-06-05_17-39-08',
 'Animal_898/2016-06-05_15-40-07',
 'Animal_898/2016-06-05_16-38-09',
 'Animal_898/2016-06-05_18-39-44',
 'Animal_898/2016-06-05_15-06-47',
 'Animal_898/2016-06-05_17-36-11',
 'Animal_898/2016-06-05_19-13-55',
 'Animal_898/2016-06-05_18-11-33',
 'Animal_898/2016-06-05_14-15-44',
 'Animal_898/2016-06-05_18-18-39',
 'Animal_898/2016-06-05_19-51-34',
 'Animal_898/2016-06-05_16-19-01',
 'Animal_898/2016-06-05_17-01-53',
 'Animal_898/2016-06-05_13-43-18',
 'Animal_898/2016-06-05_16-07-10']

PSTHErrors = ['Animal_366/2016-06-12_16-49-22',
 'Animal_366/2016-06-12_15-24-57',
 'Animal_366/2016-06-12_16-42-17',
 'Animal_897/2016-06-09_16-17-58',
 'Animal_875/2016-06-04_22-15-21',
 'Animal_870/2016-06-11_19-54-45',
 'Animal_899/2016-06-08_19-45-14',
 'Animal_899/2016-06-08_23-10-50',
 'Animal_363/2016-06-03_22-18-55',
 'Animal_868/2016-06-12_03-45-13',
 'Animal_873/2016-06-11_23-11-47',
 'Animal_871/2016-06-12_10-17-31']

PSTHErrorsLog = [IndexError('list index out of range'),
 IndexError('list index out of range'),
 IndexError('list index out of range'),
 IndexError('list index out of range'),
 IndexError('list index out of range'),
 IndexError('list index out of range'),
 KeyError('04'),
 IndexError('list index out of range'),
 IndexError('list index out of range'),
 IndexError('list index out of range'),
 IndexError('list index out of range'),
 IndexError('list index out of range')]