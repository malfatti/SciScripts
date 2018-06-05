#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: T. Malfatti <malfatti@disroot.org>
@date: 2018-02-26
@license: GNU GPLv3 <https://raw.githubusercontent.com/malfatti/SciScripts/master/LICENSE>
@homepage: https://github.com/Malfatti/SciScripts
"""

import numpy as np
import os

from DataAnalysis import DataAnalysis, Stats
from DataAnalysis.Plot import Plot
from glob import glob
from IO import Asdf, Bin, Hdf5, IO, Klusta, OpenEphys, Txt
from itertools import combinations
from klusta.kwik import KwikModel

PltParams = Plot.Set(Params=True)
from matplotlib import rcParams; rcParams.update(PltParams)
from matplotlib.gridspec import GridSpec, GridSpecFromSubplotSpec
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import pyplot as plt

#Ind = UnitRec['RI']/UnitRec['RISurr'] > 0
#Ind = (UnitRec['UnitId'] == 822) * (UnitRec['Freq'] == '8000-10000')
#      
#FR = np.array(UnitRec['FiringRate'])[Ind]
#
#for fr in FR: plt.plot(fr)
#plt.show()

#%% Detection and clustering
# Board = 'OE'
# TimeBeforeTTL = 0; TimeAfterTTL = 300; BinSize = 3

# NNAdaptor = [13, 12, 14, 11, 15, 10, 16, 9, 5, 4, 6, 3, 7, 2, 8, 1]
# CustomAdaptor = [5, 6, 7, 8, 9, 10 ,11, 12, 13, 14, 15, 16, 1, 2, 3, 4]
# Probe = {'Tip': [9, 8, 10, 7, 13, 4, 12, 5, 15, 2, 16, 1, 14, 3, 11, 6],
#          'Head': [8, 7, 6, 5, 4, 3, 2, 1, 9, 10, 11, 12, 13, 14, 15, 16],
#          'Adaptor': [13, 12, 14, 11, 15, 10, 16, 9, 5, 4, 6, 3, 7, 2, 8, 1],}
# ProbeSpacing = 50

# #Folders = glob('*UnitRec/**/*.kwd', recursive=True) + glob('*/**/*UnitRec/**/*.kwd', recursive=True)
# #Folders = ['/'.join(F.split('/')[:-2]) for F in Folders]
# #Folders = DataAnalysis.UniqueStr(Folders)

# Folders = [
#         #'DreaddsValidation/CaMKIIahM3Dn01-20150903-UnitRec/OpenEphysFiles',
#         #'Prevention/20170803-Prevention_A4-UnitRec/DatFiles',
#         'Prevention/20170815-Prevention_A5-UnitRec/KwikFiles',
#         'DreaddsValidation/CaMKIIahM3Dn02-20150921-UnitRec/OpenEphysFiles',
#         #'DreaddsValidation/CaMKIIahM4Dn04-20151012-UnitRec/OpenEphysFiles',
#         'EarBarTest/EarBarTest_01-20170216-UnitRec/KwikFiles',
#         'EarBarTest/EarBarTest_02-20170217-UnitRec/KwikFiles',
#         'EarBarTest/EarBarTest_03-20170217-UnitRec/KwikFiles',
#         'EarBarTest/EarBarTest_04-20170218-UnitRec/KwikFiles'
# ]

Log = dict(
    Done = [],
    Errors = [],
    ErrorsLog = [],
    Skipped = []
)

def ClusterizeSpks(Board, Folders, ProbeSpacing, Probe=None, Log=None):
    if not Log: Log = {K: [] for K in ['Done', 'Errors', 'ErrorsLog', 'Skipped']}
    
    for Folder in Folders:
        # try:
        if Folder in Log['Done']+Log['Errors']+Log['Skipped']: continue
        
        ExpFolder = '/'.join(Folder.split('/')[:-1]) + '/KlustaFiles'
        os.makedirs(ExpFolder, exist_ok=True)
        
        ExpName = Folder.split('/')[-2]
        Exps = glob(Folder+'/*'); Exps.sort()
        Exps = [_ for _ in Exps if '.dict' not in _
                                and '.hdf5' not in _
                                and 'BigOne' not in _] # Temp override
        ExpGroups = []
        
        while Exps:
            print('-1: Run each separately')
            print('0: All at once')
            for E, Exp in enumerate(Exps): print(str(E+1)+':', Exp)
            print(str(len(Exps)+1)+': Discard remaining experiments')
            print('\n', 'You can run all folders separately,', 
                  '\n', 'run all as one single experiment or',
                  '\n', 'group folders you want by index, comma separated.')
            
            FGroup = input('Choose folders to group: ')
            FGroup = [int(_) for _ in FGroup.split(',')]
            if len(FGroup) == 1:
                if FGroup[0] == -1: 
                    for _ in Exps: ExpGroups.append([_])
                    Exps = []; break
                
                elif FGroup[0] == 0:
                    ExpGroups.append(Exps)
                    Exps = []; break
                
                elif FGroup[0] == len(Exps)+1:
                    Exps = []; break
                
                else: 
                    FGroup = [Exps[FGroup[0]-1]]
                    ExpGroups.append(FGroup); Exps.remove(FGroup[0])
            
            else: 
                FGroup = [Exps[F-1] for F in FGroup]
                ExpGroups.append(FGroup)
                for F in FGroup: Exps.remove(F)
        
        # if not ExpGroups: print('No exps to run. Skipping folder...')
        
        for Es, Exps in enumerate(ExpGroups):
            for E, Exp in enumerate(Exps):
                Data, Rate = IO.DataLoader(Exp, 'Bits')
                if len(Data.keys()) == 1: Proc = list(Data.keys())[0]
                else: Proc = Hdf5.GetProc(Data, Board)
                
                # Keys = list(Data[Proc]['data'].keys())
                # if Data[Proc]['data'][Keys[0]].shape[1] < 16: 
                #     Log['Skipped'].append(Exp); continue
                
                for R, Rec in Data[Proc].items():
                    DataFile = ''.join([ExpName, '_Exp', "{0:02d}".format(int(Es)), 
                                       '-', "{0:02d}".format(int(E)), '_Rec', 
                                       "{0:02d}".format(int(R))])
                    RawDataInfo = {'Rate': Rate[Proc]}
                    Bin.Write(DataFile+'.dat', ExpFolder, Rec, RawDataInfo)
            
            FilesPrefix = ExpFolder+'/'+ExpName+'_Exp' + "{0:02d}".format(int(Es))
            RawDataInfo = Txt.DictRead(ExpFolder+'/'+DataFile+'-Info.dict')
            raw_data_files = glob(os.getcwd()+'/'+ FilesPrefix + '*.dat')
            raw_data_files.sort()
            
            DataInfo = glob('/'.join(Folder.split('/')[:-1]) + '/*.dict')[0]
            DataInfo = Txt.DictRead(DataInfo)
            if 'Probe' in DataInfo.keys():
                if DataInfo['Probe']['Remapped']:
                    PrbFile = os.environ['SCRIPTSPATH'] + \
                              '/Python3/Klusta/A16-'+str(ProbeSpacing)+'.prb'
                else:
                    PrbFile = FilesPrefix+'.prb'
                    Map = DataAnalysis.RemapCh(DataInfo['Probe']['Probe'], 
                                               DataInfo['Probe']['Adaptor'])
                    Klusta.PrbWrite(PrbFile, Map, ProbeSpacing)
            else:
                print('No probe info saved. Skip remap...')
                PrbFile = os.environ['SCRIPTSPATH'] + '/Python3/Klusta/A16-'+str(ProbeSpacing)+'.prb'
            
            
            Klusta.PrmWrite(FilesPrefix+'.prm', ExpName+'_Exp' + "{0:02d}".format(int(Es)), 
                            PrbFile, raw_data_files, RawDataInfo['Rate'],
                            RawDataInfo['Shape'][1], RawDataInfo['DType'])
            
            Klusta.Run(ExpName+'_Exp' + "{0:02d}".format(int(Es))+'.prm', 
                       os.getcwd()+'/'+ExpFolder, Overwrite=True)
        
        Log['Done'].append(Folder)
        
        # except Exception as e:
        #     Log['Errors'].append(Folder)
        #     Log['ErrorsLog'].append(e)
        #     print(e)


#%% After detection and clustering
# change DVCoord to real DVCoord
def FiringRateCalc(SpkClusters, SpkIds, SpkRecs, SpkSamples, Rate, RecLen, Rec, Offset=0):
    FR = ['0' for _ in range(len(SpkIds))]
    
    for I, Id in enumerate(SpkIds):
        SpksId = np.where((SpkClusters == Id) & (SpkRecs == Rec))[0]
        Spks = (SpkSamples[SpksId]-Offset)/Rate
        FR[I] = np.array([len(Spks[(Spks >= Sec) * (Spks < Sec+1)]) for Sec in range(RecLen)], dtype=int)
    
    return(FR)


def HistCalc(Spks, TTLs, Rate, HistX, Offset=0):
    Hist = np.array([])
    for TTL in TTLs:
        Firing = ((Spks-Offset)*1000/Rate) - (TTL*1000/Rate)
        Firing = Firing[(Firing >= HistX[0]) * (Firing < HistX[-1])]
        
        Hist = np.concatenate((Hist, Firing)); del(Firing)
    
    HistY, HistX = np.histogram(Hist, HistX)
    return(HistX, HistY)


def LineHistCalc(Spks, TTLs, Rate, HistX, Offset=0, Output='all'):
    """ Output: 'all', 'mean', 'norm'"""
    
    Hist = np.zeros((len(HistX)-1 ,len(TTLs)), dtype='int32')
    for T, TTL in enumerate(TTLs):
        Firing = ((Spks-Offset)*1000/Rate) - (TTL*1000/Rate)
        Firing = Firing[(Firing >= HistX[0]) * (Firing < HistX[-1])]
        
        Hist[:,T] = np.histogram(Firing, HistX)[0]; del(Firing)
    
    BinSize = HistX[1] - HistX[0]
    if Output.lower() == 'all': return(Hist)
    elif Output.lower() == 'mean': return(Hist.mean(axis=1))
    elif Output.lower() == 'norm': return(Hist.sum(axis=1)/(len(TTLs)*BinSize))
    else: print('Output should be "all", "mean" or "norm".'); return(None)


def RasterCalc(Spks, TTLs, Rate, RasterX, Offset=0):
    Raster = np.zeros((len(RasterX)-1, len(TTLs)), dtype='float32')
    
    for T, TTL in enumerate(TTLs):
        Firing = ((Spks-Offset)*1000/Rate) - (TTL*1000/Rate)
        Firing = Firing[(Firing >= RasterX[0]) * (Firing < RasterX[-1])]
        
        Raster[:,T] = np.histogram(Firing, RasterX)[0]; del(Firing)
    
    Raster[Raster == 0] = np.nan
    
    return(Raster)


def UnitResp(HistX, HistY, TimeWindow=20):
    HistX = HistX[:-1]
    RI = HistY[(HistX > 0) * (HistX <= TimeWindow)].mean() / \
         HistY[(HistX >=-TimeWindow) * (HistX < 0)].mean()
    
    if RI != RI: RI = 0
    RI = (RI - 1) * 100
    
    return(RI)


## Plots
def BarAutolabel(Ax, Bars, Color='k', Position='bottom'):
    """
    Modified from http://matplotlib.org/examples/api/barchart_demo.html
    
    Attach a text label above each bar displaying its height.
    """
    # if Position == 'bottom': Space = 1.05
    # elif Position == 'top': Space = 0.95
    
    for Bar in Bars:
        Height = Bar.get_height()
        Ax.text(Bar.get_x() + Bar.get_width()/2.,
                Height,# * Space,
                str(Height),
                color=Color,
                ha='center',
                va=Position)
    
    return(Ax)


def AllFiringRate(FR, Ax=None, Return=False, Save=True, Show=True):
    if not Ax: Ax = plt.axes()
    
    for F, Fr in enumerate(FR): Ax.plot(np.array(Fr), lw=2)
    
    if Return: return(Ax)
    
    if Save: pass
    
    if Show: plt.show()
    else: plt.close()
    return(None)


def LinePSTH(Hist, HistX, StD=False, Ax=None, Return=False, Save=True, Show=True):
    if not Ax: Ax = plt.axes()
    
    HistX = HistX[:-1]
    Mean = Hist.mean(axis=1)
    SEM = Hist.std(axis=1) / (Hist.shape[0]**0.5)
    Mins = Mean-SEM; Mins[Mins<0] = 0
    
    BLMean = Hist[HistX<0,:].mean()
    BLStD = Hist[HistX<0,:].std()# * 2
    BLSEM = Hist[HistX<0,:].std() / (Hist[HistX<0,:].shape[0]**0.5)
    BLMins = BLMean-BLSEM
    if BLMins < 0: BLMins = 0
    
    Ax.fill_between(HistX, Mean+SEM, Mins, color='r', alpha=0.3, label='SEM')
    # Ax.axhspan(BLMins, BLMean+BLSEM, color='k', alpha=0.2, lw=0)
    Ax.plot([HistX[0], HistX[-1]], [BLMean]*2, 'b--')
    if StD: Ax.plot([HistX[0], HistX[-1]], [BLStD]*2)
    Ax.plot(HistX, Mean, 'k', lw=2)
    
    if Return: return(Ax)
    
    if Save: pass
    
    if Show: plt.show()
    else: plt.close()
    return(None)


def RasterPlot(Raster, RasterX, Marker='|', Ax=None, Return=False, Save=True, Show=True):
    if not Ax: Ax = plt.axes()
    
    for R in range(Raster.shape[1]): Ax.scatter(RasterX[:-1], Raster[:,R]*(R+1), c='k', marker=Marker, s=5)
    
    if Return: return(Ax)
    if Save: pass
    
    if Show: plt.show()
    else: plt.close()
    return(None)


def RI_NaClVsCNO(RIArray, AxArgs={}, Ax=None, Return=False, FigFile='RI_StimType', Ext=['svg'], Save=True, Show=True):
    Fig = plt.figure(figsize=(4, 3), dpi=250)
    if not Ax: Ax = plt.axes()
    
    Ax.scatter(RIArray, range(len(RIArray)), c='k')
    
    if Return: return(Ax)
    
    if Save:
        Plot.Set(Ax=Ax, AxArgs=AxArgs)
        FigTitle = FigFile.split('/')[-1]
        Plot.Set(Fig=Fig, FigTitle=FigTitle, HideControls=True)
        FigPath = '/'.join(FigFile.split('/')[:-1])
        os.makedirs(FigPath, exist_ok=True)
        for E in Ext: Fig.savefig(FigFile+'.'+E, format=E, dpi=300)
    
    if Show: plt.show()
    else: plt.close()
    return(None)


def RI_StimType(RIArray, StimTypeArray, AxArgs={}, Ax=None, Return=False, FigFile='RI_StimType', Ext=['svg'], Save=True, Show=True):
    Fig = plt.figure(figsize=(4, 3), dpi=250)
    if not Ax: Ax = plt.axes()
    
    Colors = ['r', 'k']
    Stims = np.unique(StimTypeArray)
    
    for S, Stim in enumerate(Stims):
        Ax.scatter(RIArray[StimTypeArray == Stim], 
                   range(len(RIArray[StimTypeArray == Stim])), 
                   c=Colors[S])
    
    if Return: return(Ax)
    
    if Save:
        Plot.Set(Ax=Ax, AxArgs=AxArgs)
        FigTitle = FigFile.split('/')[-1]
        Plot.Set(Fig=Fig, FigTitle=FigTitle, HideControls=True)
        FigPath = '/'.join(FigFile.split('/')[:-1])
        os.makedirs(FigPath, exist_ok=True)
        for E in Ext: Fig.savefig(FigFile+'.'+E, format=E, dpi=300)
    
    if Show: plt.show()
    else: plt.close()
    return(None)


def RI_StimType_Freq(RIArray, StimTypeArray, FreqArray, AxArgs={}, Ax=None, Return=False, FigFile='RI_StimType_Freq', Ext=['svg'], Save=True, Show=True):
    Fig = plt.figure(figsize=(5, 5), dpi=250)
    if not Ax: Ax = Axes3D(Fig)
    
    Colors = ['r', 'k']
    Stims = np.unique(StimTypeArray)
    
    for S, Stim in enumerate(Stims):
        Ind = StimTypeArray == Stim
        Ax.scatter(range(len(RIArray[Ind])), 
                   RIArray[Ind], 
                   FreqArray[Ind], c=Colors[S])
    
    Ax.set_xlabel('Units')
    Ax.set_ylabel('FR Change')
    Ax.set_zlabel('MaxFreq')
    if Return: return(Ax)
    
    if Save:
        Plot.Set(Ax=Ax, AxArgs=AxArgs)
        Ax.autoscale(enable=True, axis='both', tight=True)
        # Ax.autoscale_view(tight=True, scalex=True, scaley=True, scalez=True)
        FigTitle = FigFile.split('/')[-1]
        
        Fig.suptitle(FigTitle)
        Fig.canvas.toolbar.pack_forget()
        for Obj in Fig.findobj(): Obj.set_clip_on(False)
        Fig.patch.set_visible(False)
        
        FigPath = '/'.join(FigFile.split('/')[:-1])
        os.makedirs(FigPath, exist_ok=True)
        for E in Ext: Fig.savefig(FigFile+'.'+E, format=E, dpi=300)
    
    if Show: plt.show()
    else: plt.close()
    return(None)


def WF_MeanAllChs(Waveforms, SpkX, RMSs, Ax=None, Return=False, Save=True, Show=True):
    if not Ax: Ax = plt.axes()
    ChNo = Waveforms.shape[2]
    BestCh = RMSs.index(max(RMSs))
    
    for S in range(ChNo, 0, -1):
        if S == ChNo-BestCh: Color = 'r'
        else: Color = 'k'
        
        Ax.plot(SpkX, np.mean(Waveforms[:, :, -S], axis=0)+(S*max(RMSs)), Color)
    
    Step = int(round(max(RMSs)))
    Ax.set_yticks(range(Step, int(round(ChNo*max(RMSs))), Step))
    if len(Ax.get_yticks()) < ChNo:
        Diff = ChNo - len(Ax.get_yticks())
        Ax.set_yticks(range(Step, int(round((ChNo+Diff)*max(RMSs))), Step))
    Ax.set_yticklabels(range(ChNo, 0, -1))
    
    if Return: return(Ax)
    
    if Save: pass
    
    if Show: plt.show()
    else: plt.close()
    return(None)


def WF_BestCh(Waveforms, SpkX, RMSs, SpksToPlot=100, Ax=None, Return=False, Save=True, Show=True):
    if not Ax: Ax = plt.axes()
    
    BestCh = RMSs.index(max(RMSs))
    Spks = np.arange(Waveforms.shape[0]); np.random.shuffle(Spks)
    Spks = Spks[:SpksToPlot]
    
    for Spk in Spks:  Ax.plot(SpkX, Waveforms[Spk, :, BestCh], 'r')
    Ax.plot(SpkX, np.mean(Waveforms[:, :, BestCh], axis=0), 'k')
    
    if Return: return(Ax)
    
    if Save: pass
    
    if Show: plt.show()
    else: plt.close()
    return(None)


def WF_Raster_PSTH(Waveforms, Raster, RasterX, Hist, HistX, RMSs, SpksToPlot, Rate, FigFile='WF_Raster_PSTH', Ext=['svg'], Save=True, Show=True):
    Fig = plt.figure(figsize=(7, 2.7), dpi=250)
    Grid = GridSpec(1, 3, width_ratios=[1, 3, 3])
    SubGrid = GridSpecFromSubplotSpec(2, 1, subplot_spec=Grid[-1])
    Axes = [plt.subplot(G) for G in [Grid[0], Grid[1], SubGrid[0], SubGrid[1]]]
    
    SpkX = np.arange(Waveforms.shape[1])*1000/Rate
    
    Axes[0] = WF_MeanAllChs(Waveforms, SpkX, RMSs, Axes[0], Return=True)
    Axes[1] = WF_BestCh(Waveforms, SpkX, RMSs, SpksToPlot, Ax=Axes[1], Return=True)
    Axes[2] = RasterPlot(Raster, RasterX, Marker='.', Ax=Axes[2], Return=True)
    Axes[3] = LinePSTH(Hist, HistX, StD=False, Ax=Axes[3], Return=True)
    
    YLabel = ['Channels', 'Voltage [µv]', 'Trials', 'Mean number of spikes']
    for A in range(len(Axes)):
        AxArgs = {'xlabel': 'Time [ms]', 'ylabel': YLabel[A]}
        Plot.Set(Ax=Axes[A], AxArgs=AxArgs)
    
    
    Axes[0].set_xticks(np.linspace(round(SpkX[0],1), round(SpkX[-1],1), 3))
    Axes[0].autoscale(enable=True, axis='y', tight=True)
    Axes[0].spines['left'].set_bounds(Axes[0].get_yticks()[0], Axes[0].get_yticks()[-1])
    
    Axes[1].locator_params(axis='x', nbins=4)
    
    Axes[2].spines['bottom'].set_visible(False)
    Axes[2].get_xaxis().set_visible(False)
    Axes[2].locator_params(axis='y', nbins=4)
    
    Axes[3].set_yticks(np.linspace(Axes[3].get_yticks()[0], Axes[3].get_yticks()[-1], 4))
    
    FigTitle = FigFile.split('/')[-1]
    Plot.Set(Fig=Fig, FigTitle=FigTitle, HideControls=True)
    
    if Save:
        FigPath = '/'.join(FigFile.split('/')[:-1])
        os.makedirs(FigPath, exist_ok=True)
        FigName = FigFile.split(' ')[0]
        for E in Ext: Fig.savefig(FigName+'.'+E, format=E, dpi=300)
    
    print(FigTitle, '|', Hist.sum(), 'Spks in PSTH, max of', Hist.sum(axis=1).max())
    if Show: plt.show()
    else: plt.close()
    
    return(None)


def Raster_Full(SpkClusters, UnitsId, SpkRecs, SpkSamples, Rate, Rec, Offset=0, PulseDur=None, TTLs=[], Ax=None, Return=False, FigFile='RasterFull', Ext=['svg'], Save=True, Show=True):
    if not Ax: Ax = plt.axes()
    
    if PulseDur and TTLs.size:
        for TTL in TTLs: 
            Ax.axvspan(TTL/Rate, (TTL/Rate)+PulseDur, 
                       color='k', alpha=0.3, lw=0)
    
    for I, Id in enumerate(UnitsId):
        Ids = np.where((SpkClusters == Id) & (SpkRecs == Rec))[0]
        Spks = (SpkSamples[Ids]-Offset)/Rate
        Ax.plot(Spks, np.ones(len(Spks))+I, 'ko')
    
    if Save:
        FigTitle = FigFile.split('/')[-1]
        Plot.Set(Fig=plt.figure(), FigTitle=FigTitle, HideControls=True)
        FigPath = '/'.join(FigFile.split('/')[:-1])
        os.makedirs(FigPath, exist_ok=True)
        FigName = FigFile.split(' ')[0]
        for E in Ext: plt.savefig(FigName+'.'+E, format=E, dpi=300)
    
    if Show: plt.show()
    else: plt.close()
    
    return(None)


def BoxPlots(Data, Names, LinesAmpF=2, AxArgs={}, Ax=None, Return=None, FigFile='Boxplot', Ext=['svg'], Save=True, Show=True):
    Fig = plt.figure(figsize=(4, 3), dpi=250)
    if not Ax: Ax = plt.axes()
    if type(Data) is not type(np.array([])): Data = np.array(Data).T
    
    BoxPlot = Ax.boxplot(Data, showmeans=True)
    
    for K in ['boxes', 'whiskers', 'caps', 'medians', 'fliers']:
        for I in range(Data.shape[1]): 
            BoxPlot[K][I].set(color='k')
    
    LineY = np.amax(Data) + ((np.amax(Data)-np.amin(Data))*0.1)
    
    TPairs = list(combinations(range(Data.shape[1]), 2))
    
    for TP in TPairs:
        SS = Data[:, TP[0]]; DD = Data[:, TP[1]]
        if np.mean(SS) > np.mean(DD): SS, DD = DD, SS
        
        p = Stats.RTTest(SS, DD, Alt='two.sided')
        p = round(p['p.value']*len(TPairs), 4)
        
        if p < 0.05: 
            LineY = LineY+(TP[1]*LinesAmpF)
            Plot.SignificanceBar([TP[0]+1, TP[1]+1], [LineY, LineY], str(p), Ax)
    
    Ax.set_xticks(list(range(1,Data.shape[1]+1))); Ax.set_xticklabels(Names)
    
    if Return: return(Ax)
    
    if Save:
        Plot.Set(Ax=Ax, AxArgs=AxArgs)
        FigTitle = FigFile.split('/')[-1]
        Plot.Set(Fig=Fig, FigTitle=FigTitle, HideControls=True)
        FigPath = '/'.join(FigFile.split('/')[:-1])
        os.makedirs(FigPath, exist_ok=True)
        for E in Ext: Fig.savefig(FigFile+'.'+E, format=E, dpi=300)
    
    if Show: plt.show()
    else: plt.close()
    return(None)


def PercBars(CellChanges, FigFile='./PercBars', Ext=['svg'], Save=True, Show=True):
    PercDec = CellChanges['Sound_CNOVsSound_NaCl-PercDec']
    PercDecdB = CellChanges['Sound_CNOVsSound_NaCl-PercDecdB']
    PercDecFreq = CellChanges['Sound_CNOVsSound_NaCl-PercDecFreq']
    PercDecInc = CellChanges['Sound_CNOVsSound_NaCl-PercDecInc']
    PercDecDec = CellChanges['Sound_CNOVsSound_NaCl-PercDecDec']
    
    PercInc = CellChanges['Sound_CNOVsSound_NaCl-PercInc']
    PercIncdB = CellChanges['Sound_CNOVsSound_NaCl-PercIncdB']
    PercIncFreq = CellChanges['Sound_CNOVsSound_NaCl-PercIncFreq']
    PercIncInc = CellChanges['Sound_CNOVsSound_NaCl-PercIncInc']
    PercIncDec = CellChanges['Sound_CNOVsSound_NaCl-PercIncDec']
    
    Fig, Ax = plt.subplots(1, 1, figsize=(7, 3), dpi=250)#, gridspec_kw={'width_ratios':[1, 3]})
    
    Bars = [[], [], []]
    Bars[0] = Ax.bar([0.6, 3.6, 7.6], [PercDec, PercDecInc, PercIncInc], color='gray', label='Decreased')
    Bars[1] = Ax.bar([0.6, 3.6, 7.6], [PercInc, PercDecDec, PercIncDec], bottom=[PercDec, PercDecInc, PercIncInc], label='Increased')
    
    Bars[2] = Ax.bar([1.6, 2.6, 5.6, 6.6], [PercDecdB, PercDecFreq, PercIncdB, PercIncFreq])
    
    for B, Bar in enumerate(Bars): Ax = BarAutolabel(Ax, Bar, 'w', 'top')
    
    # Ax.xticks(np.arange(8)+0.6, ('G1', 'G2', 'G3', 'G4', 'G5'))
    Ax.legend()
    Plot.Set(Ax=Ax, AxArgs={})
    FigTitle = FigFile.split('/')[-1]
    Plot.Set(Fig=Fig, FigTitle=FigTitle, HideControls=True)
    
    if Save:
        FigPath = '/'.join(FigFile.split('/')[:-1])
        os.makedirs(FigPath, exist_ok=True)
        FigName = FigFile.split(' ')[0]
        for E in Ext: Fig.savefig(FigName+'.'+E, format=E, dpi=300)
    
    if Show: plt.show()
    else: plt.close()
    
    return(None)


def FreqAnddBCurve(AllCells, FigPath='./FreqAnddBCurve', Ext=['svg'], Save=True, Show=True):
    for m, M in enumerate(AllCells):
        IntFreq = np.zeros(M['Freqs'].shape, dtype=np.int16)
        for F, Freq in enumerate(M['Freqs']):
            IntFreq[F] = sum([float(_) for _ in Freq[F].split('-')])/2
        
        for U in range(len(M['UnitId'])):
            plt.rc('font',size=8)
            Fig, Axes = plt.subplots(2,1)
            
            for d in range(M['dBCurve'].shape[0]):
                F = [3,4,0,1,2]
                Axes[0].plot(M['dBCurve'][d,F,U], label=M['Intensities'][d], lw=2)
            
            for f in [3,4,0,1,2]:    
                Axes[1].plot(M['dBCurve'][:,f,U], label=IntFreq[f], lw=2)
            
            AxArgs = {'xticks': range(len(IntFreq)),
                      'xticklabels': IntFreq[F]}
            Plot.Set(Ax=Axes[0], AxArgs=AxArgs)
            
            AxArgs = {'xticks': range(len(M['Intensities'])),
                      'xticklabels': M['Intensities']}
            Plot.Set(Ax=Axes[1], AxArgs=AxArgs)
            
            for Ax in Axes:
                # Box = Ax.get_position()
                # Ax.set_position([Box.x0, Box.y0, Box.width * 0.9, Box.height])
                Ax.legend(loc='center left', bbox_to_anchor=(1.05, 0.5),
                          prop={'size':6})
            
            FigTitle = ''.join(['Unit',
                                "{0:02d}".format(m), 
                                'x', "{0:04d}".format(M['UnitId'][U]),
                                '-', M['StimType'][U]])
            
            FigFile = '/'.join([FigPath, FigTitle])
            print(FigTitle)
            
            if Save:
                Plot.Set(Fig=Fig, FigTitle=FigTitle, HideControls=True, Tight=False)
                Fig.subplots_adjust(hspace=0.5, right=0.8)
                os.makedirs(FigPath, exist_ok=True)
                for E in Ext: Fig.savefig(FigFile+'.'+E, format=E, dpi=300)
            
            if Show: plt.show()
            else: plt.close()
    
    return(None)


def ScatterMean(Data, Names, LinesAmpF=2, Spread=0.2, LogY=False, FigFile='./ScatterDecInc', Ext=['svg'], Save=True, Show=True):
    # if type(Data) is not np.ndarray: Data = np.array(Data).T
    
    Fig, Ax = plt.subplots(1,1)
    
    if type(Data) == np.ndarray: BoxNo = Data.shape[1]
    else: BoxNo = len(Data)
    
    for P in range(BoxNo):
        if type(Data) == np.ndarray: Box = Data[:,P]
        else: Box = Data[P]
                          
        X = np.random.uniform(P+1-Spread, P+1+Spread, len(Box))
        Error = [0, np.std(Box)/len(Box)**0.5, 0]
        
        if LogY: Ax.semilogy(X, Box, 'ko')
        else: Ax.plot(X, Box, 'ko')
        
        Ax.errorbar([P+1-Spread, P+1, P+1+Spread], [np.mean(Box)]*3, Error, lw=3, elinewidth=1, capsize=10, color='k')
    
    if type(Data) == np.ndarray: 
        Margin = ((np.amax(Data)-np.amin(Data))*0.1)
        LineY = np.amax(Data) + Margin
    else:
        Margin = (max([max(m) for m in Data]) - min([min(m) for m in Data]))*0.1
        LineY = max([max(m) for m in Data]) + Margin
    
    TPairs = list(combinations(range(BoxNo), 2))
    
    for TP in TPairs:
        if type(Data) == np.ndarray: SS, DD = Data[:,TP[0]], Data[:,TP[1]]
        else: SS, DD = Data[TP[0]], Data[TP[1]]
        
        if np.mean(SS) > np.mean(DD): SS, DD = DD, SS
        if len(SS) != len(DD):
            Diff = abs(len(SS) - len(DD))
            Diff = np.full([Diff], np.nan)
            
            if len(SS) > len(DD): DD = np.hstack((DD, Diff))
            else: SS = np.hstack((SS, Diff))
            
            
        p = Stats.RTTest(SS, DD, Alt='two.sided')
        p = round(p['p.value']*len(TPairs), 4)
        
        if p < 0.05: 
            LineY = LineY+(TP[1]*LinesAmpF)
            Plot.SignificanceBar([TP[0]+1, TP[1]+1], [LineY, LineY], str(p), Ax)
        
        # if p < 0.05: 
        #     if p < 0.001: p = 'p < 0.001'
        #     else: p = 'p = ' + str(round(p, 3))
            
        #     LineY = LineY+(TP[1]*LinesAmpF)
        #     Plot.SignificanceBar([TP[0]+1, TP[1]+1], [LineY, LineY], p, Ax)
    
    AxArgs = {'ylim': [-Margin, LineY + LineY*0.05], 
              'xlim': [0, BoxNo+1],
              'ylabel': 'Firing rate change'}
    Plot.Set(Ax=Ax, AxArgs=AxArgs)
    Ax.set_xticks(list(range(1,BoxNo+1))); Ax.set_xticklabels(Names)
    
    FigTitle = FigFile.split('/')[-1]
    Plot.Set(Fig=Fig, FigTitle=FigTitle, HideControls=True)
    
    if Save:
        FigPath = '/'.join(FigFile.split('/')[:-1])
        os.makedirs(FigPath, exist_ok=True)
        FigName = FigFile.split(' ')[0]
        for E in Ext: Fig.savefig(FigName+'.'+E, format=E, dpi=300)
    
    if Show: plt.show()
    else: plt.close()
    
    return(None)


## Level 1-2
def GetAllUnits(Folder, TTLCh, ProbeChSpacing, HistX, Clusters, AnalogTTLs=True, AnalysisPath='.', SpksToPlot=100, Ext=['svg'], Save=True, Show=False):
    UnitRec = {}
    DataFolder = '/'.join(Folder.split('/')[:-1]) + '/KwikFiles'
    InfoFile = glob('/'.join(Folder.split('/')[:-1]) + '/*.dict')[0]
    FigPath = Folder.split('/')[0] + '/Figs/' + Folder.split('/')[1]
    
    DataInfo = Txt.DictRead(InfoFile)
    
    KlustaRecs = []; ExpFolders = sorted(glob(DataFolder+'/*'))
    for f, F in enumerate(ExpFolders):
        # Data = IO.DataLoader(F, 'Bits')[0]
        # Proc = list(Data.keys())[0]
        # Recs = sorted(Data[Proc].keys())
        ExpInfo = DataInfo['ExpInfo']["{0:02d}".format(f)]
        Recs = OpenEphys.GetRecs(F)
        Recs = Recs[list(Recs.keys())[0]]
        
        Freq = ExpInfo['Hz']
        DV = ExpInfo['DVCoord']
        StimType = '_'.join(ExpInfo['StimType'])
        
        for R in Recs: KlustaRecs.append([F, Freq+'Hz', str(DataInfo['Intensities'][int(R)])+'dBSPL', 
                                          StimType, DV, R])
        # del(Data)
    
    Rate = np.array(Clusters.sample_rate)
    Offsets = Clusters.all_traces.offsets
    RecLen = 50
    
    Good = Clusters.cluster_groups
    Good = [Id for Id,Key in Good.items() if Key == 'good']
    
    ChNo = len(Clusters.channels)
    # SpkLen = Clusters.n_samples_waveforms
    RasterX = np.arange(HistX[0], HistX[-1], 1000/Rate)
    
    UnitRec = {}
    for K in ['UnitId', 'DV']: UnitRec[K] = np.zeros((len(Good)*len(Clusters.recordings)), dtype=np.int16)
    for K in ['dB', 'RI']: UnitRec[K] = np.zeros((len(Good)*len(Clusters.recordings)), dtype=np.float32)
    for K in ['Freq', 'StimType']: UnitRec[K] = ['0' for _ in range(len(Good)*len(Clusters.recordings))]
    
    UnitRec['RMSs'] = np.zeros((ChNo, len(Good)*len(Clusters.recordings)), dtype=np.float32)
    
    UnitRec['PSTHX'] = np.zeros((len(HistX),len(Good)*len(Clusters.recordings)), dtype=np.float32)
    UnitRec['PSTH'] = np.zeros((1,1,1))
    
    UnitRec['FiringRate'] = np.zeros((RecLen,len(Good)*len(Clusters.recordings)), dtype=np.float32)
    UnitRec['Spks'] = [0 for _ in range(len(Good)*len(Clusters.recordings))]
    
    for rec, Rec in enumerate(Clusters.recordings):
        print('')
        print('Running folder', Folder, 'Rec', rec, 'of', len(Clusters.recordings))
        RecFolder = KlustaRecs[Rec][0]
        Freq = KlustaRecs[Rec][1]
        dB = KlustaRecs[Rec][2]
        StimType = KlustaRecs[Rec][3]
        R = KlustaRecs[Rec][5]
        
        if AnalogTTLs:
            TTLs = IO.DataLoader(RecFolder, 'Bits')[0]
            Proc = list(TTLs.keys())[0]
            TTLs = TTLs[Proc][KlustaRecs[Rec][-1]][:, TTLCh-1]
            TTLs = DataAnalysis.QuantifyTTLsPerRec(AnalogTTLs, TTLs)
            
        else:
            TTLs = DataAnalysis.QuantifyTTLsPerRec(AnalogTTLs, EventsDict='Events', 
                                                   TTLCh=TTLCh, Proc=Proc, Rec=R)
        
        if not UnitRec['PSTH'].max():
            UnitRec['PSTH'] = np.zeros((len(HistX)-1, 
                                         len(TTLs), 
                                         len(Good)*len(Clusters.recordings)), 
                                        dtype=np.float32)
        
        for I, Id in enumerate(Good):
            Ind = I + (len(Good)*rec)
            DV = KlustaRecs[Rec][4]
            SpksId = (Clusters.spike_clusters == Id) & \
                     (Clusters.spike_recordings == Rec)
            
            Waveforms = Clusters.all_waveforms[SpksId]
            ChNo = Waveforms.shape[2]
            
            RMSs = [(np.nanmean((np.nanmean(Waveforms[:, :, Ch], axis=0))**2))**0.5
                    for Ch in range(ChNo)]
            
            if RMSs[0] != RMSs[0]: print('Spk contains NaNs. Skipping...'); continue
            
            BestCh = RMSs.index(max(RMSs))
            DV = str(int(DV) - ((ChNo-BestCh) * ProbeChSpacing))
            
            Hist = LineHistCalc(Clusters.spike_samples[SpksId], TTLs, Rate, HistX, Offsets[Rec])
            Raster = RasterCalc(Clusters.spike_samples[SpksId], TTLs, Rate, RasterX, Offsets[Rec])
            
            if Hist.max() == 0: print('No Spks in PSTH. Skipping...'); continue
            
            RI = UnitResp(HistX, Hist)
            
            UnitRec['UnitId'][Ind] = Id
            UnitRec['StimType'][Ind] = StimType
            UnitRec['Freq'][Ind] = Freq[:-2]
            UnitRec['dB'][Ind] = float(dB[:-5])
            UnitRec['DV'][Ind] = int(DV)
            UnitRec['RI'][Ind] = RI
            UnitRec['Spks'][Ind] = Waveforms
            
            UnitRec['RMSs'][:, Ind] = RMSs
            UnitRec['PSTHX'][:,Ind] = HistX
            UnitRec['PSTH'][:,:len(TTLs),Ind] = Hist
            
            UnitInfo = '_'.join([
                           'Unit'+"{0:04d}".format(Id), 
                           StimType, 
                           str(DV)+'DV', 
                           Freq, 
                           dB
                        ]) + ' RI:' + str(round(RI,4))
            
            FigFile = FigPath+'/'+FigPath.split('/')[-1]+'-'+UnitInfo
            WF_Raster_PSTH(Waveforms, Raster, RasterX, Hist, HistX, RMSs, SpksToPlot, Rate, FigFile, Ext, Save, Show)
        
        
        FR = FiringRateCalc(Clusters.spike_clusters, Good, Clusters.spike_recordings, Clusters.spike_samples, Rate, RecLen, Rec, Offsets[Rec])
        for F, Fr in enumerate(FR):
            Ind = F + (len(Good)*rec)
            UnitRec['FiringRate'][:,Ind] = Fr
        
        del(FR)
    
    # Cleaning
    for K in ['Freq', 'StimType']: UnitRec[K] = np.array(UnitRec[K])
    
    Keys = ['UnitId', 'DV', 'dB', 'RI', 'Freq', 'StimType']
    ValidIds = UnitRec['UnitId'] != 0
    
    for Key in Keys: UnitRec[Key] = UnitRec[Key][ValidIds]
    for Key in ['RMSs', 'PSTHX', 'FiringRate']: UnitRec[Key] = UnitRec[Key][:,ValidIds]
    UnitRec['PSTH'] = UnitRec['PSTH'][:, :, ValidIds]
    UnitRec['Spks'] = [UnitSpks for S, UnitSpks in enumerate(UnitRec['Spks']) 
                                if ValidIds[S]]
    
    # Hdf5.DataWrite(UnitRec, '/'+Folder.split('/')[1]+'/Units', AnalysisFile, Overwrite=True)
    Asdf.Write(UnitRec, AnalysisPath, AnalysisPath + '_' +Folder.split('/')[0]+'_AllUnits.asdf')
    return(UnitRec)


def GetCellChanges(Cells, Stims, StimRef):
    CellChanges = {}
    
    for Stim in Stims:
        Ind = (Cells['StimType'] == Stim)
        CellChanges[Stim] = Cells['MaxRI'][Ind]
        CellChanges[Stim+'-Freqs'] = Cells['MaxFreq'][Ind]
        CellChanges[Stim+'-dBs'] = Cells['MaxdB'][Ind]
        CellChanges[Stim+'-Dec'] = CellChanges[Stim] < 0
    
    StimPairs = list(combinations(Stims,2))
    
    for Pair in StimPairs:
        if StimRef not in Pair: continue
        if Pair[-1] != StimRef: Pair = (Pair[1], Pair[0])
        
        Key = 'Vs'.join(Pair)
        
        CellChanges[Key] = ((CellChanges[Pair[0]]/100)+1)/((CellChanges[Pair[1]]/100)+1)
        CellChanges[Key] =  (CellChanges[Key]-1)*100
        
        Dec = CellChanges[Key] < 0
        
        CellChanges[Key+'-PercInc'] = round(len(~Dec[~Dec == True])*100/len(CellChanges[Key]), 2)
        CellChanges[Key+'-PercDec'] = round(len(Dec[Dec == True])*100/len(CellChanges[Key]), 2)
        
        ChangedFreq = CellChanges[Pair[0]+'-Freqs'] != CellChanges[Pair[1]+'-Freqs']
        ChangeddB = CellChanges[Pair[0]+'-dBs'] != CellChanges[Pair[1]+'-dBs']
        
        CellChanges[Key+'-PercIncFreq'] = round(len(CellChanges[Key][~Dec * ChangedFreq])*100/len(CellChanges[Key][~Dec]), 2)
        CellChanges[Key+'-PercIncdB'] = round(len(CellChanges[Key][~Dec * ChangeddB])*100/len(CellChanges[Key][~Dec]), 2)
        CellChanges[Key+'-PercIncInc'] = round(len(~Dec[~Dec * ~CellChanges[StimRef+'-Dec']])*100/len(~Dec[~Dec]), 2)
        CellChanges[Key+'-PercIncDec'] = round(len(~Dec[~Dec * CellChanges[StimRef+'-Dec']])*100/len(~Dec[~Dec]), 2)
        
        CellChanges[Key+'-PercDecFreq'] = round(len(CellChanges[Key][Dec * ChangedFreq])*100/len(CellChanges[Key][Dec]), 2)
        CellChanges[Key+'-PercDecdB'] = round(len(CellChanges[Key][Dec * (ChangeddB)])*100/len(CellChanges[Key][Dec]), 2)
        CellChanges[Key+'-PercDecInc'] = round(len(Dec[Dec * ~CellChanges[StimRef+'-Dec']])*100/len(Dec[Dec]), 2)
        CellChanges[Key+'-PercDecDec'] = round(len(Dec[Dec * CellChanges[StimRef+'-Dec']])*100/len(Dec[Dec]), 2)
        
        CellChanges[Key+'-Dec'] = Dec
    
    return(CellChanges)
    


def GetUnitsParameters(Folder, UnitRec, UnitsId, Stims, Freqs, Intensities, AnalysisPath):
    Cells = {}
    Cells['UnitId'] = np.zeros((Stims.size * UnitsId.size), dtype=UnitRec['UnitId'].dtype)
    Cells['StimType'] = np.zeros((Stims.size * UnitsId.size), dtype=UnitRec['StimType'].dtype)
    Cells['MaxRI'] = np.zeros((Stims.size * UnitsId.size), dtype=UnitRec['RI'].dtype)
    Cells['MaxFreq'] = np.zeros((Stims.size * UnitsId.size), dtype=UnitRec['Freq'].dtype)
    Cells['MaxdB'] = np.zeros((Stims.size * UnitsId.size), dtype=UnitRec['dB'].dtype)
    Cells['dBCurve'] = np.zeros((Intensities.size, Freqs.size, Stims.size * UnitsId.size), dtype=UnitRec['dB'].dtype)
    Cells['HzCurve'] = np.zeros((Freqs.size, Intensities.size, Stims.size * UnitsId.size), dtype=UnitRec['dB'].dtype)
    Cells['Intensities'] = Intensities
    Cells['Freqs'] = Freqs
    
    InfUnits = UnitRec['RI'] == np.inf
    Invalid = []
    
    for I, Id in enumerate(UnitsId):
        UnitId = UnitRec['UnitId'] == Id
        
        for S, Stim in enumerate(Stims):
            Ind = S + (len(Stims)*I)
            
            UnitStim = UnitRec['StimType'] == Stim
            ThisUnit = UnitId * UnitStim * ~InfUnits
            
            if not True in ThisUnit:# or np.unique(UnitRec['RI'][ThisUnit])[0] == -100:
                Invalid.append(Id)
                print("No responses found for unit", Id, 'under', Stim, 'stimulation.')
                # for Key in ['MaxFreq', 'MaxdB']:
                #     Cells[Key][Ind] = np.nan
                
                # Cells['MaxRI'][Ind] = np.nan
                # Cells['dBCurve'][:, :, Ind] = np.nan
                break
            
            MaxRI = max(UnitRec['RI'][ThisUnit])
            # MaxRI = max(UnitRec['RI'][ThisUnit * (-100 < UnitRec['RI']) * (UnitRec['RI'] < RIStd)])
            # MinRI = min(UnitRec['RI'][ThisUnit * (-RIStd < UnitRec['RI']) * (UnitRec['RI'] < RIStd)])
            # if abs(MinRI) > abs(MaxRI): MaxRI = MinRI
            MaxRI = np.where(ThisUnit * (UnitRec['RI'] == MaxRI))
            # MinRI = min(UnitRec['RI'][ThisUnit])
            # MinRI = np.where(ThisUnit * (UnitRec['RI'] == MinRI))
            
            for F, Freq in enumerate(Freqs):
                for d, dB in enumerate(Intensities):
                    ThisDV = np.unique(UnitRec['DV'][
                                (UnitRec['Freq'] == Freq) * 
                                (UnitRec['dB'] == dB) * 
                                ThisUnit
                             ])
                    
                    if not ThisDV.size: continue
                    
                    ThisDV = ThisDV[(np.abs(ThisDV - UnitRec['DV'][MaxRI][0])).argmin()]
                    
                    FinalInd = ThisUnit * \
                               (UnitRec['Freq'] == Freq) * \
                               (UnitRec['dB'] == dB) * \
                               (UnitRec['DV'] == ThisDV)
                    
                    Cells['dBCurve'][d, F, Ind] = UnitRec['RI'][FinalInd]
            
            dBFMax = Cells['dBCurve'][:, :, Ind].max()
            # dBFMin = Cells['dBCurve'][:, :, Ind].min()
            # if abs(dBFMin) > abs(dBFMax): dBFMax = dBFMin
            d, F = np.where(Cells['dBCurve'][:, :, Ind] == dBFMax)
            
            Cells['UnitId'][Ind] = Id
            Cells['StimType'][Ind] = Stim
            Cells['MaxRI'][Ind] = UnitRec['RI'][MaxRI][0]
            Cells['MaxFreq'][Ind] = Freqs[F][0]
            Cells['MaxdB'][Ind] = Intensities[d][0]
            # Cells['MaxFreq'][Ind] = UnitRec['Freq'][MaxRI][0]
            # Cells['MaxdB'][Ind] = UnitRec['dB'][MaxRI][0]
    
    
    ## Cleaning
    Keys = ['UnitId', 'StimType', 'MaxRI', 'MaxFreq', 'MaxdB']
    
    if Invalid: 
        Invalid = [Cells['UnitId'] == I for I in Invalid]
        Invalid = np.sum(Invalid, axis=0, dtype=bool)
        ValidIds = (Cells['UnitId'] != 0) * ~Invalid
    else:
        ValidIds = (Cells['UnitId'] != 0)
    
    for Key in Keys: Cells[Key] = Cells[Key][ValidIds]
    Cells['dBCurve'] = Cells['dBCurve'][:, :, ValidIds]
    
    Asdf.Write(Cells, AnalysisPath, AnalysisPath + '_' +Folder.split('/')[0]+'_Cells.asdf')
    return(Cells)


def Units(Group, TimeBeforeTTL, TimeAfterTTL, BinSize, TTLCh, ProbeChSpacing, AnalogTTLs=True, SpksToPlot=100, Ext=['svg'], Save=True, Show=False):
    HistX = np.arange(-TimeBeforeTTL, TimeAfterTTL, BinSize)
    
    Folders = sorted(glob('**/KlustaFiles', recursive=True))
    Folders = [_ for _ in Folders if Group in _]
    
    for Folder in Folders:
        AnalysisPath = '/'+Folder.split('/')[1]+'/Units'
        AnalysisPath = AnalysisPath[1:].split('/')[0]
        FigPath = Folder.split('/')[0] + '/Figs/' + Folder.split('/')[1]
        KwikFile = glob(Folder+'/*.kwik')[0]
        
        Clusters = KwikModel(KwikFile)
        UnitRec = GetAllUnits(Folder, TTLCh, ProbeChSpacing, HistX, Clusters, AnalogTTLs, AnalysisPath,  SpksToPlot, Ext, Save, Show)
#        UnitRec = Asdf.Load('/', AnalysisPath+ '_' +Folder.split('/')[0]+'_AllUnits.asdf')
        UnitsId = np.unique(UnitRec['UnitId'])
        Stims = np.unique(UnitRec['StimType'])
        Freqs = np.unique(UnitRec['Freq'])
        Intensities = np.unique(UnitRec['dB'])
        Cells = GetUnitsParameters(Folder, UnitRec, UnitsId, Stims, Freqs, Intensities, AnalysisPath)
        
        IntFreq = np.zeros(Cells['MaxFreq'].shape, dtype=np.int16)
        for F, Freq in enumerate(Cells['MaxFreq']):
            IntFreq[F] = sum([float(_) for _ in Cells['MaxFreq'][F].split('-')])/2
        
        CellChanges = GetCellChanges(Cells, Stims, Stims[1])
        
        FreqsNaCl = IntFreq[Cells['StimType'] == [Stims[1]]]
        FreqsCNO = IntFreq[Cells['StimType'] == [Stims[0]]]
        ChangeInFreq = ((FreqsCNO/FreqsNaCl)-1)*100
        Dec = CellChanges['Vs'.join(Stims)+'-Dec']
        
        ## Plots
        UnitInfo = 'PercBars-'+'-'.join(Stims)
        FigFile = FigPath+'/'+FigPath.split('/')[-1]+'-'+UnitInfo
        PercBars(CellChanges, FigFile, Ext, Save, Show)
        
        DecInc = CellChanges['Sound_CNOVsSound_NaCl']
        UnitInfo = 'ScatterDecInc-'+'-'.join(Stims)
        FigFile = FigPath+'/'+FigPath.split('/')[-1]+'-'+UnitInfo
        ScatterMean([(DecInc[DecInc<0]/100)+1, (DecInc[DecInc>0]/100)+1], 
                    ('Decreased', 'Increased'), 
                    FigFile=FigFile, Ext=['svg'], Save=Save, Show=Show)
        
        UnitInfo = 'ScatterChangedFreq-'+'-'.join(Stims)
        FigFile = FigPath+'/'+FigPath.split('/')[-1]+'-'+UnitInfo
        ScatterMean([ChangeInFreq[Dec], ChangeInFreq[~Dec]], 
                    ('Decreased', 'Increased'), 
                    FigFile=FigFile, Ext=['svg'], Save=Save, Show=Show)
        
        UnitInfo = 'RI-'+'-'.join(Stims)
        FigFile = FigPath+'/'+FigPath.split('/')[-1]+'-'+UnitInfo
        RI_StimType(Cells['MaxRI'], Cells['StimType'], FigFile=FigFile, Ext=Ext, Save=Save, Show=Show)
        RI_StimType_Freq(Cells['MaxRI'], Cells['StimType'], IntFreq, FigFile=FigFile, Ext=Ext, Save=Save, Show=Show)
        
        UnitInfo = 'Boxplot-'+'-'.join(Stims)
        FigFile = FigPath+'/'+FigPath.split('/')[-1]+'-'+UnitInfo
        AxArgs = {'xlabel': '', 'ylabel': 'Evoked firing rate change'}
        BoxPlots(np.vstack((CellChanges[Stims[1]], CellChanges[Stims[0]])).T, ['NaCl', 'CNO'], AxArgs=AxArgs, FigFile=FigFile, Ext=Ext, Save=Save, Show=Show)
        
        UnitInfo = 'NaClVsCNO-'+'-'.join(Stims)
        FigFile = FigPath+'/'+FigPath.split('/')[-1]+'-'+UnitInfo
        
        RI_NaClVsCNO(CellChanges[Stims[0]+'Vs'+Stims[1]], FigFile=FigFile, Ext=Ext, Save=Save, Show=Show)
        
        # Raster_Full(Clusters.spike_clusters, Good, Clusters.spike_recordings, 
        #             Clusters.spike_samples, Rate, StimType, Freq, dB, 
        #             Offsets[Rec], DataInfo['SoundPulseDur'], TTLs, FigPath, Ext, 
        #             Save, Show)
    
        # AllFiringRate(FR, Save=Save, Show=Show)
        
    
    AllCells = [Asdf.Load('/', Folder.split('/')[-2] + '_' + Folder.split('/')[0] + '_Cells.asdf')
                for Folder in Folders]
    
    FreqAnddBCurve(AllCells, Ext=Ext, Save=Save, Show=Show)
    
    TotalU = [np.unique(C['UnitId']) for C in AllCells]
    
    for C, Cell in enumerate(AllCells):
        
        if C == 0: NewUnitsId = np.arange(len(TotalU[C]))
        else: NewUnitsId = np.arange(len(TotalU[C])) + AllCells[C-1]['UnitId'][-1]+1
        
        for I, Id in enumerate(Cell['UnitId']):
            Ind = np.where(TotalU[C] == Id)[0][0]
            AllCells[C]['UnitId'][I] = NewUnitsId[Ind]
    
    Merge = {}
    for Key in AllCells[0].keys():
        if Key not in ['Intensities', 'Freqs']:
            Merge[Key] = np.concatenate([Cell[Key] for Cell in AllCells], axis=-1)
    
    Stims = np.unique(Merge['StimType'])
    CellChangesMerge = GetCellChanges(Merge, Stims, Stims[-1])
    
    UnitInfo = 'PercBarsAll-'+'-'.join(Stims)
    FigFile = FigPath+'/'+FigPath.split('/')[-1]+'-'+UnitInfo
    PercBars(CellChangesMerge, FigFile, Ext, Save, Show)
    
    DecInc = CellChanges['Sound_CNOVsSound_NaCl']
    UnitInfo = 'ScatterDecIncAll-'+'-'.join(Stims)
    FigFile = FigPath+'/'+FigPath.split('/')[-1]+'-'+UnitInfo
    ScatterMean([(DecInc[DecInc<0]/100)+1, (DecInc[DecInc>0]/100)+1], 
                ('Decreased', 'Increased'), 
                FigFile=FigFile, Ext=['svg'], Save=Save, Show=Show)
    
    UnitInfo = 'ScatterChangedFreqAll-'+'-'.join(Stims)
    FigFile = FigPath+'/'+FigPath.split('/')[-1]+'-'+UnitInfo
    ScatterMean([ChangeInFreq[Dec], ChangeInFreq[~Dec]], 
                ('Decreased', 'Increased'), 
                FigFile=FigFile, Ext=['svg'], Save=Save, Show=Show)
    
    UnitInfo = 'BoxplotAll-'+'-'.join(Stims)
    FigFile = FigPath+'/'+FigPath.split('/')[-1]+'-'+UnitInfo
    AxArgs = {'xlabel': '', 'ylabel': 'Evoked firing rate change'}
    BoxPlots(np.vstack((CellChangesMerge[Stims[1]], CellChangesMerge[Stims[0]])).T, ['NaCl', 'CNO'], AxArgs=AxArgs, FigFile=FigFile, Ext=Ext, Save=Save, Show=Show)
    
    UnitInfo = 'NaClVsCNOAll-'+'-'.join(Stims)
    FigFile = FigPath+'/'+FigPath.split('/')[-1]+'-'+UnitInfo
    RI_NaClVsCNO(CellChangesMerge[Stims[0]+'Vs'+Stims[1]], FigFile=FigFile, Ext=Ext, Save=Save, Show=Show)
        
        
