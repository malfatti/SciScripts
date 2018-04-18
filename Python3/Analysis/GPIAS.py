#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: T. Malfatti <malfatti@disroot.org>
@year: 2016
@license: GNU GPLv3 <https://raw.githubusercontent.com/malfatti/SciScripts/master/LICENSE>
@homepage: https://github.com/Malfatti/SciScripts
"""
#%% Import
#import GPIAZon
import numpy as np
import os

from DataAnalysis import GPIAS
from DataAnalysis.Plot import GPIAS as PlotGPIAS
from glob import glob
from IO import Hdf5, OpenEphys, Txt


#%% Batch
Group = 'ToBeAssigned'
AnalysisFile = Group + '/' + Group + '-Analysis.hdf5'

GPIASTimeBeforeTTL = 200   # in ms
GPIASTimeAfterTTL = 200    # in ms
FilterFreq = [100, 300]     # frequency for filter
FilterOrder = 3       # butter order
Filter = None
# Filter = 'butter'
Stim = 'Sound'

Ext=['svg']; Save = False; Show = True

Exps = sorted(glob(Group+'/2*IAS'))#[2:]
Exps = [Exps[-1]] # Just the last folder

for Exp in Exps:
#    Exp = Group + '/20170721-Prevention-GPIAS'
    Folders = sorted(glob(Exp + '/' + Exp.split('/')[-1][:4] + '-*'))
    Files = sorted(glob(Exp + '/' + Exp.split('/')[-1][:4] + '*dict'))
    
    # NaCl
    #del(Paths[3], Paths[2], Paths[0]); del(Files[3], Files[2], Files[0])
    # SSal
    #del(Paths[-2:], Paths[5], Paths[0]); del(Files[-2:], Files[5], Files[0])
#    20170521-Prevention_B-GPIAS
#    del(Folders[:2]); del(Files[:2])
    
    StimExps = []
    for F, Folder in enumerate(Folders):
        DataInfo = Txt.DictRead(Files[F])
        
        if Stim in DataInfo['Animal']['StimType']: StimExps.append(Folder)
    
#    if not StimExps: continue
    
    for F, Folder in enumerate(StimExps):
        if Folder == 'Recovery/20160702-Recovery-GPIAS/2016-07-02_13-05-52_CaMKIIahM4Dn09': continue
        RecFolder = Folder.split('/')[-1]
        
        Data, Rate = OpenEphys.DataLoader(Folder, AnalogTTLs=True, Unit='uV')
        if len(Data.keys()) == 1: Proc = list(Data.keys())[0]
        
        if Files[F][-4:] == 'hdf5': 
            DataInfo = Hdf5.DictLoad('/DataInfo', Files[F])
            DataInfo['PiezoCh'] = [1]; DataInfo['TTLCh'] = 0
            for Path in ['Freqs', 'FreqOrder', 'FreqSlot']:
                DataInfo[Path] = Hdf5.DatasetLoad('/DataInfo/'+Path, Files[F])
            
        else: DataInfo = Txt.DictRead(Files[F])
        
        # Check channels | File = Files[F]
#        from DataAnalysis.Plot import Plot
#        Chs = [Data[Proc]['0'][:,Ch] for Ch in range(Data[Proc]['0'].shape[1])]
#        Plot.RawCh(Chs, Lines=len(Chs), Cols=1, Save=False)
        
        # Test recs
#        SOAB = GPIAS.CheckGPIASRecs(Data[Proc], [65000, 100000])
#        if SOAB: 
#            print(); print(Folder)
#            print(); print(SOAB); print()
#        else: print(); print(Folder, 'clear.'); print()
        
    #    for Rec in Data[Proc].keys():
    #        BitVolts = 10000/(2**16)
    #        Data[Proc][Rec] = Data[Proc][Rec] * BitVolts
        
    #    for Path in ['Freqs', 'FreqOrder', 'FreqSlot']:
    #        DataInfo[Path] = Hdf5F.LoadDataset('/DataInfo/'+Path, Files[F])
    #    
    #    DataInfo['FreqOrder'][-3:][:,1] = -2
    #    DataInfo['PiezoCh'] = [int(DataInfo['PiezoCh'])]
    #    DataInfo['TTLCh'] = int(DataInfo['TTLCh'])
        
    #     DataInfo['PiezoCh'] = [3]; DataInfo['TTLCh'] = 1
        ExpStim = '_'.join(DataInfo['Animal']['StimType'])
        AnalysisKey = Files[F][:-5].split('/')[-1] + '-' + ExpStim + '-' + Group + '_' + 'GPIAS'
        FigPrefix = AnalysisKey.replace('/', '_')
        FigName = '/'.join([Group, 'Figs', FigPrefix+'_Traces'])
        
        # if Save:
        #     os.makedirs('/'.join(FigName.split('/')[:-1]), exist_ok=True)
        
        GPIASRec, XValues = GPIAS.Analysis(
                             Data[Proc], DataInfo, Rate[Proc], AnalysisFile, 
                             AnalysisKey, GPIASTimeBeforeTTL, GPIASTimeAfterTTL, 
                             FilterFreq, FilterOrder, Filter, Return=True, Overwrite=True)
        
#        GPIASData = Hdf5.DataLoad(AnalysisKey, AnalysisFile)[0]
#        GPIASRec, XValues = GPIASData['GPIAS'], GPIASData['XValues']
        
        PlotGPIAS.Traces(GPIASRec, XValues, DataInfo['Audio']['SoundLoudPulseDur'], 
                         FigName, Ext, Save, Show)
        
        del(GPIASRec, XValues)


#%% Single
GPIASTimeBeforeTTL = 200   # in ms
GPIASTimeAfterTTL = 200    # in ms
FilterFreq = [100, 300]     # frequency for filter
FilterOrder = 3       # butter order
# Filter = 'butter'
Filter = None
Stim = 'Sound'

Ext=['svg']; Save = False; Show = True

Folder = '/home/cerebro/Malfatti/Data/2018-04-04_11-21-10_A1'
InfoFile = '/home/cerebro/Data/20180404112049-A1-GPIAS.dict'
AnalysisFile = 'Test.hdf5'
AnalysisKey = 'Test'
FigPrefix = 'Test'
FigName = 'Test'
RecFolder = Folder.split('/')[-1]

Data, Rate = OpenEphys.DataLoader(Folder, AnalogTTLs=True, Unit='uV')
if len(Data.keys()) == 1: Proc = list(Data.keys())[0]

DataInfo = Txt.DictRead(InfoFile)
# DataInfo['DAqs']['PiezoCh'] = [3]
ExpStim = '_'.join(DataInfo['Animal']['StimType'])

GPIASRec, XValues = GPIAS.Analysis(
                     Data[Proc], DataInfo, Rate[Proc], AnalysisFile, 
                     AnalysisKey, GPIASTimeBeforeTTL, GPIASTimeAfterTTL, 
                     FilterFreq, FilterOrder, Filter, Return=True, Overwrite=True)

PlotGPIAS.Traces(GPIASRec, XValues, DataInfo['Audio']['SoundLoudPulseDur'], 
                 FigName, Ext, Save, Show)


#%% Olds
Animals = ['CaMKIIahM4Dn08', 'CaMKIIahM4Dn09']
ExpList = ['Before', 'After1', 'After2', 'After3', 'NaCl', 'CNO']

Group = {Animal: {Exp: {} for Exp in ExpList} for Animal in Animals}



for Animal in Animals:
    AnalysisFile = Animal + '/' + Animal + '-Analysis.hdf5'
    Exps = Hdf5.GetGroupKeys('/', AnalysisFile)
    Exps = [E for E in Exps if 'GPIAS' in E]
    
    for Exp in Exps:
        pass
    

#%% Group
#AnalysisFile = 'GPIAZon/GPIAZon-Analysis.hdf5'
#Groups = ['GPIAZon_NaCl', 'GPIAZon_SSal']
#ExpList = ['NaCl', 'SSal']#, 'Atr']
#Exps = {'GPIAZon_NaCl': ['NaCl', 'SSal', 'Atr'],
#        'GPIAZon_SSal': ['SSal', 'Atr', 'NaCl']}

Group = 'Recovery'
ExpList = ['BeforeANT', 'AfterANTNaCl', 'AfterANTCNO']
Animals = ['CaMKIIahM4Dn06', 'CaMKIIahM4Dn08', 'CaMKIIahM4Dn09']
#Animals = ['Prevention_A3', 'Prevention_A4', 'Prevention_A5']

Exps = glob(Group + '/*' + Group + '*GPIAS'); Exps.sort()
Exps = [Exps[E] for E in [0,1,4,6]]# # Recovery override
#Exps = [Exps[E] for E in [1,3,4]] # Prevention override

AnalysisFile = Group + '/' + Group + '-Analysis.hdf5'
AnalysisFolder = Group + '/Analysis'
#Dicts = sorted([E for A in Animals for E in glob(Group+'/2*IAS/*dict') if A in E])

Save = False; Invalid = False
DiffThr = 60; InvalidThr = 20

IndexPerExp = GetMAF(Group, Animals, Exps, ExpList, AnalysisFolder, DiffThr, Invalid, InvalidThr)
PVals = GetPValues(Group, Animals, Exps, ExpList, AnalysisFolder, DiffThr, Invalid, InvalidThr)
Index_Exp_BP(IndexPerExp, ExpList, PVals, Invalid, Save=Save)

IndexPerExp = GPIAS.GroupData.GetMAF(Group, Animals, Exps, ExpList, AnalysisFolder, DiffThr, Invalid, InvalidThr)
GPIAS.GroupData.Index_Exp_BP(IndexPerExp, ExpList)


#%% MatFiles
from IO import Asdf, Mat

AnalysisFolder = 'Recovery/Analysis'
Folders = sorted(glob('Recovery/2*IAS/*00-00-00*'))
InfoFiles = sorted(glob('Recovery/2*IAS/*.mat'))

GPIASTimeBeforeTTL = 200   # in ms
GPIASTimeAfterTTL = 200    # in ms
FilterFreq = [100, 300]     # frequency for filter
FilterOrder = 3       # butter order
Filter = 'butter'

Mat.GPIASAnalysis(Folders, InfoFiles, AnalysisFolder, GPIASTimeBeforeTTL, 
                  GPIASTimeAfterTTL, FilterFreq, FilterOrder, Filter)

for F, Folder in enumerate(Folders):
    AnalysisKey = InfoFiles[F].split('/')[1].split('-')
    AnalysisKey[0] = AnalysisKey[0]+'000000'
    AnalysisKey[1] = Folder.split('_')[-1]
    AnalysisKey = '-'.join(AnalysisKey) + '-Sound-Recovery_' + 'GPIAS'
    
#    Data = Hdf5.DataLoad(AnalysisKey, AnalysisFile)[0]
    Data = Asdf.Load('/', AnalysisFolder+'/'+AnalysisKey+'.asdf')
    FigName = 'Recovery/Figs/'+AnalysisKey.replace('/', '_')+'_Traces'
    
    PlotGPIAS.Traces(Data['GPIAS'], Data['XValues'], 0.05, FigName, Save=True, Show=False)


#%% convert hdf5 to dict
Files = glob('Prevention/*GPIAS/2*.hdf5'); Files.sort()

from DataAnalysis.Plot import Plot
from IO import OpenEphys
    
for File in Files:
    StimInfo, DataInfo = Hdf5.DataLoad('/DataInfo', File)
    DataInfo.update(StimInfo)
    
    Folder = sorted(glob('/'.join(File.split('/')[:-1])+'/2*'+DataInfo['AnimalName']))[0]
    Data, Rate = OpenEphys.DataLoader(Folder, AnalogTTLs=True, Unit='Bits')
    if len(Data.keys()) == 1: Proc = list(Data.keys())[0]
    
    Chs = [Data[Proc]['0'][:,Ch] for Ch in range(Data[Proc]['0'].shape[1])]
    Plot.RawCh(Chs, Lines=len(Chs), Cols=1, Save=False)
    
    PiezoCh = input('PiezoCh: '); PiezoCh = [int(PiezoCh)]
    TTLCh = input('TTLCh: '); TTLCh = int(TTLCh)
    
#    DataDict = Txt.DictRead(File[:-4]+'dict')
#    PiezoCh, TTLCh = DataDict['PiezoCh'], DataDict['TTLCh']
    
    DataInfo.update({'PiezoCh': PiezoCh, 'TTLCh': TTLCh, 'FileName': File[:-4]+'dict'})
    Txt.DictWrite(File[:-4]+'dict', DataInfo)


#%% Add StimType info to dicts
Files = glob('Recovery/*GPIAS/2*.dict'); Files.sort()

for File in Files:
    DataInfo = Txt.DictRead(File[:-4]+'dict')
    
    DataInfo['ExpInfo'] = {'StimType': ['Sound']}
    for Key in ['FreqOrder', 'FreqSlot', 'Freqs']:
        DataInfo['ExpInfo'][Key] = DataInfo[Key]
        del(DataInfo[Key])
    
    Txt.DictWrite(File[:-4]+'dict', DataInfo)


#%% Batch broken
Animal = 'GPIAZon'
AnalysisFile = Animal + '/' + Animal + '-Analysis.hdf5'

GPIASTimeBeforeTTL = 200   # in ms
GPIASTimeAfterTTL = 200    # in ms
FilterFreq = [100, 300]     # frequency for filter
FilterOrder = 3       # butter order

Paths = [
    'GPIAZon/GPIAZon_NaCl/2017-04-11_14-18-12_GPIAZon_NaCln01',
    'GPIAZon/GPIAZon_NaCl/2017-04-11_16-01-38_GPIAZon_NaCln03',
    'GPIAZon/GPIAZon_NaCl/2017-04-11_16-53-06_GPIAZon_NaCln04',
    'GPIAZon/GPIAZon_SSal/2017-04-13_12-48-26_GPIAZon_SSaln01',
    'GPIAZon/GPIAZon_SSal/2017-04-15_14-28-03_GPIAZon_SSaln04',
    'GPIAZon/GPIAZon_SSal/2017-04-21_14-47-12_GPIAZon_SSaln04',
    'GPIAZon/GPIAZon_SSal/2017-04-24_13-54-34_GPIAZon_SSaln05'
    ]

Files = [
    'GPIAZon/GPIAZon_NaCl/20170411141734-GPIAZon_NaCln01-GPIAS.hdf5',
    'GPIAZon/GPIAZon_NaCl/20170411160111-GPIAZon_NaCln03-GPIAS.hdf5',
    'GPIAZon/GPIAZon_NaCl/20170411165239-GPIAZon_NaCln04-GPIAS.hdf5',
    'GPIAZon/GPIAZon_SSal/20170413124748-GPIAZon_SSaln01-GPIAS.hdf5',
    'GPIAZon/GPIAZon_SSal/20170415142722-GPIAZon_SSaln04-GPIAS.hdf5',
    'GPIAZon/GPIAZon_SSal/20170421144530-GPIAZon_SSaln04-GPIAS.hdf5',
    'GPIAZon/GPIAZon_SSal/20170424135406-GPIAZon_SSaln05-GPIAS.hdf5'
    ]

ToDelete = [
    [['59'], [60]],
    [['110'], []],
    [['15'], []],
    [['51', '46', '36', '54', '38', '58', '23', '63', '70', '48', '34'], [35, 36, 37, 38, 39, 42, 43, 44, 47, 48, 49, 50, 51, 52, 53, 54]],
    [['13', '14', '16', '2', '20', '33', '37', '49', '50', '54', '55', '56', '70', '71', '80', '87'], []],
    [['49'], []],
    [['53'], []]
    ]

for Ind, DataPath in enumerate(Paths):
    Exp = DataPath.split('/')[1]
    RecFolder = DataPath.split('/')[-1]
    AnalysisKey = Exp + '/' + RecFolder.split('/')[-1]
    FigName = '/'.join([Animal, Exp, 'Figs', RecFolder+'-GPIAS'])
    FigPath = '/'.join(FigName.split('/')[:-1])
    makedirs(FigPath, exist_ok=True)
    
    Data = Hdf5F.LoadOEKwik(DataPath, AnalogTTLs=True, Unit='Bits')[0]
    DataInfo = Hdf5F.LoadDict('/DataInfo', Files[Ind])
    Proc = Hdf5F.GetProc(Data, 'OE')
    Rate = Data[Proc]['info']['0']['sample_rate']
    
    for Rec in Data[Proc]['data'].keys():
        BitVolts = 10000/(2**16)
        Data[Proc]['data'][Rec] = Data[Proc]['data'][Rec] * BitVolts
    
    DataInfo['PiezoCh'] = [int(DataInfo['PiezoCh'])]
    DataInfo['TTLCh'] = int(DataInfo['TTLCh'])
    
    for Path in ['Freqs', 'FreqOrder', 'FreqSlot']:
        DataInfo[Path] = Hdf5F.LoadDataset('/DataInfo/'+Path, Files[Ind])
    
    DataInfo['FreqOrder'][-3:][:,1] = -2
    
    DataInfo['FreqOrder'] = np.delete(DataInfo['FreqOrder'], ToDelete[Ind][1], 0)
    for Key in Data[Proc].keys(): 
        for r in ToDelete[Ind][0]: del(Data[Proc][Key][r])
    
    GPIASRec, XValues = GPIAS.Analysis(
                         Data[Proc]['data'], DataInfo, Rate, AnalysisFile, 
                         AnalysisKey, GPIASTimeBeforeTTL, GPIASTimeAfterTTL, 
                         FilterFreq, FilterOrder, Return=True)
    
    
    Plot.GPIAS.Traces(GPIASRec, XValues, DataInfo['SoundLoudPulseDur'], 
                            FigName, Save=True, Visible=True)
    
    del(GPIASRec, XValues)


# #%% Plot single freq

# from DataAnalysis.Plot import Plot

# Params = Plot.Set(Params=True)
# from matplotlib import rcParams; rcParams.update(Params)
# from matplotlib import pyplot as plt

# SoundPulseDur = DataInfo['Audio']['SoundLoudPulseDur']

# Ind1 = list(XValues).index(0)
# Ind2 = list(XValues).index(int(SoundPulseDur*1000))

# GPIASData = GPIASRec
# Freq = list(GPIASData['Trace'].keys())[0]

# SubTitle = Freq + ' Hz' + ' Index = ' + str(round(GPIASData['Index'][Freq]['GPIASIndex'], 4))
# LineNoGapLabel = 'No Gap'; LineGapLabel = 'Gap'
# SpanLabel = 'Sound Pulse'
# XLabel = 'time [ms]'; YLabel = 'voltage [mV]'

# Fig, Axes = plt.subplots(1, 1, figsize=(7, 12), sharex=True)
# Axes.axvspan(XValues[Ind1], XValues[Ind2], color='k', alpha=0.5, lw=0, label=SpanLabel)
# Axes.plot(XValues, GPIASData['Trace'][Freq]['NoGap'], color='r', label=LineNoGapLabel, lw=2)
# Axes.plot(XValues, GPIASData['Trace'][Freq]['Gap'], color='b', label=LineGapLabel, lw=2)

# AxArgs = {'title': SubTitle, 'ylabel': YLabel, 'xlabel': XLabel}
# Plot.Set(Ax=Axes, AxArgs=AxArgs)
# Axes.legend(loc='center left', bbox_to_anchor=(1.05, 0.5), prop={'size':6})

# plt.show()

# #%% Individual

# Animal = 'Prevention'
# Exp = 'PreventionA'
# RecFolder = '2017-05-21_11-18-08_Prevention_A2'
# ExpFile = '20170521111251-Prevention_A2-GPIAS.hdf5'

# GPIASTimeBeforeTTL = 200   # in ms
# GPIASTimeAfterTTL = 200    # in ms
# FilterFreq = [100, 300]     # frequency for filter
# FilterOrder = 3       # filter order
# Filter = 'butter'
# PiezoCh = [8]
# TTLCh = 6

# DataPath = Animal + '/' + Exp + '/' + RecFolder
# AnalysisFile = Animal + '/' + Animal + '-Analysis.hdf5'
# AnalysisKey = Exp + '/' + RecFolder.split('/')[-1]
# FigName = '/'.join([Animal, Exp, 'Figs', RecFolder+'-GPIAS'])
# FigPath = '/'.join(FigName.split('/')[:-1])
# makedirs(FigPath, exist_ok=True)

# Data = Hdf5F.LoadOEKwik(DataPath, AnalogTTLs=True, Unit='Bits')[0]
# DataInfo = Hdf5F.LoadDict('/DataInfo', Animal + '/' + Exp + '/' + ExpFile)
# Proc = Hdf5F.GetProc(Data, 'OE')
# Rate = Data[Proc]['info']['0']['sample_rate']

# for Rec in Data[Proc]['data'].keys():
#     BitVolts = 10000/(2**16)
#     Data[Proc]['data'][Rec] = Data[Proc]['data'][Rec] * BitVolts

# DataInfo['PiezoCh'] = [3]; DataInfo['TTLCh'] = 1
# #DataInfo['PiezoCh'] = PiezoCh
# #DataInfo['TTLCh'] = TTLCh
# #
# #for Path in ['Freqs', 'FreqOrder', 'FreqSlot']:
# #    DataInfo[Path] = Hdf5F.LoadDataset('/DataInfo/'+Path, Animal + '/' + Exp + '/' + ExpFile)
# #
# #DataInfo['FreqOrder'][-3:][:,1] = -2

# ## Fix stupid breaks in recs
# #import numpy as np
# #SOAB = DataAnalysis.GPIAS.CheckGPIASRecs(Data[Proc]['data'], [65000, 100000]); SOAB.sort()
# #Params = {'backend': 'TkAgg'}
# #from matplotlib import rcParams; rcParams.update(Params)
# #from matplotlib import pyplot as plt
# #R = {}
# #for r in ['1', '3', '4', '5'] + \
# #         [str(_) for _ in range(10, 23)] + \
# #         [str(_) for _ in range(32, 39)] + \
# #         [str(_) for _ in range(47, 58)] + \
# #         [str(_) for _ in range(68, 73)] + \
# #         [str(_) for _ in range(78, 82)] + \
# #         [str(_) for _ in range(85, 90)]: 
# #    print(r)
# #    R[r] = DataAnalysis.PSD(Data[Proc]['data'][r][:,2], Rate)
# #    plt.figure(); plt.semilogy(R[r][0], R[r][1])
# #    plt.figure(); plt.plot(Data[Proc]['data'][r][:,2])
# #    plt.show()

# # 20170411141734-GPIAZon_NaCln01-GPIAS.hdf5
# #DataInfo['FreqOrder'] = np.delete(DataInfo['FreqOrder'], [60], 0)
# #for Key in Data[Proc].keys(): del(Data[Proc][Key]['59'])

# # 20170411160111-GPIAZon_NaCln03-GPIAS.hdf5
# #for Key in Data[Proc].keys(): del(Data[Proc][Key]['110'])

# # 20170411165239-GPIAZon_NaCln04-GPIAS.hdf5
# #for Key in Data[Proc].keys(): del(Data[Proc][Key]['15'])

# # 20170413124748-GPIAZon_SSaln01-GPIAS.hdf5
# #SOAB = ['51', '46', '36', '54', '38', '58', '23', '63', '70', '48', '34']
# #ToDel = [_ for _ in range(35,40)] + \
# #        [_ for _ in range(42,45)] + \
# #        [_ for _ in range(47,55)]
# #DataInfo['FreqOrder'] = np.delete(DataInfo['FreqOrder'], ToDel, 0)
# #for Key in Data[Proc].keys(): 
# #    for r in SOAB: del(Data[Proc][Key][r])

# # 20170415142722-GPIAZon_SSaln04-GPIAS.hdf5
# #SOAB = ['13', '14', '16', '2', '20', '33', '37', '49', '50', '54', '55', '56', '70', '71', '80', '87']
# #for Key in Data[Proc].keys(): 
# #    for r in SOAB: del(Data[Proc][Key][r])

# # 20170521102604-Prevention_B1-GPIAS.hdf5
# #for r in ['50', '51', '90']: del(Data[Proc][r])
# #DataInfo['ExpInfo']['FreqOrder'] = np.concatenate((
# #        DataInfo['ExpInfo']['FreqOrder'][:51], 
# #        np.array([[-99, -99]]), 
# #        DataInfo['ExpInfo']['FreqOrder'][51:] ))

# # 20170521111251-Prevention_B2-GPIAS.hdf5
# #del(Data[Proc]['85'])

# ## Run Analysis
# GPIASRec, XValues = GPIAS.Analysis(
#                      Data[Proc]['data'], DataInfo, Rate, AnalysisFile, 
#                      AnalysisKey, GPIASTimeBeforeTTL, GPIASTimeAfterTTL, 
#                      FilterFreq, FilterOrder, Filter, Return=True)


# Plot.GPIAS.Traces(GPIASRec, XValues, DataInfo['SoundLoudPulseDur'], 
#                         FigName, Save=True, Visible=True)