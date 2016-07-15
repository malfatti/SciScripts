# -*- coding: utf-8 -*-
"""
    Copyright (C) 2015  T. Malfatti
    
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

import h5py
import Kwik
import numpy as np
from datetime import datetime
from numbers import Number


def CheckGroup(FileName, Group):
    with h5py.File(FileName, 'r') as F:
        if Group in F.keys():
            print(Group + ' already exists.')
            print('Running this will erase previous analysis. Be careful!')
            Ans = input('Continue? [y/N] ')
            if Ans in ['y', 'Y', 'yes', 'Yes', 'YES']:
                return(True)
            else:
                return(False)
        else:
            return(True)


def LoadExpPerStim(StimType, DirList, FileName):
    with h5py.File(FileName, 'r') as F:
        Exps = [DirList[int(Exp)] for Exp in F['ExpInfo'].keys() 
                if np.string_(StimType) in F['ExpInfo'][Exp].attrs['StimType']]
            
    return(Exps)


def ExpExpInfo(RecFolder, DirList, FileName):
    ExpInfo = {}
    with h5py.File(FileName, 'r') as F:
#        Key = str(DirList.index(RecFolder))
        Key = "{0:02d}".format(DirList.index(RecFolder))
        ExpInfo['DVCoord'] = F['ExpInfo'][Key].attrs['DVCoord']
        ExpInfo['Hz'] = F['ExpInfo'][Key].attrs['Hz']
    
    return(ExpInfo)


def GetExpKeys(ExpStr, OpenedFile):
    Keys = [Key for Key in OpenedFile.keys() if ExpStr in Key]; Keys.sort()
    if len(Keys) > 1:
        print('Choose dataset to load:')
        for Ind, Key in enumerate(Keys):
            print(str(Ind), '=' , Key)
        Key = input(': ')
        Key = Keys[int(Key)]
    else:
        Key = Keys[0]
    
    return(Key)


def GPIASDataInfo(FileName):
    DataInfo = {}
    with h5py.File(FileName, 'r') as F:
        for Key, Value in F['DataInfo'].items():
            DataInfo['SoundBackgroundAmpF'] = {}
            DataInfo['SoundPulseAmpF'] = {}
            for aKey, aValue in F['DataInfo']['SoundPulseAmpF'].items():
                DataInfo['SoundPulseAmpF'][aKey] = aValue[:]
            
            for cKey, cValue in F['DataInfo']['SoundBackgroundAmpF'].items():
                DataInfo['SoundBackgroundAmpF'][cKey] = cValue[:]
        
        for bKey, bValue in F['DataInfo'].attrs.items():
            if isinstance(bValue, Number):
                DataInfo[bKey] = float(bValue)
            else:
                DataInfo[bKey] = bValue
        
        DataInfo['Freqs'] = F['DataInfo']['Freqs'][:]
        DataInfo['FreqOrder'] = F['DataInfo']['FreqOrder'][:]
        DataInfo['FreqSlot'] = F['DataInfo']['FreqSlot'][:]
        
        return(DataInfo)


def LoadABRs(FileName):
    with h5py.File(FileName, 'r') as F:
        Key = GetExpKeys('ABRs', F)
        
        ABRs = [0]*len(F[Key])
        for Freq in range(len(F[Key])):
            ABRs[Freq] = [0]*len(F[Key][str(Freq)])
            
            for AmpF in range(len(F[Key][str(Freq)])):
                ABRs[Freq][AmpF] = {}
                
                for DV in F[Key][str(Freq)][str(AmpF)].keys():
                    ABRs[Freq][AmpF][DV] = [0]*len(F[Key][str(Freq)][str(AmpF)][DV])
                    
                    for Trial in range(len(F[Key][str(Freq)][str(AmpF)][DV])):
                        ABRs[Freq][AmpF][DV][Trial] = F[Key][str(Freq)][str(AmpF)][DV][str(Trial)][:]
        
        XValues = F[Key].attrs['XValues'][:]
        
    return(ABRs, XValues)


def LoadClusters(FileName):
    with h5py.File(FileName, 'r') as F:
        Key = GetExpKeys('SpkClusters', F)
        
        Clusters = {}
        for RKey in Clusters.keys():
            Clusters[RKey] = {}
            
            for CKey in Clusters[RKey].keys():
                Path = '/'+Key+'/'+RKey+'/'+CKey
                
                Clusters[RKey][CKey] = {}
                Clusters[RKey][CKey]['ClusterClass'] = F[Path]['ClusterClass'][:]
                Clusters[RKey][CKey]['Timestamps'] = F[Path]['Timestamps'][:]
                Clusters[RKey][CKey]['Spikes'] = F[Path]['Spikes'][:]
                
                Clusters[RKey][CKey]['Info'] = {}
                Clusters[RKey][CKey]['Info']['Parameters'] = F[Path]['Spikes'].attrs['Parameters'][:]
#                Clusters[RKey][CKey]['Info']['InSpk'] = F[Path]['Spikes'].attrs['InSpk'][:]
    
    return(Clusters)


def LoadDict(Path, FileName, Attrs=True):
    Dict = {}
    with h5py.File(FileName, 'r') as F:
        if Attrs:
            for Key, Value in F[Path].attrs.items():
                if isinstance(Value, Number): Dict[Key] = float(Value)
                else: Dict[Key] = Value
        
        else:
            
            for Key, Value in F[Path].items(): Dict[Key] = Value[:]
    
    return(Dict)


def LoadDataset(Path, FileName):
    with h5py.File(FileName, 'r') as F: Dataset = F[Path][:]
    return(Dataset)


def LoadGPIAS(FileName):
    with h5py.File(FileName, 'r') as F:
        Key = GetExpKeys('GPIAS', F)
        XValues = F[Key].attrs['XValues']
        
        GPIAS = [[0] for _ in range(len(F[Key]))]
        for Freq in range(len(GPIAS)):
            GPIAS[Freq] = {}
            GPIAS[Freq]['NoGap'] = F[Key][str(Freq)]['NoGap'][:]
            GPIAS[Freq]['Gap'] = F[Key][str(Freq)]['Gap'][:]
    
    return(GPIAS, XValues)


def LoadOEKwik(RecFolder, AnalogTTLs):
    if AnalogTTLs: Raw, _, Spks, Files = Kwik.load_all_files(RecFolder)
    else: Raw, Events, Spks, Files = Kwik.load_all_files(RecFolder)
    
    print('Check if files are ok...')
    if 'Raw' not in locals():
        print('.kwd file is corrupted. Skipping dataset...')
        return(None)
    
    if not AnalogTTLs:
        if 'Events' not in locals():
            print('.kwe/.kwik file is corrupted. Skipping dataset...')
            return(None)
    
    print('Data from ', RecFolder, ' loaded.')
    
    if AnalogTTLs: return(Raw, Spks, Files)
    else: return(Raw, Events, Spks, Files)


def LoadTTLsLatency(FileName):
    Test


def LoadUnits(FileName):
    with h5py.File(FileName, 'r') as F:
        Key = GetExpKeys('Units', F)
        XValues = F[Key].attrs['XValues'][:]
        
        Units = {}
        for SKey in F[Key].keys():
            Units[SKey] = {}
            
            for FKey in F[Key][SKey].keys():
                Units[SKey][FKey] = {}
                
                for RKey in F[Key][SKey][FKey].keys():
                    Units[SKey][FKey][RKey] = {}
                    
                    for Ch in F[Key][SKey][FKey][RKey].keys():
                        Units[SKey][FKey][RKey][Ch] = {}
                        Path = '/'+Key+'/'+SKey+'/'+FKey+'/'+RKey+'/'+Ch
                        
                        SpkWFNo = len(list(F[Path]['SpkWF'].keys()))
                        Units[SKey][FKey][RKey][Ch]['SpkWF'] = [[] for _ in range(SpkWFNo)]
                        for DKey in F[Path]['SpkWF'].keys():
                            Units[SKey][FKey][RKey][Ch]['SpkWF'][int(DKey)] = F[Path]['SpkWF'][DKey][:]
                        Units[SKey][FKey][RKey][Ch]['NoOfSpks'] = F[Path]['NoOfSpks'][:]
    
    return(Units, XValues)


def SoundCalibration(SBAmpFsFile, SoundBoard, Key):
    with h5py.File(SBAmpFsFile, 'r') as h5: 
        Var = h5[SoundBoard][Key][0]
    return(Var)


def SoundMeasurement(FileName, Var='SoundIntensity'):
    DataInfo = {}; SoundIntensity = {}
    with h5py.File(FileName, 'r') as h5:
        if Var in ['DataInfo', 'SoundRec', 'SoundIntensity']:
            for Key,Val in h5['SoundRec'].attrs.items():
                DataInfo[Key] = Val
            
            if Var == 'SoundRec':
                SoundRec = [0]*len(DataInfo['NoiseFrequency'])
                
                for Freq in range(len(SoundRec)):
                    SoundRec[Freq] = [0]*len(DataInfo['SoundAmpF'])
                    Key = str(DataInfo['NoiseFrequency'][Freq][0]) + '-' + \
                          str(DataInfo['NoiseFrequency'][Freq][1])
                    
                    for AmpF in range(len(SoundRec[Freq])):
                        aKey = str(DataInfo['SoundAmpF'][AmpF])
                        if aKey == '0.0': aKey = '0'
                        
                        SoundRec[Freq][AmpF] = list(h5['SoundRec'][Key][aKey])
                
                return(SoundRec)
            
            elif Var == 'SoundIntensity':
                for Key in h5['SoundIntensity'].keys():
                    SoundIntensity[Key] = {}
                    
                    for AmpF in h5['SoundIntensity'][Key].keys():
                        SoundIntensity[Key][AmpF] = h5['SoundIntensity'][Key][AmpF][0]
                
                return(SoundIntensity)
            
            else:
                return(DataInfo)
        
        else:
            print('Supported variables: DataInfo, SoundRec, SoundIntensity.')


def WriteABRs(ABRs, XValues, FileName):
    print('Writing data to', FileName+'... ', end='')
    Now = datetime.now().strftime("%Y%m%d%H%M%S")
    Group = 'ABRs-' + Now
    with h5py.File(FileName) as F:
        F.create_group(Group)
        F[Group].attrs['XValues'] = XValues
            
        for Freq in range(len(ABRs)):
            for AmpF in range(len(ABRs[Freq])):
                for DV in ABRs[Freq][AmpF].keys():
                    Path = str(Freq) + '/' + str(AmpF) + '/' + DV
                    
                    F[Group].create_group(Path)
                    del(Path)
                    
                    for Trial in range(len(ABRs[Freq][AmpF][DV])):
                        F[Group][str(Freq)][str(AmpF)][DV][str(Trial)] = \
                            ABRs[Freq][AmpF][DV][Trial][:]
    
    print('Done.')
    return(None)


def WriteClusters(Clusters, FileName):
    print('Writing data to', FileName+'... ', end='')
    Now = datetime.now().strftime("%Y%m%d%H%M%S")
    Group = 'SpkClusters-' + Now
    with h5py.File(FileName) as F:
        for RKey in Clusters.keys():
            for CKey in Clusters[RKey].keys():
                Path = '/'+Group+'/'+RKey+'/'+CKey
                if Path not in F: F.create_group(Path)
                
                F[Path]['ClusterClass'] = Clusters[RKey][CKey]['ClusterClass']
                F[Path]['Timestamps'] = Clusters[RKey][CKey]['Timestamps']
                F[Path]['Spikes'] = Clusters[RKey][CKey]['Spikes']
                
                F[Path]['Spikes'].attrs['Parameters'] = Clusters[RKey][CKey]['Info']['Parameters']
#                F[Path]['Spikes'].attrs['InSpk'] = Clusters[RKey][CKey]['Info']['InSpk']
    
    print('Done.')
    return(None)


def WriteDict(Dict, Path, FileName, Attrs=True):
    print('Writing dictionary at', Path+'... ', end='')
    with h5py.File(FileName) as F:
        if Path not in F: F.create_group(Path)
        
        if Attrs:
            for Key, Value in Dict.items():
                F[Path].attrs[Key] = Value
        else:
            for Key, Value in Dict.items():
                F[Path][Key] = Value
    
    print('Done.')
    return(None)


def WriteExpInfo(StimType, DVCoord, Freq, FileName):
    print('Writing ExpInfo to', FileName+'... ', end='')
    with h5py.File(FileName) as F:
        if '/ExpInfo' not in F: F.create_group('ExpInfo')
        
        Key = "{0:02d}".format(len(list(F['ExpInfo'])))
        F['ExpInfo'].create_group(Key)
        
        F['ExpInfo'][Key].attrs['StimType'] = [np.string_(StimType)]
        F['ExpInfo'][Key].attrs['DVCoord'] = DVCoord
        F['ExpInfo'][Key].attrs['Hz'] = Freq
    
    print('Done.')
    return(None)


def WriteGPIAS(GPIAS, RecFolder, XValues, FileName):
    print('Writing data to', FileName+'... ', end='')
    Now = datetime.now().strftime("%Y%m%d%H%M%S")
    Group = 'GPIAS-' + RecFolder[10:] + '-' + Now
    with h5py.File(FileName) as F:
        F.create_group(Group)
        F[Group].attrs['XValues'] = XValues
        
        for Freq in range(len(GPIAS)):
            F[Group].create_group(str(Freq))
            
            F[Group][str(Freq)]['NoGap'] = GPIAS[Freq]['NoGap']
            F[Group][str(Freq)]['Gap'] = GPIAS[Freq]['Gap']
    
    print('Done.')
    return(None)


def WriteTTLslatency(TTLsLatency, XValues, FileName):
    print('Writing data to', FileName+'... ', end='')
    Now = datetime.now().strftime("%Y%m%d%H%M%S")
    Group = 'TTLsLatency-' + Now
    with h5py.File(FileName) as F:
        for Rec in TTLsLatency.keys():
            Path = '/' + Group + '/' + Rec; F.create_group(Path)
            
            for Key, Value in TTLsLatency[Rec].items():
                F[Path][Key] = Value
                
        F[Group]['XValues'] = XValues
    
    print('Done.')
    return(None)


def WriteUnits(Units, XValues, FileName):
    print('Writing data to', FileName+'... ', end='')
    Now = datetime.now().strftime("%Y%m%d%H%M%S")
    Group = 'Units-' + Now
    Thrash = []
    with h5py.File(FileName) as F:
        for SKey in Units.keys():
            for FKey in Units[SKey].keys():
                for RKey in Units[SKey][FKey].keys():
                    for Ch in Units[SKey][FKey][RKey].keys():
                        Path = '/'+Group+'/'+SKey+'/'+FKey+'/'+RKey+'/'+Ch
                        if Path not in F: F.create_group(Path)
                        if Path+'/SpkWF' not in F: F.create_group(Path+'/SpkWF')
                        
                        ClusterNo = len(Units[SKey][FKey][RKey][Ch]['SpkWF'])
                        for Cluster in range(ClusterNo):
                            for Key, Spk in enumerate(Units[SKey][FKey][RKey][Ch]['SpkWF'][Cluster]):
                                if not len(Spk):
                                    Thrash.append([Path, Cluster, Key])
                                    del(Units[SKey][FKey][RKey][Ch]['SpkWF'][Cluster][Key])
                                    continue
                                F[Path]['SpkWF'][str(Cluster)][str(Key)] = Spk
                        
                        F[Path]['NoOfSpks'] = Units[SKey][FKey][RKey][Ch]['NoOfSpks']
        
        F[Group].attrs['XValues'] = XValues
    
    if Thrash: 
        for _ in Thrash: print(_)
    
    return(None)