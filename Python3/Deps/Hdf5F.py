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
                if F['ExpInfo'][Exp].attrs['StimType'][0].decode() 
                == StimType]
            
    return(Exps)


def ExpExpInfo(RecFolder, DirList, FileName):
    ExpInfo = {}
    with h5py.File(FileName, 'r') as F:
#        Key = str(DirList.index(RecFolder))
        Key = "{0:02d}".format(DirList.index(RecFolder))
        ExpInfo['DVCoord'] = F['ExpInfo'][Key].attrs['DVCoord']
        ExpInfo['Hz'] = F['ExpInfo'][Key].attrs['Hz']
    
    return(ExpInfo)


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
        Keys = [Key for Key in F.keys() if 'ABRs' in Key]; Keys.sort()
        if len(Keys) > 1:
            print('Choose dataset to load:')
            for Ind, Key in enumerate(Keys):
                print(str(Ind), '=' , Key)
            Key = input(': ')
            Key = Keys[int(Key)]
        else:
            Key = Keys[0]
        
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
        Keys = [Key for Key in F.keys() if 'GPIAS' in Key]; Keys.sort()
        if len(Keys) > 1:
            print('Choose dataset to load:')
            for Ind, Key in enumerate(Keys):
                print(str(Ind), '=' , Key)
            Key = input(': ')
            Key = Keys[int(Key)]
        else: Key = Keys[0]
        
        XValues = F[Key].attrs['XValues']
        
        GPIAS = [[0] for _ in range(len(F[Key]))]
        for Freq in range(len(GPIAS)):
            GPIAS[Freq] = {}
            GPIAS[Freq]['NoGap'] = F[Key][str(Freq)]['NoGap'][:]
            GPIAS[Freq]['Gap'] = F[Key][str(Freq)]['Gap'][:]
    
    return(GPIAS, XValues)


def LoadOEFiles(RecFolder, AnalogTTLs):
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

