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
from datetime import datetime
from numbers import Number


def CheckGroup(FileName, Group):
    with h5py.File(FileName) as F:
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


def ExpDataInfo(FileName, DirList, StimType, Var='DataInfo'):
    DataInfo = {}
    with h5py.File(FileName) as F:
        for Key, Value in F['DataInfo'].items():
            DataInfo['SoundAmpF'] = {}
            for aKey, aValue in F['DataInfo']['SoundAmpF'].items():
                DataInfo['SoundAmpF'][aKey] = aValue[:]
        
        for bKey, bValue in F['DataInfo'].attrs.items():
            if isinstance(bValue, Number):
                DataInfo[bKey] = float(bValue)
            else:
                DataInfo[bKey] = bValue
        
        if Var == 'DataInfo':
            return(DataInfo)
        else:
            Exps = [DirList[int(Exp)] for Exp in F['ExpInfo'].keys() 
                    if F['ExpInfo'][Exp].attrs['StimType'][0].decode() 
                    == StimType]
            
            return(DataInfo, Exps)


def ExpExpInfo(FileName, RecFolder, DirList):
    ExpInfo = {}
    with h5py.File(FileName) as F:
#        Key = str(DirList.index(RecFolder))
        Key = "{0:02d}".format(DirList.index(RecFolder))
        ExpInfo['DVCoord'] = F['ExpInfo'][Key].attrs['DVCoord']
        ExpInfo['Hz'] = F['ExpInfo'][Key].attrs['Hz']
    
    return(ExpInfo)


def GPIASDataInfo(FileName):
    DataInfo = {}
    with h5py.File(FileName) as F:
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
    with h5py.File(FileName) as F:
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


def LoadDict(FileName, GroupName):
    Dict = {}
    with h5py.File(FileName) as F:
        for Key, Value in F[GroupName].attrs.items():
            if isinstance(Value, Number):
                Dict[Key] = float(Value)
            else:
                Dict[Key] = Value
    
    return(Dict)


def SoundCalibration(SBAmpFsFile, SoundBoard, Key):
    with h5py.File(SBAmpFsFile) as h5: 
        Var = h5[SoundBoard][Key][0]
    return(Var)


def SoundMeasurement(FileName, Var='SoundIntensity'):
    DataInfo = {}; SoundIntensity = {}
    with h5py.File(FileName) as h5:
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


def WriteABRs(FileName, ABRs, XValues):
    print('Saving data to ' + FileName)
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
    
    return(None)


