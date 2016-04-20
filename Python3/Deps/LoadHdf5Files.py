# -*- coding: utf-8 -*-
"""
Created on Tue Apr 19 08:59:41 2016

@author: cerebro
"""

import h5py
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


def ExpDataInfo(FileName, DirList, StimType):
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
        
        Exps = [DirList[int(Exp)] for Exp in F['ExpInfo'].keys() 
                if F['ExpInfo'][Exp].attrs['StimType'][0].decode() == StimType]
    
    return(DataInfo, Exps)


def ExpExpInfo(FileName, DirList, RecFolder):
    ExpInfo = {}
    with h5py.File(FileName) as F:
        Key = str(DirList.index(RecFolder))
        ExpInfo['DVCoord'] = F['ExpInfo'][Key].attrs['DVCoord']
        ExpInfo['Hz'] = F['ExpInfo'][Key].attrs['Hz']
    
    return(ExpInfo)


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


