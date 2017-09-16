#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Jul  8 20:04:14 2017

@author: malfatti
"""
import OpenEphys, SettingsXML
import numpy as np

from DataAnalysis.DataAnalysis import BitsToVolts, UniqueStr
from glob import glob
from IO import Hdf5


## Level 0
def ApplyChannelMap(Data, ChannelMap):
    print('Retrieving channels according to ChannelMap... ', end='')
    ChannelMap = [_-1 for _ in ChannelMap]
    for R, Rec in Data.items():
        if Rec.shape[1] < len(ChannelMap): continue
        
        Chs = sorted(ChannelMap)
        Data[R][:, Chs] = Data[R][:, ChannelMap]
    
    return(Data)


def DataTouV(Data, RecChs):
    print('Converting to uV... ', end='')
    Data = {R: Rec.astype('float32') for R, Rec in Data.items()}
    Data = BitsToVolts(Data, RecChs)
    
    return(Data)


## Level 1
def DatLoad(Folder, Unit='uV', ChannelMap=[]):
    Files = glob(Folder+'/*.dat'); Files.sort()
    RecChs = SettingsXML.GetRecChs(Folder+'/settings.xml')[0]
    
    Data = {Proc: {} for Proc in RecChs.keys()}
    Rate = {Proc: [] for Proc in RecChs.keys()}
    for File in Files:
        Exp, Proc, Rec = File.split('/')[-1][10:-4].split('_')
        with open(File, 'rb') as F: Raw = F.read()
        
        Data[Proc][Rec] = np.fromstring(Raw, 'int16')
        ChNo = len(RecChs[Proc])
        SamplesPerCh = Data[Proc][Rec].shape[0]//ChNo
        
        Data[Proc][Rec] = Data[Proc][Rec].reshape((SamplesPerCh, ChNo))        
        # Still to be parsed. Assuming 30000.
        Rate[Proc].append(np.array(30000))
    
    for Proc in Data.keys():
        if Unit.lower() == 'uv': Data[Proc] = DataTouV(Data[Proc], RecChs[Proc])
        if ChannelMap: Data[Proc] = ApplyChannelMap(Data[Proc], ChannelMap)
        if len(np.unique(Rate[Proc])) == 1: Rate[Proc] = Rate[Proc][0]
    
    return(Data, Rate)


def KwikLoad(Folder, Unit='uV', ChannelMap=[]):
    Kwds = glob(Folder+'/*.kwd')
    XMLFile = glob(Folder+'/setting*.xml')[0]
    if Unit.lower() == 'uv': RecChs = SettingsXML.GetRecChs(XMLFile)[0]
    
    for Kwd in Kwds:
        Proc = Kwd[-11:-8]
        
        Data = {}; Rate = {}
        Data[Proc], Attrs = Hdf5.DataLoad('/', Kwd)
        Data[Proc] = {R: Rec['data'] for R, Rec in Data[Proc]['recordings'].items()}
        Rate[Proc] = [np.array(Rec['sample_rate']) for Rec in Attrs['recordings'].values()]
        if len(np.unique(Rate[Proc])) == 1: Rate[Proc] = Rate[Proc][0]
        
        if Unit.lower() == 'uv': Data[Proc] = DataTouV(Data[Proc], RecChs[Proc])
        if ChannelMap: Data[Proc] = ApplyChannelMap(Data[Proc], ChannelMap)
    
    print('Done.')
    return(Data, Rate)


def OpenEphysLoad(Folder, Unit='uV', ChannelMap=[]):
    OEs = glob(Folder+'/*continuous')
    XMLFile = glob(Folder+'/setting*.xml')[0]
    
    Chs = [_.split('/')[-1] for _ in OEs]
    Procs = UniqueStr([_[:3] for _ in Chs])
    
    Data = {_: {} for _ in Procs}
    
    if UniqueStr([len(_.split('.')[0].split('_')) for _ in Chs])[0] == 2:
        Recs = ['0']
    else:
        Recs = sorted(UniqueStr([_.split('.')[0].split('_')[-1] for _ in Chs]))
    
    Chs = {Proc: [_ for _ in Chs if _[:3] == Proc] for Proc in Procs}
    
    for Proc in Procs:
        Chs[Proc].sort(key=lambda x: [int(''.join(x.split('CH')[1:]).split('.')[0])])
        Data[Proc] = {Rec: [] for Rec in Recs}
        
        for Rec in Recs:
            Data[Proc][Rec] = OpenEphys.loadFolderToArray(Folder, channels=Chs[Proc], source=Proc)


## Level 2
def DataLoader(Folder, AnalogTTLs=True, Unit='uV', ChannelMap=[]):
    FilesExt = [F[-3:] for F in glob(Folder+'/*.*')]
    
    if 'kwd' in FilesExt:
        if AnalogTTLs:  Data, Rate = KwikLoad(Folder, Unit, ChannelMap)
#        else:  Data, Events, _, Files = Hdf5.OEKwikLoad(Folder, AnalogTTLs)
        
        return(Data, Rate)
    
    elif 'dat' in FilesExt:
        Data, Rate = DatLoad(Folder, Unit, ChannelMap)
        return(Data, Rate)
    
    elif 'ous' in FilesExt:
        print('To be implemented.'); return(None)
    
    else: print('Data format not supported.'); return(None)


