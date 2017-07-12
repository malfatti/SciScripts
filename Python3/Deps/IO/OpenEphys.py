#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Jul  8 20:04:14 2017

@author: malfatti
"""
import SettingsXML
import numpy as np

from DataAnalysis.DataAnalysis import BitsToVolts
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
    if Unit.lower() == 'uv': RecChs = SettingsXML.GetRecChs(Folder+'/settings.xml')[0]
    
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


## Level 2
def DataLoader(Folder, AnalogTTLs=True):
    FilesExt = [F[-3:] for F in glob(Folder+'/*.*')]
    
    if 'kwd' in FilesExt:
        if AnalogTTLs:  Data, Rate = KwikLoad(Folder)
#        else:  Data, Events, _, Files = Hdf5.OEKwikLoad(Folder, AnalogTTLs)
        
        return(Data, Rate)
    
    elif 'dat' in FilesExt:
        Data, Rate = DatLoad(Folder)
        return(Data, Rate)
    
    elif 'ous' in FilesExt:
        print('To be implemented.'); return(None)
    
    else: print('Data format not supported.'); return(None)


