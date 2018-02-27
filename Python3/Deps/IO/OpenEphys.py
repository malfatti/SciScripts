#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Jul  8 20:04:14 2017

@author: malfatti
"""
import OpenEphys, SettingsXML
import numpy as np

from DataAnalysis import DataAnalysis
from glob import glob
from IO import Hdf5
from IO import Intan
from IO.Txt import DictPrint


## Level 0
def ApplyChannelMap(Data, ChannelMap):
    print('Retrieving channels according to ChannelMap... ', end='')
    ChannelMap = [_-1 for _ in ChannelMap]
    for R, Rec in Data.items():
        if Rec.shape[1] < len(ChannelMap): continue
        
        Chs = sorted(ChannelMap)
        Data[R][:, Chs] = Data[R][:, ChannelMap]
    
    return(Data)


def ChooseProcs(XMLFile, Procs):
    ProcList = SettingsXML.GetRecChs(XMLFile)[1]
    ProcList = {Id: Name for Id, Name in ProcList.items() if Id in Procs}
    
    print(DictPrint(ProcList))
    Procs = input('Which Procs should be kept (comma separated) ? ')
    Procs = [_ for _ in Procs.split(',')]
    
    return(Procs)


def DataTouV(Data, RecChs):
    print('Converting to uV... ', end='')
    Data = {R: Rec.astype('float32') for R, Rec in Data.items()}
    Data = DataAnalysis.BitsToVolts(Data, RecChs)
    
    return(Data)


def EventsLoad(Folder):
    Files = glob(Folder+'/*.events'); Files.sort()
    if len(Files) > 1: print('Multiple sessions not supported yet.'); return(None)
#    for File in Files:
#        Session = File.split('.')[0].split('_')[-1]
        
    EventsDict = OpenEphys.loadEvents(Folder+'/'+Files[0])
    return(EventsDict)


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
    if Unit.lower() == 'uv': 
        XMLFile = glob(Folder+'/setting*.xml')[0]
        RecChs = SettingsXML.GetRecChs(XMLFile)[0]
    
    for Kwd in Kwds:
        Proc = Kwd[-11:-8]
        
        Data = {}; Rate = {}
        Data[Proc], Attrs = Hdf5.DataLoad('/recordings', Kwd)
        Data[Proc] = {R: Rec['data'] for R, Rec in Data[Proc].items()}
        Rate[Proc] = [np.array(Rec['sample_rate']) for Rec in Attrs.values()]
        if len(np.unique(Rate[Proc])) == 1: Rate[Proc] = Rate[Proc][0]
        
        if Unit.lower() == 'uv': Data[Proc] = DataTouV(Data[Proc], RecChs[Proc])
        if ChannelMap: Data[Proc] = ApplyChannelMap(Data[Proc], ChannelMap)
    
    print('Done.')
    return(Data, Rate)


def OpenEphysLoad(Folder, Unit='uV', ChannelMap=[]):
    OEs = glob(Folder+'/*continuous')
    XMLFile = glob(Folder+'/setting*.xml')[0]
    
    Chs = [_.split('/')[-1] for _ in OEs]
    Procs = DataAnalysis.UniqueStr([_[:3] for _ in Chs])
    
    if len(Procs) > 1: Procs = ChooseProcs(XMLFile, Procs)
    
    Data = {_: {} for _ in Procs}; Rate = {_: {} for _ in Procs}
    Chs = {Proc: [_ for _ in Chs if _[:3] == Proc] for Proc in Procs}
    
    for P, Proc in Chs.items():
        Type = Chs[P][0].split('_')[-1].split('.')[0][:-1]
        print(Type)
        Chs[P] = sorted(Proc, key=lambda x: int(x.split('_'+Type)[1].split('_')[0].split('.')[0]))
    
    for Proc in Data.keys():
        ACh = Chs[Proc][0].split('.')[0]
        OEData = OpenEphys.loadFolder(Folder, source=Proc)
        Rate[Proc] = int(OEData[ACh]['header']['sampleRate'])
        
        Recs = np.unique(OEData[ACh]['recordingNumber'])
        BlockSize = int(OEData[ACh]['header']['blockLength'])
        for Rec in Recs:
            R = str(int(Rec))
            RecInd = np.where(OEData[ACh]['recordingNumber'].repeat(BlockSize) == Rec)
            Data[Proc][R] = [OEData[_.split('.')[0]]['data'][RecInd] for _ in Chs[Proc]]
            Data[Proc][R] = np.array(Data[Proc][R]).T
    
        if Unit.lower() == 'uv': 
            ChsInfo = [OEData[_.split('.')[0]]['header']['bitVolts'] for _ in Chs[Proc]]
            ChsInfo = {str(Ch): {'gain': BitVolt} for Ch, BitVolt in enumerate(ChsInfo)}
            Data[Proc] = DataTouV(Data[Proc], ChsInfo)
        
        if ChannelMap: Data[Proc] = ApplyChannelMap(Data[Proc], ChannelMap)
    
    return(Data, Rate)


## Level 2
def GetRecs(Folder):
    FilesExt = [F[-3:] for F in glob(Folder+'/*.*')]
    
    if 'kwd' in FilesExt:
        Kwds = glob(Folder+'/*.kwd')
        
        for Kwd in Kwds:
            Proc = Kwd[-11:-8]
            
            Recs = {}
            Recs[Proc] = Hdf5.DataLoad('/recordings', Kwd)[0]
            Recs[Proc] = [R for R in Recs[Proc].keys()]
        
    elif 'dat' in FilesExt:
        Files = glob(Folder+'/*.dat'); Files.sort()
        RecChs = SettingsXML.GetRecChs(Folder+'/settings.xml')[0]
        
        Recs = {Proc: [] for Proc in RecChs.keys()}
        
        for File in Files:
            _, Proc, Rec = File.split('/')[-1][10:-4].split('_')
            Recs[Proc].append(Rec)
        
    elif 'ous' in FilesExt:
        OEs = glob(Folder+'/*continuous')
        Chs = [_.split('/')[-1] for _ in OEs]
        Procs = DataAnalysis.UniqueStr([_[:3] for _ in Chs])
        Recs = {}
        
        for Proc in Procs:
            ACh = Chs[Proc][0].split('.')[0]
            OEData = OpenEphys.loadFolder(Folder, source=Proc)
            R = np.unique(OEData[ACh]['recordingNumber'])
            Recs[Proc] = str(int(R))
    
    return(Recs)


def DataLoader(Folder, Unit='uV', ChannelMap=[], AnalogTTLs=True):
    FilesExt = [F[-3:] for F in glob(Folder+'/*.*')]
    
    if 'kwd' in FilesExt: Data, Rate = KwikLoad(Folder, Unit, ChannelMap)
    elif 'dat' in FilesExt: Data, Rate = DatLoad(Folder, Unit, ChannelMap)
    elif 'ous' in FilesExt: Data, Rate = OpenEphysLoad(Folder, Unit, ChannelMap)
    else: print('Data format not supported.'); return(None)
    
    if not AnalogTTLs:
        if 'kwd' in FilesExt: 
            Kwds = glob(Folder+'/*.events')
            if len(Kwds) > 1: print('Multiple sessions not supported yet.'); return(None)
            
            EventsDict = 'ToBeContinued'
        else:
            EventsDict = EventsLoad(Folder)
        
        return(Data, Rate, EventsDict)
    else:
        return(Data, Rate)
        

