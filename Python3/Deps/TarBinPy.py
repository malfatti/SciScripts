#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jan 18 10:28:45 2017

@author: cerebro
"""
#%%
from glob import glob

import json
import numpy as np
import os
import shutil
import tarfile

### Level 0 ###
def DictRead(DictFile):
    with open(DictFile, 'r') as File: Dict = json.load(File)
    return(Dict)


def DictWrite(DictFile, Dict):
    with open(DictFile, 'w') as File: json.dump(Dict, File)
    return(None)


def TarWrite(TarFile, FileList):
    with tarfile.open(TarFile, 'a') as F:
        FilesInTar = F.getnames()
        FilesToExclude = []
        
        for File in FileList:
            if File in FilesInTar:
                FilesToExclude.append(File)
            else:
                F.add(File)
        
        if FilesToExclude:
            


### Level 1 ###
def BinRead(DataFile, Path, ChList=[]):
    """ Read flat interleaved binary data and return it as a numpy array. Data 
        will be represented as Data[Channels, Samples]. The ChList argument 
        allows to choose which channels will be loaded. If empty, all channels 
        will be loaded.
        
        This function needs a text file containing a dictionary with the data 
        info. The minimal information needed is Info['DataShape'] and 
        Info['DataDType']. """
    
    # Read Data and Info files
    with open(Path + DataFile, 'rb') as File: Raw = File.read()
    InfoFile = ''.join(DataFile.split('.')[:-1]) + '-Info.dict'
    Info = DictRead(Path + InfoFile)
    
    # Convert bytes to linear array
    RawData = np.fromstring(Raw, Info['DataDType'])
    
    # Reshape data according to info in the Info dictionary
    if ChList:
        Data = np.zeros((len(ChList), Info['DataShape'][1]))
        Data = np.array(Data, Info['DataDType'])
        
        for Ind, Ch in enumerate(ChList):
            Data[Ind, :] = RawData[range(Ch-1,RawData.size,Info['DataShape'][0])]
    else:
        Data = np.zeros((Info['DataShape'][0], Info['DataShape'][1]))
        Data = np.array(Data, Info['DataDType'])
        
        for Ch in range(Info['DataShape'][0]):
            Data[Ch, :] = RawData[range(Ch,RawData.size,Info['DataShape'][0])]
    
    return(Data, Info)


def BinWrite(DataFile, Path, Data, Info={}):
    """ Write numpy array to flat interleaved binary file. Data will be 
        represented as ch1s1-ch2s1-ch3s1-...-chNs1-ch1s2-ch2s2...chNsN.
        
        Also, write a text file containing data info for data loading. """
    
    # Get info and generate path
    Info['DataShape'] = Data.shape
    Info['DataDType'] = str(Data.dtype)
    InfoFile = ''.join(DataFile.split('.')[:-1]) + '-Info.dict'
    
    if Path[-1] != '/': Path = Path + '/'
    os.makedirs(Path, exist_ok=True)
    
    # Write text info file
    DictWrite(Path + InfoFile, Info)
    
    # Write interleaved data
    with open(Path + DataFile, 'wb') as File: File.write(Data.T.tobytes())
    
    return(None)


### Level 2 ###
def DataWrite(DataFile, TarFile, DataPath, Data, Info):
    BinWrite(DataFile, DataPath, Data, Info)
    
    FileList = glob(DataPath + '**/*.*', recursive=True); FileList.sort()
    
    TarWrite(TarFile, FileList)
    DataDir = DataPath.split('/')[0]
    shutil.rmtree(DataDir)

    
#%% Test Write/Read
#DataFile = 'Data.flat'; InfoFile = 'Info.dict'
#Info = {'ChNo': 27,
#        'Format': '<f'}
##ChList = [1, 3, 10, 27]
#ChList = []
#
#Data = np.random.randn(Info['ChNo'], 30000)
#Data = np.array(Data, Info['Format'])
##Data = Data.T
#
#BinWrite(DataFile, InfoFile, Data, Info)
#DataTest, InfoTest = BinRead(DataFile, InfoFile, ChList)

#%% Test WriteTar
Exp = 'SC-20170121-132017'; Conn = 'I-M-M-S'; Setup = 'UnitRec'
Path = Exp + '/' + Conn + '/' + Setup + '/'
TarFile = Exp + '.tar'
DataFile = Exp +'.flat'

Info = {'ChNo': 27, 'Format': '<f'}
Data = np.random.randn(Info['ChNo'], 30000)
Data = np.array(Data, Info['Format'])

DataWrite(DataFile, TarFile, Path, Data, Info)

F = tarfile.open(TarFile, 'r')
F.getmembers()
