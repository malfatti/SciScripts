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
import tarfile

def DataRead(DataFile, InfoFile, ChList=[]):
    """ Read flat interleaved binary data and return it as a numpy array. Data 
        will be represented as Data[Channels, Samples]. The ChList argument 
        allows to choose which channels will be loaded. If empty, all channels 
        will be loaded.
        
        This function needs a text file containing a dictionary with the data 
        info. The minimal information needed is Info['DataShape'] and 
        Info['DataDType']. """
    
    # Read Data and Info files
    with open(DataFile, 'rb') as File: Raw = File.read()
    Info = DictRead(InfoFile)
    
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


def DataWrite(DataFile, InfoFile, Data, Info={}):
    """ Write numpy array to flat interleaved binary file. Data will be 
        represented as ch1s1-ch2s1-ch3s1-...-chNs1-ch1s2-ch2s2...chNsN.
        
        Also, write a text file containing data info for data loading. """
    
    # Get info
    Info['DataShape'] = Data.shape
    Info['DataDType'] = str(Data.dtype)
    
    # Write text info file
    DictWrite(InfoFile, Info)
    
    # Write interleaved data
    with open(DataFile, 'wb') as File: File.write(Data.T.tobytes())
    
    return(None)


def DictRead(DictFile):
    with open(DictFile, 'r') as File: Dict = json.load(File)
    return(Dict)


def DictWrite(DictFile, Dict):
    with open(DictFile, 'w') as File: json.dump(Dict, File)
    return(None)

    
#%% Test Write/Read
DataFile = 'Data.flat'; InfoFile = 'Info.dict'
Info = {'ChNo': 27,
        'Format': '<f'}
#ChList = [1, 3, 10, 27]
ChList = []

Data = np.random.randn(Info['ChNo'], 30000)
Data = np.array(Data, Info['Format'])
#Data = Data.T

DataWrite(DataFile, InfoFile, Data, Info)
DataTest, InfoTest = DataRead(DataFile, InfoFile, ChList)

#%% Test WriteTar
Path = 'Exp/Conn/Setup/'
os.makedirs(Path, exist_ok=True)

TarFile = 'Exp.tar'
DataFile = 'Data.flat'; InfoFile = 'Info.dict'
DataFile = Path + DataFile; InfoFile = Path + InfoFile

Info = {'ChNo': 27, 'Format': '<f'}
Data = np.random.randn(Info['ChNo'], 30000)
Data = np.array(Data, Info['Format'])

DataWrite(DataFile, InfoFile, Data, Info)
FileList = glob(Path + '*', recursive=True); FileList.sort()

with tarfile.open(TarFile, 'a') as F:
    for File in FileList:
        F.add(File)
