#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jun 12 15:42:06 2017

@author: malfatti
"""
import numpy as np
import os

from IO.Txt import DictPrint, DictRead, DictWrite


def Read(DataFile, Path, ChList=[]):
    """ Read flat interleaved binary data and return it as a numpy array. Data 
        will be represented as Data[Channels, Samples]. The ChList argument 
        allows to choose which channels will be loaded. If empty, all channels 
        will be loaded.
        
        This function needs a text file containing a dictionary with the data 
        info. The minimal information needed is Info['Shape'] and 
        Info['DType']. """
    
    # Read Data and Info files
    with open(Path + '/' + DataFile, 'rb') as File: Raw = File.read()
    InfoFile = ''.join(DataFile.split('.')[:-1]) + '-Info.dict'
    Info = DictRead(Path + '/' + InfoFile)
    
    # Convert bytes to linear array
    RawData = np.fromstring(Raw, Info['DType'])
    
    # Reshape data according to info in the Info dictionary
    if ChList:
        Data = np.zeros((Info['Shape'][0], len(ChList)), dtype=Info['DType'])
        
        for Ind, Ch in enumerate(ChList):
            Data[:, Ind] = RawData[range(Ch-1,RawData.size,Info['Shape'][1])]
    else:
        Data = np.zeros((Info['Shape'][0], Info['Shape'][1]), 
                        dtype=Info['DType'])
        
        for Ch in range(Info['Shape'][1]):
            Data[:, Ch] = RawData[range(Ch,RawData.size,Info['Shape'][1])]
    
    return(Data, Info)


def Write(DataFile, Path, Data, Info={}):
    """ Write numpy array to flat interleaved binary file. Data will be 
        represented as ch1s1-ch2s1-ch3s1-...-chNs1-ch1s2-ch2s2...chNsN.
        
        Also, write a text file containing data info for data loading. """
    
    # Get info and generate path
    Info['Shape'] = Data.shape
    Info['DType'] = str(Data.dtype)
    InfoFile = ''.join(DataFile.split('.')[:-1]) + '-Info.dict'
    
    os.makedirs(Path, exist_ok=True)
    
    # Write text info file
    DictWrite(Path + '/' + InfoFile, Info)
    
    # Write interleaved data
    with open(Path + '/' + DataFile, 'wb') as File: File.write(Data.tobytes())
    
    return(None)

