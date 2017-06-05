#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jan 18 10:28:45 2017

@author: cerebro
"""
#%%
#from datetime import datetime
from ast import literal_eval
from glob import glob

#import json
import numpy as np
import os
import shutil
import tarfile

### Level 0 ###
def DictRead(File):
    Dict = literal_eval(open(File).read())
    return(Dict)


def DictWrite(File, Str):
    with open(File, 'w') as F: F.write(Str)
    return(None)


def PrintDict(value, htchar='    ', itemchar=' ', breaklineat='auto', lfchar='\n', indent=0):
    ''' Modified from y.petremann's code.
        Added options to set item separator for list or tuple and to set a number
        of items per line, or yet, to calculate items per line so it will not 
        have more than 80 chars per line.
        Source: https://stackoverflow.com/a/26209900 '''
    
    nlch = lfchar + htchar * (indent + 1)
    if type(value) is dict:
        items = [
            nlch + repr(key) + ': ' + PrintDict(value[key], htchar, itemchar, breaklineat, lfchar, indent + 1)
            for key in value
        ]
        
        return '{%s}' % (','.join(items) + lfchar + htchar * indent)
    
    elif type(value) is list:
        items = [
            itemchar + PrintDict(item, htchar, itemchar, breaklineat, lfchar, indent + 1)
            for item in value
        ]
        
        if breaklineat == 'auto':
           bl = int((80 - (len(htchar)*(indent + 1)))/
                (int((sum([len(i)+4 for i in items])-2)/len(items))))
        else: bl = breaklineat
       
        if len(items) > bl:
            for i in list(range(bl, len(items), bl)):
                items[i] = lfchar + htchar*(indent+1) + '  ' + items[i]
        
        return '[%s]' % (','.join(items))
    
    elif type(value) is tuple:
        items = [
            itemchar + PrintDict(item, htchar, itemchar, breaklineat, lfchar, indent + 1)
            for item in value
        ]
        
        if breaklineat == 'auto':
           bl = (80 - (len(htchar)*(indent + 1)))// \
                ((sum([len(i)+4 for i in items])-2)//len(items))
        else: bl = breaklineat
       
        if len(items) > bl:
            for i in list(range(bl, len(items), bl)):
                items[i] = lfchar + htchar*(indent+1) + '  ' + items[i]
        
        return '(%s)' % (','.join(items))
    
    else:
        return repr(value)


def TarWrite(TarFile, FileList):
    with tarfile.open(TarFile, 'a') as F:
        for File in FileList: F.add(File)
    
    return(None)


#def TarWriteNew(TarFile, FileList):
#    with tarfile.open(TarFile, 'a') as F:
#        FilesInTar = F.getnames()
#        FilesToExclude = []
#        FilesToAdd = []
#        
#        for File in FileList:
#            if File in FilesInTar: FilesToExclude.append(File)
#            else: FilesToAdd.append(File)
#    
#    
#    if FilesToExclude:
#        Now = '_' + datetime.now().strftime("%Y%m%d%H%M%S")
#        
#        with tarfile.open(TarFile + Now, 'a') as NewF:
#            with tarfile.open(TarFile, 'r') as F:
#                for File in FileList: NewF.add(File)
#                
#                for File in FilesInTar:
#                    if File not in FileList:
#                        FileExtr = F.extractfile(File)
#                        NewF.add(FileExtr)
#    else:
#        for File in FileList:
#            F.add(File)


### Level 1 ###
def BinRead(DataFile, Path, ChList=[]):
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


def BinWrite(DataFile, Path, Data, Info={}):
    """ Write numpy array to flat interleaved binary file. Data will be 
        represented as ch1s1-ch2s1-ch3s1-...-chNs1-ch1s2-ch2s2...chNsN.
        
        Also, write a text file containing data info for data loading. """
    
    # Get info and generate path
    Info['Shape'] = Data.shape
    Info['DType'] = str(Data.dtype)
    InfoFile = ''.join(DataFile.split('.')[:-1]) + '-Info.dict'
    
    os.makedirs(Path, exist_ok=True)
    
    # Write text info file
    DictWrite(Path + '/' + InfoFile, PrintDict(Info))
    
    # Write interleaved data
    with open(Path + '/' + DataFile, 'wb') as File: File.write(Data.tobytes())
    
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
#Exp = 'SC-20170121-132017'; Conn = 'I-M-M-S'; Setup = 'UnitRec'
#Path = Exp + '/' + Conn + '/' + Setup + '/'
#TarFile = Exp + '.tar'
#DataFile = Exp +'.flat'
#
#Info = {'ChNo': 27, 'Format': '<f'}
#Data = np.random.randn(Info['ChNo'], 30000)
#Data = np.array(Data, Info['Format'])
#
#DataWrite(DataFile, TarFile, Path, Data, Info)
#
##F = tarfile.open(TarFile, 'r')
##F.getmembers()
