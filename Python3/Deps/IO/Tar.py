#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jan 18 10:28:45 2017

@author: cerebro
"""
#%%
#from datetime import datetime
#from glob import glob
#from IO.Txt import DictRead, DictWrite

#import json
#import numpy as np
#import os
#import shutil
import tarfile

### Level 0 ###
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

### Level 2 ###
#def DataWrite(DataFile, TarFile, DataPath, Data, Info):
#    BinWrite(DataFile, DataPath, Data, Info)
#    
#    FileList = glob(DataPath + '**/*.*', recursive=True); FileList.sort()
#    
#    TarWrite(TarFile, FileList)
#    DataDir = DataPath.split('/')[0]
#    shutil.rmtree(DataDir)

    
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
