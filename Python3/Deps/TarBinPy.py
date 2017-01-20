#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jan 18 10:28:45 2017

@author: cerebro
"""
#%%
import json
import numpy as np

def DataWrite(DataFile, InfoFile, Data, Info):
    if Info['Transpose']: 
        with open(DataFile, 'wb') as File: File.write(Data.T.tobytes())
    else:
        with open(DataFile, 'wb') as File: File.write(Data.tobytes())
    
    Info['DataLineNo'] = Data.shape[0]
    Info['DataColumnNo'] = Data.shape[1]
    
#    with open(FileName, 'w') as File: File.write(str(Info))
    with open(InfoFile, 'w') as File: json.dump(Info, File)
    
    return(None)


def DataRead(DataFile, InfoFile, ChList=[]):
    with open(DataFile, 'rb') as File: Raw = File.read()
    with open(InfoFile, 'r') as File: Info = json.load(File)
    
    RawData = np.fromstring(Raw, Info['Format'])
    
    if ChList:
#        DataR = np.zeros((len(ChList), (RawData.size//Info['ChNo'])))
        DataR = np.zeros((len(ChList), Info['DataColumnNo']))
        DataR = np.array(DataR, Info['Format'])
        
        for Ind, Ch in enumerate(ChList):
            DataR[Ind, :] = RawData[range(Ch-1,RawData.size,Info['DataLineNo'])]
    else:
#        Data = np.zeros((Info['ChNo'], (RawData.size//Info['ChNo'])))
        DataR = np.zeros((Info['DataLineNo'], Info['DataColumnNo']))
        DataR = np.array(DataR, Info['Format'])
        
        for Ch in range(Info['DataLineNo']):
            DataR[Ch, :] = RawData[range(Ch,RawData.size,Info['DataLineNo'])]
    
    return(Data)


FileName = 'numpy2.flat'
Info = {'ChNo': 27,
        'Format': '<f'}
ChList = [1, 3, 10, 27]

Data = np.random.randn(Info['ChNo'], 5*60*30000)
Data = np.array(Data, Info['Format'])

DataWrite(FileName, Data)
DataTest = DataRead(FileName, Info, ChList)
