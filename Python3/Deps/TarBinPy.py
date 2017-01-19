#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jan 18 10:28:45 2017

@author: cerebro
"""
#%%
import struct
import numpy as np

def DataWrite(Data):
    with open('numpy2.dat', 'wb') as file: file.write(Data.tobytes())
#        for i in range(len(Data)):
#            Data.tofile(stream)

def struct_approach(Data):
    with open('structnp.dat', 'wb') as stream:
        for Sample in range(Data.shape[1]):
            for Ch in range(Data.shape[0]):
                s = struct.pack('<f', Data[Ch, Sample])
                stream.write(s)

def structread():
    ChNo = 27
    with open('numpy2.dat', 'rb') as stream: Raw = stream.read()
    
    RawData = np.fromstring(Raw, '<f')
    
    DataRead = np.zeros((ChNo, (RawData.size//ChNo)))
    DataRead = np.array(DataRead, 'Float32')
    
    for Ch in range(ChNo):
        DataRead[Ch, :] = RawData[range(Ch,RawData.size,ChNo)]


Data = np.random.randn(27, 5*60*30000)
Data = np.array(Data, 'Float32')
# https://gist.github.com/deanmalmgren/fd1714799dc5b5643b87#file-write_profiler-py