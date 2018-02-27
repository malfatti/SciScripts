#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: T. Malfatti <malfatti@disroot.org>
@date: 2016-11-25
@license: GNU GPLv3 <https://raw.githubusercontent.com/malfatti/SciScripts/master/LICENSE>
@homepage: https://github.com/Malfatti/SciScripts
"""
import Intan, os
import numpy as np

from IO import Hdf5

## Level 0
def Load(File, ChannelMap=[]):
    Data = Intan.ReadData(File)
    
    try: len(Data['aux'][0])
    except TypeError: Data['aux'] = Data['aux'].reshape((Data['aux'].shape[0], 1))
    
    if ChannelMap: 
        ChannelMap = [_-1 for _ in ChannelMap]
        Data['analog'][:, (range(16))] = Data['analog'][:, ChannelMap]
    
    #XValues = Data['t']
    Data = np.hstack((Data['analog'], Data['aux']))
    
    return(Data)#, XValues)


## Level 1
def Write2Kwik(File):
    PathList = File.split('/')
    KPath = '/'.join(PathList[:-2]) + '/KwikFiles/' + PathList[-1][:-4]
    KFile = KPath + '/experiment1_100.raw.kwd'
    os.makedirs(KPath, exist_ok=True)
    
    Data = {}; Attrs = {}
    Data['data'] = Load(File)
    Attrs['sample_rate'] = np.array([25000])
    
    Hdf5.DictWrite(Data, '/recordings/0', KFile, Attrs=False)
    Hdf5.DictWrite(Attrs, '/recordings/0', KFile)
    
    
    