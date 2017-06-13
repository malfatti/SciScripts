#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Nov 25 13:29:22 2016

@author: malfatti
"""
import Intan
import numpy as np


def IntLoad(File, ChannelMap=[]):
    Data = Intan.ReadData(File)
    
    try: len(Data['aux'][0])
    except TypeError: Data['aux'] = Data['aux'].reshape((Data['aux'].shape[0], 1))
    
    if ChannelMap: 
        ChannelMap = [_-1 for _ in ChannelMap]
        Data['analog'][:, (range(16))] = Data['analog'][:, ChannelMap]
    
    XValues = Data['t']
    Data = np.hstack((Data['analog'], Data['aux']))
    
    return(Data, XValues)

