#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Apr 20 14:04:13 2018

@author: malfatti
"""
#%%
from IO import Hdf5
from IO.Txt import DictFlat
import numpy as np
import scipy.io as sio

File = '/home/malfatti/NotSynced/Downloads/Animal_868-2016-06-12_03-15-29.hdf5'

Data = Hdf5.DataLoad('/', File)[0]
Data = DictFlat(Data)

KeysToDel = []; KeysToAdd = {}
for K, V in Data.items():
    if 'Spks' in K:
        SpksKey = '_'.join(K.split('_')[:-1])
        SpkNo = int(K.split('_')[-1])
        SpksNo = len([k for k in Data.keys() if SpksKey in k])
        
        if SpksKey not in KeysToAdd: KeysToAdd[SpksKey] = np.zeros((V.shape[0], SpksNo))
        KeysToAdd[SpksKey][:,SpkNo] = V
        
        KeysToDel.append(K)

for K in KeysToDel: del(Data[K])
Data = {**Data, **KeysToAdd}

sio.savemat(File[:-4]+'mat', Data)
