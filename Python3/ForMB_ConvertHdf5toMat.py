#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Apr 20 14:04:13 2018

@author: malfatti
"""
#%%
from glob import glob
from IO import Hdf5
from IO.Txt import DictFlat
import numpy as np
import scipy.io as sio

Files = glob('/home/cerebro/Martina/**/Anim*.hdf5', recursive=True)
OutPath = '/home/cerebro/Martina/Mat'

for File in Files:
    FileName = OutPath + '/' + File.split('/')[-1][:-4]+'mat'
    if FileName in glob(OutPath+'/*'): 
        print('Skipping', File, '...')
        continue
    
    print('Loading', File, '...')
    Data = Hdf5.DataLoad('/', File)[0]
    Data = DictFlat(Data)
    
    Keys = {K: K.replace('-', '_') for K in Data.keys()}
    for Old,New in Keys.items(): Data[New] = Data.pop(Old)
    
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
    
    sio.savemat(FileName, Data)
    print('Written to', FileName, '.')