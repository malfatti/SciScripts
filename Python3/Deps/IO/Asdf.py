#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Oct  7 09:47:01 2017

@author: malfatti
"""

import os
import numpy as np
from asdf import AsdfFile

## Level 0
def FitData(Data, Path, File):
    if Path[0] != '/': Path = '/' + Path
    
    if Path == '/': 
        if not os.path.isfile(File):
            return(Data)
        else:
            Tmp = Load('/', File)
            Data = {**Tmp, **Data}
            return(Data)
    else:
        PathBlocks = Path[1:].split('/')
        Paths = ['["'+ '"]["'.join(PathBlocks[:b+1]) + '"]' for b in range(len(PathBlocks))]
        
        Tmp = {}
        for P in Paths: exec('Tmp' + P + '={}')
        exec('Tmp' + Paths[-1] + '=Data')
        Data = {**Tmp}
        
        if not os.path.isfile(File):
            return(Data)
        else:
            Tmp = Load(Path, File)


def Write(Data, Path, File):
#    Data = FitData(Data)
    
    with AsdfFile(Data) as F: F.write_to(File)
    
    return(None)


def ItemPop(I):
    if type(I) == list:
        if True in ['NDArrayType' in str(type(_)) for _ in I]:
            I = [ItemPop(_) for _ in I]
            return(I)
        else:
            return(I)
    
    elif type(I) == dict:
        I = {Key: ItemPop(Val) for Key, Val in I.items()}
        return(I)
        
    elif type(I) in [str, float, int]: return(I)
        
    elif 'NDArrayType' in str(type(I)): 
        I = np.array(I, dtype=I.dtype)
        return(I)
    
    else:
        print('Type', type(I), 'not understood.')
        return(None)


## Level 1
def Load(Path, File):
#    Dict = {}; F = AsdfFile.open(File, mode='r')
    
    with AsdfFile.open(File, mode='r') as F:
        Dict = {Key: ItemPop(F.tree.get(Key)) 
                for Key in F.tree.keys() if 'asdf' not in Key}
        
#        while F.tree.keys():
#            p = F.tree.popitem()
#            if 'asdf' in p[0]: continue
#            
#            V = ItemPop(p[1])
#            Dict[p[0]] = V
        
    F.close()
    return(Dict)

