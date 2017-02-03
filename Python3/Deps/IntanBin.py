#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Nov 25 13:29:22 2016

@author: malfatti
"""

import Hdf5F
import Intan
import numpy as np
from glob import glob
from scipy import io

def IntLoad(File):
    Data = Intan.ReadData(File)
    
    try: len(Data['aux'][0])
    except TypeError: Data['aux'] = Data['aux'].reshape((Data['aux'].shape[0], 1))
    
    XValues = Data['t']
    Data = np.hstack((Data['analog'], Data['aux']))
    
    return(Data, XValues)


#Path = 'SpkCh/'
#File = 'Test.int'
#
#ClusterList = glob(Path+'times_*'); ClusterList.sort()
#
#Clusters = {'00': {}}
#for File in ClusterList:
#    Ch = File[-6:-4]
#    
#    if Ch[0] == '_': Ch = 'Ch' + "{0:02d}".format(int(Ch[1]))
#    else: Ch = 'Ch' + Ch
#    
#    ClusterFile = io.loadmat(File)
#    Clusters['00'][Ch] = {}
#    Clusters['00'][Ch]['ClusterClass'] = ClusterFile['cluster_class'][:, 0]
#    Clusters['00'][Ch]['Timestamps'] = ClusterFile['cluster_class'][:, 1]
#    Clusters['00'][Ch]['Spikes'] = ClusterFile['spikes'][:]
#
#Hdf5F.WriteClusters(Clusters, './SpkClusters.hdf5')