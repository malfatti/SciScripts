#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jun 12 15:29:28 2017

@author: malfatti
"""
from IO.IO import RunProcess
from scipy.io import loadmat


def LoadClusters(FileList):
    Clusters = {}
    for File in FileList:
        Ch = File[-8:-4]
        
        ClusterFile = loadmat(File)
        Clusters[Ch] = {}
        Clusters[Ch]['Id'] = ClusterFile['cluster_class'][:, 0]
        Clusters[Ch]['Timestamps'] = ClusterFile['cluster_class'][:, 1]
        Clusters[Ch]['Spikes'] = ClusterFile['spikes'][:]
    
    return(Clusters)


def Run(Rate, Path):
    print('Clustering spikes...')
    
#    MLab = '/home/cerebro/Software/Programs/MatLabR2014a/bin/matlab'
    MLab = '/home/malfatti/Software/Programs/MatLabR2015a/bin/matlab'
    CmdCd = 'cd ' + Path + '; '
    CmdCluster = 'try, par.sr=' + str(int(Rate)) + '; ' + \
                 "Get_spikes('Files.txt', 'parallel', true, 'par', par);" + \
                 "Files = dir('./*spikes.mat'); Files = {Files.name}; " + \
                 "Do_clustering(Files, 'make_plots', false); end; quit"
    RunProcess([MLab, '-nosplash', '-r', CmdCd+CmdCluster])
    
    print('Done clustering.')
    return(None)

