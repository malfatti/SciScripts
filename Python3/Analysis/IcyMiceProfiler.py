#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: malfatti

Icy Mice Profiler analysis
"""

#%% Icy XMLs to Dict
import DataAnalysis
import numpy as np

from glob import glob
from itertools import combinations
from scipy import io
from xml.etree import ElementTree

Params = DataAnalysis.Plot.Set(Params=True)
from matplotlib import rcParams; rcParams.update(Params)
from matplotlib import pyplot as plt

Files = glob('/home/malfatti/Barbara/Analysis/*/*/*allEventsExported.xml'); Files.sort()
CoordFiles = ['.'.join(F.split('.')[:-2]) for F in Files]
FPS = 23

Recs = {
    # G3
    'T0000002-H264': ['G3', 'MSFS'],
    'T0000003-H264': ['G3', 'MDFD'],
    'T0000005-H264': ['G3', 'MDFS'],
    'T0000006-H264': ['G3', 'MSFD'],
    
    # G4
    'fc2_save_2017-02-21-193715-0000': ['G4', 'MDFS',],
    'fc2_save_2017-02-21-201357-0000': ['G4', 'MSFD',],
    'fc2_save_2017-02-21-220718-0000': ['G4', 'MSFS',],
    'fc2_save_2017-02-21-224114-0000': ['G4', 'MDFD',],
    
    # G5
    'fc2_save_2017-02-22-183437-0000': ['G5', 'MDFD',],
    'fc2_save_2017-02-22-190844-0000': ['G5', 'MSFS',],
    'fc2_save_2017-02-22-211421-0000': ['G5', 'MDFS',],
    'fc2_save_2017-02-22-214754-0000': ['G5', 'MSFD',],
    
    # G6
    'fc2_save_2017-04-15-180049-0000': ['G6', 'MSFS',],
    'fc2_save_2017-04-15-183537-0000': ['G6', 'MDFD',],
    'fc2_save_2017-04-15-204341-0000': ['G6', 'MSFD',],
    'fc2_save_2017-04-15-211652-0000': ['G6', 'MDFS',],
    
    # G7
    'fc2_save_2017-04-15-190919-0000': ['G7', 'MSFS',],
    'fc2_save_2017-04-15-194240-0000': ['G7', 'MDFD',],
    'fc2_save_2017-04-15-215017-0000': ['G7', 'MDFS',],
    'fc2_save_2017-04-15-225311-0000': ['G7', 'MSFD',],
}

def GetDataFromXMLs(Files, CoordFiles, Recs, Save=True):
    Data = {}
    for F, File in enumerate(Files):
        Rec = File.split('/')[-1].split('.')[0]
        Group, Rec = Recs[Rec]
        
        if Group not in Data.keys(): Data[Group] = {}
        Data[Group][Rec] = {}
        
        Tree = ElementTree.parse(File); CoordTree = ElementTree.parse(CoordFiles[F])
        Root = Tree.getroot(); CoordRoot = CoordTree.getroot()
        
        Data[Group][Rec]['Behaviours'] = {
            ''.join([b.capitalize() for b in Behav.get('type').split('_')]): 
                np.array([[Event.get('startFrame'), Event.get('endFrame')] for Event in Behav], dtype='int16')
            for Behav in Root}
        
        Data[Group][Rec]['t'] = np.array([F.get('t') for F in CoordRoot[2][0]], dtype='float32')
        Data[Group][Rec]['Coords'] = {}
        for Animal in CoordRoot[2]:
            if Animal.tag not in Data[Group][Rec]['Coords']:
                Data[Group][Rec]['Coords'][Animal.tag] = {}
            
            Keys = list(Animal[0].keys())
            Data[Group][Rec]['Coords'][Animal.tag] = {Key: np.array([F.get(Key) for F in Animal], dtype='float32')
                                                      for Key in Keys if Key != 't'}
    
    if Save:
        io.savemat('AllData-Behaviours_Coords.mat', {'Data': Data}, 
                   long_field_names=True, oned_as='column')
    
    return(Data)


Data = GetDataFromXMLs(Files, CoordFiles, Recs, Save=True)

#%% Plots

Groups = list(Data.keys()); Groups.sort(); #Groups = ['G3', 'G4', 'G5', 'G7']
Pairs = ['MSFS', 'MSFD', 'MDFS', 'MDFD']
Behs = list(Data[Groups[0]][Pairs[0]]['Behaviours'].keys())

for B in Behs:
    Fig, Ax = plt.subplots(1, 2, figsize=(10,5))
    Y = np.zeros((len(Groups), len(Pairs), 2))
    
    for P, Pair in enumerate(Pairs):
        for G, Group in enumerate(Groups):
            Y[G,P,0] = Data[Group][Pair]['Behaviours'][B].shape[0]
            Y[G,P,1] = sum([E[1]-E[0] for E in Data[Group][Pair]['Behaviours'][B]])/FPS
        
        Ax[0].plot([P+1]*len(Groups), Y[:,P,0], 'ko')
        Ax[1].plot([P+1]*len(Groups), Y[:,P,1], 'ko')
    
    Nb = Ax[0].boxplot(Y[:,:,0], showmeans=True)
    Dur = Ax[1].boxplot(Y[:,:,1], showmeans=True)
    
    for K in ['boxes', 'whiskers', 'caps', 'medians', 'fliers']:
        for I in range(len(Pairs)): 
            Nb[K][I].set(color='k')
            Dur[K][I].set(color='k')
    
    LineYNb = np.amax(Y[:,:,0]) + ((np.amax(Y[:,:,0])-np.amin(Y[:,:,0]))*0.1)
    LineYDur = np.amax(Y[:,:,1]) + ((np.amax(Y[:,:,1])-np.amin(Y[:,:,1]))*0.1)
    
    TPairs = list(combinations(range(len(Y[0,:,0])), 2))
    for TP in TPairs:
        SSNb = Y[:,TP[0],0]; DDNb = Y[:,TP[1],0]
        SSDur = Y[:,TP[0],1]; DDDur = Y[:,TP[1],1]
        if np.mean(SSNb) > np.mean(DDNb): SSNb, DDNb = DDNb, SSNb
        if np.mean(SSDur) > np.mean(DDDur): SSDur, DDDur = DDDur, SSDur
        
        pNb = DataAnalysis.Stats.RTTest(SSNb, DDNb, Confidence=1-(0.05/len(TPairs)))
        pDur = DataAnalysis.Stats.RTTest(SSDur, DDDur, Confidence=1-(0.05/len(TPairs)))
        pNb = pNb['p.value']; pDur = pDur['p.value']
        
        if pNb < 0.05: 
            LineYNb = LineYNb+(TP[1]*2)
            DataAnalysis.Plot.SignificanceBar([TP[0]+1, TP[1]+1], [LineYNb, LineYNb], str(pNb), Ax[0])
        if pDur < 0.05: 
            LineYDur = LineYDur+(TP[1]*2)
            DataAnalysis.Plot.SignificanceBar([TP[0]+1, TP[1]+1], [LineYDur, LineYDur], str(pDur), Ax[1])
    
    DataAnalysis.Plot.Set(AxesObj=Ax[0], Axes=True)
    DataAnalysis.Plot.Set(AxesObj=Ax[1], Axes=True)
    Ax[0].set_ylim([0, LineYNb+(LineYNb*0.08)])
    Ax[1].set_ylim([0, LineYDur+(LineYDur*0.08)])
    Ax[0].spines['left'].set_position(('outward', 5))
    Ax[0].spines['bottom'].set_position(('outward', 5))
    Ax[0].set_ylabel('Number of events'); Ax[0].set_xlabel('Pairs')
    Ax[1].spines['left'].set_position(('outward', 5))
    Ax[1].spines['bottom'].set_position(('outward', 5))
    Ax[1].set_ylabel('Duration of behaviour [s]'); Ax[1].set_xlabel('Pairs')
    Ax[0].set_xticks([1, 2, 3, 4]); Ax[0].set_xticklabels(Pairs)
    Ax[1].set_xticks([1, 2, 3, 4]); Ax[1].set_xticklabels(Pairs)
    Fig.suptitle(B)
    Fig.savefig('G32G7-'+B+'.pdf', format='pdf')
    Fig.savefig('G32G7-'+B+'.eps', format='eps')
plt.show()

