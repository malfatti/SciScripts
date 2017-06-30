#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: malfatti

Rat Behaviour analysis
"""
#%% Plot Vocs and Behaviours
import numpy as np
import h5py

from DataAnalysis import Stats
from DataAnalysis.DataAnalysis import Spectrogram, PSD
from DataAnalysis.Plot import Plot
from IO import Txt
from glob import glob
from itertools import combinations
from scipy import io
from xml.etree import ElementTree

Params = Plot.Set(Params=True)
from matplotlib import rcParams; rcParams.update(Params)
from matplotlib import pyplot as plt


def BoxPlots(Data, FigTitle, XLabel, YLabel, Names, FigName, LinesAmpF=2, Ext=['pdf', 'eps'], Save=True):
    if type(Data) is not type(np.array([])): Data = np.array(Data).T
    
    Fig, Ax = plt.subplots(1,1)
    BoxPlot = Ax.boxplot(Data, showmeans=True)
    
    for K in ['boxes', 'whiskers', 'caps', 'medians', 'fliers']:
        for I in range(Data.shape[1]): 
            BoxPlot[K][I].set(color='k')
            BoxPlot[K][I].set(color='k')
    
    LineY = np.amax(Data) + ((np.amax(Data)-np.amin(Data))*0.1)
    
    TPairs = list(combinations(range(Data.shape[1]), 2))
    if FigName.split('-')[1] == 'ST': TPairs = [(0,1), (2,3)]
    else: TPairs = [(0,3), (1,2)]
    
    for TP in TPairs:
        SS = Data[:, TP[0]]; DD = Data[:, TP[1]]
        if np.mean(SS) > np.mean(DD): SS, DD = DD, SS
        
        p = Stats.RTTest(SS, DD)
        p = p['p.value']*len(TPairs)
        
        if p < 0.05: 
            LineY = LineY+(TP[1]*LinesAmpF)
            Plot.SignificanceBar([TP[0]+1, TP[1]+1], [LineY, LineY], str(p), Ax)
    
    Plot.Set(AxesObj=Ax, Axes=True)
    Ax.spines['left'].set_position(('outward', 5))
    Ax.spines['bottom'].set_position(('outward', 5))
    Ax.set_ylabel(YLabel); Ax.set_xlabel(XLabel)
    Ax.set_xticks(list(range(1,Data.shape[1]+1))); Ax.set_xticklabels(Names)
    Fig.suptitle(FigTitle)
    if Save:
        for Ex in Ext:
            Fig.savefig(FigName+'.'+Ex, format=Ex)
    
    return(Fig, Ax)

def ScatterMean(Data, FigTitle, XLabel, YLabel, Names, FigName, LinesAmpF=2, Spread=0.2, Ext=['pdf', 'svg'], Save=True):
    if type(Data) is not type(np.array([])): Data = np.array(Data).T
    
    Fig, Ax = plt.subplots(1,1)
    
    for P in range(Data.shape[1]):
        X = np.random.uniform(P+1-Spread, P+1+Spread, len(Data[:,P]))
        Error = [0, np.std(Data[:,P])/len(Data[:,P])**0.5, 0]
        
        Ax.plot(X, Data[:,P], 'ko')
        Ax.errorbar([P+1-Spread, P+1, P+1+Spread], [np.mean(Data[:,P])]*3, Error, lw=3, elinewidth=1, capsize=10, color='k')
    
    Margin = ((np.amax(Data)-np.amin(Data))*0.1)
    LineY = np.amax(Data) + Margin
    
    TPairs = list(combinations(range(Data.shape[1]), 2))
    if FigName.split('-')[1] == 'ST': TPairs = [(0,1), (2,3)]
    else: TPairs = [(0,3), (1,2)]
    
    for TP in TPairs:
        SS = Data[:, TP[0]]; DD = Data[:, TP[1]]
        if np.mean(SS) > np.mean(DD): SS, DD = DD, SS
        
        p = Stats.RTTest(SS, DD)
        p = p['p.value']*len(TPairs)
        print(FigTitle, TP, round(p, 3))
        
        if p < 0.05: 
            if p < 0.001: p = 'p < 0.001'
            else: p = 'p = ' + str(round(p, 3))
            
            LineY = LineY+(TP[1]*LinesAmpF)
#            Plot.SignificanceBar([TP[0]+1, TP[1]+1], [LineY, LineY], p, Ax)
    
    Plot.Set(AxesObj=Ax, Axes=True)
    Ax.set_ylim([-Margin, LineY + LineY*0.05]); Ax.set_xlim([0, Data.shape[1]+1])
    Ax.spines['left'].set_position(('outward', 5))
    Ax.spines['bottom'].set_position(('outward', 5))
    Ax.set_ylabel(YLabel); Ax.set_xlabel(XLabel)
    Ax.set_xticks(list(range(1,Data.shape[1]+1))); Ax.set_xticklabels(Names)
    Fig.suptitle(FigTitle)
    if Save:
        for Ex in Ext:
            Fig.savefig(FigName+'.'+Ex, format=Ex)
    
    
    return(Fig, Ax)


def STSameGender(Data, FigTitle, XLabel, YLabel, Names, FigName, LinesAmpF, Ext=['pdf', 'svg'], Save=True):
    if type(Data) is not type(np.array([])): Data = np.array(Data).T
    
    Fig, Ax = plt.subplots(1,1)
    
    Margin = ((np.amax(Data)-np.amin(Data))*0.1)
    LineY = np.amax(Data) + Margin
    TPairs = [(0,1), (2,3)]
    
    for TP in TPairs:
        for G in range(Data.shape[0]):
            Ax.plot([TP[0]+1, TP[1]+1], [Data[G,TP[0]], Data[G,TP[1]]], 'ko-')
    
        SS = Data[:, TP[0]]; DD = Data[:, TP[1]]
        if np.mean(SS) > np.mean(DD): SS, DD = DD, SS
        
        p = Stats.RTTest(SS, DD)
        p = p['p.value']*len(TPairs)
        print(FigTitle, TP, round(p, 3))
        
        if p < 0.05: 
            if p < 0.001: p = 'p < 0.001'
            else: p = 'p = ' + str(round(p, 3))
            
            LineY = LineY+(TP[1]*LinesAmpF)
#            Plot.SignificanceBar([TP[0]+1, TP[1]+1], [LineY, LineY], p, Ax)
    
    Plot.Set(AxesObj=Ax, Axes=True)
    Ax.set_ylim([-Margin, LineY + LineY*0.05]); Ax.set_xlim([0, Data.shape[1]+1])
    Ax.spines['left'].set_position(('outward', 5))
    Ax.spines['bottom'].set_position(('outward', 5))
    Ax.set_ylabel(YLabel); Ax.set_xlabel(XLabel)
    Ax.set_xticks(list(range(1,Data.shape[1]+1))); Ax.set_xticklabels(Names)
    Fig.suptitle(FigTitle)
    if Save:
        for Ex in Ext:
            Fig.savefig(FigName+'.'+Ex, format=Ex)
    
    
    return(Fig, Ax)


def USVSort(File):
    Vocs = io.loadmat(File, squeeze_me=True, struct_as_record=False)
    Rate = getattr(Vocs['USVinfo'], 'soundFs')
    VocsSorted = {'22kHz': [], '50kHz': [], 'Skipped': []}
    
    for V, Voc in enumerate(Vocs['USVclip'][0]):
        if getattr(Vocs['USV'][0][V], 'q') == 'd': 
            VocsSorted['Skipped'].append(V)
            print('Skipped No', V); continue
        
        Break = False
        
        F, Pxx = PSD(Voc.sound, Rate)
        PxxSp = Pxx[:]; Pxx[F<11000] = 0
        Thr = Pxx.mean() + (3*Pxx.std())
        
        print('Voc No', str(V)+', peak at', str(F[Pxx.argmax()]/1000)+'kHz')
        if Pxx.max() > Thr:
            if F[Pxx.argmax()] > 30000: VocsSorted['50kHz'].append(V)
            else: VocsSorted['22kHz'].append(V)
        else:
            ShowVoc = True
            SF, ST, Sxx = Spectrogram(Voc.sound, Rate, 10000)
            
            while ShowVoc:
                Fig, Axes = plt.subplots(2, 1)
                Axes[0].plot(F, PxxSp, lw=2); Axes[0].plot([F[0], F[-1]], [Thr, Thr], lw=2)
                Plot.Spectrogram(Axes[1], ST, SF, Sxx, HighFreqThr=100000)
                plt.show()
                
                Ans = None
                while Ans not in range(5):
                    print('0) 22kHz')
                    print('1) 50kHz')
                    print('2) Show voc.', V, 'again')
                    print('3) Skip voc.', V)
                    print('4) Stop sorting')
                    Ans = input(': '); Ans = int(Ans)
                
                ShowVoc = False
                if Ans == 0: VocsSorted['22kHz'].append(V)
                elif Ans == 1: VocsSorted['50kHz'].append(V)
                elif Ans == 2: ShowVoc = True
                elif Ans == 3: VocsSorted['Skipped'].append(V)
                elif Ans == 4: Break = True
        
        if Break: break
    
    return(VocsSorted)


#%% Vocs
#Fields = Vocs['Vocs']._fieldnames
VocsFile = 'G32G7-Vocs.mat'
Vocs = io.loadmat(VocsFile, squeeze_me=True, struct_as_record=False)

STNames = ['STMS', 'STMD', 'STFS', 'STFD']
SINames = ['SIMSFS', 'SIMSFD', 'SIMDFS', 'SIMDFD']
OFNames = ['OFMSFS', 'OFMSFD', 'OFMDFS', 'OFMDFD']
Names = [STNames, SINames, OFNames]
ST = [getattr(Vocs['Vocs'], Pair) for Pair in STNames]
SI = [getattr(Vocs['Vocs'], Pair) for Pair in SINames]
OF = [getattr(Vocs['Vocs'], Pair) for Pair in OFNames]

Exps = [(ST,'ShortTrack'), (SI,'SexualInteraction'), (OF,'OpenField')]
for E, Exp in enumerate(Exps):
    FigName = 'G32G7-'+Names[E][0][:2]+'-Vocs'
    XLabels = [N[2:] for N in Names[E]]
    
#    Fig, Ax = BoxPlots(Exp[0], Exp[1], 'Pairs', 'Mean number of USVs per 20min', 
#                       XLabels, FigName, 20, Save=False)
    Fig, Ax = ScatterMean(Exp[0], Exp[1], 'Pairs', 'Mean number of USVs per 10min', 
                              XLabels, FigName, 20, Save=True)
        
    if Names[E][0][:2] == 'ST':
        Fig, Ax = STSameGender(Exp[0], Exp[1], 'Pairs', 'Mean number of USVs per 10min', 
                              XLabels, FigName, 20, Save=True)
plt.show()


#%% Separate 22kHz and 50kHz vocs
Files = glob('/run/media/malfatti/Malfatti1TB3/Bkp/Barbara/**/T000*.mat', recursive=True); Files.sort()
PB = [F for F in Files if 'back' in F or 'G3-' in F]
for F in PB: del(Files[Files.index(F)])

Recs = Txt.DictRead('./Data/Recs.txt')
VocsSorted = {}; Vocs = {}

Groups = ['G3', 'G4', 'G5', 'G6', 'G7']
Pairs = ['MSFS', 'MSFD', 'MDFS', 'MDFD']
Classes = ['22kHz', '50kHz']

for File in Files:
    Rec = '/'.join(File.split('/')[-2:])
    if Rec not in Recs: continue
    
    Group, Pair = Recs[Rec]
    Setup = File.split('/')[-2].split('-')[2]
    
    if Group not in ['G4', 'G5', 'G6', 'G7']: Group = 'G3'
    if Setup[-2:] == 'on': Setup = 'SI'
    elif Setup[-2:] == 'ld': Setup = 'OF'
    elif Setup[-4:] == 'back': Setup = 'OFPB'
    else: 
        print('Setup is now', Setup)
        Setup = input('What setup is this? ')
    
    G = Groups.index(Group); P = Pairs.index(Pair)
    
    if Setup not in VocsSorted: VocsSorted[Setup] = {}
    if Group not in VocsSorted[Setup]: VocsSorted[Setup][Group] = {}
    
    VocsSorted[Setup][Group][Pair] = USVSort(File)
    
    if Setup not in Vocs: Vocs[Setup] = {}
    for Class in Classes:
        if Class not in Vocs[Setup]: 
            Vocs[Setup][Class] = np.zeros((len(Pairs), len(Groups)))
        
        Vocs[Setup][Class][P,G] = len(VocsSorted[Setup][Group][Pair][Class])
    
io.savemat('VocsSorted.mat', {'VocsSorted': VocsSorted},
#                                             'Vocs': Vocs}, 
            long_field_names=True, oned_as='column')

with h5py.File('VocsSorted.hdf5') as F:
    for S, Setup in Vocs.items():
        F.create_group('/Vocs/'+S)
        
        for C, Class in Setup.items():
            F['Vocs'][S][C] = Class
        
    for S, Setup in VocsSorted.items():    
        for G, Group in Setup.items():
            for P, Pair in Group.items():
                F.create_group('/VocsSorted/'+S+'/'+G+'/'+P)
                
                for C, Class in Pair.items():
                    F['VocsSorted'][S][G][P][C] = Class

with h5py.File('VocsSorted.hdf5', 'r') as F:
    a = {}
    for S, Setup in F['Vocs'].items():
        a[S] = {}
        for C, Class in Setup.items():
            a[S][C] = Class[:]
    
    b = {}
    for S, Setup in F['VocsSorted'].items():
        b[S] = {}
        for G, Group in Setup.items():
            b[S][G] = {}
            for P, Pair in Group.items():
                b[S][G][P] = {}
                
                for C, Class in Pair.items():
                    b[S][G][P][C] = Class[:]


#%% Speaker
import matplotlib.pyplot as plt

def on_click(event):
    if event.inaxes is not None:
        print(event.xdata, event.ydata)
    else:
        print('Clicked ouside axes bounds but inside plot window')

fig, ax = plt.subplots()
fig.canvas.callbacks.connect('button_press_event', on_click)
plt.show()


#%% Behaviours
BehavioursFile = 'G32G7-Behaviours.mat'
Behaviours = io.loadmat(BehavioursFile, squeeze_me=True, struct_as_record=False)

Keys = list(Behaviours.keys())
for K in Keys:
    if K[0] == '_': del(Behaviours[K])

Groups = Behaviours['Cols']; del(Behaviours['Cols'])
Pairs = Behaviours['Lines']; del(Behaviours['Lines'])

for B, Behaviour in Behaviours.items():
    if B[:2] == 'OF': 
        if B[:4] == 'OFPB': FigName = 'G32G7-OFPB-'+B
        else: FigName = 'G32G7-OF-'+B
    else: FigName = 'G32G7-SI-'+B
    
    if B == 'Darts': B = 'Dart/Hop'
#    Fig, Ax = BoxPlots(Behaviour.T.tolist(), B, 'Pairs', 'Number of events', Pairs, FigName, 2, Save=True)
    Fig, Ax = ScatterMean(Behaviour, B, 'Pairs', 'Number of events per 20min', Pairs, FigName, 200, Save=True)
plt.show()
    

#%% Icy Mice Profiler analysis

## XMLs to Dict
Recs = Txt.DictRead('./Data/Recs.txt')
Files = glob('./Data/*/*allEventsExported.xml'); Files.sort()
#Files = glob('./InvertedXMLs/*allEventsExported.xml'); Files.sort()
#Files = glob('./G*/*allEventsExported.xml'); Files.sort()

CoordFiles = ['.'.join(F.split('.')[:-2]) for F in Files]
FPS = 23

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


Data = GetDataFromXMLs(Files, CoordFiles, Recs, Save=False)

## Plots

Groups = list(Data.keys()); Groups.sort(); #Groups = ['G3', 'G4', 'G5', 'G7']
Pairs = ['MSFS', 'MSFD', 'MDFS', 'MDFD']
Behs = list(Data[Groups[0]][Pairs[0]]['Behaviours'].keys())
Spread = 0.2

for B in Behs:
    Fig, Ax = plt.subplots(1, 2, figsize=(15,5))
    Y = np.zeros((len(Groups), len(Pairs), 2))
    
    for P, Pair in enumerate(Pairs):
        for G, Group in enumerate(Groups):
            Y[G,P,0] = Data[Group][Pair]['Behaviours'][B].shape[0]
            Y[G,P,1] = sum([E[1]-E[0] for E in Data[Group][Pair]['Behaviours'][B]])/FPS
        
        XNb = np.random.uniform(P+1-Spread, P+1+Spread, len(Y[:,P,0]))
        XDur = np.random.uniform(P+1-Spread, P+1+Spread, len(Y[:,P,1]))
        ErrorNb = [0, np.std(Y[:,P,0])/len(Y[:,P,0])**0.5, 0]
        ErrorDur = [0, np.std(Y[:,P,1])/len(Y[:,P,1])**0.5, 0]
        
        Ax[0].plot(XNb, Y[:,P,0], 'ko'); Ax[1].plot(XDur, Y[:,P,1], 'ko')
        Ax[0].errorbar([P+1-Spread, P+1, P+1+Spread], [np.mean(Y[:,P,0])]*3, ErrorNb, lw=3, elinewidth=1, capsize=10, color='k')
        Ax[1].errorbar([P+1-Spread, P+1, P+1+Spread], [np.mean(Y[:,P,1])]*3, ErrorDur, lw=3, elinewidth=1, capsize=10, color='k')
    
    MarginNb = ((np.amax(Y[:,:,0])-np.amin(Y[:,:,0]))*0.1)
    MarginDur = ((np.amax(Y[:,:,1])-np.amin(Y[:,:,1]))*0.1)
    LineYNb = np.amax(Y[:,:,0]) + MarginNb
    LineYDur = np.amax(Y[:,:,1]) + MarginDur
    
    TPairs = list(combinations(range(len(Y[0,:,0])), 2))
    TPairs = [(0,3), (1,2)]
    for TP in TPairs:
        SSNb = Y[:,TP[0],0]; DDNb = Y[:,TP[1],0]
        SSDur = Y[:,TP[0],1]; DDDur = Y[:,TP[1],1]
        if np.mean(SSNb) > np.mean(DDNb): SSNb, DDNb = DDNb, SSNb
        if np.mean(SSDur) > np.mean(DDDur): SSDur, DDDur = DDDur, SSDur
        
        pNb = Stats.RTTest(SSNb, DDNb)
        pDur = Stats.RTTest(SSDur, DDDur)
        pNb = pNb['p.value']*len(TPairs); pDur = pDur['p.value']*len(TPairs)
        
        if pNb < 0.05: 
            if pNb < 0.001: pNb = 'p < 0.001'
            else: pNb = 'p = ' + str(round(pNb, 3))
            
            LineYNb = LineYNb+(TP[1]*2)
#            Plot.SignificanceBar([TP[0]+1, TP[1]+1], [LineYNb, LineYNb], pNb, Ax[0])
        if pDur < 0.05: 
            if pDur < 0.001: pNb = 'p < 0.001'
            else: pDur = 'p = ' + str(round(pDur, 3))
            
            LineYDur = LineYDur+(TP[1]*2)
#            Plot.SignificanceBar([TP[0]+1, TP[1]+1], [LineYDur, LineYDur], pDur, Ax[1])
    
    Plot.Set(AxesObj=Ax[0], Axes=True)
    Plot.Set(AxesObj=Ax[1], Axes=True)
    
    Ax[0].set_ylim([-MarginNb, LineYNb+(LineYNb*0.08)])
    Ax[1].set_ylim([-MarginDur, LineYDur+(LineYDur*0.08)])
    Ax[0].set_xlim([0, len(Pairs)+1]); Ax[1].set_xlim([0, len(Pairs)+1])
    
    Ax[0].spines['left'].set_position(('outward', 5))
    Ax[0].spines['bottom'].set_position(('outward', 5))
    Ax[0].set_ylabel('Number of events per 20min'); Ax[0].set_xlabel('Pairs')
    
    Ax[1].spines['left'].set_position(('outward', 5))
    Ax[1].spines['bottom'].set_position(('outward', 5))
    Ax[1].set_ylabel('Duration of behaviour [s]'); Ax[1].set_xlabel('Pairs')
    
    Ax[0].set_xticks([1, 2, 3, 4]); Ax[0].set_xticklabels(Pairs)
    Ax[1].set_xticks([1, 2, 3, 4]); Ax[1].set_xticklabels(Pairs)
    
    Fig.suptitle(B)
    Fig.savefig('G32G7-OF-'+B+'.pdf', format='pdf')
    Fig.savefig('G32G7-OF-'+B+'.svg', format='svg')
#    Fig.suptitle('BFollowA')
#    Fig.savefig('G32G7-OF-'+'BFollowA'+'.pdf', format='pdf')
#    Fig.savefig('G32G7-OF-'+'BFollowA'+'.svg', format='svg')
plt.show()

