#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: malfatti

Rat Behaviour analysis
"""
#%% Plot Vocs and Behaviours
import numpy as np
import h5py

from copy import deepcopy
from DataAnalysis import Stats
from DataAnalysis.DataAnalysis import Spectrogram, PSD
from DataAnalysis.Plot import Plot
from IO import Txt
from glob import glob
from itertools import combinations
from PIL import Image
from scipy import io
from xml.etree import ElementTree

Params = Plot.Set(Params=True)
from matplotlib import rcParams; rcParams.update(Params)
from matplotlib import pyplot as plt


def ClickCallback(Event):
    """ Taken and modified from Joe Kington's code.
        Source at https://stackoverflow.com/a/27566966 """
    
    global X, Y
    while X == Y == None:
        print('Click on the speaker position.')
        if Event.inaxes is not None:
            print(Event.xdata, Event.ydata)
            X, Y = Event.inaxes.transData.inverted().transform((Event.x, Event.y))
        else:
            print('Clicked coordinate is ouside axes limits')
    
    plt.close()
    return(None)


def GetDataFromXMLs(Files, CoordFiles, Recs, Behs=True, Coords=True, Save=True, SaveFile=''):
    Data = {}
    for F, File in enumerate(Files):
        Rec = File.split('/')[-1].split('.')[0]
        Group, Rec = Recs[Rec]
        
        if Group not in Data.keys(): Data[Group] = {}
        Data[Group][Rec] = {}
        
        Tree = ElementTree.parse(File); CoordTree = ElementTree.parse(CoordFiles[F])
        Root = Tree.getroot(); CoordRoot = CoordTree.getroot()
        
        if Behs:
            Data[Group][Rec]['Behaviours'] = {
                ''.join([b.capitalize() for b in Behav.get('type').split('_')]): 
                    np.array([[Event.get('startFrame'), Event.get('endFrame')] for Event in Behav], dtype='int16')
                for Behav in Root}
        
        if Coords:
            Data[Group][Rec]['Coords'] = {}
            for Animal in CoordRoot[2]:
                if Animal.tag not in Data[Group][Rec]['Coords']:
                    Data[Group][Rec]['Coords'][Animal.tag] = {}
                
                Keys = list(Animal[0].keys())
                Data[Group][Rec]['Coords'][Animal.tag] = {Key: np.array([F.get(Key) for F in Animal], dtype='float32')
                                                          for Key in Keys if Key != 't'}
        
        Data[Group][Rec]['t'] = np.array([F.get('t') for F in CoordRoot[2][0]], dtype='float32')
    
    if Save:
        io.savemat(SaveFile, {'Data': Data}, 
                   long_field_names=True, oned_as='column')
    
    return(Data)


def USVSort(File):
    Vocs = io.loadmat(File, squeeze_me=True, struct_as_record=False)
    Rate = getattr(Vocs['USVinfo'], 'soundFs')
    VocsSorted = {'22kHz': [], '50kHz': [], 'Skipped': []}
    
    if type(Vocs['USVclip'][0]) is not io.matlab.mio5_params.mat_struct:
        Vocs['USVclip'] = Vocs['USVclip'][0]
        Vocs['USV'] = Vocs['USV'][0]
    
    for V, Voc in enumerate(Vocs['USVclip']):
        if getattr(Vocs['USV'][V], 'q') == 'd': 
            VocsSorted['Skipped'].append(V)
            print('Skipped No', V); continue
        
        Break = False
        
        F, Pxx = PSD(Voc.sound, Rate)
        PxxSp = Pxx[:]; Pxx[F<12000] = 0
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


def USVSortedWrite(FileName, Vocs, VocsSorted):
    with h5py.File(FileName) as F:
        for S, Setup in Vocs.items():
            Path = '/Vocs/'+S
            if Path not in F: F.create_group(Path)
            
            for C, Class in Setup.items():
                F['Vocs'][S][C] = Class
            
        for S, Setup in VocsSorted.items():    
            for G, Group in Setup.items():
                for P, Pair in Group.items():
                    Path = '/VocsSorted/'+S+'/'+G+'/'+P
                    if Path not in F: F.create_group(Path)
                    
                    for C, Class in Pair.items():
                        F['VocsSorted'][S][G][P][C] = Class
    
    return(None)


def USVSortedRead(FileName):
    with h5py.File(FileName, 'r') as F:
        Vocs = {}
        for S, Setup in F['Vocs'].items():
            Vocs[S] = {}
            for C, Class in Setup.items():
                Vocs[S][C] = Class[:]
        
        VocsSorted = {}
        for S, Setup in F['VocsSorted'].items():
            VocsSorted[S] = {}
            for G, Group in Setup.items():
                VocsSorted[S][G] = {}
                for P, Pair in Group.items():
                    VocsSorted[S][G][P] = {}
                    
                    for C, Class in Pair.items():
                        VocsSorted[S][G][P][C] = Class[:]
    
    return(Vocs, VocsSorted)


def DictFixKeys(Dict, Copy=True):
    if Copy: DictCopy = deepcopy(Dict)
    else: DictCopy = Dict
    
    for K,Key in DictCopy.items():
        if type(Key) is dict: DictCopy[K] = DictFixKeys(Key)
        else: 
            if K[0].isdigit(): 
                Kk = 'x'+K
                DictCopy[Kk] = Key
                del(DictCopy[K])
    
    return(DictCopy)


def DictUndoFixKeys(Dict, Copy=True):
    if Copy: DictCopy = deepcopy(Dict)
    else: DictCopy = Dict
    
    for K,Key in DictCopy.items():
        if type(Key) is dict: DictCopy[K] = DictUndoFixKeys(Key)
        else: 
            if K[0] == 'x': 
                Kk = K[1:]
                DictCopy[Kk] = Key
                del(DictCopy[K])
    
    return(DictCopy)


def RecursiveDict(Dict, Copy=True):
    if Copy: DictCopy = deepcopy(Dict)
    else: DictCopy = Dict
    
    for K,Key in DictCopy.items():
        if type(Key) is dict: DictCopy[K] = RecursiveDict(Key)
        elif type(Key) is list: pass#DictCopy[K] = np.array(Key, dtype='float32')
        elif type(Key) is np.ndarray: pass#DictCopy[K] = Key.tolist()
    
    return(DictCopy)


## Plots

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
    
#    TPairs = list(combinations(range(Data.shape[1]), 2))
#    if FigName.split('-')[1] == 'ST': TPairs = [(0,1), (2,3)]
#    else: TPairs = [(0,3), (1,2)]
#    
#    for TP in TPairs:
#        SS = Data[:, TP[0]]; DD = Data[:, TP[1]]
#        if np.mean(SS) > np.mean(DD): SS, DD = DD, SS
#        
#        p = Stats.RTTest(SS, DD)
#        p = p['p.value']*len(TPairs)
#        print(FigTitle, TP, round(p, 3))
#        
#        if p < 0.05: 
#            if p < 0.001: p = 'p < 0.001'
#            else: p = 'p = ' + str(round(p, 3))
#            
#            LineY = LineY+(TP[1]*LinesAmpF)
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


def ScatterMeanSorted(Data, FigTitle, XLabel, YLabel, Names, FigName, LinesAmpF=2, Spread=0.2, Ext=['pdf', 'svg'], Save=True):
    Fig, Ax = plt.subplots(1,1)
    Colors = ['m', 'y']; XTicks = []
    for P in range(len(Data)):
        I = [P*3+1, P*3+2]; XTicks.append(sum(I)/2)
        
        for C in range(len(Data[P])):
            X = np.random.uniform(I[C]-Spread, I[C]+Spread, len(Data[P][C]))
            Error = [0, np.std(Data[P][C])/len(Data[P][C])**0.5, 0]
            
            Ax.plot(X, Data[P][C], Colors[C]+'o')
            Ax.errorbar([I[C]-Spread, I[C], I[C]+Spread], [np.mean(Data[P][C])]*3, Error, lw=3, elinewidth=1, capsize=10, color='k')
    
    Margin = ((np.amax(Data)-np.amin(Data))*0.1)
    LineY = np.amax(Data) + Margin
    
    Plot.Set(AxesObj=Ax, Axes=True)
    Ax.set_ylim([-Margin, LineY + LineY*0.05]); Ax.set_xlim([0, XTicks[-1]+1.5])
    Ax.spines['left'].set_position(('outward', 5))
    Ax.spines['bottom'].set_position(('outward', 5))
    Ax.set_ylabel(YLabel); Ax.set_xlabel(XLabel)
    Ax.set_xticks(XTicks); Ax.set_xticklabels(Names)
    Ax.legend(['22kHz', '50kHz'], loc='best')
    
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
#    TPairs = [(0,1), (2,3)]
#    
#    for TP in TPairs:
#        for G in range(Data.shape[0]):
#            Ax.plot([TP[0]+1, TP[1]+1], [Data[G,TP[0]], Data[G,TP[1]]], 'ko-')
#    
#        SS = Data[:, TP[0]]; DD = Data[:, TP[1]]
#        if np.mean(SS) > np.mean(DD): SS, DD = DD, SS
#        
#        p = Stats.RTTest(SS, DD)
#        p = p['p.value']*len(TPairs)
#        print(FigTitle, TP, round(p, 3))
#        
#        if p < 0.05: 
#            if p < 0.001: p = 'p < 0.001'
#            else: p = 'p = ' + str(round(p, 3))
#            
#            LineY = LineY+(TP[1]*LinesAmpF)
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


#%% Separate 22kHz and 50kHz vocs
Files = glob('/run/media/malfatti/Malfatti1TB3/Bkp/Barbara/**/T000*.mat', recursive=True); Files.sort()
PB = [F for F in Files if 'back' in F or 'G3-OF' in F]
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
    try: Setup = File.split('/')[-2].split('-')[2]
    except IndexError: Setup = 'on'
    
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


io.savemat('Analysis/VocsSorted.mat', {'VocsSorted': DictFixKeys(VocsSorted), 
                                       'Vocs': DictFixKeys(Vocs), 
                                       'Cols': Groups, 'Lines': Pairs}, 
            long_field_names=True, oned_as='column')


#%% Vocs
#Fields = Vocs['Vocs']._fieldnames
VocsFile = 'Analysis/G32G7-Vocs.mat'
AllVocs = io.loadmat(VocsFile, squeeze_me=True, struct_as_record=False)
Vocs, VocsSorted = USVSortedRead('Analysis/VocsSorted.hdf5')

Pairs = ['MSFS', 'MSFD', 'MDFS', 'MDFD']
STNames = ['STMS', 'STMD', 'STFS', 'STFD']
SINames = ['SIMSFS', 'SIMSFD', 'SIMDFS', 'SIMDFD']
OFNames = ['OFMSFS', 'OFMSFD', 'OFMDFS', 'OFMDFD']
Names = [STNames, SINames, OFNames]
ST = [getattr(AllVocs['Vocs'], Pair) for Pair in STNames]
SI = [getattr(AllVocs['Vocs'], Pair) for Pair in SINames]
OF = [getattr(AllVocs['Vocs'], Pair) for Pair in OFNames]

Exps = [(ST,'ShortTrack'), (SI,'SexualInteraction'), (OF,'OpenField')]
for E, Exp in enumerate(Exps):
    FigName = 'G32G7-'+Names[E][0][:2]+'-Vocs'
    XLabels = [N[2:] for N in Names[E]]
    
#    Fig, Ax = BoxPlots(Exp[0], Exp[1], 'Pairs', 'Mean number of USVs per 20min', 
#                       XLabels, FigName, 20, Save=False)
    Fig, Ax = ScatterMean(Exp[0], '', 'Pairs', 'Number of USVs per 20min', 
                              XLabels, FigName, 20, Save=True)
        
    if Names[E][0][:2] == 'ST':
        Fig, Ax = ScatterMean(Exp[0], '', 'Pairs', 'Mean number of USVs per 10min', 
                              XLabels, FigName, 20, Save=True)
plt.show()

for S, Setup in Vocs.items():
    Data = [[Setup['22kHz'][A], Setup['50kHz'][A]] for A in range(len(Pairs))]
    FigName = 'G32G7-'+S+'-VocsSorted'
    
    ScatterMeanSorted(Data, '', 'Pairs', 'Number of USVs per 20min', Pairs, FigName, 20, Save=True)
plt.show()


#%% Speaker
Px2M = 0.43/128
r = 0.50

Files = glob('./Tracking/OFPB/G*/*allEventsExported.xml'); Files.sort()
CoordFiles = ['.'.join(F.split('.')[:-2]) for F in Files]
Recs = Txt.DictRead('./Data/Recs.txt')

r = round(r/Px2M)

## Find speaker coordinate
#X, Y = None, None
#Frame = Image.open('./Data/OFPB.bmp')
#
#Fig, Ax = plt.subplots(figsize=(_/200 for _ in Frame.size), dpi=150)
#Ax.imshow(Frame, vmin=0, vmax=(2**24)-1)
#Ax.tick_params(bottom='off', top='off', labelbottom='off',
#               right='off', left='off', labelleft='off')
#
#Fig.tight_layout(pad=0)
#Fig.suptitle('Click on the speaker position')
#Fig.canvas.callbacks.connect('button_press_event', ClickCallback)
#plt.show()
X, Y = (532.5, 1014.75)

Data = GetDataFromXMLs(Files, CoordFiles, Recs, Save=False)

for G, Group in Data.items():
    for P, Pair in Group.items():
        for A, Animal in Pair['Coords'].items():
            Coords = np.vstack((Animal['bodyx'], Animal['bodyy'])).T
            Events = np.zeros(Coords.shape[0])
            Events[np.where(((Coords[:,0] - X)**2 + (Coords[:,1] - Y)**2)**0.5 < r)] = 1
            if Events[0] == 1: Events[0] = 0; Events[1] = 1
            
            Events = np.diff(Events)
            Events = np.hstack((np.argwhere(Events == 1), np.argwhere(Events == -1)))
            
            Data[G][P]['Behaviours'][A[-1]+'GoToSpeaker'] = Events

Groups = list(Data.keys()); Groups.sort(); #Groups = ['G3', 'G4', 'G5', 'G7']
Pairs = ['MSFS', 'MSFD', 'MDFS', 'MDFD']
Spread = 0.2
FPS = 23

Fig, Ax = plt.subplots(1, 2, figsize=(20,5))
YA, YB = np.zeros((len(Groups), len(Pairs), 2)), np.zeros((len(Groups), len(Pairs), 2))

Colors = ['b', 'g']; XTicks = []
for P, Pair in enumerate(Pairs):
    I = [P*3+1, P*3+2]; XTicks.append(sum(I)/2); 
    
    for G, Group in enumerate(Groups):
        YA[G,P,0] = Data[Group][Pair]['Behaviours']['AGoToSpeaker'].shape[0]
        YA[G,P,1] = sum([E[1]-E[0] for E in Data[Group][Pair]['Behaviours']['AGoToSpeaker']])/FPS
        YB[G,P,0] = Data[Group][Pair]['Behaviours']['BGoToSpeaker'].shape[0]
        YB[G,P,1] = sum([E[1]-E[0] for E in Data[Group][Pair]['Behaviours']['BGoToSpeaker']])/FPS
    
    XANb = np.random.uniform(I[0]-Spread, I[0]+Spread, len(YA[:,P,0]))
    XADur = np.random.uniform(I[0]-Spread, I[0]+Spread, len(YA[:,P,1]))
    XBNb = np.random.uniform(I[1]-Spread, I[1]+Spread, len(YA[:,P,0]))
    XBDur = np.random.uniform(I[1]-Spread, I[1]+Spread, len(YA[:,P,1]))
    
    ErrorANb = [0, np.std(YA[:,P,0])/len(YA[:,P,0])**0.5, 0]
    ErrorADur = [0, np.std(YA[:,P,1])/len(YA[:,P,1])**0.5, 0]
    ErrorBNb = [0, np.std(YB[:,P,0])/len(YB[:,P,0])**0.5, 0]
    ErrorBDur = [0, np.std(YB[:,P,1])/len(YB[:,P,1])**0.5, 0]
    
    Ax[0].plot(XANb, YA[:,P,0], Colors[0]+'o')
    Ax[0].plot(XBNb, YB[:,P,0], Colors[1]+'o')
    Ax[1].plot(XADur, YA[:,P,1], Colors[0]+'o')
    Ax[1].plot(XBDur, YB[:,P,1], Colors[1]+'o')
    
    Ax[0].errorbar([I[0]-Spread, I[0], I[0]+Spread], [np.mean(YA[:,P,0])]*3, ErrorANb, lw=3, elinewidth=1, capsize=10, color='k')
    Ax[1].errorbar([I[0]-Spread, I[0], I[0]+Spread], [np.mean(YA[:,P,1])]*3, ErrorADur, lw=3, elinewidth=1, capsize=10, color='k')
    Ax[0].errorbar([I[1]-Spread, I[1], I[1]+Spread], [np.mean(YB[:,P,0])]*3, ErrorBNb, lw=3, elinewidth=1, capsize=10, color='k')
    Ax[1].errorbar([I[1]-Spread, I[1], I[1]+Spread], [np.mean(YB[:,P,1])]*3, ErrorBDur, lw=3, elinewidth=1, capsize=10, color='k')

MarginNb = max([((np.amax(YA[:,:,0])-np.amin(YA[:,:,0]))*0.1), ((np.amax(YB[:,:,0])-np.amin(YB[:,:,0]))*0.1)])
MarginDur = max([((np.amax(YA[:,:,1])-np.amin(YA[:,:,1]))*0.1), ((np.amax(YB[:,:,1])-np.amin(YB[:,:,1]))*0.1)])
LineYNb = max([np.amax(YA[:,:,0]), np.amax(YB[:,:,0])]) + MarginNb
LineYDur = max([np.amax(YA[:,:,1]), np.amax(YB[:,:,1])]) + MarginDur

Plot.Set(AxesObj=Ax[0], Axes=True)
Plot.Set(AxesObj=Ax[1], Axes=True)
for i in range(2): Ax[i].legend(['Male', 'Female'], loc='best')

Ax[0].set_ylim([-MarginNb, LineYNb+(LineYNb*0.08)])
Ax[1].set_ylim([-MarginDur, LineYDur+(LineYDur*0.08)])
Ax[0].set_xlim([0, XTicks[-1]+1.5]); Ax[1].set_xlim([0, XTicks[-1]+1.5])

Ax[0].spines['left'].set_position(('outward', 5))
Ax[0].spines['bottom'].set_position(('outward', 5))
Ax[0].set_ylabel('Number of events per 20min'); Ax[0].set_xlabel('Pairs')

Ax[1].spines['left'].set_position(('outward', 5))
Ax[1].spines['bottom'].set_position(('outward', 5))
Ax[1].set_ylabel('Duration of behaviour [s]'); Ax[1].set_xlabel('Pairs')

for i in range(2): Ax[i].set_xticks(XTicks); Ax[i].set_xticklabels(Pairs)

#print('G32G7-OF-'+B)
Fig.savefig('G32G7-OFPB-AAndBGoToSpeaker.pdf', format='pdf')
Fig.savefig('G32G7-OFPB-AAndBGoToSpeaker.svg', format='svg')
#    Fig.savefig('G32G7-OF-'+'BFollowA'+'.pdf', format='pdf')
#    Fig.savefig('G32G7-OF-'+'BFollowA'+'.svg', format='svg')
plt.show()


#%% Behaviours
BehavioursFile = 'Analysis/G32G7-BehavioursPB.mat'
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
    Fig, Ax = ScatterMean(Behaviour, '', 'Pairs', 'Number of events per 20min', Pairs, FigName, 200, Save=True)
plt.show()
    

#%% Icy Mice Profiler analysis

## XMLs to Dict
Recs = Txt.DictRead('./Data/Recs.txt')
#Files = glob('./Tracking/OF/*/*allEventsExported.xml'); Files.sort()
#Files = glob('./Tracking/iOF/*allEventsExported.xml'); Files.sort()
#Files = glob('./Tracking/OFPB/G*/*allEventsExported.xml'); Files.sort()
Files = glob('./Tracking/iOFPB/**/*allEventsExported.xml', recursive=True); Files.sort()
CoordFiles = ['.'.join(F.split('.')[:-2]) for F in Files]

Data = GetDataFromXMLs(Files, CoordFiles, Recs, Coords=True, Save=False)

#for G, Group in Data.items():
#    for P, Pair in Group.items():
#        Data[G][P]['Behaviours']['BFollowA'] = iData[G][P]['Behaviours']['AFollowB']

## Plots
Groups = list(Data.keys()); Groups.sort(); #Groups = ['G3', 'G4', 'G5', 'G7']
Pairs = ['MSFS', 'MSFD', 'MDFS', 'MDFD']
Behs = list(Data[Groups[0]][Pairs[0]]['Behaviours'].keys())
Spread = 0.2
FPS = 23

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
    
#    TPairs = list(combinations(range(len(Y[0,:,0])), 2))
#    TPairs = [(0,3), (1,2)]
#    for TP in TPairs:
#        SSNb = Y[:,TP[0],0]; DDNb = Y[:,TP[1],0]
#        SSDur = Y[:,TP[0],1]; DDDur = Y[:,TP[1],1]
#        if np.mean(SSNb) > np.mean(DDNb): SSNb, DDNb = DDNb, SSNb
#        if np.mean(SSDur) > np.mean(DDDur): SSDur, DDDur = DDDur, SSDur
#        
#        pNb = Stats.RTTest(SSNb, DDNb)
#        pDur = Stats.RTTest(SSDur, DDDur)
#        pNb = pNb['p.value']*len(TPairs); pDur = pDur['p.value']*len(TPairs)
#        
#        if pNb < 0.05: 
#            if pNb < 0.001: pNb = 'p < 0.001'
#            else: pNb = 'p = ' + str(round(pNb, 3))
#            
#            LineYNb = LineYNb+(TP[1]*2)
##            Plot.SignificanceBar([TP[0]+1, TP[1]+1], [LineYNb, LineYNb], pNb, Ax[0])
#        if pDur < 0.05: 
#            if pDur < 0.001: pNb = 'p < 0.001'
#            else: pDur = 'p = ' + str(round(pDur, 3))
#            
#            LineYDur = LineYDur+(TP[1]*2)
##            Plot.SignificanceBar([TP[0]+1, TP[1]+1], [LineYDur, LineYDur], pDur, Ax[1])
    
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
    
#    Fig.savefig('G32G7-OFPB-'+B+'.pdf', format='pdf')
#    Fig.savefig('G32G7-OFPB-'+B+'.svg', format='svg')
    Fig.savefig('G32G7-OFPB-'+'BFollowA'+'.pdf', format='pdf')
    Fig.savefig('G32G7-OFPB-'+'BFollowA'+'.svg', format='svg')
plt.show()

