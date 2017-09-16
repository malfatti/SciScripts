#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
GPIAZon group analysis
"""
import numpy as np
import os

from DataAnalysis import Stats
from DataAnalysis.Plot import Plot
from IO import Hdf5
from itertools import combinations

Params = Plot.Set(Params=True)
from matplotlib import rcParams; rcParams.update(Params)
from matplotlib import pyplot as plt


Colors = ['k', 'r', 'b', 'm', 'g', '#ffa500', '#00b2b2']

## Level 0
def ClearPairs(Pairs):
    ToDelete = []; EmptyPairs = []
    for P, Pair in Pairs.items():
        for F, Freq in Pair.items():
            if Freq['p.value'] > 0.05: ToDelete.append([P, F])
            if Freq['p.value'] != Freq['p.value']: ToDelete.append([P, F])
    
    for KeyPair in ToDelete: del(Pairs[KeyPair[0]][KeyPair[1]])
    
    for P, Pair in Pairs.items(): 
        if len(Pair) == 0: EmptyPairs.append(P)
    
    for KeyPair in EmptyPairs: del(Pairs[KeyPair])
    
    return(Pairs)


def GetIndex(Group, Animals, Exps, ExpList, AnalysisFile):
    Index = {}
    for A, Animal in enumerate(Animals):
        Index[Animal] = {}
        
        Done = []
        for Exp in Exps:
            AnalysisKey = [AK for AK in Hdf5.GetGroupKeys('/', AnalysisFile) 
                           if Animal in AK 
                           and Exp.split('/')[-1][:8] in AK]
            
            if AnalysisKey: AnalysisKey = AnalysisKey[0]
            else: continue
            
#            print(Animal, Exp)
            Ei = len(Done); Done.append(Exp)
            GPIAS = Hdf5.DataLoad(AnalysisKey, AnalysisFile)[0]
            Index[Animal][ExpList[Ei]] = {}
#            Freqs = list(GPIAS['Index'].keys())
#            Freqs.sort(key=lambda x: [int(y) for y in x.split('-')])
            
            for F, Freq in GPIAS[list(GPIAS.keys())[0]]['GPIAS']['Index'].items(): 
                if F == '9000-11000': continue
                Index[Animal][ExpList[Ei]][F] = abs(Freq['GPIASIndex'])
                if Index[Animal][ExpList[Ei]][F] > 1:
                    Index[Animal][ExpList[Ei]][F] = 1
            
            del(GPIAS)
    
    return(Index)


#STD = []
#for G, Group in enumerate(Groups):
#    Animals = [Group.split('_')[-1] + 'n0' + str(_) for _ in range(1,6)]
#    
#    for A, Animal in enumerate(Animals):
#        Freqs = list(Index[Group][Animal]['NaCl'].keys())
#        Freqs.sort(key=lambda x: [int(y) for y in x.split('-')])
#        Freqs = [Freqs[0]] + Freqs[2:]# + [Freqs[1]]
#        
#        STD.append([abs(Index[Group][Animal]['NaCl'][Freq]) for Freq in Freqs])
#
#MuSigma = []
#for F in range(len(Freqs)):
#    List = [A[F] for A in STD]
#    MuSigma.append([np.mean(List), np.std(List)])
#
#
#Sim = np.random.normal(MuSigma[0][0], MuSigma[0][1], 10)
#_, SimBins, _ = plt.hist(Sim, 10, normed=True)
#Gauss = 1/(MuSigma[0][1] * (2 * np.pi)**0.5) * np.exp(-(SimBins - MuSigma[0][0])**2 / (2 * MuSigma[0][1]**2) )
#plt.plot(SimBins, Gauss, linewidth=2, color='r'); plt.show()


## Level 1
def GetFreqPerExp(Group, Animals, Exps, ExpList, AnalysisFile):
    Index = GetIndex(Group, Animals, Exps, ExpList, AnalysisFile)
    
    FreqPerExp = {}
    for A, Animal in Index.items():
        for E, Exp in Animal.items():
            if E not in FreqPerExp: FreqPerExp[E] = {}
            
            for F, Freq in Exp.items():
                if F not in FreqPerExp[E]: FreqPerExp[E][F] = []
                
                FreqPerExp[E][F].append(abs(Freq))
    
#    for E, Exp in FreqPerExp.items():
#        for F, Freq in Exp.items():
#            if len(Freq) < len(Index.keys()): 
#                Miss = len(Index.keys()) - len(Freq)
#                FreqPerExp[E][F] = FreqPerExp[E][F] + [np.mean(FreqPerExp[E][F])]*Miss
    
    
    return(FreqPerExp)


def GetFreqPerAnimal(Group, Animals, Exps, ExpList, AnalysisFile, Invalid=False, InvalidThr=0.3):
    Index = GetIndex(Group, Animals, Exps, ExpList, AnalysisFile)
    FreqPerAnimal = {}
    
    for A, Animal in Index.items():
        if A not in FreqPerAnimal: FreqPerAnimal[A] = []
        
        for E, Exp in enumerate(ExpList):
            if Exp not in Animal: continue
            
            Freqs = list(Animal[Exp].keys())
            Freqs.sort(key=lambda x: [int(y) for y in x.split('-')])
#            Freqs = [Freqs[0]] + Freqs[2:]# + [Freqs[1]]
            
            Y = [abs(Animal[Exp][Freq]) for Freq in Freqs]
            FreqPerAnimal[A].append(np.array(Y))
        
        if not Invalid: 
            for F, Freq in enumerate(FreqPerAnimal[A]):
                inv = np.where((Freq < InvalidThr)+
                               (1 < Freq))[0]
    #            Diff[A][0] = np.delete(Diff[A][0], inv)
    #            Diff[A][1] = np.delete(Diff[A][1], inv)
                FreqPerAnimal[A][F][inv] = -1
    
    return(FreqPerAnimal)


## Level 2
def GetValid(Group, Animals, Exps, ExpList, AnalysisFile, Invalid=True, InvalidThr=0.3):
    Index = GetIndex(Group, Animals, Exps, ExpList, AnalysisFile)
    Valid = GetFreqPerAnimal(Group, Animals, Exps, ExpList, AnalysisFile, Invalid, InvalidThr)
    
    for A, Animal in Index.items():
        for E, Exp in enumerate(ExpList):
            if Exp not in Animal: continue
            
            Freqs = list(Animal[Exp].keys())
            Freqs.sort(key=lambda x: [int(y) for y in x.split('-')])
#            Freqs = [Freqs[0]] + Freqs[2:] + [Freqs[1]]
            
            if max(Valid[A][E]) < InvalidThr: del(Valid[A][E])
            else:
                F = np.where((Valid[A][E] > InvalidThr)*(1 > Valid[A][E]))[0]
                freqs = [Freqs[f] for f in F]
                Valid[A][E] = [freqs, Valid[A][E][F]]
                
                if len(Valid[A][E]) == 0: del(Valid[A][E])
    
    return(Valid)


def GetMeans(Group, Animals, Exps, ExpList, AnalysisFile):
    FreqPerExp = GetFreqPerExp(Group, Animals, Exps, ExpList, AnalysisFile)
    
    SEMs = {}; Means = {}
    for E, Exp in FreqPerExp.items():
        SEMs[E] = {Freq: np.std(Val)/len(Val) for Freq, Val in Exp.items()}
        Means[E] = {Freq: np.nanmean(Val) for Freq, Val in Exp.items()}
    
    return(Means, SEMs)


def GetDiff(Group, Animals, Exps, ExpList, AnalysisFile, DiffThr=0.6, Invalid=False, InvalidThr=0.3):
    Index = GetIndex(Group, Animals, Exps, ExpList, AnalysisFile)
    FreqPerAnimal = GetFreqPerAnimal(Group, Animals, Exps, ExpList, AnalysisFile, Invalid, InvalidThr)
    FreqPerExp = GetFreqPerExp(Group, Animals, Exps, ExpList, AnalysisFile)
    
    PairList = list(combinations(FreqPerExp.keys(), 2))
    Diff = {}
    
    for A, Animal in FreqPerAnimal.items():
        Diff[A] = {}
        
        for Pair in PairList:
            PKey = '_'.join(Pair)
            E0 = ExpList.index(Pair[0]); E1 = ExpList.index(Pair[1])
            
            if ExpList[E0] not in Index[A] or ExpList[E1] not in Index[A]:
                continue
            
            Freqs = list(Index[A][ExpList[E0]].keys())
            Freqs.sort(key=lambda x: [int(y) for y in x.split('-')])
    #        Freqs = [Freqs[0]] + Freqs[2:] + [Freqs[1]]
            
            if len(FreqPerAnimal[A][E0]) == 0: continue
            
            Ratio = FreqPerAnimal[A][E1]/FreqPerAnimal[A][E0]
            Ratio[Ratio>1] = 1 - (1/Ratio[Ratio>1])
            Ratio[Ratio<=1] = 1 - Ratio[Ratio<=1]
            
            if max(Ratio) > DiffThr:
                F = np.where(Ratio > DiffThr)[0]
                freqs = [Freqs[f] for f in F]
                
                Diff[A][PKey] = [freqs, Ratio[F]]
                
                inv = np.where(Ratio[F] > 1)
                Diff[A][PKey][0] = np.delete(Diff[A][PKey][0], inv).tolist()
                Diff[A][PKey][1] = np.delete(Diff[A][PKey][1], inv)
                
                if len(Diff[A][PKey][0]) == 0: del(Diff[A][PKey])
            del(Ratio)
        
        if len(Diff[A]) == 0: del(Diff[A])
    
    return(Diff)


def GetMAF(Group, Animals, Exps, ExpList, AnalysisFile, DiffThr=0.5, Invalid=True, InvalidThr=0.1):
    Index = GetIndex(Group, Animals, Exps, ExpList, AnalysisFile)
    Diff = GetDiff(Group, Animals, Exps, ExpList, AnalysisFile, DiffThr, Invalid, InvalidThr)
    FreqPerAnimal = GetFreqPerAnimal(Group, Animals, Exps, ExpList, AnalysisFile, Invalid, InvalidThr)
    #MAF = np.zeros((len(Animals), len(ExpList)))
    MAF = [[] for _ in ExpList]
    
    for A, Animal in Diff.items():
        Freqs = list(Index[A][ExpList[0]].keys())
        Freqs.sort(key=lambda x: [int(y) for y in x.split('-')])
        
        BestFreq = np.argmax(Animal[ExpList[0]+'_'+ExpList[1]][1])
        BestFreq = Animal[ExpList[0]+'_'+ExpList[1]][0][BestFreq]
        BestFreq = Freqs.index(BestFreq)
        
#        An = Animals.index(A)
        for E in range(len(ExpList)): 
            if E < len(FreqPerAnimal[A]):
                MAF[E].append(FreqPerAnimal[A][E][BestFreq])
        
    return(MAF)


def GetPValues(Group, Animals, Exps, ExpList, AnalysisFile, DiffThr=0.5, Invalid=True, InvalidThr=0.1):
    IndexPerExp = GetMAF(Group, Animals, Exps, ExpList, AnalysisFile, DiffThr, Invalid, InvalidThr)
    PVals = {}; PairList = list(combinations(ExpList, 2))
    
    for Pair in PairList:
        PKey = '_'.join(Pair)
        E0 = ExpList.index(Pair[0]); E1 = ExpList.index(Pair[1])
        
        CL = 1 - 0.05
        DataA = IndexPerExp[E0]; DataB = IndexPerExp[E1]
#        if Group == 'Recovery' and PKey == 'BeforeANT_AfterANTCNO': 
#            DataB.append(np.mean(DataB)) # Recovery override
        
        if np.mean(DataA) > np.mean(DataB): DataA, DataB = DataB, DataA
        
        PVals[PKey] = Stats.RTTest(DataA, DataB, Confidence=CL)
#        PVals[PKey]['p.value'] *= len(PairList)
    
    return(PVals)


def Index_Exp_BP(Data, ExpList, PVals, Invalid=False, Show=True, Save=False, FigName=None, Ext=['svg']):
    Fig, Ax = plt.subplots(1, 1, figsize=(len(Data)*3,4))
    BoxPlot = Ax.boxplot(Data, showmeans=True)
    
    for I in range(len(Data)):
        BoxPlot['boxes'][I].set(label=ExpList[I])
        
    for K in ['boxes', 'whiskers', 'caps', 'medians', 'fliers']:
        for I in range(len(Data)): BoxPlot[K][I].set(color=Colors[I])
    
    Plot.Set(AxesObj=Ax, Axes=True)
    
#    Ax.legend(loc='best')
    Ax.set_ylabel('GPIASIndex')
    Ax.set_xticklabels(ExpList)
    
    Ax.spines['left'].set_position(('outward', 5))
    Ax.spines['bottom'].set_position(('outward', 5))
    
    Max = [max(Ax.get_xticks()), max(Ax.get_yticks())]
    Min = [min(Ax.get_xticks()), min(Ax.get_yticks())]
    
    for P, Pair in PVals.items():
        Text = str(round(Pair['p.value'], 4))
        P0 = P.split('_')[0]; P1 = P.split('_')[1]
        E0 = ExpList.index(P0); E1 = ExpList.index(P1)
        
        Plot.SignificanceBar([E0+1, E1+1], [Max[1]]*2, Text, Ax, TicksDir='down')
        Max = [max(Ax.get_xticks()), max(Ax.get_yticks())]
    
    Ax.spines['bottom'].set_bounds(Min[0], Max[0])
    Ax.spines['left'].set_bounds(Min[1], Max[1])
    
    Ax.set_ylim(Min[1], Max[1])
    
    if not FigName: FigName = './GPIASIndexBoxPlot'
    if Invalid: FigName = FigName + 'Invalid'
    
    if Save: 
        for E in Ext: Fig.savefig(FigName+'.'+E, format=E)
    
    if Show: plt.show()
    return(None)


#class Plot():
#    def Index_Freq_Exp_Group(Index, Groups, ExpList, Save=False):
#        Fig, Axes = plt.subplots(len(Groups), 1, sharex=True, figsize=(12,3*len(Groups)))
#        for G, Group in enumerate(Groups):
#            Animals = list(Index[Group].keys()); Animals.sort()
#            
#            for A, Animal in enumerate(Animals):
#                Freqs = list(Index[Group][Animal]['NaCl'].keys())
#                Freqs.sort(key=lambda x: [int(y) for y in x.split('-')])
#                Freqs = [Freqs[0]] + Freqs[2:] + [Freqs[1]]
#                
#                for F, Freq in enumerate(Freqs):
#                    Y = [abs(Index[Group][Animal][Exp][Freq]) for Exp in ExpList]
#                    X = np.arange(len(ExpList))
#                    
#                    if A == 0:
#                        Axes[G].plot(X+(len(ExpList)*F), Y, Colors[F]+'o-', label=Freq)
#                    else: 
#                        Axes[G].plot(X+(len(ExpList)*F), Y, Colors[F]+'o-')
#            
#            Axes[G].legend(bbox_to_anchor=(0.95, 0.55), loc='lower left', borderaxespad=0)
#            Plot.Set(AxesObj=Axes[G], Axes=True)
#            Axes[G].set_title(Group)
#            Axes[G].set_xticks(np.arange(len(Freqs))*3+1); Axes[G].set_xticklabels(Freqs)
#            Axes[G].set_ylabel('Mean GPIAS index')
#        
#        FigName = FigPath + '/GPIAZon-GPIASIndexMeanPerFreqPerExpPerAnimal.svg'
#        if Save: Fig.savefig(FigName, format='svg')
#        plt.show()
#
#    def AnimalNo_Freq_Exp_Thr_Bar(MeansFull, ExpList, Thrs=[0.1, 0.2, 0.3], Save=False):
#        YMax=0; Wid = (1/len(ExpList)) * 0.8
#        Fig, Axes = plt.subplots(len(Thrs), 1, sharex=True, figsize=(12,3.5*len(Thrs)))
#        for E, Exp in enumerate(ExpList):
#            Freqs = list(MeansFull[Exp].keys())
#            Freqs.sort(key=lambda x: [int(y) for y in x.split('-')])
#            Freqs = [Freqs[0]] + Freqs[2:] + [Freqs[1]]
#            
#            X = np.arange(len(Freqs))
#            Y = [MeansFull[Exp][Freq] for Freq in Freqs]; Y = np.array(Y)
#            
#            for T, Thr in enumerate(Thrs):
#                y = [len(a[a>Thr]) for a in Y]
#                if max(y) > YMax: YMax = max(y)
#                
#                Axes[T].bar(X+(E*Wid), y, width=Wid, color=Colors[E], alpha=0.4, label=Exp)
#        
#        for T, Thr in enumerate(Thrs):
#            Axes[T].legend(bbox_to_anchor=(1, 0.5), loc='lower left', borderaxespad=0)
#            Plot.Set(AxesObj=Axes[T], Axes=True)
#            Axes[T].set_ylim(0, YMax)
#            Axes[T].set_title(str(100*Thr)+'% decrease')
#            Axes[T].set_xticks(np.arange(len(Freqs))+0.3); Axes[T].set_xticklabels(Freqs)
#            Axes[T].spines['left'].set_position(('outward', 5))
#            Axes[T].spines['bottom'].set_position(('outward', 5))
#            Axes[T].set_ylabel('No of animals')
#        
#        FigName = FigPath + '/GPIAZon-GPIASIndexPerAnimalPerThreshold.svg'
#        if Save: Fig.savefig(FigName, format='svg')
#        plt.show()
#
#    def Index_Exp_BP(NaCl, SSal, ExpList, YMax=0.4, Invalid=False, Save=False):
#        Fig, Ax = plt.subplots(1, 1, figsize=(6,4))
#        BoxPlot = Ax.boxplot([NaCl, SSal], showmeans=True)
#        
#        for I in [0,1]: BoxPlot['boxes'][I].set(label=ExpList[I])
#        for K in ['boxes', 'whiskers', 'caps', 'medians', 'fliers']:
#            BoxPlot[K][0].set(color='r')
#            BoxPlot[K][1].set(color='k')
#        
#        Text = Stats.RTTest(SSal, NaCl, Confidence=1-0.05)
#        Text = str(round(Text['p.value'], 4))
#        Plot.SignificanceBar([1, 2], [YMax]*2, Text, Ax, TicksDir='down')
#        
#        Ax.legend(loc='best')
#        Ax.set_ylabel('GPIASIndex')
#        Ax.set_ylim(0, YMax)
#        Ax.set_xticklabels(ExpList)
#        Plot.Set(AxesObj=Ax, Axes=True)
#        Ax.spines['left'].set_position(('outward', 5))
#        Ax.spines['bottom'].set_position(('outward', 5))
#        
#        if Invalid: FigName = FigPath + '/GPIAZon-GPIASIndexDecreaseBoxPlotInvalid.svg'
#        else: FigName = FigPath + '/GPIAZon-GPIASIndexDecreaseBoxPlot.svg'
#        
#        if Save: Fig.savefig(FigName, format='svg')
#        plt.show()
#
#    def Index_Freq_Exp_Group_BP(MeansV, Pairs, Groups, ExpList, Save=False):
#        Fig, Axes = plt.subplots(len(Groups), 1, sharex=True, figsize=(14,4*len(Groups)))
#        for G, Group in enumerate(Groups):
#            Wid = (1/len(ExpList)) * 0.8
#            
#            for E, Exp in enumerate(ExpList):
#                Freqs = list(MeansV[Group][Exp].keys())
#                Freqs.sort(key=lambda x: [int(y) for y in x.split('-')])
#                Freqs = [Freqs[0]] + Freqs[2:] + [Freqs[1]]
#                
#                X = np.arange(len(Freqs))
#                Y = [MeansV[Group][Exp][Freq] for Freq in Freqs]
#        #        Error = [SEMs[Group][Exp][Freq] for Freq in Freqs]
#                
#                BoxPlot = Axes[G].boxplot(Y, positions=X+(E*Wid)+Wid/2, 
#                                          widths=[0.15]*len(X), showmeans=True)
#                
#                for Box in BoxPlot['boxes']: Box.set(color=Colors[E])
#                for Whisker in BoxPlot['whiskers']: Whisker.set(color=Colors[E])
#                for Cap in BoxPlot['caps']: Cap.set(color=Colors[E])
#                for Median in BoxPlot['medians']: Median.set(color=Colors[E])
#                for Flier in BoxPlot['fliers']: Flier.set(color=Colors[E])
#                Box.set(label=Exp)
#                    
#            
#            Axes[G].legend(bbox_to_anchor=(1, 0.5), loc='lower left', borderaxespad=0)
#            Plot.Set(AxesObj=Axes[G], Axes=True)
#            Axes[G].set_title(Group)
#            Axes[G].set_xticks(np.arange(len(Freqs))+0.3); Axes[G].set_xticklabels(Freqs)
#            Axes[G].spines['left'].set_position(('outward', 5))
#            Axes[G].spines['bottom'].set_position(('outward', 5))
#            Axes[G].set_ylabel('Mean GPIAS index')
#        
#        for G, Group in enumerate(Pairs.keys()):
#            for Pair in Pairs[Group].keys():
#                P = Pair.split('_')
#                
#                for FInd, Freq in enumerate(Freqs):
#                    if Freq not in Pairs[Group][Pair]: continue
#                    
#                    X = [FInd + (ExpList.index(P[_])*Wid) + (Wid/2) for _ in [0,1]]
#                    Y = max([max(MeansV[Group][P[_]][Freq]) for _ in [0,1]])
#                    AmpF = np.random.choice(np.arange(2, 4.25, 0.25))
#                    Y = Y + Y/AmpF; Y = [Y, Y]
#                    
#                    Text = str(round(Pairs[Group][Pair][Freq]['p.value'], 4))
#                    Plot.SignificanceBar(X, Y, Text, Axes[G], TicksDir='down')
#        
#        FigName = FigPath + '/GPIAZon-GPIASIndexMeanPerFreqPerExpBoxPlot.svg'
#        if Save: Fig.savefig(FigName, format='svg')
#        plt.show()
#
#    def Index_Freq_Exp_BP(MeansFull, PairsFull, ExpList, Save=False):
#        Wid = (1/len(ExpList)) * 0.8
#        Fig, Axes = plt.subplots(1, 1, sharex=True, figsize=(14,4))
#        for E, Exp in enumerate(ExpList):
#            Freqs = list(MeansFull[Exp].keys())
#            Freqs.sort(key=lambda x: [int(y) for y in x.split('-')])
#            Freqs = [Freqs[0]] + Freqs[2:] + [Freqs[1]]
#            
#            X = np.arange(len(Freqs))
#            Y = [MeansFull[Exp][Freq] for Freq in Freqs]
#            
#            BoxPlot = Axes.boxplot(Y, positions=X+(E*Wid)+Wid/2, 
#                                      widths=[0.15]*len(X), showmeans=True)
#            
#            for Box in BoxPlot['boxes']: Box.set(color=Colors[E])
#            for Whisker in BoxPlot['whiskers']: Whisker.set(color=Colors[E])
#            for Cap in BoxPlot['caps']: Cap.set(color=Colors[E])
#            for Median in BoxPlot['medians']: Median.set(color=Colors[E])
#            for Flier in BoxPlot['fliers']: Flier.set(color=Colors[E])
#            Box.set(label=Exp)
#        
#        for Pair in PairsFull.keys():
#            P = Pair.split('_')
#            
#            for FInd, Freq in enumerate(Freqs):
#                if Freq not in PairsFull[Pair]: continue
#                
#                X = [FInd + (ExpList.index(P[_])*Wid) + (Wid/2) for _ in [0,1]]
#                Y = max([max(MeansFull[P[_]][Freq]) for _ in [0,1]])
#                AmpF = np.random.choice(np.arange(2, 4.25, 0.25))
#                Y = Y + Y/AmpF; Y = [Y, Y]
#                
#                Text = str(round(PairsFull[Pair][Freq]['p.value'], 4))
#                Plot.SignificanceBar(X, Y, Text, Axes, TicksDir='down')
#        
#        Axes.legend(bbox_to_anchor=(1, 0.5), loc='lower left', borderaxespad=0)
#        Plot.Set(AxesObj=Axes, Axes=True)
#        Axes.set_xticks(np.arange(len(Freqs))+0.3); Axes.set_xticklabels(Freqs)
#        Axes.spines['left'].set_position(('outward', 5))
#        Axes.spines['bottom'].set_position(('outward', 5))
#        Axes.set_ylabel('Mean GPIAS index')
#        
#        FigName = FigPath + '/GPIAZon-GPIASIndexMeanPerFreqPerExpFullBoxPlot.svg'
#        if Save: Fig.savefig(FigName, format='svg')
#        plt.show()
#
#    def Index_Freq_Exp_Group_Bar(Means, SEMs, Pairs, Groups, ExpList, Save=False):
#        Wid = (1/len(ExpList)) * 0.8
#        Fig, Axes = plt.subplots(len(Groups), 1, sharex=True, figsize=(10,3*len(Groups)))
#        for G, Group in enumerate(Groups):
#            for E, Exp in enumerate(ExpList):
#                Freqs = list(Means[Group][Exp].keys())
#                Freqs.sort(key=lambda x: [int(y) for y in x.split('-')])
#                Freqs = [Freqs[0]] + Freqs[2:] + [Freqs[1]]
#                
#                X = np.arange(len(Freqs))
#                Y = [Means[Group][Exp][Freq] for Freq in Freqs]
#                Error = [SEMs[Group][Exp][Freq] for Freq in Freqs]
#                
#                Axes[G].bar(X+(E*Wid), Y, width=Wid, color=Colors[E], alpha=0.4, label=Exp)
#                Axes[G].errorbar(X+(E*Wid)+(Wid/2), Y, Error, color='k', fmt='.')
#            
#            Axes[G].legend(bbox_to_anchor=(1, 0.5), loc='lower left', borderaxespad=0)
#            Plot.Set(AxesObj=Axes[G], Axes=True)
#            Axes[G].set_title(Group)
#            Axes[G].set_xticks(np.arange(len(Freqs))+0.3); Axes[G].set_xticklabels(Freqs)
#            Axes[G].set_ylabel('Mean GPIAS index')
#        
#        for G, Group in enumerate(Pairs.keys()):
#            for Pair in Pairs[Group].keys():
#                P = Pair.split('_')
#                
#                for FInd, Freq in enumerate(Freqs):
#                    if Freq not in Pairs[Group][Pair]: continue
#                    
#                    X = [FInd + (ExpList.index(P[_])*Wid) + (Wid/2) for _ in [0,1]]
#                    Y = max([Means[Group][P[_]][Freq] + SEMs[Group][P[_]][Freq] for _ in [0,1]])
#                    AmpF = np.random.choice(np.arange(2, 4.25, 0.25))
#                    Y = Y + Y/AmpF; Y = [Y, Y]
#                    
#                    Text = str(round(Pairs[Group][Pair][Freq]['p.value'], 4))
#                    Plot.SignificanceBar(X, Y, Text, Axes[G], TicksDir='down')
#        
#        FigName = FigPath + '/GPIAZon-GPIASIndexMeanPerFreqPerExp.svg'
#        if Save: Fig.savefig(FigName, format='svg')
#        plt.show()
#
#    def Index_Freq_Exp_Group_Animal(Index, Groups, ExpList, Save=False):
#        Fig, Axes = plt.subplots(5, len(Groups), sharex=True, figsize=(5*len(Groups),2*5))
#        Colors = ['r', 'g', 'b', 'm', 'k', '#ffa500', '#00b2b2']
#        YMax = 0
#        for G, Group in enumerate(Groups):
#            Animals = [Group.split('_')[-1] + 'n0' + str(_) for _ in range(1,6)]
#            
#            for A, Animal in enumerate(Animals):
#                for E, Exp in enumerate(ExpList):
#                    if Exp not in Index[Group][Animal]: continue
#                    
#                    Freqs = list(Index[Group][Animal][Exp].keys())
#                    Freqs.sort(key=lambda x: [int(y) for y in x.split('-')])
#                    Freqs = [Freqs[0]] + Freqs[2:]# + [Freqs[1]]
#                    
#                    Y = [abs(Index[Group][Animal][Exp][Freq]) for Freq in Freqs]
#                    
#                    if max(Y) > YMax: YMax = max(Y)
#                    
#                    Axes[A][G].plot(Y, Colors[E]+'o-', label=Exp)
#                
#                Axes[A][G].legend(loc='best')
#                Plot.Set(AxesObj=Axes[A][G], Axes=True)
#                Axes[A][G].set_title(Animal)
#                Axes[A][G].set_ylabel('Mean GPIAS index')
#        
#        for G, Group in enumerate(Groups):
#            Animals = [Group.split('_')[-1] + 'n0' + str(_) for _ in range(1,6)]
#            for A, Animal in enumerate(Animals):
#                Axes[A][G].set_ylim(0, YMax)
#        
#        plt.show()    
#
#    def RawTraces(Freq, Keys, AnalysisFile, Save=False):
#        Freq = '8000-16000'
#        Keys = ['GPIAZon_NaCl/2017-04-11_16-01-38_GPIAZon_NaCln03', 
#                'GPIAZon_NaCl/2017-04-19_18-28-39_GPIAZon_NaCln03']
#        Fig, Axes = plt.subplots(2, 1, sharex=True, figsize=(6,3*2))
#        for K, Key in enumerate(Keys):
#            GPIAS, XValues = Hdf5.GPIASLoad(AnalysisFile, Key)
#            Ind1 = list(XValues).index(0)
#            Ind2 = list(XValues).index(int(0.05*1000))
#        
#            Axes[K].plot(XValues, GPIAS['Trace'][Freq]['NoGap'], 'r', label='NoGap')
#            Axes[K].plot(XValues, GPIAS['Trace'][Freq]['Gap'], 'b', label='Gap')
#            Axes[K].axvspan(XValues[Ind1], XValues[Ind2], color='k', alpha=0.5, lw=0, label='Sound pulse')
#            
#            Axes[K].legend(loc='best')
#            Plot.Set(AxesObj=Axes[K], Axes=True)
#            Axes[K].set_ylim(-5, 5)
#            Axes[K].spines['left'].set_position(('outward', 5))
#            Axes[K].spines['bottom'].set_position(('outward', 5))
#            Axes[K].set_ylabel('Mean GPIAS index')
#        
#        FigName = FigPath + '/GPIAZon-GPIASTraces.svg'
#        if Save: Fig.savefig(FigName, format='svg')
#        plt.show()
#
#
