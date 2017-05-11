#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
GPIAZon group analysis
"""
## Import
import DataAnalysis, Hdf5F
import numpy as np
from glob import glob
from itertools import combinations
from os import makedirs

Params = DataAnalysis.Plot.Set(Params=True)
from matplotlib import rcParams; rcParams.update(Params)
from matplotlib import pyplot as plt


## Analysis
FigPath = 'GPIAZon/Figs'; makedirs(FigPath, exist_ok=True)
Colors = ['k', 'r', 'b', 'm', 'g', '#ffa500', '#00b2b2']    

def GetIndex(Groups, Exps, AnalysisFile):
    Index = {}
    for Group in Groups:
        Animals = [Group.split('_')[-1] + 'n0' + str(_) for _ in range(1,6)]
        Index[Group] = {}
        
        for Animal in Animals:
            Paths = glob('GPIAZon/' + Group + '/*' + Animal); Paths.sort()
            Index[Group][Animal] = {}
            
            for P,Path in enumerate(Paths):
                AnalysisKey = Group + '/' + Path.split('/')[-1]
                GPIAS = Hdf5F.LoadGPIAS(AnalysisFile, AnalysisKey)[0]
                Index[Group][Animal][Exps[Group][P]] = {}
    #            Freqs = list(GPIAS['Index'].keys())
    #            Freqs.sort(key=lambda x: [int(y) for y in x.split('-')])
                
                for Freq in GPIAS['Index'].keys(): 
                    if Freq == '9000-11000': continue
                    Index[Group][Animal][Exps[Group][P]][Freq] = GPIAS['Index'][Freq]['GPIASIndex']
                
                del(GPIAS)
    
    return(Index)

def GetMeansV(Index, Groups, ExpList):
    MeansV = {}; MeansFull = {}
    for Group in Groups:
        Animals = [Group.split('_')[-1] + 'n0' + str(_) for _ in range(1,6)]
        MeansV[Group] = {}
        
        for Animal in Animals:
            for Exp in ExpList:
                if Exp not in MeansV[Group]: MeansV[Group][Exp] = {}
                if Exp not in MeansFull: MeansFull[Exp] = {}
                
                for Freq in Index[Group][Animal][Exp].keys():
                    if Freq not in MeansV[Group][Exp]: MeansV[Group][Exp][Freq] = []
                    if Freq not in MeansFull[Exp]: MeansFull[Exp][Freq] = []
                    
                    MeansV[Group][Exp][Freq].append(abs(Index[Group][Animal][Exp][Freq]))
                    MeansFull[Exp][Freq].append(abs(Index[Group][Animal][Exp][Freq]))
    
    for Group in Groups:
        for Exp in ExpList:
            for F, Freq in MeansV[Group][Exp].items():
                if len(Freq) == 3: 
                    MeansV[Group][Exp][F] = MeansV[Group][Exp][F] + [np.mean(MeansV[Group][Exp][F])]*2
                
                if len(MeansFull[Exp][F]) == 8:
                    MeansFull[Exp][F] = MeansFull[Exp][F] + [np.mean(MeansFull[Exp][F])]*2
    
    return(MeansV, MeansFull)

def GetDiff(Index, Groups, ExpList, DiffThr=0.6, Invalid=False, InvalidThr=0.1):
    Diff = {}
    for G, Group in enumerate(Groups):
        Animals = [Group.split('_')[-1] + 'n0' + str(_) for _ in range(1,6)]
        
        for A, Animal in enumerate(Animals):
            if Animal not in Diff: Diff[Animal] = []
            
            for E, Exp in enumerate(ExpList):
                if Exp not in Index[Group][Animal]: continue
                
                Freqs = list(Index[Group][Animal][Exp].keys())
                Freqs.sort(key=lambda x: [int(y) for y in x.split('-')])
                Freqs = [Freqs[0]] + Freqs[2:]# + [Freqs[1]]
                
                Y = [abs(Index[Group][Animal][Exp][Freq]) for Freq in Freqs]
                Diff[Animal].append(np.array(Y))
            
            if not Invalid: 
                inv = np.where(Diff[Animal][0] < InvalidThr)[0]
    #            Diff[Animal][0] = np.delete(Diff[Animal][0], inv)
    #            Diff[Animal][1] = np.delete(Diff[Animal][1], inv)
                Diff[Animal][0][inv] = -1
            
            if len(Diff[Animal][0]) == 0: del(Diff[Animal]); continue
            
            diff = 1 - (Diff[Animal][1]/Diff[Animal][0])
            if max(diff) < DiffThr: del(Diff[Animal])
            else:
                F = np.where(diff > DiffThr)[0]
                freqs = [Freqs[f] for f in F]
                Diff[Animal] = [freqs, diff[F]]
                
                inv = np.where(diff[F] > 1)
                Diff[Animal][0] = np.delete(Diff[Animal][0], inv).tolist()
                Diff[Animal][1] = np.delete(Diff[Animal][1], inv)
                
                if len(Diff[Animal][0]) == 0: del(Diff[Animal])
            del(diff)
    
    return(Diff)

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

def GetMeans(MeansV, Groups, ExpList):
    SEMs = {}; Means = {}
    for Group in Groups:
        SEMs[Group] = {}; Means[Group] = {}
        for Exp in ExpList:
            SEMs[Group][Exp] = {Freq: np.std(Val)/len(Val) for Freq, Val in MeansV[Group][Exp].items()}
            Means[Group][Exp] = {Freq: np.nanmean(Val) for Freq, Val in MeansV[Group][Exp].items()}
    
    return(Means, SEMs)

def GetIndFreqMean(Index, Diff, Groups):
    NaCl = []; SSal = []
    for G, Group in enumerate(Groups):
        Animals = [Group.split('_')[-1] + 'n0' + str(_) for _ in range(1,6)]
        
        for A, Animal in enumerate(Animals):
            if Animal not in  Diff: continue
            
            F = np.where(Diff[Animal][1] == max(Diff[Animal][1]))[0][0]
            Freq = Diff[Animal][0][F]
            NaCl.append(abs(Index[Group][Animal]['NaCl'][Freq]))
            SSal.append(abs(Index[Group][Animal]['SSal'][Freq]))
    
    return(NaCl, SSal)

def GetPairs(MeansV, MeansFull, Groups, ExpList):
    Pairs = {}; PairsFull = {}; PairList = list(combinations(ExpList, 2))
    for Group in Groups:
        Pairs[Group] = {}
        
        for Pair in PairList:
            PKey = '_'.join(Pair)
            Pairs[Group][PKey] = {}
            PairsFull[PKey] = {}
            
            for Freq in MeansV[Group][Pair[0]].keys():
                CL = 1 - (0.05/len(PairList))
                DataA = MeansV[Group][Pair[0]][Freq][:]
                DataB = MeansV[Group][Pair[1]][Freq][:]
                if np.mean(DataA) > np.mean(DataB): DataA, DataB = DataB, DataA
                
                DataC = MeansFull[Pair[0]][Freq][:]
                DataD = MeansFull[Pair[1]][Freq][:]
                if np.mean(DataC) > np.mean(DataD): DataC, DataD = DataD, DataC
                
                Pairs[Group][PKey][Freq] = DataAnalysis.Stats.RTTest(DataA, DataB, Confidence=CL)
                PairsFull[PKey][Freq] = DataAnalysis.Stats.RTTest(DataC, DataD, Confidence=CL)
    
    return(Pairs, PairsFull)

def ClearPairs(Pairs, PairsFull):
    ToDelete = []; ToDeleteFull = []
    for Group in Pairs.keys():
        for Pair in Pairs[Group].keys():
            for Freq in Pairs[Group][Pair].keys():
                if Pairs[Group][Pair][Freq]['p.value'] > 0.05:
                    ToDelete.append([Group, Pair, Freq])
                
                if Pairs[Group][Pair][Freq]['p.value'] != Pairs[Group][Pair][Freq]['p.value']:
                    ToDelete.append([Group, Pair, Freq])
    
    for Pair in PairsFull.keys():
        for Freq in PairsFull[Pair].keys():
            if PairsFull[Pair][Freq]['p.value'] > 0.05:
                ToDeleteFull.append([Pair, Freq])
            
            if PairsFull[Pair][Freq]['p.value'] != PairsFull[Pair][Freq]['p.value']:
                ToDeleteFull.append([Pair, Freq])
    
    for KeyPair in ToDelete: del(Pairs[KeyPair[0]][KeyPair[1]][KeyPair[2]])
    for KeyPair in ToDeleteFull: del(PairsFull[KeyPair[0]][KeyPair[1]])
    
    EmptyPairs = []; EmptyPairsFull = []
    for Group in Pairs.keys():
        for Pair in Pairs[Group].keys():
            if len(Pairs[Group][Pair]) == 0: EmptyPairs.append([Group, Pair])
    
    for Pair in PairsFull.keys():
        if len(PairsFull[Pair]) == 0: EmptyPairsFull.append([Pair])
    
    for KeyPair in EmptyPairs: del(Pairs[KeyPair[0]][KeyPair[1]])
    for KeyPair in EmptyPairsFull: del(PairsFull[KeyPair[0]])
    
    return(Pairs, PairsFull)

class Plot():
    def Index_Freq_Exp_Group(Index, Groups, ExpList, Save=False):
        Fig, Axes = plt.subplots(len(Groups), 1, sharex=True, figsize=(12,3*len(Groups)))
        for G, Group in enumerate(Groups):
            Animals = list(Index[Group].keys()); Animals.sort()
            
            for A, Animal in enumerate(Animals):
                Freqs = list(Index[Group][Animal]['NaCl'].keys())
                Freqs.sort(key=lambda x: [int(y) for y in x.split('-')])
                Freqs = [Freqs[0]] + Freqs[2:] + [Freqs[1]]
                
                for F, Freq in enumerate(Freqs):
                    Y = [abs(Index[Group][Animal][Exp][Freq]) for Exp in ExpList]
                    X = np.arange(len(ExpList))
                    
                    if A == 0:
                        Axes[G].plot(X+(len(ExpList)*F), Y, Colors[F]+'o-', label=Freq)
                    else: 
                        Axes[G].plot(X+(len(ExpList)*F), Y, Colors[F]+'o-')
            
            Axes[G].legend(bbox_to_anchor=(0.95, 0.55), loc='lower left', borderaxespad=0)
            DataAnalysis.Plot.Set(AxesObj=Axes[G], Axes=True)
            Axes[G].set_title(Group)
            Axes[G].set_xticks(np.arange(len(Freqs))*3+1); Axes[G].set_xticklabels(Freqs)
            Axes[G].set_ylabel('Mean GPIAS index')
        
        FigName = FigPath + '/GPIAZon-GPIASIndexMeanPerFreqPerExpPerAnimal.svg'
        if Save: Fig.savefig(FigName, format='svg')
        plt.show()

    def AnimalNo_Freq_Exp_Thr_Bar(MeansFull, ExpList, Thrs=[0.1, 0.2, 0.3], Save=False):
        YMax=0; Wid = (1/len(ExpList)) * 0.8
        Fig, Axes = plt.subplots(len(Thrs), 1, sharex=True, figsize=(12,3.5*len(Thrs)))
        for E, Exp in enumerate(ExpList):
            Freqs = list(MeansFull[Exp].keys())
            Freqs.sort(key=lambda x: [int(y) for y in x.split('-')])
            Freqs = [Freqs[0]] + Freqs[2:] + [Freqs[1]]
            
            X = np.arange(len(Freqs))
            Y = [MeansFull[Exp][Freq] for Freq in Freqs]; Y = np.array(Y)
            
            for T, Thr in enumerate(Thrs):
                y = [len(a[a>Thr]) for a in Y]
                if max(y) > YMax: YMax = max(y)
                
                Axes[T].bar(X+(E*Wid), y, width=Wid, color=Colors[E], alpha=0.4, label=Exp)
        
        for T, Thr in enumerate(Thrs):
            Axes[T].legend(bbox_to_anchor=(1, 0.5), loc='lower left', borderaxespad=0)
            DataAnalysis.Plot.Set(AxesObj=Axes[T], Axes=True)
            Axes[T].set_ylim(0, YMax)
            Axes[T].set_title(str(100*Thr)+'% decrease')
            Axes[T].set_xticks(np.arange(len(Freqs))+0.3); Axes[T].set_xticklabels(Freqs)
            Axes[T].spines['left'].set_position(('outward', 5))
            Axes[T].spines['bottom'].set_position(('outward', 5))
            Axes[T].set_ylabel('No of animals')
        
        FigName = FigPath + '/GPIAZon-GPIASIndexPerAnimalPerThreshold.svg'
        if Save: Fig.savefig(FigName, format='svg')
        plt.show()

    def Index_Exp_BP(NaCl, SSal, ExpList, YMax=0.4, Invalid=False, Save=False):
        Fig, Ax = plt.subplots(1, 1, figsize=(6,4))
        BoxPlot = Ax.boxplot([NaCl, SSal], showmeans=True)
        
        for I in [0,1]: BoxPlot['boxes'][I].set(label=ExpList[I])
        for K in ['boxes', 'whiskers', 'caps', 'medians', 'fliers']:
            BoxPlot[K][0].set(color='r')
            BoxPlot[K][1].set(color='k')
        
        Text = DataAnalysis.Stats.RTTest(SSal, NaCl, Confidence=1-0.05)
        Text = str(round(Text['p.value'], 4))
        DataAnalysis.Plot.SignificanceBar([1, 2], [YMax]*2, Text, Ax, TicksDir='down')
        
        Ax.legend(loc='best')
        Ax.set_ylabel('GPIASIndex')
        Ax.set_ylim(0, YMax)
        Ax.set_xticklabels(ExpList)
        DataAnalysis.Plot.Set(AxesObj=Ax, Axes=True)
        Ax.spines['left'].set_position(('outward', 5))
        Ax.spines['bottom'].set_position(('outward', 5))
        
        if Invalid: FigName = FigPath + '/GPIAZon-GPIASIndexDecreaseBoxPlotInvalid.svg'
        else: FigName = FigPath + '/GPIAZon-GPIASIndexDecreaseBoxPlot.svg'
        
        if Save: Fig.savefig(FigName, format='svg')
        plt.show()

    def Index_Freq_Exp_Group_BP(MeansV, Pairs, Groups, ExpList, Save=False):
        Fig, Axes = plt.subplots(len(Groups), 1, sharex=True, figsize=(14,4*len(Groups)))
        for G, Group in enumerate(Groups):
            Wid = (1/len(ExpList)) * 0.8
            
            for E, Exp in enumerate(ExpList):
                Freqs = list(MeansV[Group][Exp].keys())
                Freqs.sort(key=lambda x: [int(y) for y in x.split('-')])
                Freqs = [Freqs[0]] + Freqs[2:] + [Freqs[1]]
                
                X = np.arange(len(Freqs))
                Y = [MeansV[Group][Exp][Freq] for Freq in Freqs]
        #        Error = [SEMs[Group][Exp][Freq] for Freq in Freqs]
                
                BoxPlot = Axes[G].boxplot(Y, positions=X+(E*Wid)+Wid/2, 
                                          widths=[0.15]*len(X), showmeans=True)
                
                for Box in BoxPlot['boxes']: Box.set(color=Colors[E])
                for Whisker in BoxPlot['whiskers']: Whisker.set(color=Colors[E])
                for Cap in BoxPlot['caps']: Cap.set(color=Colors[E])
                for Median in BoxPlot['medians']: Median.set(color=Colors[E])
                for Flier in BoxPlot['fliers']: Flier.set(color=Colors[E])
                Box.set(label=Exp)
                    
            
            Axes[G].legend(bbox_to_anchor=(1, 0.5), loc='lower left', borderaxespad=0)
            DataAnalysis.Plot.Set(AxesObj=Axes[G], Axes=True)
            Axes[G].set_title(Group)
            Axes[G].set_xticks(np.arange(len(Freqs))+0.3); Axes[G].set_xticklabels(Freqs)
            Axes[G].spines['left'].set_position(('outward', 5))
            Axes[G].spines['bottom'].set_position(('outward', 5))
            Axes[G].set_ylabel('Mean GPIAS index')
        
        for G, Group in enumerate(Pairs.keys()):
            for Pair in Pairs[Group].keys():
                P = Pair.split('_')
                
                for FInd, Freq in enumerate(Freqs):
                    if Freq not in Pairs[Group][Pair]: continue
                    
                    X = [FInd + (ExpList.index(P[_])*Wid) + (Wid/2) for _ in [0,1]]
                    Y = max([max(MeansV[Group][P[_]][Freq]) for _ in [0,1]])
                    AmpF = np.random.choice(np.arange(2, 4.25, 0.25))
                    Y = Y + Y/AmpF; Y = [Y, Y]
                    
                    Text = str(round(Pairs[Group][Pair][Freq]['p.value'], 4))
                    DataAnalysis.Plot.SignificanceBar(X, Y, Text, Axes[G], TicksDir='down')
        
        FigName = FigPath + '/GPIAZon-GPIASIndexMeanPerFreqPerExpBoxPlot.svg'
        if Save: Fig.savefig(FigName, format='svg')
        plt.show()

    def Index_Freq_Exp_BP(MeansFull, PairsFull, ExpList, Save=False):
        Wid = (1/len(ExpList)) * 0.8
        Fig, Axes = plt.subplots(1, 1, sharex=True, figsize=(14,4))
        for E, Exp in enumerate(ExpList):
            Freqs = list(MeansFull[Exp].keys())
            Freqs.sort(key=lambda x: [int(y) for y in x.split('-')])
            Freqs = [Freqs[0]] + Freqs[2:] + [Freqs[1]]
            
            X = np.arange(len(Freqs))
            Y = [MeansFull[Exp][Freq] for Freq in Freqs]
            
            BoxPlot = Axes.boxplot(Y, positions=X+(E*Wid)+Wid/2, 
                                      widths=[0.15]*len(X), showmeans=True)
            
            for Box in BoxPlot['boxes']: Box.set(color=Colors[E])
            for Whisker in BoxPlot['whiskers']: Whisker.set(color=Colors[E])
            for Cap in BoxPlot['caps']: Cap.set(color=Colors[E])
            for Median in BoxPlot['medians']: Median.set(color=Colors[E])
            for Flier in BoxPlot['fliers']: Flier.set(color=Colors[E])
            Box.set(label=Exp)
        
        for Pair in PairsFull.keys():
            P = Pair.split('_')
            
            for FInd, Freq in enumerate(Freqs):
                if Freq not in PairsFull[Pair]: continue
                
                X = [FInd + (ExpList.index(P[_])*Wid) + (Wid/2) for _ in [0,1]]
                Y = max([max(MeansFull[P[_]][Freq]) for _ in [0,1]])
                AmpF = np.random.choice(np.arange(2, 4.25, 0.25))
                Y = Y + Y/AmpF; Y = [Y, Y]
                
                Text = str(round(PairsFull[Pair][Freq]['p.value'], 4))
                DataAnalysis.Plot.SignificanceBar(X, Y, Text, Axes, TicksDir='down')
        
        Axes.legend(bbox_to_anchor=(1, 0.5), loc='lower left', borderaxespad=0)
        DataAnalysis.Plot.Set(AxesObj=Axes, Axes=True)
        Axes.set_xticks(np.arange(len(Freqs))+0.3); Axes.set_xticklabels(Freqs)
        Axes.spines['left'].set_position(('outward', 5))
        Axes.spines['bottom'].set_position(('outward', 5))
        Axes.set_ylabel('Mean GPIAS index')
        
        FigName = FigPath + '/GPIAZon-GPIASIndexMeanPerFreqPerExpFullBoxPlot.svg'
        if Save: Fig.savefig(FigName, format='svg')
        plt.show()

    def Index_Freq_Exp_Group_Bar(Means, SEMs, Pairs, Groups, ExpList, Save=False):
        Wid = (1/len(ExpList)) * 0.8
        Fig, Axes = plt.subplots(len(Groups), 1, sharex=True, figsize=(10,3*len(Groups)))
        for G, Group in enumerate(Groups):
            for E, Exp in enumerate(ExpList):
                Freqs = list(Means[Group][Exp].keys())
                Freqs.sort(key=lambda x: [int(y) for y in x.split('-')])
                Freqs = [Freqs[0]] + Freqs[2:] + [Freqs[1]]
                
                X = np.arange(len(Freqs))
                Y = [Means[Group][Exp][Freq] for Freq in Freqs]
                Error = [SEMs[Group][Exp][Freq] for Freq in Freqs]
                
                Axes[G].bar(X+(E*Wid), Y, width=Wid, color=Colors[E], alpha=0.4, label=Exp)
                Axes[G].errorbar(X+(E*Wid)+(Wid/2), Y, Error, color='k', fmt='.')
            
            Axes[G].legend(bbox_to_anchor=(1, 0.5), loc='lower left', borderaxespad=0)
            DataAnalysis.Plot.Set(AxesObj=Axes[G], Axes=True)
            Axes[G].set_title(Group)
            Axes[G].set_xticks(np.arange(len(Freqs))+0.3); Axes[G].set_xticklabels(Freqs)
            Axes[G].set_ylabel('Mean GPIAS index')
        
        for G, Group in enumerate(Pairs.keys()):
            for Pair in Pairs[Group].keys():
                P = Pair.split('_')
                
                for FInd, Freq in enumerate(Freqs):
                    if Freq not in Pairs[Group][Pair]: continue
                    
                    X = [FInd + (ExpList.index(P[_])*Wid) + (Wid/2) for _ in [0,1]]
                    Y = max([Means[Group][P[_]][Freq] + SEMs[Group][P[_]][Freq] for _ in [0,1]])
                    AmpF = np.random.choice(np.arange(2, 4.25, 0.25))
                    Y = Y + Y/AmpF; Y = [Y, Y]
                    
                    Text = str(round(Pairs[Group][Pair][Freq]['p.value'], 4))
                    DataAnalysis.Plot.SignificanceBar(X, Y, Text, Axes[G], TicksDir='down')
        
        FigName = FigPath + '/GPIAZon-GPIASIndexMeanPerFreqPerExp.svg'
        if Save: Fig.savefig(FigName, format='svg')
        plt.show()

    def Index_Freq_Exp_Group_Animal(Index, Groups, ExpList, Save=False):
        Fig, Axes = plt.subplots(5, len(Groups), sharex=True, figsize=(5*len(Groups),2*5))
        Colors = ['r', 'g', 'b', 'm', 'k', '#ffa500', '#00b2b2']
        YMax = 0
        for G, Group in enumerate(Groups):
            Animals = [Group.split('_')[-1] + 'n0' + str(_) for _ in range(1,6)]
            
            for A, Animal in enumerate(Animals):
                for E, Exp in enumerate(ExpList):
                    if Exp not in Index[Group][Animal]: continue
                    
                    Freqs = list(Index[Group][Animal][Exp].keys())
                    Freqs.sort(key=lambda x: [int(y) for y in x.split('-')])
                    Freqs = [Freqs[0]] + Freqs[2:]# + [Freqs[1]]
                    
                    Y = [abs(Index[Group][Animal][Exp][Freq]) for Freq in Freqs]
                    
                    if max(Y) > YMax: YMax = max(Y)
                    
                    Axes[A][G].plot(Y, Colors[E]+'o-', label=Exp)
                
                Axes[A][G].legend(loc='best')
                DataAnalysis.Plot.Set(AxesObj=Axes[A][G], Axes=True)
                Axes[A][G].set_title(Animal)
                Axes[A][G].set_ylabel('Mean GPIAS index')
        
        for G, Group in enumerate(Groups):
            Animals = [Group.split('_')[-1] + 'n0' + str(_) for _ in range(1,6)]
            for A, Animal in enumerate(Animals):
                Axes[A][G].set_ylim(0, YMax)
        
        plt.show()    

    def RawTraces(Freq, Keys, AnalysisFile, Save=False):
        Freq = '8000-16000'
        Keys = ['GPIAZon_NaCl/2017-04-11_16-01-38_GPIAZon_NaCln03', 
                'GPIAZon_NaCl/2017-04-19_18-28-39_GPIAZon_NaCln03']
        Fig, Axes = plt.subplots(2, 1, sharex=True, figsize=(6,3*2))
        for K, Key in enumerate(Keys):
            GPIAS, XValues = Hdf5F.LoadGPIAS(AnalysisFile, Key)
            Ind1 = list(XValues).index(0)
            Ind2 = list(XValues).index(int(0.05*1000))
        
            Axes[K].plot(XValues, GPIAS['Trace'][Freq]['NoGap'], 'r', label='NoGap')
            Axes[K].plot(XValues, GPIAS['Trace'][Freq]['Gap'], 'b', label='Gap')
            Axes[K].axvspan(XValues[Ind1], XValues[Ind2], color='k', alpha=0.5, lw=0, label='Sound pulse')
            
            Axes[K].legend(loc='best')
            DataAnalysis.Plot.Set(AxesObj=Axes[K], Axes=True)
            Axes[K].set_ylim(-5, 5)
            Axes[K].spines['left'].set_position(('outward', 5))
            Axes[K].spines['bottom'].set_position(('outward', 5))
            Axes[K].set_ylabel('Mean GPIAS index')
        
        FigName = FigPath + '/GPIAZon-GPIASTraces.svg'
        if Save: Fig.savefig(FigName, format='svg')
        plt.show()


