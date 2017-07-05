#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jun 12 14:12:37 2017

@author: malfatti
"""
from rpy2 import robjects as RObj
from rpy2.robjects import packages as RPkg

## Level 0
def RCheckPackage(Packages):
    RPacksToInstall = [Pack for Pack in Packages 
                       if not RPkg.isinstalled(Pack)]
    if len(RPacksToInstall) > 0:
        print(str(RPacksToInstall), 'not installed. Install now?')
        Ans = input('[y/N]: ')
        
        if Ans.lower() in ['y', 'yes']:
            from rpy2.robjects.vectors import StrVector as RStrVector
            
            RUtils = RPkg.importr('utils')
            RUtils.chooseCRANmirror(ind=1)
            
            RUtils.install_packages(RStrVector(RPacksToInstall))
        
        else: print('Aborted.')
    
    else: print('Packages', str(Packages), 'installed.')
    
    return(None)


def AdjustNaNs(Array):
    NaN = RObj.NA_Real
    
    for I, A in enumerate(Array):
        if A != A: Array[I] = NaN
        
    return(Array)


## Level 1
def RAnOVa(DataA, DataB):
    Raov = RObj.r['aov']
    
    Data = DataA + DataB
    Trtmnt = ['A']*len(DataA) + ['B']*len(DataB)
    Formula = RObj.Formula('Data ~ Trtmnt'); FEnv = Formula.environment
    FEnv['Data'] = RObj.FloatVector(Data)
    FEnv['Trtmnt'] = RObj.StrVector(Trtmnt)
    
    Results = Raov(Formula)
    print(RObj.r['summary'](Results))
    
    return(Results)


def RPwrAnOVa(GroupNo=RObj.NULL, SampleSize=RObj.NULL, Power=RObj.NULL, 
           SigLevel=RObj.NULL, EffectSize=RObj.NULL):
    RCheckPackage(['pwr']); Rpwr = RPkg.importr('pwr')
    
    Results = Rpwr.pwr_anova_test(k=GroupNo, power=Power, sig_level=SigLevel, 
                                  f=EffectSize, n=SampleSize)
    
    print('Running', Results.rx('method')[0][0] + '... ', end='')
    AnOVaResults = {}
    for Key, Value in {'k': 'GroupNo', 'n': 'SampleSize', 'f': 'EffectSize', 
                       'power':'Power', 'sig.level': 'SigLevel'}.items():
        AnOVaResults[Value] = Results.rx(Key)[0][0]
    
    print('Done.')
    return(AnOVaResults)


def RTTest(DataA, DataB, Paired=True, Alt='less', Confidence=0.95):
    Rttest = RObj.r['t.test']
    
    DataA = AdjustNaNs(DataA); DataB = AdjustNaNs(DataB)
    
    Results = Rttest(RObj.FloatVector(DataA), RObj.FloatVector(DataB), 
                     paired=Paired, var_equal=False, alternative=Alt, 
                     conf_level=RObj.FloatVector([Confidence]), 
                     na_action=RObj.r['na.omit'])
    
    print('Calculating', Results.rx('method')[0][0] + '... ', end='')
    TTestResults = {}; Names = list(Results.names)
    for Name in Names:
        TTestResults[Name] = Results.rx(Name)[0][0]
    
    print('Done.')
    return(TTestResults)

