#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Treadmill analysis
@author: malfatti
"""
#%% Import
import DataAnalysis, Hdf5F
import numpy as np

from glob import glob
from os import makedirs


#%% Calculate spectrograms data
Delta, Theta = [2, 4], [7, 10]
SensorCh = 17
Diameter = 0.6; PeaksPerCycle = 12
Lowpass = 5; FilterOrder = 2

#cd ~/Martina/Treadmill/Data/201703/
AnalysisFile = '201703-Treadmill.hdf5'
Paths = glob('2017-03-08*'); Paths.sort()

PeakDist = (np.pi*Diameter)/PeaksPerCycle

Done, Errors, ErrorsLog = [], [], []
for Path in Paths:
    Data = Hdf5F.LoadOEKwik(Path)[0]
    Proc = Hdf5F.GetProc(Data, 'OE')
    
#    Recs = list(Data[Proc]['data'].keys())
    for Rec in Data[Proc]['data'].keys():
        print('Processing Rec', "{0:02d}".format(int(Rec)), '...')
        RecData = Data[Proc]['data'][Rec]
        Rate = Data[Proc]['info'][Rec]['sample_rate']
        
        try:
            Treadmill = DataAnalysis.Treadmill.Analysis(RecData, Rate, SensorCh, 
                                                        PeakDist, Lowpass, 
                                                        FilterOrder, Theta, 
                                                        Delta)
            
            AnalysisKey = Path + '/' + "{0:02d}".format(int(Rec))
            Hdf5F.WriteTreadmill(Treadmill, AnalysisKey, AnalysisFile)
            del(Treadmill)
            Done.append([Path, Rec])
        
        except Exception as E:
            Errors.append([Path, Rec])
            ErrorsLog.append(E)
            print(''); print('!!!==========')
            print(E)
            print('!!!=========='); print(''); 


#%% Plot
makedirs('Figs', exist_ok=True)
AnalysisFile = '201703-Treadmill.hdf5'
Paths = glob('2017-03-08*'); Paths.sort()
Exps = ['DMSO', 'CNO']
TDIndexes = {}
for Path in Paths:
    Treadmill = Hdf5F.LoadTreadmill(Path, AnalysisFile)
#    Recs = list(Treadmill.keys()).sort()
    Animal = Path.split('_')[-2:]
    Animal, Exp = Animal[1], Animal[0]
    FigName = 'Figs/' + Path; Ext = 'png'
    
    if Animal not in TDIndexes: TDIndexes[Animal] = {}
    if Exp not in TDIndexes[Animal]: TDIndexes[Animal][Exp] = {}
    BestCh = [0,0,0]
    for R, Rec in Treadmill.items():
        Max = max([Ch['TDIndex'] for C, Ch in Rec.items() if C != 'V'])
        TDIndexes[Animal][Exp][R] = [[C, Ch['TDIndex']] 
                                     for C, Ch in Rec.items()
                                     if C != 'V'
                                     if Ch['TDIndex'] == Max]
        
        if Max > BestCh[2]: BestCh = [R] + TDIndexes[Animal][Exp][R][0]
#        DataAnalysis.Plot.Treadmill_AllChs(Rec, FigName+'-AllChs-'+R+'.'+Ext, Ext, Save=True, Visible=False)
    
    TDIndexes[Animal][Exp]['BestCh'] = BestCh
    
    DataAnalysis.Plot.Treadmill_AllPerCh(Treadmill[BestCh[0]][BestCh[1]]['T'], 
                                         Treadmill[BestCh[0]][BestCh[1]]['F'], 
                                         Treadmill[BestCh[0]][BestCh[1]]['Sxx'], 
                                         Treadmill[BestCh[0]][BestCh[1]]['SxxMaxs'], 
                                         Treadmill[BestCh[0]][BestCh[1]]['SxxPerV'], 
                                         Treadmill[BestCh[0]][BestCh[1]]['VMeans'], 
                                         Treadmill[BestCh[0]][BestCh[1]]['VMeansSorted'], 
                                         Treadmill[BestCh[0]][BestCh[1]]['VInds'], 
                                         FigName+'-BestCh.'+Ext, Ext, Save=True, Visible=False)

Ext = 'png'
for A, Animal in TDIndexes.items():
    R = Animal['DMSO']['BestCh'][0]; C = Animal['DMSO']['BestCh'][1]
    F = 'Figs/'+A+'-DMSOxCNO-'+R+'-'+C+'.'+Ext
    
    Paths = glob('2017-03-08*'+A); Paths.sort()
    Ch1 = Hdf5F.LoadDict(Paths[0]+'/'+R+'/'+C, AnalysisFile, False)
    Ch2 = Hdf5F.LoadDict(Paths[1]+'/'+R+'/'+C, AnalysisFile, False)
    
    DataAnalysis.Plot.Treadmill_ChPair(Ch1, Ch2, Exps[0], Exps[1], F, Ext, Save=True, Visible=False)

