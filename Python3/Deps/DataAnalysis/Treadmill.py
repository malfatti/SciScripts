#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jun 12 14:12:37 2017

@author: malfatti
"""
import numpy as np

from DataAnalysis.DataAnalysis import FilterSignal, FindPeaks, QuantifyTTLsPerRec, Spectrogram
from scipy.interpolate import interp1d

def Analysis(Data, Rate, SensorCh, PeakDist, Lowpass=5, FilterOrder=2, Theta=[7, 10], Delta=[2, 4]):
    SensorData = Data[:,SensorCh-1]*-1
    SensorData = FilterSignal(SensorData, Rate, [Lowpass], FilterOrder, 
                              'butter', 'lowpass')
    
    Peaks = QuantifyTTLsPerRec(True, SensorData)
    
    V = np.zeros(len(SensorData))
    for P in range(1, len(Peaks)):
        Samples = Peaks[P] - Peaks[P-1]; Time = Samples/Rate
        Speed = PeakDist/Time; V[Peaks[P-1]:Peaks[P]] = [Speed]*Samples
    
    VInd, VPeaks = FindPeaks(V, Rate*3)
    f = interp1d(VInd, VPeaks, fill_value=0.0, bounds_error=False)
    V = f(np.arange(len(V))); V[V!=V] = 0.0
    
    Treadmill = {}
#        Treadmill = {'V': V}
    for C in range(Data.shape[1]-1):
        Ch = "{0:02d}".format(C+1); Treadmill[Ch] = {}
        print('Processing Ch', Ch, '...')
        
        F, T, Sxx = Spectrogram(Data[:,C], Rate, max(Theta))
        
        VMeans = [np.mean(V[int(T[t]*Rate):int(T[t+1]*Rate)]) 
                  for t in range(len(T)-1)] + [0.0]
        VMeansSorted = sorted(VMeans)
        VInds = [VMeansSorted.index(v) for v in VMeans]
        SxxPerV = Sxx[:,VInds]
        
        ThetaMaxs = [max(Sxx[(F>=min(Theta))*(F<max(Theta)),t]) for t in range(len(T))]
        DeltaMaxs = [max(Sxx[(F>=min(Delta))*(F<max(Delta)),t]) for t in range(len(T))]
        SxxMaxs = np.array(ThetaMaxs)/np.array(DeltaMaxs)
        
        Start = np.where(T>Peaks[0]/Rate)[0][0]
        End = np.where(T<Peaks[-1]/Rate)[0][-1]
        TDIndex = max(ThetaMaxs[Start:End])/max(DeltaMaxs[Start:End])
#            SxxSR = 1/(T[1]-T[0]); SxxLowPass = SxxSR*5/100
#            SxxMeans = DataAnalysis.FilterSignal(SxxMeans, SxxSR, 
#                                                  [SxxLowPass], 1, 'butter', 
#                                                    'lowpass')
        
        Treadmill[Ch] = {
            'F': F, 'T': T, 'Sxx': Sxx, 'SxxMaxs': SxxMaxs, 
            'SxxPerV': SxxPerV, 'TDIndex': TDIndex, 'VMeans': VMeans, 
            'VMeansSorted': VMeansSorted, 'VInds':VInds
        }
    
    print('Done.')
    return(Treadmill)

