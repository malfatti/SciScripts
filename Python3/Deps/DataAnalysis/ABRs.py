#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
    Copyright (C) 2017  T. Malfatti
    
    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.
    
    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.
    
    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.

"""
import numpy as np

from DataAnalysis.DataAnalysis import FilterSignal, QuantifyTTLsPerRec, SliceData
from itertools import tee
from scipy import signal


## Level 0
def ABRAnalysis(FileName, ABRCh=[1], ABRTTLCh=1, ABRTimeBeforeTTL=0, 
                ABRTimeAfterTTL=12, FilterFreq=[300, 3000], FilterOrder=4, 
                StimType='Sound', AnalogTTLs=False, Board='OE', Type='kwik',
                Override={}):
    print('Load DataInfo...')
    DirList = glob('KwikFiles/*'); DirList.sort()
    DataInfo = Hdf5F.LoadDict('/DataInfo', FileName)
    
    AnalysisFile = './' + DataInfo['AnimalName'] + '-Analysis.hdf5'
    Now = datetime.now().strftime("%Y%m%d%H%M%S")
    Here = os.getcwd().split(sep='/')[-1]
    Group = Here + '-ABRs_' + Now
        
    for Stim in StimType:
        if Override != {}: 
            if 'Stim' in Override.keys(): Stim = Override['Stim']
        
        Exps = Hdf5F.LoadExpPerStim(Stim, DirList, FileName)
        
        for RecFolder in Exps:
            if AnalogTTLs: 
                Raw, _, Files = Hdf5F.LoadOEKwik(RecFolder, AnalogTTLs)
            else: 
                Raw, Events, _, Files = Hdf5F.LoadOEKwik(RecFolder, AnalogTTLs)
            
            ExpInfo = Hdf5F.ExpExpInfo(RecFolder, DirList, FileName)
            OEProc, RHAProc, ABRProc = Hdf5F.GetProc(Raw, Board)
            
            if AnalogTTLs: Raw = Hdf5F.GetRecKeys(Raw, [0], AnalogTTLs)
            else:
                Raw, EventRec = Hdf5F.GetRecKeys(Raw, Events, AnalogTTLs)
#                TTLsPerRec = GetTTLInfo(Events, EventRec, ABRTTLCh)
            
            Rate = Raw['100']['info']['0']['sample_rate']
            ABRs, Info = ABRCalc(Raw, Rate, ExpInfo)
            Hdf5F.WriteABR(ABRs, Info['XValues'], Group, Info['Path'], AnalysisFile)
        

#FileName, ABRCh=[1], ABRTTLCh=1, ABRTimeBeforeTTL=0, 
#                ABRTimeAfterTTL=12, FilterFreq=[300, 3000], FilterOrder=4, 
#                StimType='Sound', AnalogTTLs=False, Board='OE', Type='kwik',
#                Override={}
def ABRCalc(Data, Rate, ExpInfo, DataInfo, Stim, TimeBeforeTTL=3, 
            TimeAfterTTL=12, AnalogTTLs=True, ABRCh=5, TTLCh=17, 
            FilterFreq=[300, 3000], FilterOrder=5, TTLsPerRec=[]):
    ABRs = {}; Info = {}
    
    NoOfSamplesBefore = TimeBeforeTTL*int(Rate*10**-3)
    NoOfSamplesAfter = TimeAfterTTL*int(Rate*10**-3)
    NoOfSamples = NoOfSamplesBefore + NoOfSamplesAfter
    
    Info['Frequency'] = ExpInfo['Hz']
    
    Info['XValues'] = (range(-NoOfSamplesBefore, 
                             NoOfSamples-NoOfSamplesBefore)/Rate)*10**3
    
    for R, Rec in Data.items():
        print('Slicing and filtering ABRs Rec ', str(Rec), '...')
        
        if len(Rec) < 50*Rate:
            print('Rec', R, 'is broken!!!')
            continue
        
        if AnalogTTLs:
            TTLs = QuantifyTTLsPerRec(AnalogTTLs, Rec[:,TTLCh-1])
            ABR = SliceData(Rec[:, ABRCh-1], TTLs, NoOfSamplesBefore, 
                            NoOfSamplesAfter, NoOfSamples, AnalogTTLs)
#        else:
#            RawTime, TTLs = QuantifyTTLsPerRec(Data, Rec, AnalogTTLs, 
#                                               TTLsPerRec=TTLsPerRec)
#            ABR = SliceData(Data[Rec][:,ABRCh-1], TTLs,NoOfSamplesBefore, 
#                            NoOfSamplesAfter, NoOfSamples, AnalogTTLs, RawTime)
        
        for TTL in range(len(TTLs)):
            ABR[TTL] = FilterSignal(ABR[TTL], Rate, [min(FilterFreq)], 
                                    FilterOrder, 'butter', 'highpass')
        
        ABR = np.mean(ABR, axis=0)
        ABR = FilterSignal(ABR, Rate, [max(FilterFreq)], FilterOrder, 
                           'butter', 'lowpass')
        
        dB = str(DataInfo['Intensities'][int(R)]) + 'dB'
        ABRs[dB] = ABR[:]; del(ABR)
    
    Info['Path'] = Stim+'/'+ExpInfo['DVCoord']+'/'+Info['Frequency']
    
    return(ABRs, Info)

