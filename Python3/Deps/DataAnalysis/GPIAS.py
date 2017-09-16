#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jun 12 14:12:37 2017

@author: malfatti
"""
import numpy as np

from DataAnalysis.DataAnalysis import FilterSignal, QuantifyTTLsPerRec, SliceData
from IO import Hdf5
from scipy import signal


## Level 0
def CheckGPIASRecs(Data, SizeLimits, Plot=False):
    ToCheck = [Rec for Rec in Data.keys() 
                   if len(Data[Rec])<min(SizeLimits)
                   or len(Data[Rec])>max(SizeLimits)]
    
    if ToCheck:
        if Plot:
            Params = {'backend': 'TkAgg'}
            from matplotlib import rcParams; rcParams.update(Params)
            import matplotlib.pyplot as plt
            
            for Rec in ToCheck:
                print('Showing Rec', Rec+', size', Data[Rec].shape[0])
                plt.plot(Data[Rec])
                plt.show()
            
        return(ToCheck)
    else:
        print('All recs within expected size.')
        return(None)

def IndexCalc(Data, Keys, PulseSampleStart, SliceSize):
    Index = {}
    for Key in Keys:
        BGStart = 0; BGEnd = SliceSize
        PulseStart = PulseSampleStart; PulseEnd = PulseSampleStart + SliceSize
        
        ResRMSBG = (np.mean(Data[Key[0]][BGStart:BGEnd]**2))**0.5
        ResRMSPulse = (np.mean(Data[Key[0]][PulseStart:PulseEnd]**2))**0.5
#            ResRMS = ResRMSPulse
        if ResRMSPulse < ResRMSBG: ResRMS = ResRMSPulse
        else: ResRMS = ResRMSPulse - ResRMSBG
        
        RefRMSBG = (np.mean(Data[Key[1]][BGStart:BGEnd]**2))**0.5
        RefRMSPulse = (np.mean(Data[Key[1]][PulseStart:PulseEnd]**2))**0.5
#            RefRMS = RefRMSPulse
        if RefRMSPulse < RefRMSBG: RefRMS = RefRMSPulse
        else: RefRMS = RefRMSPulse - RefRMSBG
        
        # GPIAS index (How much Res is different from Ref)
        Index[Key[2]] = (RefRMS-ResRMS)/RefRMS
    
    return(Index)


def PreallocateDict(DataInfo, PrePostFreq):
    Dict = {}
    for Key in ['Trace', 'Index']:
        Dict[Key] = {''.join([str(Freq[0]), '-', str(Freq[1])]): {}
                     for Freq in DataInfo['NoiseFrequency']}
    
    if PrePostFreq:
        Dict['Trace'][PrePostFreq]['Pre'] = []
        Dict['Trace'][PrePostFreq]['Post'] = []
        Dict['Index'][PrePostFreq]['Pre'] = []
        Dict['Index'][PrePostFreq]['Post'] = []
    
    for Freq in Dict['Trace'].keys():
        Dict['Trace'][Freq]['NoGap'] = []; Dict['Trace'][Freq]['Gap'] = []
        Dict['Index'][Freq]['NoGap'] = []; Dict['Index'][Freq]['Gap'] = []
    
    return(Dict)


def OrganizeRecs(Dict, Data, DataInfo, AnalogTTLs, NoOfSamplesBefore, 
                 NoOfSamplesAfter, NoOfSamples):
    for R, Rec in Data.items():
        print('Slicing and filtering Rec ', R, '...')
        Freq = DataInfo['FreqOrder'][int(R)][0]; 
        Trial = DataInfo['FreqOrder'][int(R)][1];
        
        SFreq = ''.join([str(DataInfo['NoiseFrequency'][Freq][0]), '-', 
                         str(DataInfo['NoiseFrequency'][Freq][1])])
        
        if Trial == -1: STrial = 'Pre'
        elif Trial == -2: STrial = 'Post'
        elif Trial % 2 == 0: STrial = 'NoGap'
        else: STrial = 'Gap'
        
        if AnalogTTLs:
            TTLs = QuantifyTTLsPerRec(AnalogTTLs, Rec[:,DataInfo['TTLCh']-1], StdNo=2)
            if len(TTLs) > 1: TTLs = [np.argmax(Rec[:,DataInfo['TTLCh']-1])]
            print(TTLs)
            
            if not TTLs: print('No TTL detected. Skipping trial...'); continue
            
            GD = SliceData(Rec[:,DataInfo['PiezoCh'][0]-1], TTLs, 
                           NoOfSamplesBefore, NoOfSamplesAfter, NoOfSamples, 
                           AnalogTTLs)
#        else:
#            RawTime, TTLs = QuantifyTTLsPerRec(Data, R, AnalogTTLs, 
#                                               TTLsPerRec=TTLsPerRec)
#            GD = SliceData(Raw, OEProc, R, TTLs, GPIASCh, NoOfSamplesBefore, 
#                           NoOfSamplesAfter, NoOfSamples, AnalogTTLs, RawTime)
        
        Dict['Index'][SFreq][STrial].append(GD[0])
        Dict['Trace'][SFreq][STrial].append(GD[0])
    
    return(Dict)
    

## Level 1    
def Analysis(Data, DataInfo, Rate, AnalysisFile, AnalysisKey, 
             GPIASTimeBeforeTTL=100, GPIASTimeAfterTTL=150, 
             FilterFreq=[70, 400], FilterOrder=4, Filter='butter', 
             SliceSize=100, AnalogTTLs=True, Return=False):
    
    NoOfSamplesBefore = int(round((GPIASTimeBeforeTTL*Rate)*10**-3))
    NoOfSamplesAfter = int(round((GPIASTimeAfterTTL*Rate)*10**-3))
    NoOfSamples = NoOfSamplesBefore + NoOfSamplesAfter
    
    XValues = (range(-NoOfSamplesBefore, NoOfSamples-NoOfSamplesBefore)
               /Rate)*10**3
    
    PrePostFreq = DataInfo['FreqOrder'][0][0]
    PrePostFreq = '-'.join([str(DataInfo['NoiseFrequency'][PrePostFreq][0]),
                            str(DataInfo['NoiseFrequency'][PrePostFreq][1])])
    
    # Temporary override
    if AnalysisFile.split('/')[0] == 'Recovery': PrePostFreq = []
    
    GPIASData = PreallocateDict(DataInfo, PrePostFreq)
    GPIASData = OrganizeRecs(GPIASData, Data, DataInfo, AnalogTTLs,
                                   NoOfSamplesBefore, NoOfSamplesAfter, 
                                   NoOfSamples)
    
    SliceSize = int(SliceSize * (Rate/1000))
    
    for Freq in GPIASData['Index'].keys():
        for Key in GPIASData['Index'][Freq].keys():
            # Average trials for traces
            GPIASData['Trace'][Freq][Key] = np.mean(GPIASData['Trace'][Freq][Key], axis=0)
            
            # Bandpass filter
            GPIASData['Trace'][Freq][Key] = FilterSignal(GPIASData['Trace'][Freq][Key], 
                                                       Rate, FilterFreq, FilterOrder, 
                                                       Filter, 'bandpass')
            
            for Tr in range(len(GPIASData['Index'][Freq][Key])):
                # Bandpass filter
#                    TR = FilterSignal(GPIASData['Trace'][Freq][Key][Tr], Rate, 
#                                      FilterFreq, FilterOrder, Filter, 'bandpass')
                
                AE = FilterSignal(GPIASData['Index'][Freq][Key][Tr], Rate, 
                                     FilterFreq, FilterOrder, Filter, 'bandpass')
                
                # Amplitude envelope
                AE = abs(signal.hilbert(AE))
                
#                    GPIASData['Trace'][Freq][Key][Tr] = TR
                GPIASData['Index'][Freq][Key][Tr] = AE
            
#                GPIASData['Trace'][Freq][Key] = np.mean(GPIASData['Trace'][Freq][Key], axis=0)
            GPIASData['Index'][Freq][Key] = np.mean(GPIASData['Index'][Freq][Key], axis=0)
            
        # RMS
        if Freq == PrePostFreq: Keys = [['Gap', 'NoGap', 'GPIASIndex'], 
                                        ['Post', 'Pre', 'PrePost']]
        else: Keys = [['Gap', 'NoGap', 'GPIASIndex']]
        
        GPIASData['Index'][Freq] = IndexCalc(
                                       GPIASData['Index'][Freq], Keys, 
                                       NoOfSamplesBefore, SliceSize)
    
    Hdf5.DataWrite({'GPIAS': GPIASData, 'XValues': XValues}, AnalysisKey, AnalysisFile)
    
    if Return: return(GPIASData, XValues)
    else: return(None)
