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
from IO import Hdf5, OpenEphys, Txt


## Level 0
def Calc(Data, Rate, ExpInfo, DataInfo, Stim='', ABRCh=[1], TTLCh=0,
         TimeBeforeTTL=3, TimeAfterTTL=12, AnalogTTLs=True, 
         FilterFreq=[300, 3000], FilterOrder=5, TTLsPerRec=[]):
    ABRs = {}; Info = {}
    
    NoOfSamplesBefore = TimeBeforeTTL*int(Rate*10**-3)
    NoOfSamplesAfter = TimeAfterTTL*int(Rate*10**-3)
    NoOfSamples = NoOfSamplesBefore + NoOfSamplesAfter
    
    if type(ExpInfo['Hz']) == str:
        Info['Frequency'] = ExpInfo['Hz']
    else:
        Freq = DataInfo['NoiseFrequency'][ExpInfo['Hz']]
        Info['Frequency'] = '-'.join([str(_) for _ in Freq])
    
    Info['XValues'] = (range(-NoOfSamplesBefore, 
                             NoOfSamples-NoOfSamplesBefore)/Rate)*10**3
    
    if len(Data.keys()) > len(DataInfo['Intensities']):
        print('There are more recs than Intensities. Skipping...')
        return({}, {})
    
#    Broken = []; Good = []
    Done = []
    for R, Rec in Data.items():
        print('Slicing and filtering ABRs Rec ', R, '...')
        Start = np.where(Rec[:,0] != 0)[0][0]
        End = np.where(Rec[:,0] != 0)[0][-1]
        Rec = Rec[Start:End,:]
        
        if len(Rec) < 20*Rate:
            print('Rec', R, 'is broken!!!')#; Broken.append(R)
            continue
        
        ABR = np.zeros((NoOfSamples, len(ABRCh)))
        for C, Ch in enumerate(ABRCh):
            if AnalogTTLs:
                # Temporary override
                if 4.5 < np.mean(Rec[:,TTLCh-1]) < 5.5: TTLCh = 21
                print('TTL mean and max:', np.mean(Rec[:,TTLCh-1]), np.max(Rec[:,TTLCh-1]))
                TTLs = QuantifyTTLsPerRec(AnalogTTLs, Rec[:,TTLCh-1])
                print(len(TTLs), 'TTLs')
                ABRData = SliceData(Rec[:, Ch-1], TTLs, NoOfSamplesBefore, 
                                NoOfSamplesAfter, NoOfSamples, AnalogTTLs)
    #        else:
    #            RawTime, TTLs = QuantifyTTLsPerRec(Data, Rec, AnalogTTLs, 
    #                                               TTLsPerRec=TTLsPerRec)
    #            ABRData = SliceData(Data[Rec][:,Ch-1], TTLs,NoOfSamplesBefore, 
    #                            NoOfSamplesAfter, NoOfSamples, AnalogTTLs, RawTime)
            
            ABR[:,C] = np.mean(ABRData, axis=0)
            ABR[:,C] = FilterSignal(ABR[:,C], Rate, FilterFreq, FilterOrder,  
                                    'butter', 'bandpass')
        
        RealR = len(Done)
        dB = str(DataInfo['Intensities'][int(RealR)]) + 'dB'
        
#        if len(DataInfo['Intensities'])-1 < int(R): R = Broken[0]; del(Broken[0])
#        if Broken: R = str(len(Good))
#        print(R)
#         try:
#             dB = str(DataInfo['Intensities'][int(R)]) + 'dB'
#         except IndexError:
#             print('Recs and Intensities indexes do not match.')
#             print('Intensities are:')
#             for I, Int in enumerate(DataInfo['Intensities']): print(str(I) + ')' , Int)
#             print('')
#             print('Recs are:')
#             for Ri, Rk in enumerate(sorted(list(Data.keys()))): print(str(Ri) + ')' , Rk)
#             print('')
#             print('Current Rec is', R)
#             R = input('Choose intensity (using above intensity index): ')
#             dB = str(DataInfo['Intensities'][int(R)]) + 'dB'
        ABRs[dB] = ABR[:]; del(ABR, ABRData)
        Done.append(R)
#        Good.append(R)
        
#    if Broken: print("There were broken recs, intensities may be wrong!")
    Info['Path'] = Stim+'/'+ExpInfo['DVCoord']+'/'+Info['Frequency']
    
    return(ABRs, Info)


## Level 1
def Analysis(Exp, Folders, InfoFile, AnalysisFile='', TimeBeforeTTL=3, 
             TimeAfterTTL=12, FilterFreq=[300, 3000], FilterOrder=4, 
             StimType=['Sound'], AnalogTTLs=True, Proc='100'):
    InfoType = InfoFile[-4:]
    
    if InfoType == 'hdf5': DataInfo = Hdf5.DictLoad('/DataInfo', InfoFile)
    else: DataInfo = Txt.DictRead(InfoFile)
    
    if not AnalysisFile: 
        AnalysisFile = './' + DataInfo['AnimalName'] + '-Analysis.hdf5'
        
    for Stim in StimType:
        if InfoType == 'hdf5': Exps = Hdf5.ExpPerStimLoad(Stim, Folders, InfoFile)
        else: 
            Exps = [Folders[int(R)] for R in DataInfo['ExpInfo']
                    if Stim in DataInfo['ExpInfo'][R]['StimType']]
            Exps.sort()
        
        Freqs = []; Trial = 0
        for F, Folder in enumerate(Exps):
            Data, Rate = OpenEphys.DataLoader(Folder, AnalogTTLs)
            if len(Data.keys()) == 1: Proc = list(Data.keys())[0]
            
            if InfoType == 'hdf5': ExpInfo = Hdf5.ExpExpInfo(Folder, F, InfoFile); R = ''
            else: 
                R = "{0:02d}".format(Folders.index(Folder))
                ExpInfo = DataInfo['ExpInfo'][R]
            
            if R: SStim = '_'.join(DataInfo['ExpInfo'][R]['StimType'])
            else: SStim = Stim
            
            ABRs, Info = Calc(Data[Proc], Rate[Proc], ExpInfo, DataInfo, 
                              SStim, DataInfo['ABRCh'], DataInfo['TTLCh'], TimeBeforeTTL, 
                              TimeAfterTTL, AnalogTTLs, FilterFreq, 
                              FilterOrder)
            
            if not ABRs: print('Empty ABRs. Skipping...'); continue
            
            Freq = Info['Frequency']
            if Freq in Freqs: Freqs = []; Trial += 1
            AnalysisPath = '-'.join(InfoFile.split('/')[-1].split('-')[:2]) + '-ABRs'
            Path = '/'+AnalysisPath+'/ABRs'+'/'+Info['Path']+'/Trial'+str(Trial)
            XValuesPath = '/'+AnalysisPath+'/XValues'+'/'+Info['Path']+'/Trial'+str(Trial)
            print(Path)
            
            Hdf5.DataWrite(ABRs, Path, AnalysisFile, Overwrite=True)
            Hdf5.DataWrite(Info['XValues'], XValuesPath, AnalysisFile, Overwrite=True)
            Freqs.append(Freq)
    
    return(None)

