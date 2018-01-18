#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: T. Malfatti <malfatti@disroot.org>
@year: 2017
@license: GNU GPLv3 <https://raw.githubusercontent.com/malfatti/SciScripts/master/LICENSE>
@homepage: https://github.com/Malfatti/SciScripts
"""
import numpy as np
import os

from DataAnalysis import DataAnalysis
from IO import Hdf5, OpenEphys, Txt

NoTTLsFile = os.environ['NEBULAPATH']+'/Documents/PhD/Notes/Experiments/NoTTLRecs.dict'

## Level 0
def _PrintOptions():
    print('')
    print('This plot shows the first 0.5s of this recording.')
    print('What to do?')
    print('1) Decrease TTL Threshold')
    print('2) Change TTL channel')
    print('3) Decrease TTL Threshold and change TTL channel')
    print('4) Choose a startpoint for TTLs')
    print('5) Plot TTLCh and TTL threshold')
    print('6) Change ABR and TTL channels')
    print('7) Ignore PulseNo restritions')
    print('8) Abort')
    print('')
    
    return(None)


def ABRPerCh(Rec, R, Rate, ABRCh, TTLCh, Info, DataInfo, NoOfSamples, NoOfSamplesBefore, NoOfSamplesAfter, FilterFreq, FilterOrder, AnalogTTLs=True, Proc=None, Events=None):
    ABR = np.zeros((NoOfSamples, len(ABRCh)))
    for C, Ch in enumerate(ABRCh):
        if AnalogTTLs:
            # Temporary override
#                if 4.5 < np.mean(Rec[:,TTLCh-1]) < 5.5: TTLCh = 21
            if TTLCh > Rec.shape[1]:
                print('TTLCh > ChNo; replaced by -1.')
                TTLCh = 0
            
            print('TTL mean and max:', np.mean(Rec[:,TTLCh-1]), np.max(Rec[:,TTLCh-1]))
            TTLs = DataAnalysis.QuantifyTTLsPerRec(AnalogTTLs, Rec[:,TTLCh-1])
            if len(TTLs) > DataInfo['Audio']['SoundPulseNo']+4 or len(TTLs) < 100:
                TTLs = []
            
            print(len(TTLs), 'TTLs')
            if max(ABRCh) > Rec.shape[1]:
                print(''); print('Wrong ABR Channels!'); print('')
                TTLs = FixTTLs(Rec, ABRCh, TTLCh, Rate, DataInfo, Info, R, AnalogTTLs)
            
            # if type(TTLs) == dict: return(TTLs)
            
            while len(TTLs) == 0:
                TTLs = FixTTLs(Rec, ABRCh, TTLCh, Rate, DataInfo, Info, R, AnalogTTLs)
            
            if type(TTLs) == dict: return(TTLs)
        
        else:
            TTLs = DataAnalysis.QuantifyTTLsPerRec(AnalogTTLs, EventsDict=Events, TTLCh=TTLCh, 
                                      Proc=Proc, Rec=R)
        
        ABRData = DataAnalysis.SliceData(Rec[:, Ch-1], TTLs, NoOfSamplesBefore, 
                        NoOfSamplesAfter, NoOfSamples, AnalogTTLs)
        
        ABR[:,C] = np.mean(ABRData, axis=0)
        ABR[:,C] = DataAnalysis.FilterSignal(ABR[:,C], Rate, FilterFreq, FilterOrder,  
                                'butter', 'bandpass')
    
    return(ABR)


def FixTTLs(Data, DataCh, TTLCh, Rate, DataInfo, Info, Rec, AnalogTTLs):
    NoTTLRecs = Txt.DictRead(NoTTLsFile)
    F = DataInfo['InfoFile'].split('/')[-1].split('.')[0]
    
    WhatToDo = None; IgnoreRestritions = False
    if F in NoTTLRecs:
        if Info['Frequency'] in NoTTLRecs[F]:
            if Rec in NoTTLRecs[F][Info['Frequency']]:
                WhatToDo = NoTTLRecs[F][Info['Frequency']][Rec]['WhatToDo']
                Std = NoTTLRecs[F][Info['Frequency']][Rec]['Std']
                SampleStart = NoTTLRecs[F][Info['Frequency']][Rec]['SampleStart']
                
                if NoTTLRecs[F][Info['Frequency']][Rec]['TTLCh']:
                    TTLCh = NoTTLRecs[F][Info['Frequency']][Rec]['TTLCh']
                
                if WhatToDo == '6': 
                    DataCh = NoTTLRecs[F][Info['Frequency']][Rec]['DataCh']
                    ChangedDataCh = True
                
    
    if not WhatToDo:
        if 'Plot' not in globals():
            from DataAnalysis.Plot import Plot
        
        print('DataCh:', DataCh)
        print('TTLCh:', TTLCh)
        print('')
        _PrintOptions()
        
        Chs = [Data[:int(Rate/2),Ch] for Ch in range(Data.shape[1])]
        
        if DataInfo['InfoFile'].split('/')[-1].split('-')[0] in ['20170623145925', '20170623162416']:
            Chs[-1] = DataAnalysis.FilterSignal(Chs[-1], Rate, [8000, 14000])
        
        Plot.RawCh(Chs, Lines=len(Chs), Cols=1, Save=False)
        WhatToDo = input(': ')
    
        if F not in NoTTLRecs: NoTTLRecs[F] = {}
        if Info['Frequency'] not in NoTTLRecs[F]: NoTTLRecs[F][Info['Frequency']] = {}
        if Rec not in NoTTLRecs[F][Info['Frequency']]: NoTTLRecs[F][Info['Frequency']][Rec] = {}
        
        NoTTLRecs[F][Info['Frequency']][Rec]['WhatToDo'] = WhatToDo
        for K in ['SampleStart', 'Std', 'TTLCh']:
            NoTTLRecs[F][Info['Frequency']][Rec][K] = None
        
        Std = NoTTLRecs[F][Info['Frequency']][Rec]['Std']
        SampleStart = NoTTLRecs[F][Info['Frequency']][Rec]['SampleStart']
        ChangedDataCh = False
    
    
    if WhatToDo == '1':
        if not Std:
            Std = input('How many std above the mean should the threshold be? ')
            Std = int(Std)
        
        TTLs = DataAnalysis.QuantifyTTLsPerRec(AnalogTTLs, Data[:,TTLCh-1], Std)
        
        NoTTLRecs[F][Info['Frequency']][Rec]['Std'] = Std
    
    elif WhatToDo == '2':
        if not TTLCh:
            TTLCh = input('TTL channel: '); TTLCh = int(TTLCh)
            
        TTLs = DataAnalysis.QuantifyTTLsPerRec(AnalogTTLs, Data[:,TTLCh-1])
        
        NoTTLRecs[F][Info['Frequency']][Rec]['TTLCh'] = TTLCh
    
    elif WhatToDo == '3':
        if not Std:
            Std = input('How many std above the mean should the threshold be? ')
            Std = float(Std)
        
        if not TTLCh:
            TTLCh = input('TTL channel: '); TTLCh = int(TTLCh)
            TTLs = DataAnalysis.QuantifyTTLsPerRec(AnalogTTLs, Data[:,TTLCh-1], Std)
        
        NoTTLRecs[F][Info['Frequency']][Rec]['Std'] = Std
        NoTTLRecs[F][Info['Frequency']][Rec]['TTLCh'] = TTLCh
    
    elif WhatToDo == '4':
        if not SampleStart:
            SampleStart = input('Sample of the rising edge of the 1st pulse: ')
            SampleStart = int(SampleStart)
        
        TTLs = DataAnalysis.GenFakeTTLsRising(Rate, DataInfo['Audio']['SoundPulseDur'], 
                                 DataInfo['Audio']['SoundPauseBeforePulseDur'], 
                                 DataInfo['Audio']['SoundPauseAfterPulseDur'], SampleStart, 
                                 DataInfo['Audio']['SoundPulseNo'])
        
        NoTTLRecs[F][Info['Frequency']][Rec]['SampleStart'] = SampleStart
        
    elif WhatToDo == '5':
        Threshold = DataAnalysis.GetTTLThreshold(Data[:int(Rate/2),TTLCh-1], Std)
        Plot.TTLCh(Data[:,TTLCh-1], Std, Threshold)
        TTLs = []
    
    elif WhatToDo == '6':
        if not ChangedDataCh:
            NewABRCh = input('ABR channels (comma separated): ')
            NewABRCh = [int(_) for _ in NewABRCh.split(',')]
            
            TTLCh = input('TTL channel: '); TTLCh = int(TTLCh)
        else:
            NewABRCh = DataCh
        
        TTLs = DataAnalysis.QuantifyTTLsPerRec(AnalogTTLs, Data[:,TTLCh-1])
        
        NoTTLRecs[F][Info['Frequency']][Rec]['DataCh'] = NewABRCh
        NoTTLRecs[F][Info['Frequency']][Rec]['TTLCh'] = TTLCh
        
        Info = {'ABRCh': NewABRCh, 'TTLCh': TTLCh}
        
        print('ABR Channels changed. Restarting calculations...')
        print('')
        
    
    elif WhatToDo == '7':
        TTLs = DataAnalysis.QuantifyTTLsPerRec(AnalogTTLs, Data[:,TTLCh-1])
        IgnoreRestritions = True
    
    else:
        print('Aborted.')
        raise(SystemExit)
    
    if not IgnoreRestritions:
        if len(TTLs) > DataInfo['Audio']['SoundPulseNo']+4 or len(TTLs) < 100:
            TTLs = []
    
    if len(TTLs) != 0: Txt.DictWrite(NoTTLsFile, NoTTLRecs)
    print(len(TTLs), 'TTLs')
    if WhatToDo == '6': return(Info)
    else: return(TTLs)


def Calc(Data, Rate, ExpInfo, DataInfo, Stim='', ABRCh=[1], TTLCh=0,
         TimeBeforeTTL=3, TimeAfterTTL=12, AnalogTTLs=True, 
         FilterFreq=[300, 3000], FilterOrder=5, Proc=None, Events=None):
    
    ABRs = {}; Info = {}
    
    NoOfSamplesBefore = TimeBeforeTTL*int(Rate*10**-3)
    NoOfSamplesAfter = TimeAfterTTL*int(Rate*10**-3)
    NoOfSamples = NoOfSamplesBefore + NoOfSamplesAfter
    
    if type(ExpInfo['Hz']) == str: Info['Frequency'] = ExpInfo['Hz']
    else:
        Freq = DataInfo['Audio']['NoiseFrequency'][ExpInfo['Hz']]
        Info['Frequency'] = '-'.join([str(_) for _ in Freq])
    
    Info['XValues'] = (range(-NoOfSamplesBefore, 
                             NoOfSamples-NoOfSamplesBefore)/Rate)*10**3
    
    if len(Data.keys()) > len(DataInfo['Audio']['Intensities']):
        print('There are more recs than Intensities. Skipping...')
        return({}, {})
    
    Done = []
    for R, Rec in Data.items():
        if len(Rec) < 20*Rate: print('Rec', R, 'is broken!!!'); continue
        
        print('Slicing and filtering ABRs Rec ', R, '...')
        Start = np.where(Rec[:,0] != 0)[0]
        if Start.size: Start = Start[0]
        else: print('Rec', R, 'is broken!!!'); continue
        
        End = np.where(Rec[:,0] != 0)[0][-1]
        Rec = Rec[Start:End,:]
        
        ABR = {'ABRCh': ABRCh, 'TTLCh': TTLCh}
        while type(ABR) == dict:
            ABR = ABRPerCh(Rec, R, Rate, ABR['ABRCh'], ABR['TTLCh'], Info, DataInfo, NoOfSamples, NoOfSamplesBefore, NoOfSamplesAfter, FilterFreq, FilterOrder, AnalogTTLs)
        
        RealR = len(Done)
        dB = str(DataInfo['Audio']['Intensities'][int(RealR)]) + 'dB'
        
#        if len(DataInfo['Audio']['Intensities'])-1 < int(R): R = Broken[0]; del(Broken[0])
#        if Broken: R = str(len(Good))
#        print(R)
#         try:
#             dB = str(DataInfo['Audio']['Intensities'][int(R)]) + 'dB'
#         except IndexError:
#             print('Recs and Intensities indexes do not match.')
#             print('Intensities are:')
#             for I, Int in enumerate(DataInfo['Audio']['Intensities']): print(str(I) + ')' , Int)
#             print('')
#             print('Recs are:')
#             for Ri, Rk in enumerate(sorted(list(Data.keys()))): print(str(Ri) + ')' , Rk)
#             print('')
#             print('Current Rec is', R)
#             R = input('Choose intensity (using above intensity index): ')
#             dB = str(DataInfo['Audio']['Intensities'][int(R)]) + 'dB'
        ABRs[dB] = ABR[:]#; del(ABR, ABRData)
        Done.append(R)
#        Good.append(R)
        
#    if Broken: print("There were broken recs, intensities may be wrong!")
    Info['Path'] = Stim+'/'+ExpInfo['DV']+'/'+Info['Frequency']
    
    print('')
    return(ABRs, Info)


## Level 1
def Analysis(Exp, Folders, InfoFile, AnalysisFile='', TimeBeforeTTL=3, 
             TimeAfterTTL=12, FilterFreq=[300, 3000], FilterOrder=4, 
             StimType=['Sound'], AnalogTTLs=True, Proc='100'):
    InfoType = InfoFile[-4:]
    
    # ABRDict = {K}
    
    if InfoType == 'hdf5': DataInfo = Hdf5.DictLoad('/DataInfo', InfoFile)
    else: DataInfo = Txt.DictRead(InfoFile)
    
    if not AnalysisFile: 
        AnalysisFile = './' + DataInfo['Animal']['AnimalName'] + '-Analysis.hdf5'
        
    for Stim in StimType:
        if InfoType == 'hdf5': Exps = Hdf5.ExpPerStimLoad(Stim, Folders, InfoFile)
        else: 
            Exps = [Folders[int(R)] for R in DataInfo['ExpInfo']
                    if Stim in DataInfo['ExpInfo'][R]['StimType']]
            Exps.sort()
        
        Freqs = []; Trial = 0
        for F, Folder in enumerate(Exps):
            print(''); print(Folder); print('')
            Data, Rate = OpenEphys.DataLoader(Folder, 'uV', [], AnalogTTLs)
            if len(Data.keys()) == 1: Proc = list(Data.keys())[0]
            
            if InfoType == 'hdf5': ExpInfo = Hdf5.ExpExpInfo(Folder, F, InfoFile); R = ''
            else: 
                R = "{0:02d}".format(Folders.index(Folder))
                ExpInfo = DataInfo['ExpInfo'][R]
            
            if R: SStim = '_'.join(DataInfo['ExpInfo'][R]['StimType'])
            else: SStim = Stim
            
            Info = {'Calc': None, 
                    'ABRCh': DataInfo['DAqs']['ABRCh'], 
                    'TTLCh': DataInfo['DAqs']['TTLCh']}
            
            while 'Calc' in Info:
                ABRs, Info = Calc(Data[Proc], Rate[Proc], ExpInfo, DataInfo, 
                                  SStim, Info['ABRCh'], Info['TTLCh'], TimeBeforeTTL, 
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
            
            # DV = Info['Path'].split('/')[1]
            # Asdf.Write({'ABRs': ABRs, 'XValues': XValues}, AnalysisPath, AnalysisPath + '_' +Folder.split('/')[0]+'_ABRs.asdf')
            Freqs.append(Freq)
            print('')
    
    return(None)

