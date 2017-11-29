#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Nov 23 11:36:25 2017

@author: malfatti
"""


import numpy as np
import sounddevice as SD
import random

from datetime import datetime
from IO import Arduino, SigGen, Txt


def AudioSet(Rate, BlockSize, Channels):
    SD.default.device = 'system'
    SD.default.samplerate = Rate
    SD.default.blocksize = BlockSize
    SD.default.channels = Channels
    Stim = SD.OutputStream(dtype='float32')
    
    return(Stim)
    

def InfoWrite(AnimalName, StimType, BGIntensity, PulseIntensity, 
              NoiseFrequency, SoundBGDur, SoundGapDur, SoundBGPrePulseDur, 
              SoundLoudPulseDur, SoundBGAfterPulseDur, SoundBetweenStimDur, 
              NoOfTrials, SoundSystem, Setup, SoundBGAmpF, SoundPulseAmpF, SoundCh, 
              TTLCh, PiezoCh, Rate, BlockSize, Channels, BaudRate, InfoFile):
    
    CalibrationFile = SigGen.CalibrationFile
    
    DataInfo = {'InfoFile': InfoFile, 'Animal': {}, 'DAqs': {}, 'Audio': {}}
    for K in ['AnimalName', 'StimType']: DataInfo['Animal'][K] = locals()[K]
    
    for K in ['SoundCh', 'TTLCh', 'PiezoCh', 'BaudRate']: 
        DataInfo['DAqs'][K] = locals()[K]
    
    for K in ['BGIntensity', 'PulseIntensity', 'NoiseFrequency', 'SoundBGDur', 
              'SoundGapDur', 'SoundBGPrePulseDur', 'SoundLoudPulseDur', 
              'SoundBGAfterPulseDur', 'SoundBetweenStimDur', 'NoOfTrials', 
              'CalibrationFile', 'SoundSystem', 'Setup', 'Rate', 'BlockSize', 
              'Channels']:
        DataInfo['Audio'][K] = locals()[K]
    
    DataInfo['ExpInfo'] = {}
    # DataInfo['Audio']['SoundBGAmpF'] = {K: Key.tolist() for K, Key in SoundBGAmpF.items()}
    # DataInfo['Audio']['SoundPulseAmpF'] = {K: Key.tolist() for K, Key in SoundPulseAmpF.items()}
    DataInfo['Audio']['SoundBGAmpF'] = SoundBGAmpF
    DataInfo['Audio']['SoundPulseAmpF'] = SoundPulseAmpF
    
    Txt.DictWrite(InfoFile, DataInfo)
    
    return(DataInfo)


def Play(Sound, Stim, ArduinoObj, NoiseFrequency, SoundBetweenStimDur, NoOfTrials, SoundBGAmpF, SoundPulseAmpF, Rate, DataInfo):
    print('Preallocating memory and pseudo-randomizing the experiment...')
    FreqsStr = ['-'.join([str(a) for a in b]) for b in NoiseFrequency]
    TrialsStr = ['NoGap', 'Gap']
    Freqs = [In for In, El in enumerate(NoiseFrequency)]*NoOfTrials
    np.random.shuffle(Freqs)
    
    FreqSlot = [[0] for _ in range(len(Freqs)*2)]
    for FE in range(len(Freqs)):
        FreqSlot[FE*2] = (Freqs[0:FE+1].count(Freqs[FE])-1)*2
        FreqSlot[FE*2+1] = (Freqs[0:FE+1].count(Freqs[FE])-1)*2+1
    
    FreqOrder = [[0]]
    Rec = -1
    
    print("Running...")
    Stim.start()
    # Play the Pre-trials
    for Pre in range(3):
        Rec += 1
        RealFreq = FreqsStr[-1]
        FreqOrder[len(FreqOrder)-1] = [-1, -1]
        FreqOrder.append([0])
        ABGKey = str(SoundBGAmpF[RealFreq][0])
        APulseKey = str(SoundPulseAmpF[RealFreq][0])
        SBSDur = random.randrange(SoundBetweenStimDur[0], SoundBetweenStimDur[1])
        
        print('Playing ', RealFreq, ' Pre-trial', 'Rec', Rec)
        Stim.write(Sound['BetweenStim'][RealFreq][ABGKey][:SBSDur*Rate, :])        
        ArduinoObj.write(b'd')
        Stim.write(Sound['BG'][RealFreq][ABGKey])
        Stim.write(Sound['Gap']['NoGap'][RealFreq][ABGKey])
        Stim.write(Sound['BGPrePulse'][RealFreq][ABGKey])
        Stim.write(Sound['LoudPulse'][RealFreq][APulseKey])
        Stim.write(Sound['BGAfterPulse'][RealFreq][ABGKey])
        ArduinoObj.write(b'w')
    
    # Play the test trials
    for Hz in range(len(Freqs)):
        Trials = [0, 1]
        random.shuffle(Trials)
        print(str(Hz+1), 'of', str(len(Freqs)))
        
        for Trial in Trials:
            Rec += 1
            RealFreq = FreqsStr[Freqs[Hz]]; RealTrial = FreqSlot[Hz*2+Trial]
            FreqOrder[len(FreqOrder)-1] = [Freqs[Hz], RealTrial]
            FreqOrder.append([0])
            ABGKey = str(SoundBGAmpF[RealFreq][0])
            APulseKey = str(SoundPulseAmpF[RealFreq][0])
            SBSDur = random.randrange(SoundBetweenStimDur[0], SoundBetweenStimDur[1])
            
            print('Playing ', RealFreq, ' trial ', TrialsStr[Trial], 'Rec', Rec)
            Stim.write(Sound['BetweenStim'][RealFreq][ABGKey][:SBSDur*Rate, :])        
            ArduinoObj.write(b'd')
            Stim.write(Sound['BG'][RealFreq][ABGKey])
            Stim.write(Sound['Gap'][TrialsStr[Trial]][RealFreq][ABGKey])
            Stim.write(Sound['BGPrePulse'][RealFreq][ABGKey])
            Stim.write(Sound['LoudPulse'][RealFreq][APulseKey])
            Stim.write(Sound['BGAfterPulse'][RealFreq][ABGKey])
            ArduinoObj.write(b'w')
    
    # Play the Post-trials
    for Pre in range(3):
        Rec += 1
        RealFreq = FreqsStr[-1]
        FreqOrder[len(FreqOrder)-1] = [-1, -2]
        FreqOrder.append([0])
        ABGKey = str(SoundBGAmpF[RealFreq][0])
        APulseKey = str(SoundPulseAmpF[RealFreq][0])
        SBSDur = random.randrange(SoundBetweenStimDur[0], SoundBetweenStimDur[1])
        
        print('Playing ', RealFreq, ' Post-trial', 'Rec', Rec)
        Stim.write(Sound['BetweenStim'][RealFreq][ABGKey][:SBSDur*Rate, :])        
        ArduinoObj.write(b'd')
        Stim.write(Sound['BG'][RealFreq][ABGKey])
        Stim.write(Sound['Gap']['NoGap'][RealFreq][ABGKey])
        Stim.write(Sound['BGPrePulse'][RealFreq][ABGKey])
        Stim.write(Sound['LoudPulse'][RealFreq][APulseKey])
        Stim.write(Sound['BGAfterPulse'][RealFreq][ABGKey])
        ArduinoObj.write(b'w')
    
    Stim.stop()
    FreqOrder.remove([0])
    
    DataInfo['FreqOrder'] = FreqOrder
    DataInfo['FreqSlot'] = FreqSlot
    DataInfo['Freqs'] = Freqs
    
    Txt.DictWrite(DataInfo['InfoFile'], DataInfo)


def Run(AnimalName, StimType, BGIntensity, PulseIntensity, 
        NoiseFrequency, SoundBGDur, SoundGapDur, 
        SoundBGPrePulseDur, SoundLoudPulseDur, 
        SoundBGAfterPulseDur, SoundBetweenStimDur, NoOfTrials, 
        SoundSystem, Setup, SoundCh, TTLCh, PiezoCh, Rate=192000, BlockSize=384, 
        Channels=2, BaudRate=115200):
    
    SoundBGAmpF = SigGen.dBToAmpF(BGIntensity, SoundSystem+'/'+Setup)
    SoundPulseAmpF = SigGen.dBToAmpF(PulseIntensity, SoundSystem+'/'+Setup)
    
    Date = datetime.now().strftime("%Y%m%d%H%M%S")
    InfoFile = '-'.join([Date, AnimalName, 'GPIAS.dict'])
    
    DataInfo = InfoWrite(AnimalName, StimType, BGIntensity, PulseIntensity, 
                         NoiseFrequency, SoundBGDur, SoundGapDur, 
                         SoundBGPrePulseDur, SoundLoudPulseDur, 
                         SoundBGAfterPulseDur, SoundBetweenStimDur, NoOfTrials, 
                         SoundSystem, Setup, SoundBGAmpF, SoundPulseAmpF, 
                         SoundCh, TTLCh, PiezoCh, Rate, BlockSize, Channels, 
                         BaudRate, InfoFile)
    
    # TTLs Amplification factor. DO NOT CHANGE unless you know what you're doing.
    DataInfo['Audio']['TTLAmpF'] = 0.4
    
    Sound = SigGen.GPIAS(**DataInfo['Audio'])
    Stim = AudioSet(Rate, BlockSize, Channels)
    ArduinoObj = Arduino.CreateObj(BaudRate)
    
    Play(Sound, Stim, ArduinoObj, DataInfo=DataInfo, **DataInfo['Audio'])
    
    print('Finished', AnimalName, 'GPIAS.')

