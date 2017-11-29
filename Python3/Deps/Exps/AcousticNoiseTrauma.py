#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Nov 23 08:29:06 2017

@author: T. Malfatti
"""
import datetime
import sounddevice as SD

from IO import SigGen, Txt


## Level 0
def AudioSet(Rate, BlockSize, Channels, Intensities, NoiseFrequency, SoundPulseDur, SoundPulseNo, System, SoundAmpF):
    # Set sound stimulation
    Sound = SigGen.SoundStim(Rate, SoundPulseDur, SoundAmpF, NoiseFrequency, 
                             0, System, TTLs=False, Map=[2,1])
    
    # Set audio objects
    SD.default.device = 'system'
    SD.default.samplerate = Rate
    SD.default.blocksize = BlockSize
    SD.default.channels = Channels
    Stim = SD.OutputStream(dtype='float32')
    
    return(Sound, Stim)


def InfoWrite(AnimalName, StimType, Rate, BlockSize, Channels, Intensities, NoiseFrequency, SoundPulseDur, SoundPulseNo, SoundAmpF, System, Setup, InfoFile):
    CalibrationFile = SigGen.CalibrationFile
    
    DataInfo = {'InfoFile': InfoFile}
    DataInfo['Animal'] = {K: locals()[K] for K in ['AnimalName', 'StimType']}
    DataInfo['Audio'] = {K: locals()[K]
                             for K in ['Rate', 'BlockSize', 'Channels', 
                                       'Intensities', 'NoiseFrequency',
                                       'SoundPulseDur', 'SoundPulseNo', 
                                       'CalibrationFile', 'System', 'Setup']}
    
    DataInfo['ExpInfo'] = {}
    DataInfo['Audio']['SoundAmpF'] = {K: Key.tolist() 
                                          for K, Key in SoundAmpF.items()}
    
    Txt.DictWrite(InfoFile, DataInfo)
    
    return(DataInfo)


def Play(Sound, Stim, Intensities, SoundPulseNo, DataInfo):
    FKey = list(Sound.keys())[0]
    AKey = list(Sound[FKey].keys())[0]
    
    Stim.start()
    print('Playing', FKey, 'at', str(Intensities[0]), 'dB')
    for Pulse in range(SoundPulseNo): Stim.write(Sound[FKey][AKey])
    Stim.stop()
    
    DataInfo['ExpInfo']['0'] = {'DVCoord': None, 
                                'StimType': DataInfo['StimType'], 'Hz': FKey}
    
    Txt.DictWrite(DataInfo['InfoFile'], DataInfo)
    
    return(None)


## Level 1
def Run(AnimalName, StimType, Intensities, NoiseFrequency, SoundPulseDur, System, Setup, Rate=192000, BlockSize=384, Channels=2):
    SoundPulseNo = round((SoundPulseDur*60)/20)
    SoundPulseDur = 20
    
    Date = datetime.now().strftime("%Y%m%d%H%M%S")
    InfoFile = '-'.join([Date, AnimalName, 'AcousticNoiseTrauma.dict'])
    
    SoundAmpF = SigGen.dBToAmpF(Intensities, System+'/'+Setup)
    DataInfo = InfoWrite(AnimalName, StimType, Rate, BlockSize, Channels, Intensities, NoiseFrequency, SoundPulseDur, SoundPulseNo, SoundAmpF, InfoFile)
    Sound, Stim = AudioSet(**DataInfo['Audio'])
    Play(Sound, Stim, Intensities, SoundPulseNo, DataInfo)
    
    print('Animal', DataInfo['AnimalName'], 'successfully traumatized :)')

