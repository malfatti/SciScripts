#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: T. Malfatti
@date: 2017-11-23
@license: GNU GPLv3 <https://raw.githubusercontent.com/malfatti/SciScripts/master/LICENSE>
@homepage: https://github.com/Malfatti/SciScripts
"""
from datetime import datetime
import sounddevice as SD

from IO import SigGen, Txt


## Level 0
def AudioSet(Rate, BlockSize, Channels, Intensities, NoiseFrequency, SoundPulseDur, SoundPulseNo, System, SoundAmpF, **Kws):
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
    
    DataInfo = {'InfoFile': InfoFile, 'Animal': {}, 'Audio': {}}
    for K in ['AnimalName', 'StimType']: DataInfo['Animal'][K] = locals()[K]
    
    for K in [
        'Rate', 'BlockSize', 'Channels', 'Intensities', 'NoiseFrequency', 
        'SoundPulseDur', 'SoundPulseNo', 'CalibrationFile', 'System', 'Setup'
    ]:
        DataInfo['Audio'][K] = locals()[K]
    
    DataInfo['Audio']['SoundAmpF'] = {K: Key for K, Key in SoundAmpF.items()}
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
    DataInfo = InfoWrite(AnimalName, StimType, Rate, BlockSize, Channels, Intensities, NoiseFrequency, SoundPulseDur, SoundPulseNo, SoundAmpF, System, Setup, InfoFile)
    Sound, Stim = AudioSet(**DataInfo['Audio'])
    Play(Sound, Stim, Intensities, SoundPulseNo, DataInfo)
    
    print('Animal', DataInfo['AnimalName'], 'successfully traumatized :)')

