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
def AudioSet(Rate, Intensities, NoiseFrequency, SoundPulseDur, SoundPulseNo, System, SoundAmpF):
    # Set sound stimulation
    Sound = SigGen.SoundStim(Rate, SoundPulseDur, SoundAmpF, NoiseFrequency, 
                             0, System, TTLs=False, Map=[2,1])
    
    # Set audio objects
    SD.default.device = 'system'
    SD.default.samplerate = Rate
    SD.default.blocksize = 384
    SD.default.channels = 2
    Stim = SD.OutputStream(dtype='float32')
    
    return(Sound, Stim)


def InfoWrite(AnimalName, StimType, Rate, Intensities, NoiseFrequency, SoundPulseDur, SoundPulseNo, SoundAmpF, InfoFile):
    DataInfo = dict((Name, eval(Name)) 
                    for Name in ['AnimalName', 'Rate', 'SoundPulseDur', 
                                 'SoundPulseNo', 'Intensities', 'NoiseFrequency', 
                                 'SigGen.CalibrationFile', 'StimType', 'InfoFile'])
    
#    Hdf5.DictWrite(DataInfo, '/DataInfo', InfoFile)
#    Hdf5.DictWrite(SoundAmpF, '/DataInfo/SoundAmpF', InfoFile)
    DataInfo['ExpInfo'] = {}
    DataInfo['SoundAmpF'] = {K: Key.tolist() for K, Key in SoundAmpF.items()}
    Txt.DictWrite(InfoFile.split('.')[0]+'.dict', DataInfo)
    
    return(DataInfo)


def Play(Sound, Stim, Intensities, SoundPulseDur, SoundPulseNo, StimType, DataInfo):
    FKey = list(Sound.keys())[0]
    AKeys = list(Sound[FKey].keys()); AKeys = sorted(AKeys, reverse=True)
    
    Stim.start()
    for AmpF, AKey in enumerate(AKeys):
    #        SS = Sound[FKey][AKey].T
    #        for Pulse in range(SoundPulseNo-1):
    #            SS = np.concatenate((SS, Sound[FKey][AKey].T))
    #    SS = np.concatenate([Sound[FKey][AKey] for _ in range(SoundPulseNo)])
        
        print('Playing', FKey, 'at', str(Intensities[AmpF]), 'dB')
        for Pulse in range(SoundPulseNo): Stim.write(Sound[FKey][AKey])
    #    del(SS)
    Stim.stop()
    
#    Hdf5.ExpInfoWrite('Sound', None, FKey, InfoFile)
    DataInfo['ExpInfo']['0'] = {'DVCoord': None, 'StimType': StimType, 'Hz': FKey}
    Txt.DictWrite(DataInfo['InfoFile'].split('.')[0]+'.dict', DataInfo)
    
    return(None)


## Level 1
def Run(AnimalName, StimType, Intensities, NoiseFrequency, SoundPulseDur, System, Setup):
    SoundPulseNo = round((SoundPulseDur*60)/20)
    SoundPulseDur = 20
    Rate = 192000
    
    Date = datetime.datetime.now()
    InfoFile = ''.join([Date.strftime("%Y%m%d%H%M%S"), '-', AnimalName, 
                        '-SoundStim.dict'])
    
    SoundAmpF = SigGen.dBToAmpF(Intensities, System+'/'+Setup)
    Sound, Stim = AudioSet(Rate, Intensities, NoiseFrequency, SoundPulseDur, SoundPulseNo, System, SoundAmpF)
    DataInfo = InfoWrite(AnimalName, StimType, Rate, Intensities, NoiseFrequency, SoundPulseDur, SoundPulseNo, SoundAmpF, InfoFile)
    Play(Sound, Stim, Intensities, SoundPulseDur, SoundPulseNo, StimType, DataInfo)
    
    print('Animal', DataInfo['AnimalName'], 'successfully traumatized :)')

