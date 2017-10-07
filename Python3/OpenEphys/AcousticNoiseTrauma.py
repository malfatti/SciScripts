#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Oct  4 11:43:58 2017

@author: malfatti
"""

#%% Acoustic trauma
import datetime, os
import numpy as np
import sounddevice as SD

from IO import Hdf5, SigGen, Txt


## Experiment parameters
AnimalName = 'Prevention_A3_4_5'

# Fill all durations in SECONDS!
SoundPulseDur = 20
SoundPulseNo = 360
Intensities = [95]
PauseBetweenIntensities=0
NoiseFrequency = [[9000, 11000]]
StimType = ['Sound']


## Hardware parameters
System = 'Jack-IntelOut-MackieIn-MackieOut-IntelIn'
Setup = 'GPIAS'
Rate = 192000

# Set sound stimulation
CalibrationFile = os.environ['DATAPATH']+'/Tests/SoundMeasurements/SoundMeasurements.hdf5'
SoundAmpF = SigGen.dBToAmpF(Intensities, CalibrationFile, System+'/'+Setup)
Sound = SigGen.SoundStim(Rate, SoundPulseDur, SoundAmpF, NoiseFrequency, 
                         0, System, TTLs=False, Map=[2,1])
Pause = np.zeros((PauseBetweenIntensities*Rate,2), dtype='float32')

# Set audio objects
SD.default.device = 'system'
SD.default.samplerate = Rate
SD.default.blocksize = 384
SD.default.channels = 2


## Write info
Date = datetime.datetime.now()
FileName = ''.join([Date.strftime("%Y%m%d%H%M%S"), '-', AnimalName, 
                    '-SoundStim.hdf5'])

DataInfo = dict((Name, eval(Name)) 
                for Name in ['AnimalName', 'Rate', 'SoundPulseDur', 
                             'SoundPulseNo', 'Intensities', 'NoiseFrequency', 
                             'CalibrationFile', 'FileName'])

Hdf5.DictWrite(DataInfo, '/DataInfo', FileName)
Hdf5.DictWrite(SoundAmpF, '/DataInfo/SoundAmpF', FileName)
DataInfo['ExpInfo'] = {}
DataInfo['SoundAmpF'] = {K: Key.tolist() for K, Key in SoundAmpF.items()}
Txt.DictWrite(FileName.split('.')[0]+'.dict', DataInfo)


#%% Run sound
FKey = list(Sound.keys())[0]
AKeys = list(Sound[FKey].keys()); AKeys = sorted(AKeys, reverse=True)
for AmpF, AKey in enumerate(AKeys):
#        SS = Sound[FKey][AKey].T
#        for Pulse in range(SoundPulseNo-1):
#            SS = np.concatenate((SS, Sound[FKey][AKey].T))
    SS = np.concatenate([Sound[FKey][AKey] for _ in range(SoundPulseNo)])
    
    print('Playing', FKey, 'at', str(Intensities[AmpF]), 'dB')
    SD.play(SS, blocking=True)
    del(SS)

Hdf5.ExpInfoWrite('Sound', None, FKey, FileName)
DataInfo['ExpInfo']['0'] = {'DVCoord': None, 'StimType': StimType, 'Hz': FKey}
Txt.DictWrite(FileName.split('.')[0]+'.dict', DataInfo)

print('Animal', AnimalName, 'successfully traumatized :)')

