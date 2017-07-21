# -*- coding: utf-8 -*-
"""
    Copyright (C) 2015  T. Malfatti
    
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

This is a script to generate pulses and send to the soundboard, and then to a 
sound amplifier and an Arduino board. Basically it generates sound pulses, 
sound square pulses (TTLs), and laser square pulses. The square pulses will be 
sent to the left channel and the sound pulses will be sent to the right 
channel. 
"""
#%% Import
import datetime, os
import numpy as np
import sounddevice as SD

from IO import Arduino, Hdf5, SigGen, Txt


#%% Set Parameters
# Order: [4, 3, 5]
AnimalName = 'Prevention_A3'
Rate = 192000
BaudRate = 115200

# Sound setup and system used
System = 'Jack-IntelOut-MackieIn-MackieOut-IntelIn'
Setup = 'UnitRec'

# TTLs Amplification factor. DO NOT CHANGE unless you know what you're doing.
TTLAmpF = 0.6

## Fill all durations in SECONDS!

## Sound
# Silence before pulse
SoundPrePauseDur = 0.004
# Pulse duration
SoundPulseDur = 0.003
# Silence after pulse
SoundPostPauseDur = 0.093
# Amount of pulses per block
SoundPulseNo = 529
# Intensities tested, in order, in dB. Supports floats :)
#Intensities = [80, 75, 70, 65, 60, 55, 50, 45, 40, 35]
Intensities = [80, 70, 60, 50, 40]
#Intensities = [80]
PauseBetweenIntensities = 10
# Noise frequency. If using one freq., keep the list in a list, [[like this]].
# USE ONLY FREQUENCY BANDS THAT WERE CALIBRATED. To check the calibrated freqs, 
# just run the cell once and then list(SoundIntensity).
NoiseFrequency = [[8000, 10000], [9000, 11000], [10000, 12000], 
                  [12000, 14000], [14000, 16000]]
#NoiseFrequency = [[8000, 10000], [12000, 14000]]


CalibrationFile = os.environ['DATAPATH']+'/Tests/SoundMeasurements/SoundMeasurements.hdf5'
Date = datetime.datetime.now()
FileName = ''.join([Date.strftime("%Y%m%d%H%M%S"), '-', AnimalName, 
                    '-SoundStim.hdf5'])

#SoundAmpF = SigGen.dBToAmpF(Intensities, CalibrationFile, System+'/'+Setup)
# Temporary override
SoundAmpF = Hdf5.DataLoad('/DataInfo/SoundAmpF', 'Prevention/20170704-Prevention_A1-ABRs/20170704102002-Prevention_A1-SoundStim.hdf5')[1]

DataInfo = dict((Name, eval(Name)) 
                for Name in ['AnimalName', 'Rate', 'BaudRate', 
                             'SoundPrePauseDur', 'SoundPulseDur', 
                             'SoundPostPauseDur', 'SoundPulseNo', 
                             'Intensities', 'NoiseFrequency', 
                             'PauseBetweenIntensities',
#                             'LaserPrePauseDur', 'LaserPulseDur', 
#                             'LaserPostPauseDur', 'LaserPulseNo', 
#                             'LaserStimBlockNo', 
#                             'LaserPauseBetweenStimBlocksDur', 
                             'CalibrationFile', 'FileName'])

Hdf5.DictWrite(DataInfo, '/DataInfo', FileName)
Hdf5.DictWrite(SoundAmpF, '/DataInfo/SoundAmpF', FileName)
DataInfo['ExpInfo'] = {}
DataInfo['SoundAmpF'] = {K: Key.tolist() for K, Key in SoundAmpF.items()}
Txt.DictWrite(FileName.split('.')[0]+'.dict', DataInfo)

ArduinoObj = Arduino.CreateObj(BaudRate)


#%% Prepare sound stimulation
Sound = SigGen.SoundStim(Rate, SoundPulseDur, SoundAmpF, NoiseFrequency, 
                         TTLAmpF, System, SoundPrePauseDur, 
                         SoundPostPauseDur)

Pause = np.zeros((PauseBetweenIntensities*Rate,2), dtype='float32')

# Set audio objects
SD.default.device = 'system'
SD.default.samplerate = Rate
SD.default.blocksize = 384
SD.default.channels = 2
Stim = SD.OutputStream(dtype='float32')
ArduinoObj.write(b'd')

#%% Run sound
DVCoord = 'Out'
StimType = ['Sound', 'CNO']
#Freq = 4
#Freq = int(Freq)

# FKeys = list(Sound.keys()); ToPrepend = []
# for FF in ['8000-10000', '9000-11000']:
#     if FF in FKeys:
#         del(FKeys[FKeys.index(FF)]); ToPrepend.append(FF)
# 
# ToPrepend.sort(); FKeys = ToPrepend + FKeys
FKeys = list(Sound.keys())
FKeys.sort(key=lambda x: [int(y) for y in x.split('-')])

# FreqOrder = ['8000-10000', '9000-11000', '12000-14000', '14000-16000', '10000-12000']
Stim.start()
while True:
    print('Remember to change folder name in OE!')
    print('Choose frequency:')
    for Ind, K in enumerate(FKeys): print(str(Ind) + ')' , K)
    print(str(len(FKeys)) + ')', 'Cancel')
    FKey = input(': ')
    
    if FKey == str(len(FKeys)): break
    
    try:
        FKey = FKeys[int(FKey)]
    except IndexError:
        print('=== Wrong Freq index. Stopping... ===')
        print('')
        break
    
    AKeys = list(Sound[FKey].keys()); AKeys = sorted(AKeys, reverse=True)
    for AmpF, AKey in enumerate(AKeys):
#        SS = Sound[FKey][AKey].T
#        for Pulse in range(SoundPulseNo-1):
#            SS = np.concatenate((SS, Sound[FKey][AKey].T))
        SS = np.concatenate([Sound[FKey][AKey] for _ in range(SoundPulseNo)])
        
        print('Playing', FKey, 'at', str(Intensities[AmpF]), 'dB')
        ArduinoObj.write(b'w')
        Stim.write(SS)
        ArduinoObj.write(b'd')
        Stim.write(Pause)
        del(SS)
    
    Hdf5.ExpInfoWrite('Sound', DVCoord, FKey, FileName)
    
    Rec = "{0:02d}".format(len(DataInfo['ExpInfo']))
    DataInfo['ExpInfo'][Rec] = {'DVCoord': DVCoord, 'StimType': StimType, 'Hz': FKey}
    Txt.DictWrite(FileName.split('.')[0]+'.dict', DataInfo)
    
    print('Played Freq', FKey, 'at', DVCoord, 'µm DV')

Stim.stop()

#%% Acoustic trauma

AnimalName = 'Prevention_A3_4_5'
Rate = 192000

# Sound setup and system used
System = 'Jack-IntelOut-MackieIn-MackieOut-IntelIn'
Setup = 'GPIAS'

## Sound
SoundPulseDur = 20
SoundPulseNo = 360
Intensities = [80]
PauseBetweenIntensities=0
NoiseFrequency = [[9000, 11000]]


CalibrationFile = os.environ['DATAPATH']+'/Tests/SoundMeasurements/SoundMeasurements.hdf5'
Date = datetime.datetime.now()
FileName = ''.join([Date.strftime("%Y%m%d%H%M%S"), '-', AnimalName, 
                    '-SoundStim.hdf5'])

SoundAmpF = SigGen.dBToAmpF(Intensities, CalibrationFile, System+'/'+Setup)

DataInfo = dict((Name, eval(Name)) 
                for Name in ['AnimalName', 'Rate', 'SoundPulseDur', 
                             'SoundPulseNo', 'Intensities', 'NoiseFrequency', 
                             'CalibrationFile', 'FileName'])

Hdf5.DictWrite(DataInfo, '/DataInfo', FileName)
Hdf5.DictWrite(SoundAmpF, '/DataInfo/SoundAmpF', FileName)


Sound = SigGen.SoundStim(Rate, SoundPulseDur, SoundAmpF, NoiseFrequency, 
                         0, System, TTLs=False, Map=[2,1])

Pause = np.zeros((PauseBetweenIntensities*Rate,2), dtype='float32')

# Set audio objects
SD.default.device = 'system'
SD.default.samplerate = Rate
SD.default.blocksize = 384
SD.default.channels = 2

DVCoord = 'Out'
#Freq = 4
#Freq = int(Freq)

FKeys = list(Sound.keys()); ToPrepend = []
for FF in ['8000-10000', '9000-11000']:
    if FF in FKeys:
        del(FKeys[FKeys.index(FF)]); ToPrepend.append(FF)

ToPrepend.sort(); FKeys = ToPrepend + FKeys

while True:
    print('Remember to change folder name in OE!')
    print('Choose frequency:')
    for Ind, K in enumerate(FKeys):
        print(str(Ind) + ')' , K)
    
    print(str(len(FKeys)) + ')', 'Cancel')
    FKey = input(': ')
    
    if FKey == str(len(FKeys)): break
    
    try:
        FKey = FKeys[int(FKey)]
    except IndexError:
        print('=== Wrong Freq index. Stopping... ===')
        print('')
        break
    
    AKeys = list(Sound[FKey].keys()); AKeys = sorted(AKeys, reverse=True)
    for AmpF, AKey in enumerate(AKeys):
#        SS = Sound[FKey][AKey].T
#        for Pulse in range(SoundPulseNo-1):
#            SS = np.concatenate((SS, Sound[FKey][AKey].T))
        SS = np.concatenate([Sound[FKey][AKey] for _ in range(SoundPulseNo)])
        
        print('Playing', FKey, 'at', str(Intensities[AmpF]), 'dB')
        SD.play(SS, blocking=True)
        del(SS)
    
    Hdf5.ExpInfoWrite('Sound', DVCoord, FKey, FileName)
    print('Played Freq', FKey, 'at', DVCoord, 'µm DV')



#%% Laser
## Silence before pulse
#LaserPrePauseDur = 0
## Pulse duration
#LaserPulseDur = 0.01
## Silence after pulse
#LaserPostPauseDur = 0.09
## Amount of pulses per block
#LaserPulseNo = 200
## Number of blocks
#LaserStimBlockNo = 1
## Duration of pause between blocks
#LaserPauseBetweenStimBlocksDur = 5

## Prepare laser stimulation
#Laser, LaserPauseBetweenStimBlocks, _ = \
#    SigGen.LaserStim(Rate, LaserPulseDur, LaserPulseNo, TTLAmpF, 
#                               CalibrationFile, SoundBoard, LaserPrePauseDur, 
#                               LaserPostPauseDur, LaserStimBlockNo, 
#                               LaserPauseBetweenStimBlocksDur)

### Run laser
#DVCoord = input('Choose DVCoord (in µm): '); 
#DVCoord = 'Out'
##
#print('Running...')
#Arduino.write(b'P')
#for OneBlock in range(LaserStimBlockNo):
#    for OnePulse in range(LaserPulseNo):
##        Arduino.write(b'b')
#        Stimulation.write(Laser)
##        Arduino.write(b'y')
#    
#    Stimulation.write(LaserPauseBetweenStimBlocks)
#Arduino.write(b'P')
#
##print('Done. Saving info...')
##lHz = 1000/round((LaserPulseDur+LaserPostPauseDur)*1000)
##with h5py.File(FileName) as h5:
##    h5.create_group(str(len(list(h5)) - 1))
##    h5[list(h5.keys())[-2]].attrs['StimType'] = [np.string_('Laser')]
##    h5[list(h5.keys())[-2]].attrs['DVCoord'] = DVCoord
##    h5[list(h5.keys())[-2]].attrs['lHz'] = lHz
##
##print('Saved.')
##print('Ran laser pulses at ' + str(lHz) + ' at ' + DVCoord + 'µm DV')
#
#

#%% Prepare sound and laser simultaneous stimulation
#SoundAndLaser, SoundAndLaserPauseBetweenStimBlocks, _ = \
#    SigGen.SoundLaserStim(Rate, SoundPulseDur, SoundPulseNo, 
#                                    SoundAmpF, NoiseFrequency, LaserPulseDur, 
#                                    LaserPulseNo, TTLAmpF, CalibrationFile, 
#                                    SoundBoard, SoundPrePauseDur, 
#                                    SoundPostPauseDur, SoundStimBlockNo, 
#                                    SoundPauseBetweenStimBlocksDur, 
#                                    LaserPrePauseDur, LaserPostPauseDur, 
#                                    LaserStimBlockNo, 
#                                    LaserPauseBetweenStimBlocksDur)

## Run sound and laser
##Hz = input('Choose Freq index: ')
##DVCoord = input('Choose DVCoord (in µm): '); 
#DVCoord = 'Out'
#Hz = 0
#Hz = int(Hz)
#
#print('Running...')
#Key = str(NoiseFrequency[Hz][0]) + '-' + str(NoiseFrequency[Hz][1])
#for AmpF in range(len(SoundAmpF[Key])):
#    Arduino.write(b'P')
#    for OnePulse in range(SoundPulseNo):
#        Stimulation.write(SoundAndLaser[Hz][AmpF])
#
#    Arduino.write(b'P')
#    Stimulation.write(SoundAndLaserPauseBetweenStimBlocks)
 
#print('Done. Saving info...')
#lHz = 1000/round((LaserPulseDur+LaserPostPauseDur)*1000)
#with h5py.File(FileName) as h5:
#    h5.create_group(str(len(list(h5)) - 1))
#    h5[list(h5.keys())[-2]].attrs['StimType'] = [np.string_('Sound'), 
#                                                 np.string_('Laser')]
#    h5[list(h5.keys())[-2]].attrs['DVCoord'] = DVCoord
#    h5[list(h5.keys())[-2]].attrs['Hz'] = Hz
#    h5[list(h5.keys())[-2]].attrs['lHz'] = lHz
#
#print('Saved.')
#print('Played Freq ' + str(Hz) + ' at ' + DVCoord + 'µm DV')
#print('Ran laser pulses at ' + str(lHz) + ' at ' + DVCoord + 'µm DV')