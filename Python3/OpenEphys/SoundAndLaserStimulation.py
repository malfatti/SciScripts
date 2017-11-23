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
sound square waves (TTLs), and laser square waves. The square waves will be 
sent to the left channel and the sound pulses will be sent to the right 
channel. 
"""
#%% Settings
import datetime
import numpy as np
import sounddevice as SD

from IO import Arduino, Hdf5, SigGen, Txt

## Experiment parameters
AnimalName = 'Prevention_A5'

ABRCh = list(range(1,17))
SoundCh = 26
TTLCh = 27
AnalogTTLs = True

# Fill all durations in SECONDS!
SoundPauseBeforePulseDur = 0.004
SoundPulseDur = 0.003
SoundPauseAfterPulseDur = 0.093
SoundPulseNo = 529
Intensities = [80, 70, 60, 50, 40]
PauseBetweenIntensities = 10
NoiseFrequency = [[8000, 10000], [9000, 11000], [10000, 12000], 
                  [12000, 14000], [14000, 16000]]


## Hardware parameters
System = 'Jack-IntelOut-MackieIn-MackieOut-IntelIn'
Setup = 'GPIAS'
Rate = 192000
BaudRate = 115200

# TTLs Amplification factor. DO NOT CHANGE unless you know what you're doing.
TTLAmpF = 0.6 # for analog TTLs only
#TTLAmpF = 6.8 # for analog and digital TTLs

# Set sound stimulation
SoundAmpF = SigGen.dBToAmpF(Intensities, System+'/'+Setup)
# Temporary override
#SoundAmpF = Hdf5.DataLoad('/DataInfo/SoundAmpF', 'Prevention/20170704-Prevention_A1-ABRs/20170704102002-Prevention_A1-SoundStim.hdf5')[1]
Sound = SigGen.SoundStim(Rate, SoundPulseDur, SoundAmpF, NoiseFrequency, 
                         TTLAmpF, System, SoundPauseBeforePulseDur, 
                         SoundPauseAfterPulseDur)

Pause = np.zeros((PauseBetweenIntensities*Rate,2), dtype='float32')

# Set audio objects
SD.default.device = 'system'
SD.default.samplerate = Rate
SD.default.blocksize = 384
SD.default.channels = 2
Stim = SD.OutputStream(dtype='float32')

# Set arduino object
ArduinoObj = Arduino.CreateObj(BaudRate)
# ArduinoObj.write(b'd')


## Write info
Date = datetime.datetime.now()
FileName = ''.join([Date.strftime("%Y%m%d%H%M%S"), '-', AnimalName, 
                    '-SoundStim.hdf5'])

DataInfo = dict((Name, eval(Name)) 
                for Name in ['AnimalName', 'Rate', 'BaudRate', 
                             'SoundPrePauseDur', 'SoundPulseDur', 
                             'SoundPostPauseDur', 'SoundPulseNo', 
                             'Intensities', 'NoiseFrequency', 
                             'PauseBetweenIntensities',
                             'SoundCh', 'TTLCh', 'ABRCh', 'AnalogTTLs'
#                             'LaserPrePauseDur', 'LaserPulseDur', 
#                             'LaserPostPauseDur', 'LaserPulseNo', 
#                             'LaserStimBlockNo', 
#                             'LaserPauseBetweenStimBlocksDur', 
                             'SigGen.CalibrationFile', 'FileName'])

Hdf5.DictWrite(DataInfo, '/DataInfo', FileName)
Hdf5.DictWrite(SoundAmpF, '/DataInfo/SoundAmpF', FileName)
DataInfo['ExpInfo'] = {}
DataInfo['SoundAmpF'] = {K: Key.tolist() for K, Key in SoundAmpF.items()}
Txt.DictWrite(FileName.split('.')[0]+'.dict', DataInfo)

#%% Run sound

## Trial info
DVCoord = '12752'
StimType = ['Sound', 'CNO']


## Trial run
#Freq = 4
#Freq = int(Freq)

FKeys = list(Sound.keys())
FKeys.sort(key=lambda x: [int(y) for y in x.split('-')])

Stim.start()
while True:
    print('Remember to change folder name in OE!')
    print('Choose frequency:')
    print('-1)', 'Baseline (No stimulus)')
    for Ind, K in enumerate(FKeys): print(str(Ind) + ')' , K)
    print(str(len(FKeys)) + ')', 'Cancel')
    FKey = input(': ')
    
    if FKey == str(len(FKeys)): break
    if FKey == str(-1):
        Hdf5.ExpInfoWrite('Sound', DVCoord, FKey, FileName)
        Rec = "{0:02d}".format(len(DataInfo['ExpInfo']))
        DataInfo['ExpInfo'][Rec] = {'DVCoord': DVCoord, 'StimType': StimType, 'Hz': 'Baseline'}
        Txt.DictWrite(FileName.split('.')[0]+'.dict', DataInfo)
        continue
    
    try:
        FKey = FKeys[int(FKey)]
    except IndexError:
        print('=== Wrong Freq index. Stopping... ===')
        print('')
        break
    
    AKeys = list(Sound[FKey].keys()); AKeys = sorted(AKeys, reverse=True)
    for AmpF, AKey in enumerate(AKeys):
        SS = np.concatenate([Sound[FKey][AKey] for _ in range(SoundPulseNo)])
        
        print('Playing', FKey, 'at', str(Intensities[AmpF]), 'dB')
        ArduinoObj.write(b'd')
        Stim.write(SS)
        ArduinoObj.write(b'w')
        Stim.write(Pause)
        del(SS)
    
    Hdf5.ExpInfoWrite('Sound', DVCoord, FKey, FileName)
    
    Rec = "{0:02d}".format(len(DataInfo['ExpInfo']))
    DataInfo['ExpInfo'][Rec] = {'DVCoord': DVCoord, 'StimType': StimType, 'Hz': FKey}
    Txt.DictWrite(FileName.split('.')[0]+'.dict', DataInfo)
    
    print('Played Freq', FKey, 'at', DVCoord, 'µm DV')

Stim.stop()


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