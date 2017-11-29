#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: T. Malfatti
@year: 2015
@licence: GNU GPLv3 <https://raw.githubusercontent.com/malfatti/SciScripts/master/LICENSE>
@homepage: https://github.com/Malfatti/SciScripts

This is a script to generate pulses and send to the soundboard, and then to a 
sound amplifier and an Arduino board. Basically it generates sound pulses, 
sound square waves (TTLs), and laser square waves. The square waves will be 
sent to the left channel and the sound pulses will be sent to the right 
channel. 
"""
#%% Settings
from datetime import datetime
import numpy as np
import sounddevice as SD

from IO import Arduino, Hdf5, SigGen, Txt


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
              TTLCh, PiezoCh, AnalogTTLs, Rate, BlockSize, Channels, BaudRate, InfoFile):
    
    CalibrationFile = SigGen.CalibrationFile
    
    DataInfo = {'InfoFile': InfoFile, 'Animal': {}, 'DAqs': {}, 'Audio': {}}
    for K in ['AnimalName', 'StimType']: DataInfo['Animal'][K] = locals()[K]
    
    for K in ['SoundCh', 'TTLCh', 'PiezoCh', 'BaudRate', 'AnalogTTLs']: 
        DataInfo['DAqs'][K] = locals()[K]
    
    for K in ['Rate', 'BaudRate', 
                             'SoundPrePauseDur', 'SoundPulseDur', 
                             'SoundPostPauseDur', 'SoundPulseNo', 
                             'Intensities', 'NoiseFrequency', 
                             'PauseBetweenIntensities',
                             'SoundCh', 'TTLCh', 'ABRCh', 'AnalogTTLs'
                            'LaserPrePauseDur', 'LaserPulseDur', 
                            'LaserPostPauseDur', 'LaserPulseNo', 
                            'LaserStimBlockNo', 
                            'LaserPauseBetweenStimBlocksDur', 
                             'SigGen.CalibrationFile', 'FileName']:
        DataInfo['Audio'][K] = locals()[K]
    
    DataInfo['ExpInfo'] = {}
    # DataInfo['Audio']['SoundBGAmpF'] = {K: Key.tolist() for K, Key in SoundBGAmpF.items()}
    # DataInfo['Audio']['SoundPulseAmpF'] = {K: Key.tolist() for K, Key in SoundPulseAmpF.items()}
    DataInfo['Audio']['SoundBGAmpF'] = SoundBGAmpF
    DataInfo['Audio']['SoundPulseAmpF'] = SoundPulseAmpF
    
    Txt.DictWrite(InfoFile, DataInfo)
    
    return(DataInfo)


def Run(AnimalName, StimType, Intensities, NoiseFrequency, SoundPulseNo, 
        SoundPauseBeforePulseDur, SoundPulseDur, SoundPauseAfterPulseDur, 
        PauseBetweenIntensities, SoundSystem, Setup, SoundCh, SoundTTLCh, 
        LaserStimBlockNo, LaserPulseNo, LaserPauseBeforePulseDur, 
        LaserPulseDur, LaserPauseAfterPulseDur, LaserPauseBetweenStimBlocksDur, 
        ABRCh=[], AnalogTTLs=True, Rate=192000, BlockSize=384, 
        Channels=2, BaudRate=115200, **Kws):
    
    SoundAmpF = SigGen.dBToAmpF(Intensities, SoundSystem+'/'+Setup)
    
    Date = datetime.now().strftime("%Y%m%d%H%M%S")
    InfoFile = '-'.join([Date, AnimalName, '_'.join(StimType)+'.dict'])
    
    DataInfo = InfoWrite(**Kws)
    # DataInfo = InfoWrite(AnimalName, StimType, BGIntensity, PulseIntensity, 
    #                      NoiseFrequency, SoundBGDur, SoundGapDur, 
    #                      SoundBGPrePulseDur, SoundLoudPulseDur, 
    #                      SoundBGAfterPulseDur, SoundBetweenStimDur, NoOfTrials, 
    #                      SoundSystem, Setup, SoundBGAmpF, SoundPulseAmpF, 
    #                      SoundCh, TTLCh, PiezoCh, Rate, BlockSize, Channels, 
    #                      BaudRate, InfoFile)
    
    # TTLs Amplification factor. DO NOT CHANGE unless you know what you're doing.
    DataInfo['Audio']['TTLAmpF'] = 0.6 # for analog TTLs only
    # DataInfo['Audio']['TTLAmpF'] = 6.8 # for analog and digital TTLs
    
    Sound = SigGen.SoundStim(**DataInfo['Audio'])
    Pause = np.zeros((PauseBetweenIntensities*Rate,2), dtype='float32')
    Stim = AudioSet(Rate, BlockSize, Channels)
    ArduinoObj = Arduino.CreateObj(BaudRate)
    
    Play(Sound, Stim, ArduinoObj, DataInfo=DataInfo, **DataInfo['Audio'])
    
    print('Finished', AnimalName, 'GPIAS.')







# TTLs Amplification factor. DO NOT CHANGE unless you know what you're doing.
TTLAmpF = 0.6 # for analog TTLs only
#TTLAmpF = 6.8 # for analog and digital TTLs

# Sound = SigGen.SoundStim(Rate, SoundPulseDur, SoundAmpF, NoiseFrequency, 
#                          TTLAmpF, System, SoundPauseBeforePulseDur, 
#                          SoundPauseAfterPulseDur)

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