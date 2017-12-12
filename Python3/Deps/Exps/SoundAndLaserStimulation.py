#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: T. Malfatti <malfatti@disroot.org>
@year: 2015
@license: GNU GPLv3 <https://raw.githubusercontent.com/malfatti/SciScripts/master/LICENSE>
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

from IO import Arduino, SigGen, Txt


def AudioSet(Rate, BlockSize, Channels):
    SD.default.device = 'system'
    SD.default.samplerate = Rate
    SD.default.blocksize = BlockSize
    SD.default.channels = Channels
    Stim = SD.OutputStream(dtype='float32')
    
    return(Stim)
    
def InfoWrite(AnimalName, StimType, Intensities, SoundAmpF, NoiseFrequency, SoundPulseNo, 
              SoundPauseBeforePulseDur, SoundPulseDur, SoundPauseAfterPulseDur, 
              PauseBetweenIntensities, System, Setup, SoundCh, TTLCh, 
              LaserStimBlockNo, LaserPulseNo, LaserPauseBeforePulseDur, 
              LaserPulseDur, LaserPauseAfterPulseDur, LaserPauseBetweenStimBlocksDur,
              LaserType, LaserDur, LaserFreq, ABRCh, AnalogTTLs, Rate, BlockSize, 
              Channels, BaudRate, InfoFile):
    
    CalibrationFile = SigGen.CalibrationFile
    
    DataInfo = {'InfoFile': InfoFile}
    for K in ['Animal', 'DAqs', 'Audio', 'Laser', 'ExpInfo']: DataInfo[K] = {}
    
    for K in ['AnimalName', 'StimType']: DataInfo['Animal'][K] = locals()[K]
    
    for K in ['SoundCh', 'TTLCh', 'ABRCh', 'BaudRate', 'AnalogTTLs']: 
        DataInfo['DAqs'][K] = locals()[K]
    
    for K in ['Rate', 'Intensities', 'NoiseFrequency', 
              'SoundPulseNo', 'SoundPauseBeforePulseDur', 'SoundPulseDur', 
              'SoundPauseAfterPulseDur', 'PauseBetweenIntensities', 
              'System', 'Setup', 'CalibrationFile']:
        DataInfo['Audio'][K] = locals()[K]
    
    for K in ['LaserStimBlockNo', 'LaserPulseNo', 'LaserPauseBeforePulseDur', 
              'LaserPulseDur', 'LaserPauseAfterPulseDur', 
              'LaserPauseBetweenStimBlocksDur', 'LaserType', 'LaserDur', 
              'LaserFreq']:
        DataInfo['Laser'][K] = locals()[K]
    
    DataInfo['Audio']['SoundAmpF'] = SoundAmpF
    
    Txt.DictWrite(InfoFile, DataInfo)
    
    return(DataInfo)


def Prepare(AnimalName, StimType, Intensities, NoiseFrequency, SoundPulseNo, 
            SoundPauseBeforePulseDur, SoundPulseDur, SoundPauseAfterPulseDur, 
            PauseBetweenIntensities, System, Setup, SoundCh, TTLCh, 
            LaserStimBlockNo, LaserPulseNo, LaserPauseBeforePulseDur, 
            LaserPulseDur, LaserPauseAfterPulseDur, 
            LaserPauseBetweenStimBlocksDur, LaserType, LaserDur, LaserFreq, 
            ABRCh=[], AnalogTTLs=True, Rate=192000, BlockSize=384, Channels=2, 
            BaudRate=115200):
    Kws = {**locals()}
    SoundAmpF = SigGen.dBToAmpF(Intensities, System+'/'+Setup)
    
    Date = datetime.now().strftime("%Y%m%d%H%M%S")
    InfoFile = '-'.join([Date, AnimalName, '_'.join(StimType)+'.dict'])
    
    DataInfo = InfoWrite(SoundAmpF=SoundAmpF, InfoFile=InfoFile, **Kws)
    # DataInfo = InfoWrite(AnimalName, StimType, BGIntensity, PulseIntensity, 
    #                      NoiseFrequency, SoundBGDur, SoundGapDur, 
    #                      SoundBGPrePulseDur, SoundLoudPulseDur, 
    #                      SoundBGAfterPulseDur, SoundBetweenStimDur, NoOfTrials, 
    #                      System, Setup, SoundBGAmpF, SoundPulseAmpF, 
    #                      SoundCh, TTLCh, PiezoCh, Rate, BlockSize, Channels, 
    #                      BaudRate, InfoFile)
    
    # TTLs Amplification factor. DO NOT CHANGE unless you know what you're doing.
    DataInfo['Audio']['TTLAmpF'] = 0.6 # for analog TTLs only
    # DataInfo['Audio']['TTLAmpF'] = 6.8 # for analog and digital TTLs
    
    
    # Stimulation = {Key: [] for Key in ['Sound', 'SoundPause', 'Laser', 
    #                                    'LaserPause', 'SoundLaser', 
    #                                    'SoundLaserPause']}
    Stimulation = {}
    Stimulation['Stim'] = AudioSet(Rate, BlockSize, Channels)
    Stimulation['ArduinoObj'] = Arduino.CreateObj(BaudRate)
    
    if 'Sound' in StimType:
        Stimulation['Sound'] = SigGen.SoundStim(**DataInfo['Audio'])
        
        Stimulation['SoundPause'] = np.zeros(
                (PauseBetweenIntensities*Rate,2), dtype='float32')
    
    if 'Laser' in StimType:
        Stimulation['Laser'] = SigGen.LaserStim(
                Rate=Rate, 
                TTLAmpF=DataInfo['Audio']['TTLAmpF'], 
                System=System, 
                **DataInfo['Laser'])
        
        Stimulation['LaserPause'] = np.zeros(
                (LaserPauseBetweenStimBlocksDur*Rate,2), dtype='float32')
    
    if 'SoundLaser' in StimType:
        Stimulation['SoundLaser'] = SigGen.SoundLaserStim(
                TTLAmpF=DataInfo['Audio']['TTLAmpF'], 
                Map=[1,2],
                **DataInfo['Audio'], 
                **DataInfo['Laser'])
        
        Stimulation['SoundLaserPause'] = np.zeros(
                (PauseBetweenIntensities*Rate,2), dtype='float32')
    
    return(Stimulation, InfoFile)


def PlaySound(Sound, Pause, Stim, ArduinoObj, InfoFile, StimType, DV='Out'):
    DataInfo = Txt.DictRead(InfoFile)
    
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
            Rec = "{0:02d}".format(len(DataInfo['ExpInfo']))
            DataInfo['ExpInfo'][Rec] = {'DV': DV, 'StimType': StimType, 'Hz': 'Baseline'}
            Txt.DictWrite(InfoFile, DataInfo)
            continue
        
        try:
            FKey = FKeys[int(FKey)]
        except IndexError:
            print('=== Wrong Freq index. Stopping... ===')
            print('')
            break
        
        AKeys = list(Sound[FKey].keys()); AKeys = sorted(AKeys, reverse=True)
        for AmpF, AKey in enumerate(AKeys):
            SS = np.concatenate([Sound[FKey][AKey] for _ in range(DataInfo['Audio']['SoundPulseNo'])])
            
            print('Playing', FKey, 'at', str(DataInfo['Audio']['Intensities'][AmpF]), 'dB')
            ArduinoObj.write(b'd')
            Stim.write(SS)
            ArduinoObj.write(b'w')
            Stim.write(Pause)
            del(SS)
        
        Rec = "{0:02d}".format(len(DataInfo['ExpInfo']))
        DataInfo['ExpInfo'][Rec] = {'DV': DV, 'StimType': StimType, 'Hz': FKey}
        Txt.DictWrite(InfoFile, DataInfo)
        
        print('Played Freq', FKey, 'at', DV, 'µm DV')
    
    Stim.stop()
    return(None)


def PlayLaser(Laser, LaserStimBlockNo, Pause, Stim, ArduinoObj, InfoFile, StimType, DV='Out'):
    DataInfo = Txt.DictRead(InfoFile)
    
    Stim.start()
    while True:
        print('What to do?')
        print('-1) Baseline (No stimulus)')
        print('0) Run stimulation')
        print('1) Cancel')
        Ans = input(': ')
        
        if Ans == '1': break
        if Ans == str(-1):
            Rec = "{0:02d}".format(len(DataInfo['ExpInfo']))
            DataInfo['ExpInfo'][Rec] = {'DV': DV, 'StimType': StimType, 'Hz': 'Baseline'}
            Txt.DictWrite(InfoFile, DataInfo)
            continue
        
        for Block in range(LaserStimBlockNo):
            # LL = np.concatenate([Laser for _ in range(DataInfo['Laser']['LaserPulseNo'])])
            print('Running laser stimulation, block', Block+1, 'of', LaserStimBlockNo)
            
            ArduinoObj.write(b'd')
            for Pulse in range(DataInfo['Laser']['LaserPulseNo']): 
                Stim.write(Laser)
            ArduinoObj.write(b'w')
            
            Stim.write(Pause)
            # del(LL)
        
        Rec = "{0:02d}".format(len(DataInfo['ExpInfo']))
        if DataInfo['Laser']['LaserType'] == 'Sin':
            DataInfo['ExpInfo'][Rec] = {'DV': DV, 'StimType': StimType, 'Hz': DataInfo['Laser']['LaserFreq']}
        else:
            DataInfo['ExpInfo'][Rec] = {'DV': DV, 'StimType': StimType, 'Hz': 'LaserPulses'}
        
        Txt.DictWrite(InfoFile, DataInfo)
        
        print('Finished laser stimulation at', DV, 'µm DV')
    
    Stim.stop()
    return(None)


def Play(Stimulation, InfoFile, StimType, DV):
    DataInfo = Txt.DictRead(InfoFile)
    
    if 'Sound' in StimType and 'Laser' in StimType:
        PlaySound(Stimulation['SoundLaser'], Stimulation['SoundLaserPause'], 
                  Stimulation['Stim'], Stimulation['ArduinoObj'], InfoFile, 
                  StimType, DV)
    
    elif 'Sound' in StimType and 'Laser' not in StimType:
        PlaySound(Stimulation['Sound'], Stimulation['SoundPause'], 
                  Stimulation['Stim'], Stimulation['ArduinoObj'], InfoFile, 
                  StimType, DV)
        
    elif 'Sound' not in StimType and 'Laser' in StimType:
        PlayLaser(Stimulation['Laser'], DataInfo['Laser']['LaserStimBlockNo'], 
                  Stimulation['LaserPause'], Stimulation['Stim'], 
                  Stimulation['ArduinoObj'], InfoFile, StimType, DV)
    
    else:
        print(""" StimType should contain 'Sound', 'Laser' or both, otherwise, 
                  this function is useless :)""")
    
    return(None)