#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: T. Malfatti <malfatti@disroot.org>
@date: 2017-01-24
@license: GNU GPLv3 <https://raw.githubusercontent.com/malfatti/SciScripts/master/LICENSE>
@homepage: https://github.com/Malfatti/SciScripts
"""

import numpy as np
import sounddevice as SD
from datetime import datetime
from IO import Arduino, SigGen, Txt

def AudioSet(Rate, BlockSize, Channels):
    SD.default.device = 'system'
    SD.default.samplerate = Rate
    SD.default.blocksize = BlockSize
    SD.default.channels = Channels
    Stim = SD.OutputStream(dtype='float32')
    
    return(Stim)
    

def InfoWrite(AnimalName, StimType, ToneIntensities, ToneAmpF, ToneFrequency, 
              TonePulseNo, TonePauseBeforePulseDur, TonePulseDur, 
              TonePauseAfterPulseDur, TonePauseBetweenIntensities, 
              ThunderIntensities, ThunderFrequency, ThunderAmpF, ThunderPulseNo, 
              ThunderPauseBeforePulseDur, ThunderPulseDur, ThunderPauseAfterPulseDur, 
              ThunderPauseBetweenIntensities, AnalogTTLs, Rate, BlockSize, 
              Channels, System, Setup, SoundCh, TTLCh, BaudRate, InfoFile):
    
    CalibrationFile = SigGen.CalibrationFile
    
    DataInfo = {'InfoFile': InfoFile}
    for K in ['Animal', 'DAqs', 'Audio', 'ExpInfo']: DataInfo[K] = {}
    
    for K in ['AnimalName', 'StimType']: DataInfo['Animal'][K] = locals()[K]
    
    for K in ['SoundCh', 'TTLCh', 'BaudRate', 'AnalogTTLs']: 
        DataInfo['DAqs'][K] = locals()[K]
    
    for K in ['ToneIntensities', 'ToneAmpF', 'ToneFrequency', 
              'TonePulseNo', 'TonePauseBeforePulseDur', 'TonePulseDur', 
              'TonePauseAfterPulseDur', 'TonePauseBetweenIntensities', 
              'ThunderIntensities', 'ThunderFrequency', 'ThunderAmpF', 'ThunderPulseNo', 
              'ThunderPauseBeforePulseDur', 'ThunderPulseDur', 'ThunderPauseAfterPulseDur', 
              'ThunderPauseBetweenIntensities', 'Rate', 'BlockSize', 
              'Channels', 'System', 'Setup']:
        DataInfo['Audio'][K] = locals()[K]
    
    Txt.DictWrite(InfoFile, DataInfo)
    
    return(DataInfo)


def Prepare(AnimalName, StimType, ToneIntensities, ToneFrequency, TonePulseNo, 
            TonePauseBeforePulseDur, TonePulseDur, TonePauseAfterPulseDur, 
            TonePauseBetweenIntensities, ThunderIntensities, ThunderFrequency, 
            ThunderPulseNo, ThunderPauseBeforePulseDur, ThunderPulseDur,
            ThunderPauseAfterPulseDur, ThunderPauseBetweenIntensities, 
            System, Setup, SoundCh, TTLCh, AnalogTTLs=True, Rate=192000, 
            BlockSize=384, Channels=2, BaudRate=115200):
    Kws = {**locals()}
    ToneAmpF = SigGen.dBToAmpF(ToneIntensities, System+'/'+Setup)
    ThunderAmpF = SigGen.dBToAmpF(ThunderIntensities, System+'/'+Setup)
    
    for F in ToneFrequency: 
        if str(F) not in ToneAmpF: ToneAmpF[str(F)] = ToneAmpF['8000-18000']
    
    for F in ThunderFrequency: 
        FKey = '-'.join([str(_) for _ in F])
        ThunderAmpF[FKey] = ThunderAmpF['8000-18000']
        # if FKey not in ThunderAmpF: 
        #     ThunderAmpF[FKey] = ThunderAmpF['8000-18000']
    
    
    Date = datetime.now().strftime("%Y%m%d%H%M%S")
    InfoFile = '-'.join([Date, AnimalName, '_'.join(StimType)+'.dict'])
    
    DataInfo = InfoWrite(ToneAmpF=ToneAmpF, ThunderAmpF=ThunderAmpF, 
                         InfoFile=InfoFile, **Kws)
    
    # TTLs Amplification factor. DO NOT CHANGE unless you know what you're doing.
    DataInfo['Audio']['TTLAmpF'] = 0.6 # for analog TTLs only
    
    Stimulation = {}
    Stimulation['Stim'] = AudioSet(Rate, BlockSize, Channels)
    Stimulation['ArduinoObj'] = Arduino.CreateObj(BaudRate)
    
    Prefix = []
    if 'Tone' in StimType: Prefix.append('Tone')
    if 'Thunder' in StimType: Prefix.append('Thunder')
    
    for P in Prefix:
        Stimulation[P] = SigGen.SoundStim(
                DataInfo['Audio']['Rate'], 
                DataInfo['Audio'][P+'PulseDur'], 
                DataInfo['Audio'][P+'AmpF'], 
                DataInfo['Audio'][P+'Frequency'], 
                DataInfo['Audio']['TTLAmpF'], 
                DataInfo['Audio']['System'], 
                DataInfo['Audio'][P+'PauseBeforePulseDur'], 
                DataInfo['Audio'][P+'PauseAfterPulseDur'], 
                True, 
                [1,2])
        
        Stimulation[P+'Pause'] = np.zeros(
                (DataInfo['Audio'][P+'PauseBetweenIntensities']*Rate,2), dtype='float32')

    return(Stimulation, InfoFile)


def PlaySound(Sound, Pause, Stim, ArduinoObj, InfoFile, StimType, DV='None', Ramp=False):
    DataInfo = Txt.DictRead(InfoFile)
    
    FKeys = list(Sound.keys())
    FKeys.sort(key=lambda x: [int(y) for y in x.split('-')])
    
    if 'Tone' in StimType: Prefix = 'Tone'
    elif 'Thunder' in StimType: Prefix = 'Thunder'
    else:
        print(""" StimType should contain 'Tone' or 'Thunder', otherwise, 
                  this function is useless :)""")
    
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
            # SS = np.concatenate([Sound[FKey][AKey] for _ in range(DataInfo['Audio'][Prefix+'PulseNo'])])
            
            if Ramp:
                NoOfSamples = int(0.003*DataInfo['Audio']['Rate'])
                Sound[FKey][AKey][:NoOfSamples,1] *= np.linspace(0, 1, NoOfSamples)
                Sound[FKey][AKey][-NoOfSamples:,1] *= np.linspace(1, 0, NoOfSamples)
                
            print('Playing', FKey, 'at', str(DataInfo['Audio'][Prefix+'Intensities'][AmpF]), 'dB')
            for Pulse in range(DataInfo['Audio'][Prefix+'PulseNo']):
                ArduinoObj.write(b'D')
                Stim.write(Sound[FKey][AKey])
                # ArduinoObj.write(b'w')
            
            Stim.write(Pause)
            
        
        Rec = "{0:02d}".format(len(DataInfo['ExpInfo']))
        DataInfo['ExpInfo'][Rec] = {'DV': DV, 'StimType': StimType, 'Hz': FKey}
        Txt.DictWrite(InfoFile, DataInfo)
        
        print('Played Freq', FKey, 'at', DV, 'µm DV')
    
    Stim.stop()
    return(None)


# def PlayThunder(Sound, Stim, ArduinoObj, DataInfo, StimType, DV=None):
#     FKey = list(Sound.keys())[0]
#     AKey = list(Sound[FKey].keys())[0]
    
#     Stim.start()
#     print('Playing Thunder...')
    
#     AKeys = list(); AKeys = sorted(AKeys, reverse=True)
#     for AKey in Sound[FKey].keys():
#         for Pulse in range(DataInfo['Audio']['SoundPulseNo']):
#             ArduinoObj.write(b'a')
#             Stim.write(Sound[FKey][AKey])
#             ArduinoObj.write(b'z')
#     Stim.stop()
    
#     Rec = "{0:02d}".format(len(DataInfo['ExpInfo']))
#     DataInfo['ExpInfo'][Rec] = {'DV': DV, 'StimType': StimType, 'Hz': FKey}
#     Txt.DictWrite(DataInfo['InfoFile'], DataInfo)
    
#     return(None)


def Play(Stimulation, InfoFile, StimType, DV, Ramp=False):
    if 'Tone' in StimType: Prefix = 'Tone'
    elif 'Thunder' in StimType: Prefix = 'Thunder'
    else:
        print(""" StimType should contain 'Tone' or 'Thunder', otherwise, 
                  this function is useless :)""")
    
    PlaySound(Stimulation[Prefix], Stimulation[Prefix+'Pause'], 
              Stimulation['Stim'], Stimulation['ArduinoObj'], InfoFile, 
              StimType, DV, Ramp)
    
    return(None)