#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: T. Malfatti <malfatti@disroot.org>
@date: 2017-06-12
@license: GNU GPLv3 <https://raw.githubusercontent.com/malfatti/SciScripts/master/LICENSE>
@homepage: https://github.com/Malfatti/SciScripts
"""
import os
import numpy as np

from IO import Hdf5
from DataAnalysis.DataAnalysis import FilterSignal
from scipy import signal


# Use one that was used in SoundBoardCalibration.py
CalibrationFile = os.environ['DATAPATH']+'/Tests/SoundMeasurements/SoundMeasurements.hdf5'
SoundTTLVal = 0.6; LaserTTLVal = 0.3

## Level 0
def ApplySoundAmpF(SoundPulseFiltered, Rate, SoundAmpF, NoiseFrequency, 
                   SBOutAmpF, SoundPauseBeforePulseDur=0, SoundPauseAfterPulseDur=0):
    print('Applying amplification factors...')
    # Preallocating memory
    SoundPauseBeforePulse = np.zeros(round(Rate * SoundPauseBeforePulseDur), dtype=np.float32)
    SoundPauseAfterPulse = np.zeros(round(Rate * SoundPauseAfterPulseDur), dtype=np.float32)
    SoundUnit = {}
    
    for FKey in SoundPulseFiltered:
        SoundUnit[FKey] = {}
        
        for AmpF in range(len(SoundAmpF[FKey])):       
            if SoundAmpF[FKey][AmpF] > 1/SBOutAmpF:
                print(SoundAmpF[FKey][AmpF], 
                      'AmpF out of range. Decreasing to', 1/SBOutAmpF, '.')
                SoundAmpF[FKey][AmpF] = 1/SBOutAmpF
            
            AKey = str(SoundAmpF[FKey][AmpF])
            SoundUnit[FKey][AKey] = np.concatenate([SoundPauseBeforePulse, 
                                                    SoundPulseFiltered[FKey],
                                                    SoundPauseAfterPulse])
            SoundUnit[FKey][AKey] = (SoundUnit[FKey][AKey]
                                     * SoundAmpF[FKey][AmpF]) * SBOutAmpF
    return(SoundUnit)


def BandpassFilterSound(SoundPulse, Rate, NoiseFrequency):
    # Preallocating memory
    SoundPulseFiltered = {}
    PulseAmp = (max(SoundPulse)-min(SoundPulse))/2
    
    print('Filtering sound: ', end='')
    for Freq in NoiseFrequency:
        FKey = '-'.join([str(_) for _ in Freq])
        print(FKey, end='...')
        
        SoundPulseFiltered[FKey] = FilterSignal(SoundPulse, Rate, Freq)
        SoundPulseFiltered[FKey] = SoundPulseFiltered[FKey].astype('float32')
        PulseFilteredAmp = (max(SoundPulseFiltered[FKey])-min(SoundPulseFiltered[FKey]))/2
        FilterAmpF = PulseAmp/PulseFilteredAmp
        SoundPulseFiltered[FKey] = SoundPulseFiltered[FKey]*FilterAmpF
        SoundPulseFiltered[FKey][-1] = 0
        
    print(end='\n')
    return(SoundPulseFiltered)


def dBToAmpF(Intensities, Path, CalibrationFile=CalibrationFile):
    print('Converting dB to AmpF...')
    SoundIntensity = Hdf5.DataLoad(Path+'/SoundIntensity', CalibrationFile)[0]
    
    SoundAmpF = {Hz: [float(min(SoundIntensity[Hz].keys(), 
                                key=lambda i: abs(SoundIntensity[Hz][i]['dB']-dB)))
                      for dB in Intensities]
                 for Hz in list(SoundIntensity)}
    
    return(SoundAmpF)     


def InterleaveChannels(Right, Left, SoundAmpF, NoiseFrequency):
    print('Interleaving channels...')
    if Right == [0]:
        List = [0]*(2*len(Left))
        for _ in range(len(Left)):
            List[_ *2] = 0
            List[_ *2+1] = Left[_]
    
    elif Left == [0]:
        List = [0]*(2*len(Right))
        for _ in range(len(Right)):
            List[_ *2] = Right[_]
            List[_ *2+1] = 0
    
    else:
        List = [0]*len(NoiseFrequency)
        
        for Freq in range(len(NoiseFrequency)):
            Key= str(NoiseFrequency[Freq][0]) + '-' + \
                 str(NoiseFrequency[Freq][1])
            List[Freq] = [0]*len(SoundAmpF[Key])
            
            for AmpF in range(len(SoundAmpF[Key])):
                List[Freq][AmpF] = [0]*(2*len(Right[Freq][AmpF]))
                
                for _ in range(len(Right[Freq][AmpF])):
                    List[Freq][AmpF][_ *2] = Right[Freq][AmpF][_]
                    List[Freq][AmpF][_ *2+1] = Left[_]
    
    return(List)


#def ListToByteArray(List, SoundAmpF, NoiseFrequency):
#    print('Converting list to byte array...')
#    if SoundAmpF == [0]:
#        ByteArray = array.array('f', List)
#        ByteArray = bytes(ByteArray)
#    
#    else:
#        ByteArray = [0]*len(NoiseFrequency)
#        
#        for Freq in range(len(NoiseFrequency)):
#            Key= str(NoiseFrequency[Freq][0]) + '-' + \
#                 str(NoiseFrequency[Freq][1])
#            ByteArray[Freq] = [0]*len(SoundAmpF[Key])
#            
#            for AmpF in range(len(SoundAmpF[Key])):
#                ByteArray[Freq][AmpF] = array.array('f', List[Freq][AmpF])
#                ByteArray[Freq][AmpF] = bytes(ByteArray[Freq][AmpF])
#    
#    return(ByteArray)


def Noise(Rate, SoundPulseDur):
    print('Generating noise...')
    Noise = np.random.uniform(-1, 1, size=round(Rate*SoundPulseDur))
    Noise[-1] = 0
    Noise = np.array(Noise, dtype=np.float32)
    
    return(Noise)


def SineWave(Rate, Freq, AmpF, Time):
    print('Generating sine wave...')
    
    if type(Freq) in [int, float]:
        ## Ensure that there will be a sample at each peak
        P = 1/(Freq*4)
        TimeShift = (P*Rate) - int(P*Rate)
        Shift = 2 * np.pi * Freq * TimeShift/Rate
        
        ## Fast way
        Pulse = np.zeros(int(Rate*(1/Freq)), dtype=np.float32)
        for s in range(int(Rate*(1/Freq))):
            Pulse[s] = np.sin((2 * np.pi * Freq * (s/Rate)) - Shift) * AmpF
        Pulse = np.tile(Pulse, int(Time/(1/Freq)))
        Left = (Time%(1/Freq)) * (1/Freq)
        if Left:
            L = np.zeros(int(Rate*Left), dtype=np.float32)
            for s in range(int(Rate*Left)):
                L[s] = np.sin((2 * np.pi * Freq * (s/Rate)) - Shift) * AmpF
            Pulse = np.concatenate((Pulse, L))
        
        ## Obvious way
        # Pulse = [np.sin((2 * np.pi * Freq * (_/Rate)) - Shift) * AmpF
        #          for _ in range(round(Rate*SoundPulseDur))]
    
    else:
        # Example:
        # Freq = [800, 1600]
        # Freq = [(Freq[1] - Freq[0])/Time, Freq[0]]
        # Time = np.linspace(0, Time, int(Rate*Time))
        
        Pulse = signal.sweep_poly(Time, Freq)
    
    Pulse[-1] = 0
    
    return(Pulse)


def SqPulse(Rate, PulseDur, TTLAmpF, TTLVal, SBOutAmpF, PauseBeforePulseDur=0, 
           PauseAfterPulseDur=0):
    Pulse = np.concatenate([
                np.zeros(PauseBeforePulseDur*Rate, dtype=np.float32),
                np.ones(PulseDur*Rate, dtype=np.float32),
                np.zeros(PauseAfterPulseDur*Rate, dtype=np.float32)])
    
    return(Pulse)


def SqWave(Rate, PulseDur, TTLAmpF, TTLVal, SBOutAmpF, PauseBeforePulseDur=0, 
           PauseAfterPulseDur=0):
    
    print('Generating Sound TTL...')
    TTLSpace = PulseDur + PauseAfterPulseDur
    if TTLSpace < 2*PulseDur:
        TTLPulse = np.concatenate([
                  np.array([TTLVal] * round(Rate*PulseDur/2), dtype=np.float32),
                  np.array([TTLVal*-1] * round(Rate*PulseDur/2), dtype=np.float32)
                  ])
    else:
        TTLPulse = np.concatenate([
                  np.array([TTLVal] * round(Rate*PulseDur), dtype=np.float32),
                  np.array([TTLVal*-1] * round(Rate*PulseDur), dtype=np.float32)
                  ])
    
    TTLPulse[-1] = 0
    
    if PauseBeforePulseDur == 0:
        if PauseAfterPulseDur == 0:
            TTLUnit = TTLPulse
        else:
            TTLPauseAfterPulse = np.zeros(round((PauseAfterPulseDur-PulseDur) * Rate), 
                            dtype=np.float32)
            TTLUnit = np.concatenate([TTLPulse, TTLPauseAfterPulse])
    else:
        TTLPauseBeforePulse = np.zeros(round(PauseBeforePulseDur * Rate), dtype=np.float32)
        if PauseAfterPulseDur == 0:
            TTLUnit = np.concatenate([TTLPauseBeforePulse, TTLPulse])
        else:
            TTLPauseAfterPulse = np.zeros(round((PauseAfterPulseDur-PulseDur) * Rate), 
                            dtype=np.float32)
            TTLUnit = np.concatenate([TTLPauseBeforePulse, TTLPulse, TTLPauseAfterPulse])
    
    TTLUnit = (TTLUnit * TTLAmpF) * SBOutAmpF
    
    return(TTLUnit)


## Level 1
def LaserStim(Rate, LaserPulseDur, LaserType, LaserDur, LaserFreq, TTLAmpF, System, LaserPauseBeforePulseDur=0, LaserPauseAfterPulseDur=0, Ch=1):
    """ if LaserType == 'Sq':
            Generate square waves in one channel that works as TTLs for laser.
            
            WARNING: The signal generated is composed of square WAVES, not pulses,
            meaning that it reaches positive AND NEGATIVE values. If your device 
            handles only positive voltage, use a diode on the input of your device.
            
            https://en.wikipedia.org/wiki/Diode
        
        elif LaserType == 'Sin':
            Generate sine waves in one channel.
            
            WARNING: The signal generated is a sine wave that reaches positive 
            AND NEGATIVE values. If your device handles only positive voltage, 
            use a DC offset circuit to shift the wave to the positive range.
            
            https://en.wikipedia.org/wiki/Voltage_divider
    """
    
    SBOutAmpF = Hdf5.DataLoad('/'+System+'/SBOutAmpF', CalibrationFile)[0]
    
    if TTLAmpF > 1/SBOutAmpF:
        print('AmpF out of range. Decreasing to', 1/SBOutAmpF, '.')
        TTLAmpF = 1/SBOutAmpF
    
    if LaserType == 'Sq':
        LaserUnit = SqWave(Rate, LaserPulseDur, TTLAmpF, LaserTTLVal, SBOutAmpF, 
                       LaserPauseBeforePulseDur, LaserPauseAfterPulseDur)
    
    elif LaserType == 'Sin':
        LaserUnit = SineWave(Rate, LaserFreq, TTLAmpF*LaserTTLVal, LaserDur)
    
    Laser = np.zeros((LaserUnit.shape[0], 2))
    Laser[:,Ch-1] = LaserUnit.T
    Laser = np.ascontiguousarray(Laser)
    
    print('Done generating laser stimulus.')
    return(Laser)


def SoundLaserStim(Rate, SoundPulseDur, SoundAmpF, NoiseFrequency, LaserPulseDur, LaserType, LaserDur, LaserFreq, TTLAmpF, System, SoundPauseBeforePulseDur=0,
                   SoundPauseAfterPulseDur=0, LaserPauseBeforePulseDur=0, LaserPauseAfterPulseDur=0, Map=[1,2]):
    """ Generate sound pulses in one channel and a mix of square waves that 
        works as TTLs for both sound and laser in the other channel.
        
        WARNING: The signal generated in the TTLs channel is composed of square 
        WAVES, not pulses, meaning that it reaches positive AND NEGATIVE values. 
        If your device  handles only positive voltage, use a diode on the 
        input of your device.
        
        https://en.wikipedia.org/wiki/Diode
    """
    
    SBOutAmpF = Hdf5.DataLoad('/'+System+'/SBOutAmpF', CalibrationFile)[0]
    
    if TTLAmpF > 1/SBOutAmpF:
        print('AmpF out of range. Decreasing to', 1/SBOutAmpF, '.')
        TTLAmpF = 1/SBOutAmpF
    
    SoundPulse = Noise(Rate, SoundPulseDur)
    print('   ', end='')
    SoundPulseFiltered = BandpassFilterSound(SoundPulse, Rate, NoiseFrequency)
    print('   ', end='')
    SoundUnit = ApplySoundAmpF(SoundPulseFiltered, Rate, SoundAmpF, 
                               NoiseFrequency, SBOutAmpF, SoundPauseBeforePulseDur, 
                               SoundPauseAfterPulseDur)
    
    SoundTTLUnit = SqWave(Rate, SoundPulseDur, TTLAmpF, SoundTTLVal, 
                          SBOutAmpF, SoundPauseBeforePulseDur, 
                          SoundPauseAfterPulseDur)
    
    if LaserType == 'Sq':
        LaserUnit = SqWave(Rate, LaserPulseDur, TTLAmpF, LaserTTLVal, SBOutAmpF, 
                       LaserPauseBeforePulseDur, LaserPauseAfterPulseDur)
    
    elif LaserType == 'Sin':
        LaserUnit = SineWave(Rate, LaserFreq, TTLAmpF*LaserTTLVal, LaserDur)
    
    SoundLaser = {}
    for FKey in SoundUnit:
        SoundLaser[FKey] = {}
        
        for AKey in SoundUnit[FKey]:
            if Map[0] == 2:
                SoundLaser[FKey][AKey] = np.vstack((SoundUnit[FKey][AKey], SoundTTLUnit+LaserUnit)).T
                SoundLaser[FKey][AKey] = np.ascontiguousarray(SoundLaser[FKey][AKey])
            else:
                SoundLaser[FKey][AKey] = np.vstack((SoundTTLUnit+LaserUnit, SoundUnit[FKey][AKey])).T
                SoundLaser[FKey][AKey] = np.ascontiguousarray(SoundLaser[FKey][AKey])
    
    print('Done generating sound and laser stimulus.')
    return(SoundLaser)


def SoundStim(Rate, SoundPulseDur, SoundAmpF, NoiseFrequency, TTLAmpF, 
              System, SoundPauseBeforePulseDur=0, SoundPauseAfterPulseDur=0, TTLs=True,
              Map=[1,2], **Kws):
    """ Generate sound pulses in one channel and TTLs in the other channel.
        
        WARNING: The signal generated in the TTLs channel is composed of square 
        WAVES, not pulses, meaning that it reaches positive AND NEGATIVE values. 
        If your device  handles only positive voltage, use a diode on the input 
        of your device.
        
        https://en.wikipedia.org/wiki/Diode
    """
    
    SBOutAmpF = Hdf5.DataLoad('/'+System+'/SBOutAmpF', CalibrationFile)[0]
    
    if TTLAmpF > 1/SBOutAmpF:
        print('AmpF out of range. Decreasing to', 1/SBOutAmpF, '.')
        TTLAmpF = 1/SBOutAmpF
    
    if type(NoiseFrequency[0]) == list:
        SoundPulse = Noise(Rate, SoundPulseDur)
        print('   ', end='')
        SoundPulseFiltered = BandpassFilterSound(SoundPulse, Rate, NoiseFrequency)
        print('   ', end='')
        SoundUnit = ApplySoundAmpF(SoundPulseFiltered, Rate, SoundAmpF, 
                                   NoiseFrequency, SBOutAmpF, SoundPauseBeforePulseDur, 
                                   SoundPauseAfterPulseDur)
    else:
        print('Generating tones... ', end='')
        SoundPulseFiltered = {}
        for Freq in NoiseFrequency:
            FKey = str(Freq)
            SoundPulseFiltered[FKey] = SineWave(Rate, Freq, 1, SoundPulseDur)
        print('Done.')
        
        SoundUnit = ApplySoundAmpF(SoundPulseFiltered, Rate, SoundAmpF, 
                                   NoiseFrequency, SBOutAmpF, SoundPauseBeforePulseDur, 
                                   SoundPauseAfterPulseDur)
        
        
    if TTLs:
        SoundTTLUnit = SqWave(Rate, SoundPulseDur, TTLAmpF, SoundTTLVal, 
                              SBOutAmpF, SoundPauseBeforePulseDur, 
                              SoundPauseAfterPulseDur)
    else:
        SoundTTLUnit = np.zeros((round(Rate*SoundPulseDur)), dtype='float32')
    
    Sound = {}
    for FKey in SoundUnit:
        Sound[FKey] = {}
        
        for AKey in SoundUnit[FKey]:
            if Map[0] == 2:
                Sound[FKey][AKey] = np.vstack((SoundUnit[FKey][AKey], SoundTTLUnit)).T
                Sound[FKey][AKey] = np.ascontiguousarray(Sound[FKey][AKey])
            else:
                Sound[FKey][AKey] = np.vstack((SoundTTLUnit, SoundUnit[FKey][AKey])).T
                Sound[FKey][AKey] = np.ascontiguousarray(Sound[FKey][AKey])
    
    print('Done generating sound stimulus.')
    return(Sound)


def SoundStim(Rate, SoundPulseDur, SoundAmpF, NoiseFrequency, TTLAmpF, 
              System, SoundPauseBeforePulseDur=0, SoundPauseAfterPulseDur=0, TTLs=True,
              Map=[1,2], **Kws):
    """ Generate sound pulses in one channel and TTLs in the other channel.
        
        WARNING: The signal generated in the TTLs channel is composed of square 
        WAVES, not pulses, meaning that it reaches positive AND NEGATIVE values. 
        If your device  handles only positive voltage, use a diode on the input 
        of your device.
        
        https://en.wikipedia.org/wiki/Diode
    """
    
    SBOutAmpF = Hdf5.DataLoad('/'+System+'/SBOutAmpF', CalibrationFile)[0]
    
    if TTLAmpF > 1/SBOutAmpF:
        print('AmpF out of range. Decreasing to', 1/SBOutAmpF, '.')
        TTLAmpF = 1/SBOutAmpF
    
    if type(NoiseFrequency[0]) == list:
        SoundPulse = Noise(Rate, SoundPulseDur)
        print('   ', end='')
        SoundPulseFiltered = BandpassFilterSound(SoundPulse, Rate, NoiseFrequency)
        print('   ', end='')
        SoundUnit = ApplySoundAmpF(SoundPulseFiltered, Rate, SoundAmpF, 
                                   NoiseFrequency, SBOutAmpF, SoundPauseBeforePulseDur, 
                                   SoundPauseAfterPulseDur)
    else:
        print('Generating tones... ', end='')
        SoundPulseFiltered = {}
        for Freq in NoiseFrequency:
            FKey = str(Freq)
            SoundPulseFiltered[FKey] = SineWave(Rate, Freq, 1, SoundPulseDur)
        print('Done.')
        
        SoundUnit = ApplySoundAmpF(SoundPulseFiltered, Rate, SoundAmpF, 
                                   NoiseFrequency, SBOutAmpF, SoundPauseBeforePulseDur, 
                                   SoundPauseAfterPulseDur)
        
        
    if TTLs:
        SoundTTLUnit = SqWave(Rate, SoundPulseDur, TTLAmpF, SoundTTLVal, 
                              SBOutAmpF, SoundPauseBeforePulseDur, 
                              SoundPauseAfterPulseDur)
    else:
        SoundTTLUnit = np.zeros((round(Rate*SoundPulseDur)), dtype='float32')
    
    Sound = {}
    for FKey in SoundUnit:
        Sound[FKey] = {}
        
        for AKey in SoundUnit[FKey]:
            if Map[0] == 2:
                Sound[FKey][AKey] = np.vstack((SoundUnit[FKey][AKey], SoundTTLUnit)).T
                Sound[FKey][AKey] = np.ascontiguousarray(Sound[FKey][AKey])
            else:
                Sound[FKey][AKey] = np.vstack((SoundTTLUnit, SoundUnit[FKey][AKey])).T
                Sound[FKey][AKey] = np.ascontiguousarray(Sound[FKey][AKey])
    
    print('Done generating sound stimulus.')
    return(Sound)


## Level 2
def GPIAS(Rate, SoundBGDur, SoundGapDur, SoundBGPrePulseDur, 
          SoundLoudPulseDur, SoundBGAfterPulseDur, 
          SoundBetweenStimDur, SoundBGAmpF, SoundPulseAmpF, TTLAmpF, 
          NoiseFrequency, System, Map=[2,1], **Kws):
    
    Sound = {}
    print('Creating SoundBG...')
    Sound['BG'] = SoundStim(Rate, SoundBGDur, SoundBGAmpF, NoiseFrequency, 
                            TTLAmpF, System, TTLs=False, Map=Map)
    
    print('Creating SoundGap...')
    Sound['Gap'] = {}
    Sound['Gap']['NoGap'] = SoundStim(Rate, SoundGapDur, SoundBGAmpF, 
                                      NoiseFrequency, TTLAmpF, System, 
                                      TTLs=False, Map=Map)
    Sound['Gap']['Gap'] = {
        FKey: {
            AKey: np.zeros(Sound['Gap']['NoGap'][FKey][AKey].shape, dtype='float32')
            for AKey in Sound['Gap']['NoGap'][FKey]
        }
        for FKey in Sound['Gap']['NoGap']
    }
    
    print('Creating SoundBGPrePulse...')
    Sound['BGPrePulse'] = SoundStim(Rate, SoundBGPrePulseDur, SoundBGAmpF, 
                                    NoiseFrequency, TTLAmpF, System, 
                                    TTLs=False, Map=Map)
    
    print('Creating SoundLoudPulse...')
    Sound['LoudPulse'] = SoundStim(Rate, SoundLoudPulseDur, SoundPulseAmpF, 
                                   NoiseFrequency, TTLAmpF, System, Map=Map)
    
    print('Creating SoundBGAfterPulse...')
    Sound['BGAfterPulse'] = SoundStim(Rate, SoundBGAfterPulseDur, SoundBGAmpF, 
                                      NoiseFrequency, TTLAmpF, System, 
                                      TTLs=False, Map=Map)
    
    print('Creating SoundBetweenStim...')
    Sound['BetweenStim'] = SoundStim(Rate, max(SoundBetweenStimDur), 
                                     SoundBGAmpF, NoiseFrequency, TTLAmpF, 
                                     System, TTLs=False, Map=Map)
    
    return(Sound)


#def GPIAS(FileList, CalibrationFile, GPIASTimeBeforeTTL, GPIASTimeAfterTTL, FilterLow, 
#          FilterHigh, FilterOrder):
#    """ Analyze GPIAS recorded using sound board input. """
#    
#    for File in FileList:
#        with shelve.open(CalibrationFile) as Shelve:
#            SoundRec = Shelve['SoundRec']
#            FakeTTLs = Shelve['FakeTTLs']
#            DataInfo = Shelve['DataInfo']
#        
#        with shelve.open(CalibrationFile) as Shelve:
#            SBInAmpF = Shelve['SBInAmpF']
#        
#        NoOfSamplesBefore = int(round((GPIASTimeBeforeTTL*DataInfo['Rate'])*10**-3))
#        NoOfSamplesAfter = int(round((GPIASTimeAfterTTL*DataInfo['Rate'])*10**-3))
#        NoOfSamples = NoOfSamplesBefore + NoOfSamplesAfter
#        XValues = (range(-NoOfSamplesBefore, 
#                         NoOfSamples-NoOfSamplesBefore)/DataInfo['Rate'])*10**3
#        
#        print('Preallocate memory...')
#        GPIAS = [[0] for _ in range(len(DataInfo['NoiseFrequency']))]
#        for Hz in range(len(DataInfo['NoiseFrequency'])):
#            GPIAS[Hz] = [[0] for _ in range(DataInfo['NoOfTrials']*2)]
#        
#        AllTTLs = [[0] for _ in range(len(DataInfo['NoiseFrequency']))]
#        for Hz in range(len(DataInfo['NoiseFrequency'])):
#            AllTTLs[Hz] = [[0] for _ in range(DataInfo['NoOfTrials']*2)]
#        
#        for Hz in range(len(SoundRec)):        
#            for Trial in range(DataInfo['NoOfTrials']*2):
#                SoundStart = ((FakeTTLs[Hz][Trial][0] * 
#                              len(SoundRec[Hz][Trial][0]))//4)-1
#                SoundEnd = ((FakeTTLs[Hz][Trial][1] * 
#                            len(SoundRec[Hz][Trial][0]))//4)-1
#                
#                Start = int(SoundStart-NoOfSamplesBefore)
#                End = int(SoundEnd+NoOfSamplesAfter)
#        
#                GPIAS[Hz][Trial] = array.array('f', b''.join(SoundRec[Hz][Trial]))
#                GPIAS[Hz][Trial] = GPIAS[Hz][Trial][Start:End]
#                GPIAS[Hz][Trial] = [_/SBInAmpF for _ in GPIAS[Hz][Trial]]
#                GPIAS[Hz][Trial] = abs(signal.hilbert(GPIAS[Hz][Trial]))
#                
#                passband = [FilterLow/(DataInfo['Rate']/2), 
#                            FilterHigh/(DataInfo['Rate']/2)]
#                f2, f1 = signal.butter(FilterOrder, passband, 'bandpass')
#                GPIAS[Hz][Trial] = signal.filtfilt(f2, f1, GPIAS[Hz][Trial], 
#                                                 padtype='odd', padlen=0)
#                
#                AllTTLs[Hz][Trial] = [SoundStart, SoundEnd]
#            
#            gData = GPIAS[Hz][:]
#            NoGapAll = [gData[_] for _ in range(len(gData)) if _%2 == 0]
#            GapAll = [gData[_] for _ in range(len(gData)) if _%2 != 0]
#            NoGapSum = list(map(sum, zip(*NoGapAll)))
#            GapSum = list(map(sum, zip(*GapAll)))
#            
#            gData = [0, 0]
#            gData[0] = [_/DataInfo['NoOfTrials'] for _ in NoGapSum]
#            gData[1] = [_/DataInfo['NoOfTrials'] for _ in GapSum]
#            gData[0] = signal.savgol_filter(gData[0], 5, 2, mode='nearest')
#            gData[1] = signal.savgol_filter(gData[1], 5, 2, mode='nearest')
#            GPIAS[Hz] = gData[:]
#            
#            
#            tData = AllTTLs[Hz][:]
#            TTLNoGapAll = [tData[_] for _ in range(len(tData)) if _%2 == 0]
#            TTLGapAll = [tData[_] for _ in range(len(tData)) if _%2 != 0]
#            TTLNoGapSum = list(map(sum, zip(*TTLNoGapAll)))
#            TTLGapSum = list(map(sum, zip(*TTLGapAll)))
#            
#            tData = [0, 0]
#            tData[0] = [int(round(_/DataInfo['NoOfTrials'])) for _ in TTLNoGapSum]
#            tData[1] = [int(round(_/DataInfo['NoOfTrials'])) for _ in TTLGapSum]
#            AllTTLs[Hz] = tData[:]
#            del(NoGapAll, GapAll, NoGapSum, GapSum, gData, tData)
#        
#        print('Saving data to ' + DataInfo['FileName'])
#        with shelve.open(DataInfo['FileName']) as Shelve:
#            Shelve['GPIAS'] = GPIAS
#            Shelve['AllTTLs'] = AllTTLs
#            Shelve['XValues'] = XValues
#        print('Done.')
#    
#    return(None)


#def PlotGPIAS(FileList):
#    Params = {'backend': 'TkAgg'}
#    from matplotlib import rcParams; rcParams.update(Params)
#    from matplotlib import pyplot as plt
#    
#    for File in FileList:
#        print('Loading data from ', File, ' ...')
#        with shelve.open(File[:-3]) as Shelve:
#            GPIAS = Shelve['GPIAS']
#            AllTTLs = Shelve['AllTTLs']
#            XValues = Shelve['XValues']
#            DataInfo = Shelve['DataInfo']
#        
#        print('Plotting...')
#        for Freq in range(len(DataInfo['NoiseFrequency'])):
#            FigTitle = str(DataInfo['NoiseFrequency'][Freq]) + '\ Hz'
#            Line0Label = 'No\ Gap'; Line1Label = 'Gap'
#            SpanLabel = 'Sound\ Pulse'
#            XLabel = 'time\ [ms]'; YLabel = 'voltage\ [mV]'
#            
#            plt.figure(Freq)
#            plt.plot(XValues, GPIAS[Freq][0], 
#                     color='r', label='$'+Line0Label+'$')
#            plt.plot(XValues, GPIAS[Freq][1], 
#                     color='b', label='$'+Line1Label+'$')
#            plt.axvspan(XValues[AllTTLs[Freq][0][0]], XValues[AllTTLs[Freq][0][1]], 
#                        color='k', alpha=0.5, lw=0, label='$'+SpanLabel+'$')
##            plt.axvspan(XValues[AllTTLs[Freq][1][0]], XValues[AllTTLs[Freq][1][1]], 
##                        color='b', alpha=0.5, lw=0, label='Sound pulse (Gap)')
#            plt.suptitle('$'+FigTitle+'$')
#            plt.ylabel('$'+YLabel+'$'); plt.xlabel('$'+XLabel+'$')
#            plt.legend(loc='lower right')
#            plt.locator_params(tight=True)
#            plt.axes().spines['right'].set_visible(False)
#            plt.axes().spines['top'].set_visible(False)
#            plt.axes().yaxis.set_ticks_position('left')
#            plt.axes().xaxis.set_ticks_position('bottom')
#            plt.savefig('Figs/' + File[:-3] + '-' + 
#                        str(DataInfo['NoiseFrequency'][Freq][0]) + '_' + 
#                        str(DataInfo['NoiseFrequency'][Freq][1]) + '.svg', 
#                        format='svg')
#        print('Done.')
#    
#    return(None)


