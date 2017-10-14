#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jun 12 16:02:43 2017

@author: malfatti
"""
import os
import numpy as np

from IO import Hdf5
from DataAnalysis.DataAnalysis import FilterSignal


# Use one that was used in SoundBoardCalibration.py
CalibrationFile = os.environ['DATAPATH']+'/Tests/SoundMeasurements/SoundMeasurements.hdf5'
SoundTTLVal = 0.6; LaserTTLVal = 0.3

## Level 0
def ApplySoundAmpF(SoundPulseFiltered, Rate, SoundAmpF, NoiseFrequency, 
                   SBOutAmpF, SoundPrePauseDur=0, SoundPostPauseDur=0):
    print('Applying amplification factors...')
    # Preallocating memory
    SoundPrePause = np.zeros(round(Rate * SoundPrePauseDur), dtype=np.float32)
    SoundPostPause = np.zeros(round(Rate * SoundPostPauseDur), dtype=np.float32)
    SoundUnit = {}
    
    for FKey in SoundPulseFiltered:
#    for Freq in range(len(NoiseFrequency)):
        SoundUnit[FKey] = {}
        
        for AmpF in range(len(SoundAmpF[FKey])):       
            if SoundAmpF[FKey][AmpF] > 1/SBOutAmpF:
                print(SoundAmpF[FKey][AmpF], 
                      'AmpF out of range. Decreasing to', 1/SBOutAmpF, '.')
                SoundAmpF[FKey][AmpF] = 1/SBOutAmpF
            
            AKey = str(SoundAmpF[FKey][AmpF])
            SoundUnit[FKey][AKey] = np.concatenate([SoundPrePause, 
                                                    SoundPulseFiltered[FKey],
                                                    SoundPostPause])
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


def dBToAmpF(Intensities, Path):
    print('Converting dB to AmpF...')
    SoundIntensity = Hdf5.SoundMeasurementLoad(CalibrationFile, Path, 
                                                'SoundIntensity')
    
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


def SineWave(Rate, Freq, AmpF, SoundPulseDur):
    print('Generating sine wave...')
    Pulse = [np.sin(2 * np.pi * Freq * (_/Rate)) * AmpF
             for _ in range(round(Rate*SoundPulseDur))]
    Pulse[-1] = 0
    Pulse = np.array(Pulse, dtype=np.float32)
    
    return(Pulse)


def SqPulse(Rate, PulseDur, TTLAmpF, TTLVal, SBOutAmpF, PrePauseDur=0, 
           PostPauseDur=0):
    Pulse = np.concatenate([
                np.zeros(PrePauseDur*Rate, dtype=np.float32),
                np.ones(PulseDur*Rate, dtype=np.float32),
                np.zeros(PostPauseDur*Rate, dtype=np.float32)])
    
    return(Pulse)


def SqWave(Rate, PulseDur, TTLAmpF, TTLVal, SBOutAmpF, PrePauseDur=0, 
           PostPauseDur=0):
    
    print('Generating Sound TTL...')
    TTLSpace = PulseDur + PostPauseDur
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
    
    if PrePauseDur == 0:
        if PostPauseDur == 0:
            TTLUnit = TTLPulse
        else:
            TTLPostPause = np.zeros(round((PostPauseDur-PulseDur) * Rate), 
                            dtype=np.float32)
            TTLUnit = np.concatenate([TTLPulse, TTLPostPause])
    else:
        TTLPrePause = np.zeros(round(PrePauseDur * Rate), dtype=np.float32)
        if PostPauseDur == 0:
            TTLUnit = np.concatenate([TTLPrePause, TTLPulse])
        else:
            TTLPostPause = np.zeros(round((PostPauseDur-PulseDur) * Rate), 
                            dtype=np.float32)
            TTLUnit = np.concatenate([TTLPrePause, TTLPulse, TTLPostPause])
    
    TTLUnit = (TTLUnit * TTLAmpF) * SBOutAmpF
    
    return(TTLUnit)


## Level 1
def LaserSinStim(Rate, Dur, Freq, SoundSystem, Ch=1):
    """ Generate square waves in one channel that works as TTLs for laser.
        
        WARNING: The signal generated is composed of square WAVES, not pulses,
        meaning that it reaches positive AND NEGATIVE values. If your device 
        handles only positive voltage, use a diode on the input of your device.
        
        https://en.wikipedia.org/wiki/Diode
        """
    
    SBOutAmpF = Hdf5.DataLoad('/'+SoundSystem+'/SBOutAmpF', CalibrationFile)[0]
    
    if TTLAmpF > 1/SBOutAmpF:
        print('AmpF out of range. Decreasing to', 1/SBOutAmpF, '.')
        TTLAmpF = 1/SBOutAmpF
    
    LaserUnit = SqWave(Rate, LaserPulseDur, TTLAmpF, LaserTTLVal, SBOutAmpF, LaserPrePauseDur, LaserPostPauseDur)
    
    Laser = np.zeros((LaserUnit.shape[0], 2))
    Laser[:,Ch-1] = LaserUnit.T
    Laser = np.ascontiguousarray(Laser)
    
    print('Done generating laser stimulus.')
    return(Laser)


def LaserSqStim(Rate, LaserPulseDur, TTLAmpF, SoundSystem, LaserPrePauseDur=0, LaserPostPauseDur=0, Ch=1):
    """ Generate square waves in one channel that works as TTLs for laser.
        
        WARNING: The signal generated is composed of square WAVES, not pulses,
        meaning that it reaches positive AND NEGATIVE values. If your device 
        handles only positive voltage, use a diode on the input of your device.
        
        https://en.wikipedia.org/wiki/Diode
        """
    
    SBOutAmpF = Hdf5.DataLoad('/'+SoundSystem+'/SBOutAmpF', CalibrationFile)[0]
    
    if TTLAmpF > 1/SBOutAmpF:
        print('AmpF out of range. Decreasing to', 1/SBOutAmpF, '.')
        TTLAmpF = 1/SBOutAmpF
    
    LaserUnit = SqWave(Rate, LaserPulseDur, TTLAmpF, LaserTTLVal, SBOutAmpF, LaserPrePauseDur, LaserPostPauseDur)
    
    Laser = np.zeros((LaserUnit.shape[0], 2))
    Laser[:,Ch-1] = LaserUnit.T
    Laser = np.ascontiguousarray(Laser)
    
    print('Done generating laser stimulus.')
    return(Laser)


def SoundLaserStim(Rate, SoundPulseDur, SoundAmpF, NoiseFrequency, LaserPulseDur, TTLAmpF, SoundSystem, SoundPrePauseDur=0,
                   SoundPostPauseDur=0, LaserPrePauseDur=0, LaserPostPauseDur=0, Map=[1,2]):
    """ Generate sound pulses in one channel and a mix of square waves that 
        works as TTLs for both sound and laser in the other channel.
        
        WARNING: The signal generated in the TTLs channel is composed of square 
        WAVES, not pulses, meaning that it reaches positive AND NEGATIVE values. 
        If your device  handles only positive voltage, use a diode on the 
        input of your device.
        
        https://en.wikipedia.org/wiki/Diode
    """
    
    SBOutAmpF = Hdf5.DataLoad('/'+SoundSystem+'/SBOutAmpF', CalibrationFile)[0]
    
    if TTLAmpF > 1/SBOutAmpF:
        print('AmpF out of range. Decreasing to', 1/SBOutAmpF, '.')
        TTLAmpF = 1/SBOutAmpF
    
    SoundPulse = Noise(Rate, SoundPulseDur)
    print('   ', end='')
    SoundPulseFiltered = BandpassFilterSound(SoundPulse, Rate, NoiseFrequency)
    print('   ', end='')
    SoundUnit = ApplySoundAmpF(SoundPulseFiltered, Rate, SoundAmpF, 
                               NoiseFrequency, SBOutAmpF, SoundPrePauseDur, 
                               SoundPostPauseDur)
    
    SoundTTLUnit = SqWave(Rate, SoundPulseDur, TTLAmpF, SoundTTLVal, 
                          SBOutAmpF, SoundPrePauseDur, 
                          SoundPostPauseDur)
    
    LaserUnit = SqWave(Rate, LaserPulseDur, TTLAmpF, LaserTTLVal, SBOutAmpF, 
                       LaserPrePauseDur, LaserPostPauseDur)
    
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
              SoundSystem, SoundPrePauseDur=0, SoundPostPauseDur=0, TTLs=True,
              Map=[1,2]):
    """ Generate sound pulses in one channel and TTLs in the other channel.
        
        WARNING: The signal generated in the TTLs channel is composed of square 
        WAVES, not pulses, meaning that it reaches positive AND NEGATIVE values. 
        If your device  handles only positive voltage, use a diode on the input 
        of your device.
        
        https://en.wikipedia.org/wiki/Diode
    """
    
    SBOutAmpF = Hdf5.DataLoad('/'+SoundSystem+'/SBOutAmpF', CalibrationFile)[0]
    
    if TTLAmpF > 1/SBOutAmpF:
        print('AmpF out of range. Decreasing to', 1/SBOutAmpF, '.')
        TTLAmpF = 1/SBOutAmpF
    
    SoundPulse = Noise(Rate, SoundPulseDur)
    print('   ', end='')
    SoundPulseFiltered = BandpassFilterSound(SoundPulse, Rate, NoiseFrequency)
    print('   ', end='')
    SoundUnit = ApplySoundAmpF(SoundPulseFiltered, Rate, SoundAmpF, 
                               NoiseFrequency, SBOutAmpF, SoundPrePauseDur, 
                               SoundPostPauseDur)
    if TTLs:
        SoundTTLUnit = SqWave(Rate, SoundPulseDur, TTLAmpF, SoundTTLVal, 
                              SBOutAmpF, SoundPrePauseDur, 
                              SoundPostPauseDur)
    else:
        SoundTTLUnit = np.zeros((SoundPulse.size), dtype='float32')
    
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
def GPIASStim(Rate, SoundBackgroundDur, SoundGapDur, SoundBackgroundPrePulseDur, 
              SoundLoudPulseDur, SoundBackgroundAfterPulseDur, 
              SoundBetweenStimDur, SoundBackgroundAmpF, SoundPulseAmpF, TTLAmpF, 
              NoiseFrequency, SoundBoard, Map=[2,1]):
    print('Creating SoundBackground...')
    SoundBackground = SoundStim(Rate, SoundBackgroundDur, SoundBackgroundAmpF, 
                                NoiseFrequency, TTLAmpF, SoundBoard, TTLs=False, 
                                Map=Map)
    
    print('Creating SoundGap...')
    SoundGap = {}
    SoundGap['NoGap'] = SoundStim(Rate, SoundGapDur, SoundBackgroundAmpF, 
                                  NoiseFrequency, TTLAmpF, SoundBoard, 
                                  TTLs=False, Map=Map)
    
    SoundGap['Gap'] = {FKey: {AKey: np.zeros(SoundGap['NoGap'][FKey][AKey].shape, 
                                             dtype='float32')
                              for AKey in SoundGap['NoGap'][FKey]} 
                       for FKey in SoundGap['NoGap']}
    
    print('Creating SoundBackgroundPrePulse...')
    SoundBackgroundPrePulse = SoundStim(Rate, SoundBackgroundPrePulseDur, 
                                        SoundBackgroundAmpF, NoiseFrequency, 
                                        TTLAmpF, SoundBoard, TTLs=False, Map=Map)
    
    print('Creating SoundLoudPulse...')
    SoundLoudPulse = SoundStim(Rate, SoundLoudPulseDur, SoundPulseAmpF, 
                               NoiseFrequency, TTLAmpF, SoundBoard, Map=Map)
    
    print('Creating SoundBackgroundAfterPulse...')
    SoundBackgroundAfterPulse = SoundStim(Rate, SoundBackgroundAfterPulseDur, 
                                          SoundBackgroundAmpF, NoiseFrequency, 
                                          TTLAmpF, SoundBoard, TTLs=False, Map=Map)
    
    print('Creating SoundBetweenStim...')
    SoundBetweenStim = SoundStim(Rate, max(SoundBetweenStimDur), 
                                 SoundBackgroundAmpF, NoiseFrequency, TTLAmpF, 
                                 SoundBoard, TTLs=False, Map=Map)
    
    return(SoundBackground, SoundGap, SoundBackgroundPrePulse, SoundLoudPulse, 
           SoundBackgroundAfterPulse, SoundBetweenStim)


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


