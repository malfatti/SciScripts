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

This is a script that defines functions allowing the use of a computer's sound
board as an analog I/O board.

"""

import array
import Hdf5F
import math
import matplotlib.animation as animation
import matplotlib.pyplot as plt
import pyaudio
import shelve
import random
from scipy import signal
import threading

SoundTTLVal = 0.6; LaserTTLVal = 0.3
SBAmpFsFile = '/home/cerebro/Malfatti/Data/Test/20160418173048-SBAmpFs.hdf5'

## Lower-level functions
def dBToAmpF(Intensities, CalibrationFile):
    print('Converting dB to AmpF...')
    SoundIntensity = Hdf5F.SoundMeasurement(CalibrationFile, 
                                                    'SoundIntensity')
    
    SoundAmpF = {Hz: [float(min(SoundIntensity[Hz].keys(), 
                                key=lambda i: abs(SoundIntensity[Hz][i]-dB))) 
                      for dB in Intensities]
                 for Hz in list(SoundIntensity)}
    
    return(SoundAmpF)


def GenNoise(Rate, SoundPulseDur):
    print('Generating noise...')
    SoundNoise = [random.random() \
                  for _ in range(round(Rate*SoundPulseDur))]
    SoundPulse = [SoundNoise[ElI]*2-1 for ElI,ElV in enumerate(SoundNoise)]
    
    return(SoundPulse)


def GenSineWave(Rate, Freq, AmpF, SoundPulseDur):
    print('Generating sine wave...')
    Pulse = [math.sin(2*math.pi*Freq*(_/Rate)) * AmpF
             for _ in range(round(Rate*SoundPulseDur))]
    Pulse[-1] = 0
    
    return(Pulse)


def BandpassFilterSound(SoundPulse, Rate, NoiseFrequency):
    # Preallocating memory
    SoundPulseFiltered = [0]*len(NoiseFrequency)
    
    print('Filtering sound: ', end='')
    for Freq in range(len(NoiseFrequency)):
#        if len(NoiseFrequency[Freq]) == 1:
#            print('Filtering sound: ', NoiseFrequency[Freq], '...')
#            passband = [(NoiseFrequency[Freq][0]-1)/(Rate/2), 
#                        (NoiseFrequency[Freq][0]+1)/(Rate/2)]
#            f2, f1 = signal.butter(4, passband, 'bandpass')
#        else:        
        print(NoiseFrequency[Freq], end='...')
        passband = [_/(Rate/2) for _ in NoiseFrequency[Freq]]
        
        f2, f1 = signal.butter(4, passband, 'bandpass')
        SoundPulseFiltered[Freq] = signal.filtfilt(f2, f1, SoundPulse, \
                                                         padtype='odd', \
                                                         padlen=0)
        SoundPulseFiltered[Freq] = SoundPulseFiltered[Freq].tolist()
        SoundPulseFiltered[Freq][-1] = 0
    print(end='\n')
    return(SoundPulseFiltered)


def ApplySoundAmpF(SoundPulseFiltered, Rate, SoundAmpF, NoiseFrequency, 
                   SBOutAmpF, SoundPrePauseDur=0, SoundPostPauseDur=0):
    print('Applying amplification factors...')
    # Preallocating memory
    SoundPrePause = [0] * round(Rate * SoundPrePauseDur)
    SoundPostPause = [0] * round(Rate * SoundPostPauseDur)
    SoundUnit = [0]*len(NoiseFrequency)
    
    for Freq in range(len(NoiseFrequency)):
        Key= str(NoiseFrequency[Freq][0]) + '-' + str(NoiseFrequency[Freq][1])
        
        SoundUnit[Freq] = [0]*len(SoundAmpF[Key])
        
        for AmpF in range(len(SoundAmpF[Key])):
            SoundUnit[Freq][AmpF] = SoundPrePause + \
                                    SoundPulseFiltered[Freq] + \
                                    SoundPostPause
            Max = max(SoundUnit[Freq][AmpF])
            SoundUnit[Freq][AmpF] = [(SoundEl * (SBOutAmpF/Max))
                                     * SoundAmpF[Key][AmpF]
                                     for SoundEl in SoundUnit[Freq][AmpF]]
    
    return(SoundUnit)


def GenTTL(Rate, PulseDur, TTLAmpF, SoundBoard, SBOutAmpF, PrePauseDur=0, 
           PostPauseDur=0):
    
    print('Generating Sound TTL...')
    
    TTLPulse = [1]*round((Rate*PulseDur)/2) + \
                    [-1]*round((Rate*PulseDur)/2)
    TTLPulse[-1] = 0

    TTLPrePause = [0] * round(PrePauseDur * Rate)
    TTLPostPause = [0] * round(PostPauseDur * Rate)
    
#    if SoundPulseDur < 0.01:
#        SoundTTLPulse = [round(SoundTTLVal/SBOutAmpF, 3)] * round(SoundPulseDur * Rate)
#    else:
#        Middle = [0]*round((SoundPulseDur-0.01) * Rate)
#        Border = [round(SoundTTLVal/SBOutAmpF, 3)] * round(0.005 * Rate)
#        SoundTTLPulse = Border + Middle + Border
    
    TTLUnit = TTLPrePause + TTLPulse + TTLPostPause
    TTLUnit = [TTLEl*TTLAmpF for TTLEl in TTLUnit]
    
    return(TTLUnit)


def InterleaveChannels(Right, Left, SoundAmpF, NoiseFrequency):
    print('Interleaving channels...')
    if len(Right) == 1:
        List = [0]*(2*len(Left))
        for _ in range(len(Left)):
            List[_ *2] = 0
            List[_ *2+1] = Left[_]
    
    elif len(Left) == 1:
        List = [0]*(2*len(Left))
        for _ in range(len(Left)):
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


def ListToByteArray(List, SoundAmpF, NoiseFrequency):
    print('Converting list to byte array...')
    if SoundAmpF == [0]:
        ByteArray = array.array('f', List)
        ByteArray = bytes(ByteArray)
    
    else:
        ByteArray = [0]*len(NoiseFrequency)
        
        for Freq in range(len(NoiseFrequency)):
            Key= str(NoiseFrequency[Freq][0]) + '-' + \
                 str(NoiseFrequency[Freq][1])
            ByteArray[Freq] = [0]*len(SoundAmpF[Key])
            
            for AmpF in range(len(SoundAmpF[Key])):
                ByteArray[Freq][AmpF] = array.array('f', List[Freq][AmpF])
                ByteArray[Freq][AmpF] = bytes(ByteArray[Freq][AmpF])
    
    return(ByteArray)


def GenPause(Rate, PauseBetweenStimBlocksDur):
    print('Generating pause...')
    PauseBetweenStimBlocksList = [0] * \
                                 round(PauseBetweenStimBlocksDur * Rate) * \
                                 2
    PauseBetweenStimBlocks = array.array('f', PauseBetweenStimBlocksList)
    PauseBetweenStimBlocks = bytes(PauseBetweenStimBlocks)
    
    return(PauseBetweenStimBlocks)


def GenAudioObj(Rate, Direction='out', Stream=False):
    print('Generating audio object...')
    if Stream:
        q = pyaudio.PyAudio()
        InOn = False
        RecStop = False
        def InCallBack(in_data, frame_count, time_info, status):
            if InOn:
                global SoundRec
                SoundRec[RealFreq][RealTrial].append(in_data)
                
            if RecStop:
                InFlag = pyaudio.paComplete
            else:
                InFlag = pyaudio.paContinue
            return(None, InFlag)
        
        Reading = q.open(format=pyaudio.paFloat32,
                         channels=1,
                         rate=Rate,
                         input=True,
                         output=False,
                         stream_callback=InCallBack)
        return(Reading, InCallBack)
    
    else:
        if Direction == 'out':
            p = pyaudio.PyAudio()
            Stimulation = p.open(format=pyaudio.paFloat32,
                            channels=2,
                            rate=Rate,
                            output=True)
            return(Stimulation)
            
        elif Direction == 'in':
            q = pyaudio.PyAudio()
            Reading = q.open(format=pyaudio.paFloat32,
                             channels=2,
                             rate=Rate,
                             input=True)
            return(Reading)
            
        elif Direction == 'inout':
            print('To be implemented :)')
            raise SystemExit(0)
        else:
            print("Choose 'in', 'out' or 'inout'")
            raise SystemExit(0)


def RunSound(Array, PauseBetweenStimBlocks, PulseNo, Stimulation, SoundAmpF, 
             NoiseFrequency, Complexity='All', StimBlockNo=1, 
             Multithread=False):
    if Multithread:
        if Complexity == 'AllPulses':
            class Play(threading.Thread):
                def run(self):
                    for OnePulse in range(PulseNo):
                        Stimulation.write(Array)
            return(Play)
        
        elif Complexity == 'All':
            class Play(threading.Thread):
                def run(self):
                    for Freq in range(len(NoiseFrequency)):
                        if len(NoiseFrequency[Freq]) == 1:
                            Key= str(NoiseFrequency[Freq][0])
                        else:
                            Key= str(NoiseFrequency[Freq][0]) + '-' \
                                 + str(NoiseFrequency[Freq][1])
                        
                        for AmpF in range(len(SoundAmpF[Key])):
                            for OneBlock in range(StimBlockNo):
                                for OnePulse in range(PulseNo):
                                    Stimulation.write(Array[Freq][AmpF])
                        
                                Stimulation.write(PauseBetweenStimBlocks)
            return(Play)
        else:
            print("Choose 'Allpulses' or 'All'.")
            raise SystemExit(0)
    else:
        if Complexity == 'AllPulses':
            if SoundAmpF == [0]:
                def Play():
                    for OnePulse in range(PulseNo):
                        Stimulation.write(Array)
                    
                return(Play)
            else:
                def Play(Freq, AmpF):
                    for OnePulse in range(PulseNo):
                        Stimulation.write(Array[Freq][AmpF])
                    
                return(Play)
        
        elif Complexity == 'AllBlocks':
            if SoundAmpF == [0]:
                def Play():
                    for OneBlock in range(StimBlockNo):
                        for OnePulse in range(PulseNo):
                            Stimulation.write(Array)
                        
                        Stimulation.write(PauseBetweenStimBlocks)
                return(Play)
            else:
                def Play(Freq, AmpF):
                    for OneBlock in range(StimBlockNo):
                        for OnePulse in range(PulseNo):
                            Stimulation.write(Array[Freq][AmpF])
                
                        Stimulation.write(PauseBetweenStimBlocks)
                return(Play)
        
        elif Complexity == 'AllAmpFs':
            def Play(Freq):
                if len(NoiseFrequency[Freq]) == 1:
                    Key= str(NoiseFrequency[Freq][0])
                else:
                    Key= str(NoiseFrequency[Freq][0]) + '-' \
                         + str(NoiseFrequency[Freq][1])
                for AmpF in range(len(SoundAmpF[Key])):
                    for OneBlock in range(StimBlockNo):
                        for OnePulse in range(PulseNo):
                            Stimulation.write(Array[Freq][AmpF])
                
                        Stimulation.write(PauseBetweenStimBlocks)
            return(Play)
        
        elif Complexity == 'All':
            def Play():
                for Freq in range(len(NoiseFrequency)):
                    if len(NoiseFrequency[Freq]) == 1:
                        Key= str(NoiseFrequency[Freq][0])
                    else:
                        Key= str(NoiseFrequency[Freq][0]) + '-' \
                             + str(NoiseFrequency[Freq][1])
                    
                    for AmpF in range(len(SoundAmpF[Key])):
                        for OneBlock in range(StimBlockNo):
                            for OnePulse in range(PulseNo):
                                Stimulation.write(Array[Freq][AmpF])
                    
                            Stimulation.write(PauseBetweenStimBlocks)
            return(Play)
        
        else:
            print("Choose 'Allpulses', 'AllBlocks', 'AllAmpFs' or 'All'.")
            raise SystemExit(0)


## Higher-level functions
def SoundStim(Rate, SoundPulseDur, SoundPulseNo, SoundAmpF, NoiseFrequency, 
              TTLAmpF, SoundBoard, Complexity='All', SoundPrePauseDur=0, 
              SoundPostPauseDur=0, SoundStimBlockNo=1, 
              SoundPauseBetweenStimBlocksDur=0):
    """ Generate sound pulses in one channel and TTLs in the other channel 
    (Check ControlArduinoWithSoundBoard.ino code)."""
    
    SBOutAmpF = Hdf5F.SoundCalibration(SBAmpFsFile, SoundBoard,
                                               'SBOutAmpF')
    
    SoundPulse = GenNoise(Rate, SoundPulseDur)
    print('   ', end='')
    SoundPulseFiltered = BandpassFilterSound(SoundPulse, Rate, NoiseFrequency)
    print('   ', end='')
    SoundUnit = ApplySoundAmpF(SoundPulseFiltered, Rate, SoundAmpF, 
                               NoiseFrequency, SBOutAmpF, SoundPrePauseDur, 
                               SoundPostPauseDur)
    SoundTTLUnit = GenTTL(Rate, SoundPulseDur, TTLAmpF, SoundBoard, SBOutAmpF, 
                          SoundPrePauseDur, SoundPostPauseDur)    
    SoundList = InterleaveChannels(SoundUnit, SoundTTLUnit, SoundAmpF, 
                                   NoiseFrequency)    
    Sound = ListToByteArray(SoundList, SoundAmpF, NoiseFrequency)
    SoundPauseBetweenStimBlocks = GenPause(Rate, SoundPauseBetweenStimBlocksDur)    
    
    Stimulation = GenAudioObj(Rate, 'out')
    PlaySound = RunSound(Sound, SoundPauseBetweenStimBlocks, SoundPulseNo, 
                         Stimulation, SoundAmpF, NoiseFrequency, Complexity, 
                         SoundStimBlockNo, Multithread=False)
    
    print('Done generating sound stimulus.')
    return(Sound, PlaySound)


def LaserStim(Rate, LaserPulseDur, LaserPulseNo, TTLAmpF, SoundBoard, 
              Complexity='AllBlocks', LaserPrePauseDur=0, LaserPostPauseDur=0, 
              LaserStimBlockNo=1, LaserPauseBetweenStimBlocksDur=0):
    """ Generate square pulses in one channel that works as TTL for laser 
    (Check ControlArduinoWithSoundBoard.ino code)."""
    
    SBOutAmpF = Hdf5F.SoundCalibration(SBAmpFsFile, SoundBoard,
                                               'SBOutAmpF')
    
    LaserUnit = GenTTL(Rate, LaserPulseDur, TTLAmpF, SoundBoard, SBOutAmpF, 
                       LaserPrePauseDur, LaserPostPauseDur)
    LaserList = InterleaveChannels([0], LaserUnit, [0], [0])
    Laser = ListToByteArray(LaserList, [0], [0])
    LaserPauseBetweenStimBlocks = GenPause(Rate, LaserPauseBetweenStimBlocksDur)    
    
    Stimulation = GenAudioObj(Rate, 'out')
    PlayLaser = RunSound(Laser, LaserPauseBetweenStimBlocks, LaserPulseNo, 
                         Stimulation, [0], [0], Complexity, LaserStimBlockNo, 
                         Multithread=False)
    
    print('Done generating laser stimulus.')
    return(Laser, PlayLaser)


def SoundLaserStim(Rate, SoundPulseDur, SoundPulseNo, SoundAmpF, 
                   NoiseFrequency, LaserPulseDur, LaserPulseNo, TTLAmpF, 
                   SoundBoard, Complexity='AllPulses', SoundPrePauseDur=0, 
                   SoundPostPauseDur=0, SoundStimBlockNo=1, 
                   SoundPauseBetweenStimBlocksDur=0, LaserPrePauseDur=0, 
                   LaserPostPauseDur=0, LaserStimBlockNo=1, 
                   LaserPauseBetweenStimBlocksDur=0):
    """ Generate sound pulses in one channel and TTLs for sound and laser in 
    the other channel (Check ControlArduinoWithSoundBoard.ino code)."""
    
    SBOutAmpF = Hdf5F.SoundCalibration(SBAmpFsFile, SoundBoard,
                                               'SBOutAmpF')
    
    LaserUnit = GenTTL(Rate, LaserPulseDur, TTLAmpF, SoundBoard, SBOutAmpF, 
                       LaserPrePauseDur, LaserPostPauseDur)
    
    SoundTTLUnit = GenTTL(Rate, SoundPulseDur, TTLAmpF, SoundBoard, SBOutAmpF, 
                          SoundPrePauseDur, SoundPostPauseDur)
    
    print('Summing sound TTL and laser units...')
    SoundTTLLaserUnit = [LaserUnit[_]+SoundTTLUnit[_] \
                         for _ in range(len(SoundTTLUnit))]
    
    SoundPulse = GenNoise(Rate, SoundPulseDur)
    print('   ', end='')
    SoundPulseFiltered = BandpassFilterSound(SoundPulse, Rate, NoiseFrequency)
    print('   ', end='')
    SoundUnit = ApplySoundAmpF(SoundPulseFiltered, Rate, SoundAmpF, 
                               NoiseFrequency, SBOutAmpF, SoundPrePauseDur, 
                               SoundPostPauseDur)
    
    SoundLaserList = InterleaveChannels(SoundUnit, SoundTTLLaserUnit, 
                                        SoundAmpF, NoiseFrequency)
    SoundLaser = ListToByteArray(SoundLaserList, SoundAmpF, NoiseFrequency)
    SoundLaserPauseBetweenStimBlocks = GenPause(Rate, 
                                                SoundPauseBetweenStimBlocksDur)
    
    Stimulation = GenAudioObj(Rate, 'out')
    PlaySoundLaser = RunSound(SoundLaser, SoundLaserPauseBetweenStimBlocks, 
                              SoundPulseNo, Stimulation, SoundAmpF, 
                              NoiseFrequency, Complexity, SoundStimBlockNo, 
                              Multithread=False)
    
    print('Done generating sound and laser stimulus.')
    return(SoundLaser, PlaySoundLaser)

def GPIAS(FileList, CalibrationFile, GPIASTimeBeforeTTL, GPIASTimeAfterTTL, FilterLow, 
          FilterHigh, FilterOrder):
    """ Analyze GPIAS recorded using sound board input. """
    
    for File in FileList:
        with shelve.open(CalibrationFile) as Shelve:
            SoundRec = Shelve['SoundRec']
            FakeTTLs = Shelve['FakeTTLs']
            DataInfo = Shelve['DataInfo']
        
        with shelve.open(CalibrationFile) as Shelve:
            SBInAmpF = Shelve['SBInAmpF']
        
        NoOfSamplesBefore = int(round((GPIASTimeBeforeTTL*DataInfo['Rate'])*10**-3))
        NoOfSamplesAfter = int(round((GPIASTimeAfterTTL*DataInfo['Rate'])*10**-3))
        NoOfSamples = NoOfSamplesBefore + NoOfSamplesAfter
        XValues = (range(-NoOfSamplesBefore, 
                         NoOfSamples-NoOfSamplesBefore)/DataInfo['Rate'])*10**3
        
        print('Preallocate memory...')
        GPIAS = [[0] for _ in range(len(DataInfo['NoiseFrequency']))]
        for Hz in range(len(DataInfo['NoiseFrequency'])):
            GPIAS[Hz] = [[0] for _ in range(DataInfo['NoOfTrials']*2)]
        
        AllTTLs = [[0] for _ in range(len(DataInfo['NoiseFrequency']))]
        for Hz in range(len(DataInfo['NoiseFrequency'])):
            AllTTLs[Hz] = [[0] for _ in range(DataInfo['NoOfTrials']*2)]
        
        for Hz in range(len(SoundRec)):        
            for Trial in range(DataInfo['NoOfTrials']*2):
                SoundStart = ((FakeTTLs[Hz][Trial][0] * 
                              len(SoundRec[Hz][Trial][0]))//4)-1
                SoundEnd = ((FakeTTLs[Hz][Trial][1] * 
                            len(SoundRec[Hz][Trial][0]))//4)-1
                
                Start = int(SoundStart-NoOfSamplesBefore)
                End = int(SoundEnd+NoOfSamplesAfter)
        
                GPIAS[Hz][Trial] = array.array('f', b''.join(SoundRec[Hz][Trial]))
                GPIAS[Hz][Trial] = GPIAS[Hz][Trial][Start:End]
                GPIAS[Hz][Trial] = [_/SBInAmpF for _ in GPIAS[Hz][Trial]]
                GPIAS[Hz][Trial] = abs(signal.hilbert(GPIAS[Hz][Trial]))
                
                passband = [FilterLow/(DataInfo['Rate']/2), 
                            FilterHigh/(DataInfo['Rate']/2)]
                f2, f1 = signal.butter(FilterOrder, passband, 'bandpass')
                GPIAS[Hz][Trial] = signal.filtfilt(f2, f1, GPIAS[Hz][Trial], 
                                                 padtype='odd', padlen=0)
                
                AllTTLs[Hz][Trial] = [SoundStart, SoundEnd]
            
            gData = GPIAS[Hz][:]
            NoGapAll = [gData[_] for _ in range(len(gData)) if _%2 == 0]
            GapAll = [gData[_] for _ in range(len(gData)) if _%2 != 0]
            NoGapSum = list(map(sum, zip(*NoGapAll)))
            GapSum = list(map(sum, zip(*GapAll)))
            
            gData = [0, 0]
            gData[0] = [_/DataInfo['NoOfTrials'] for _ in NoGapSum]
            gData[1] = [_/DataInfo['NoOfTrials'] for _ in GapSum]
            gData[0] = signal.savgol_filter(gData[0], 5, 2, mode='nearest')
            gData[1] = signal.savgol_filter(gData[1], 5, 2, mode='nearest')
            GPIAS[Hz] = gData[:]
            
            
            tData = AllTTLs[Hz][:]
            TTLNoGapAll = [tData[_] for _ in range(len(tData)) if _%2 == 0]
            TTLGapAll = [tData[_] for _ in range(len(tData)) if _%2 != 0]
            TTLNoGapSum = list(map(sum, zip(*TTLNoGapAll)))
            TTLGapSum = list(map(sum, zip(*TTLGapAll)))
            
            tData = [0, 0]
            tData[0] = [int(round(_/DataInfo['NoOfTrials'])) for _ in TTLNoGapSum]
            tData[1] = [int(round(_/DataInfo['NoOfTrials'])) for _ in TTLGapSum]
            AllTTLs[Hz] = tData[:]
            del(NoGapAll, GapAll, NoGapSum, GapSum, gData, tData)
        
        print('Saving data to ' + DataInfo['FileName'])
        with shelve.open(DataInfo['FileName']) as Shelve:
            Shelve['GPIAS'] = GPIAS
            Shelve['AllTTLs'] = AllTTLs
            Shelve['XValues'] = XValues
        print('Done.')
    
    return(None)


def PlotGPIAS(FileList):
    for File in FileList:
        print('Loading data from ', File, ' ...')
        with shelve.open(File[:-3]) as Shelve:
            GPIAS = Shelve['GPIAS']
            AllTTLs = Shelve['AllTTLs']
            XValues = Shelve['XValues']
            DataInfo = Shelve['DataInfo']
        
        print('Plotting...')
        for Freq in range(len(DataInfo['NoiseFrequency'])):
            FigTitle = str(DataInfo['NoiseFrequency'][Freq]) + '\ Hz'
            Line0Label = 'No\ Gap'; Line1Label = 'Gap'
            SpanLabel = 'Sound\ Pulse'
            XLabel = 'time\ [ms]'; YLabel = 'voltage\ [mV]'
            
            plt.figure(Freq)
            plt.plot(XValues, GPIAS[Freq][0], 
                     color='r', label='$'+Line0Label+'$')
            plt.plot(XValues, GPIAS[Freq][1], 
                     color='b', label='$'+Line1Label+'$')
            plt.axvspan(XValues[AllTTLs[Freq][0][0]], XValues[AllTTLs[Freq][0][1]], 
                        color='k', alpha=0.5, lw=0, label='$'+SpanLabel+'$')
#            plt.axvspan(XValues[AllTTLs[Freq][1][0]], XValues[AllTTLs[Freq][1][1]], 
#                        color='b', alpha=0.5, lw=0, label='Sound pulse (Gap)')
            plt.suptitle('$'+FigTitle+'$')
            plt.ylabel('$'+YLabel+'$'); plt.xlabel('$'+XLabel+'$')
            plt.legend(loc='lower right')
            plt.locator_params(tight=True)
            plt.axes().spines['right'].set_visible(False)
            plt.axes().spines['top'].set_visible(False)
            plt.axes().yaxis.set_ticks_position('left')
            plt.axes().xaxis.set_ticks_position('bottom')
            plt.savefig('Figs/' + File[:-3] + '-' + 
                        str(DataInfo['NoiseFrequency'][Freq][0]) + '_' + 
                        str(DataInfo['NoiseFrequency'][Freq][1]) + '.svg', 
                        format='svg')
        print('Done.')
    
    return(None)


def MicrOscilloscope(Rate, SoundBoard, XLim, YLim, FramesPerBuf=512):
    """ Read data from sound board input and plot it until the windows is 
    closed. """
    
    SBInAmpF = Hdf5F.SoundCalibration(SBAmpFsFile, SoundBoard,
                                              'SBInAmpF')
    
    r = pyaudio.PyAudio()
    
    Plotting = r.open(format=pyaudio.paFloat32,
                         channels=1,
                         rate=Rate,
                         input=True,
                         output=False,
                         frames_per_buffer=FramesPerBuf)
                         #stream_callback=InCallBack)
    
    Fig = plt.figure()
    Ax = plt.axes(xlim=XLim, ylim=YLim)
    Plot, = Ax.plot([float('nan')]*(Rate//10), lw=1)
    
    def AnimInit():
        Data = array.array('f', [])
        Plot.set_ydata(Data)
        return Plot,
    
    def PltUp(n):
        Data = array.array('f', Plotting.read(Rate//10))
        Data = [_/SBInAmpF for _ in Data]
        Plot.set_ydata(Data)
        return Plot,
    
    Anim = animation.FuncAnimation(Fig, PltUp, frames=FramesPerBuf, 
                                   interval=16, blit=False)
    
    return(None)


def MicrOscilloscopeRec(Rate, XLim, YLim, SoundBoard, FramesPerBuf=512):
    """ Read data from sound board input, record it to a video and plot it 
    until the windows is closed (with a delay). """
    
    SBInAmpF = Hdf5F.SoundCalibration(SBAmpFsFile, SoundBoard,
                                              'SBInAmpF')
    
    r = pyaudio.PyAudio()
    
    Plotting = r.open(format=pyaudio.paFloat32,
                         channels=1,
                         rate=Rate,
                         input=True,
                         output=False,
                         frames_per_buffer=FramesPerBuf)
                         #stream_callback=InCallBack)
    
    Writers = animation.writers['ffmpeg']
    Writer = Writers(fps=15, metadata=dict(artist='Me'))
    
    Fig = plt.figure()
    Ax = plt.axes(xlim=XLim, ylim=YLim)
    Plot, = Ax.plot([float('nan')]*(Rate//10), lw=1)
    
    def AnimInit():
        Data = array.array('f', [])
        Plot.set_ydata(Data)
        return Plot,
    
    def PltUp(n):
        Data = array.array('f', Plotting.read(Rate//10))
        Data = [_/SBInAmpF for _ in Data]
        Plot.set_ydata(Data)
        return Plot,
    
    Anim = animation.FuncAnimation(Fig, PltUp, frames=FramesPerBuf, interval=16, blit=False)
    Anim.save('MicrOscilloscope.mp4', writer=Writer)
    
    return(None)


def SoundCalOut(Rate, Freq, WaveDur):
    """ Generate a sine wave from 1 to -1 """
    SoundPulseNo = WaveDur/0.1
    
    Pulse = GenSineWave(Rate, Freq, 1, 0.1)
    List = InterleaveChannels(Pulse, [0], [0], [0])
    Pulse = ListToByteArray(List, [0], [0])
    
    Stimulation = GenAudioObj(Rate, 'out')
    Play = RunSound(Pulse, 0, SoundPulseNo, Stimulation, [0], [0], 
                    'AllPulses', 1, Multithread=False)
    
    print('Playing...')
    Play()
    
    return(None)


def SoundCalIn(Rate, Freq, WaveDur, SBOutAmpF):
    """ Generate sine wave (1V to -1V) and read 1s of it. """
    SoundPulseNo = WaveDur/0.1
    
    Pulse = GenSineWave(Rate, Freq, SBOutAmpF, 0.1)
    List = InterleaveChannels(Pulse, [0], [0], [0])
    Pulse = ListToByteArray(List, [0], [0])
    
    Stimulation = GenAudioObj(Rate, 'out')
    Reading = GenAudioObj(Rate, 'in')
    Play = RunSound(Pulse, 0, SoundPulseNo, Stimulation, [0], [0], 
                    'AllPulses', 1, Multithread=True)
    
    print('Playing...')
    Play().start()
    Data = Reading.read(Rate)
    Data = array.array('f', Data)
    
    return(Data)


def SoundMeasurementOut(Rate, SoundPulseDur, SoundPulseNo, SoundAmpF, 
                        NoiseFrequency, TTLAmpF, SBOutAmpF):
    """ This will generate the sound pulses and the sound output objects.
    Remember to create input objects and call them accordingly. """
    
    print('Generating Sound TTL...')
    SoundTTLPulse = [SoundTTLVal] * round(SoundPulseDur * Rate)
    SoundTTLPulse[-1] = 0
    
    SoundTTLUnit = [round(_/SBOutAmpF, 3)*TTLAmpF for _ in SoundTTLPulse]
    
    print('Generating sound pulse...')         
    SoundNoise = [random.random() 
                  for _ in range(round(Rate*SoundPulseDur))]
    SoundPulse = [_*2-1 for _ in SoundNoise]
    
    # Preallocating memory
    SoundPulseFiltered = [0]*len(NoiseFrequency)
    SoundUnit = [0]*len(NoiseFrequency)
    SoundList = [0]*len(NoiseFrequency)
    Sound = [0]*len(NoiseFrequency)
    SoundRec = [0]*len(NoiseFrequency)
    
    for Freq in range(len(NoiseFrequency)):        
        print('Filtering sound: ', NoiseFrequency[Freq], '...')
        passband = [_/(Rate/2) for _ in NoiseFrequency[Freq]]
        f2, f1 = signal.butter(4, passband, 'bandpass')
        SoundPulseFiltered[Freq] = signal.filtfilt(f2, f1, SoundPulse, 
                                                         padtype='odd', 
                                                         padlen=0)
        SoundPulseFiltered[Freq] = SoundPulseFiltered[Freq].tolist()
        SoundPulseFiltered[Freq][-1] = 0
        
        # Preallocating memory
        SoundUnit[Freq] = [0]*len(SoundAmpF)
        SoundList[Freq] = [0]*len(SoundAmpF)
        Sound[Freq] = [0]*len(SoundAmpF)
        SoundRec[Freq] = [[] for _ in range(len(SoundAmpF))]
        
        print('Applying amplification factors...')
        for AmpF in range(len(SoundAmpF)):
            SoundUnit[Freq][AmpF] = SoundPulseFiltered[Freq]
                                    
            Max = max(SoundUnit[Freq][AmpF])
            SoundUnit[Freq][AmpF] = [(SoundEl * (SBOutAmpF/Max))
                                     * SoundAmpF[AmpF]
                                     for SoundEl in SoundUnit[Freq][AmpF]]
            
            # Preallocating memory
            SoundList[Freq][AmpF] = [0]*(2*len(SoundUnit[Freq][AmpF]))
            Sound[Freq][AmpF] = [0]*(2*len(SoundUnit[Freq][AmpF]))
            
            for _ in range(len(SoundUnit[Freq][AmpF])):
                SoundList[Freq][AmpF][_ *2] = SoundUnit[Freq][AmpF][_]
                SoundList[Freq][AmpF][_ *2+1] = SoundTTLUnit[_]
            
            Sound[Freq][AmpF] = array.array('f')
            Sound[Freq][AmpF].fromlist(SoundList[Freq][AmpF])
            Sound[Freq][AmpF] = bytes(Sound[Freq][AmpF])
    
    print('Generating sound objects...')
    p = pyaudio.PyAudio()
    Stimulation = p.open(format=pyaudio.paFloat32,
                         channels=2,
                         rate=Rate,
                         #frames_per_buffer = len(Sound[0][0]),
                         input=False,
                         output=True)

    return(Sound, SoundRec, Stimulation)