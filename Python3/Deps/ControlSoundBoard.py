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
import math
import matplotlib.animation as animation
import matplotlib.pyplot as plt
import pyaudio
import shelve
import random
from scipy import signal
import threading

SoundTTLVal = 0.6; LaserTTLVal = 0.3

def GenSound(Rate, SoundPulseDur, SoundPulseNo, SoundAmpF, NoiseFrequency, 
             TTLAmpF, CalibrationFile, SoundPrePauseDur=0, SoundPostPauseDur=0, 
             SoundStimBlockNo=1, SoundPauseBetweenStimBlocksDur=0):
    """ Generate sound pulses in one channel and TTLs in the other channel 
    (Check ControlArduinoWithSoundBoard.ino code)."""
    
    
    with shelve.open(CalibrationFile) as Shelve:
        SBOutAmpF = Shelve['SBOutAmpF']

    print('Generating Sound TTL...')
    SoundTTLPrePause = [0] * round(SoundPrePauseDur * Rate)
    SoundTTLPostPause = [0] * round(SoundPostPauseDur * Rate)
    if SoundPulseDur < 0.01:
        SoundTTLPulse = [round(SoundTTLVal/SBOutAmpF, 3)] * round(SoundPulseDur * Rate)
    else:
        Middle = [0]*round((SoundPulseDur-0.01) * Rate)
        Border = [round(SoundTTLVal/SBOutAmpF, 3)] * round(0.005 * Rate)
        SoundTTLPulse = Border + Middle + Border
    
    SoundTTLPulse[-1] = 0
    
    SoundTTLUnit = SoundTTLPrePause + SoundTTLPulse + SoundTTLPostPause
    SoundTTLUnit = [SoundTTLEl*TTLAmpF for SoundTTLEl in SoundTTLUnit]
    
    print('Generating sound pulse...')
    SoundPrePause = [0] * round(Rate * SoundPrePauseDur)
    SoundPostPause = [0] * round(Rate * SoundPostPauseDur)

    # Preallocating memory
    SoundPulseFiltered = [0]*len(NoiseFrequency)
    SoundUnit = [0]*len(NoiseFrequency)
    SoundList = [0]*len(NoiseFrequency)
    Sound = [0]*len(NoiseFrequency)
    
    SoundNoise = [random.random() \
                  for _ in range(round(Rate*SoundPulseDur))]
    SoundPulse = [SoundNoise[ElI]*2-1 for ElI,ElV in enumerate(SoundNoise)]
    
    for Freq in range(len(NoiseFrequency)):
        Key= str(NoiseFrequency[Freq][0]) + '-' + str(NoiseFrequency[Freq][1])
        
        print('Filtering sound: ', NoiseFrequency[Freq], '...')
        passband = [NoiseFrequency[Freq][i]/(Rate/2) \
                    for i,j in enumerate(NoiseFrequency[Freq])]
        f2, f1 = signal.butter(4, passband, 'bandpass')
        SoundPulseFiltered[Freq] = signal.filtfilt(f2, f1, SoundPulse, \
                                                         padtype='odd', \
                                                         padlen=0)
        SoundPulseFiltered[Freq] = SoundPulseFiltered[Freq].tolist()
        SoundPulseFiltered[Freq][-1] = 0
        
        # Preallocating memory
        SoundUnit[Freq] = [0]*len(SoundAmpF[Key])
        SoundList[Freq] = [0]*len(SoundAmpF[Key])
        Sound[Freq] = [0]*len(SoundAmpF[Key])
        
        print('Applying amplification factors and interleaving channels...')
        for AmpF in range(len(SoundAmpF[Key])):
            SoundUnit[Freq][AmpF] = SoundPrePause + \
                                    SoundPulseFiltered[Freq] + \
                                    SoundPostPause
            Max = max(SoundUnit[Freq][AmpF])
            SoundUnit[Freq][AmpF] = [(SoundEl * (SBOutAmpF/Max))
                                     * SoundAmpF[Key][AmpF]
                                     for SoundEl in SoundUnit[Freq][AmpF]]
            
            # Preallocating memory
            SoundList[Freq][AmpF] = [0]*(2*len(SoundUnit[Freq][AmpF]))
            Sound[Freq][AmpF] = [0]*(2*len(SoundUnit[Freq][AmpF]))  
            
            for _ in range(len(SoundUnit[Freq][AmpF])):
                SoundList[Freq][AmpF][_ *2] = SoundUnit[Freq][AmpF][_]
                SoundList[Freq][AmpF][_ *2+1] = SoundTTLUnit[_]
            
            Sound[Freq][AmpF] = array.array('f', SoundList[Freq][AmpF])
            Sound[Freq][AmpF] = bytes(Sound[Freq][AmpF])
    
    print('Generating pause...')
    SoundPauseBetweenStimBlocksList = \
                         [0] * round(SoundPauseBetweenStimBlocksDur * Rate) * 2
    SoundPauseBetweenStimBlocks = array.array('f')
    SoundPauseBetweenStimBlocks.fromlist(SoundPauseBetweenStimBlocksList)
    SoundPauseBetweenStimBlocks = bytes(SoundPauseBetweenStimBlocks)
    
    print('Generating sound objects...')
    p = pyaudio.PyAudio()
    Stimulation = p.open(format=pyaudio.paFloat32,
                    channels=2,
                    rate=Rate,
                    output=True)
    
    class StartSound(threading.Thread):
        def run(self):
            for Freq in range(len(NoiseFrequency)):
                Key= str(NoiseFrequency[Freq][0]) + '-' \
                     + str(NoiseFrequency[Freq][1])
                
                for AmpF in range(len(SoundAmpF[Key])):
                    for OneBlock in range(SoundStimBlockNo):
                        for OnePulse in range(SoundPulseNo):
                            Stimulation.write(Sound[Freq][AmpF])
                
                        Stimulation.write(SoundPauseBetweenStimBlocks)
    
    
    print('Done generating sound stimulus.')
    return(Sound, SoundPauseBetweenStimBlocks, StartSound)


def GenLaser(Rate, LaserPulseDur, LaserPulseNo, TTLAmpF, CalibrationFile, 
             LaserPrePauseDur=0, LaserPostPauseDur=0, LaserStimBlockNo=1, 
             LaserPauseBetweenStimBlocksDur=0):
    """ Generate square pulses in one channel that works as TTL for laser 
    (Check ControlArduinoWithSoundBoard.ino code)."""
    
    with shelve.open(CalibrationFile) as Shelve:
        SBOutAmpF = Shelve['SBOutAmpF']
    
    print('Generating laser pulse...')
    LaserPulse = [round(LaserTTLVal/SBOutAmpF, 3)] * round(LaserPulseDur * Rate)
    LaserPulse[-1] = 0
    
    LaserPrePause = [0] * round(LaserPrePauseDur * Rate)
    LaserPostPause = [0] * round(LaserPostPauseDur * Rate)
    LaserUnit = LaserPrePause + LaserPulse + LaserPostPause
    LaserUnit = [LaserEl*TTLAmpF for LaserEl in LaserUnit]
    
    print('Interleaving channels...')
    LaserList = [0]*(2*len(LaserUnit))
    for _ in range(len(LaserUnit)):
        LaserList[_ *2] = 0
        LaserList[_ *2+1] = LaserUnit[_]
    
    Laser = array.array('f')
    Laser.fromlist(LaserList)
    Laser = bytes(Laser)
    
    print('Generating pause...')
    LaserPauseBetweenStimBlocksList = \
                        [0] * round(LaserPauseBetweenStimBlocksDur * Rate)  * 2
    LaserPauseBetweenStimBlocks = array.array('f')
    LaserPauseBetweenStimBlocks.fromlist(LaserPauseBetweenStimBlocksList)
    LaserPauseBetweenStimBlocks = bytes(LaserPauseBetweenStimBlocks)
    
    print('Generating sound objects...')
    p = pyaudio.PyAudio()
    Stimulation = p.open(format=pyaudio.paFloat32,
                    channels=2,
                    rate=Rate,
                    output=True)
    
    class StartLaser(threading.Thread):
        def run(self):
            for OneBlock in range(LaserStimBlockNo):
                for OnePulse in range(LaserPulseNo):
                    Stimulation.write(Laser)
                
                Stimulation.write(LaserPauseBetweenStimBlocks)
    
    
    print('Done generating laser stimulus.')
    return(Laser, LaserPauseBetweenStimBlocks, StartLaser)


def GenSoundLaser(Rate, SoundPulseDur, SoundPulseNo, SoundAmpF, NoiseFrequency, 
                  LaserPulseDur, LaserPulseNo, TTLAmpF, CalibrationFile, 
                  SoundPrePauseDur=0, SoundPostPauseDur=0, SoundStimBlockNo=1, 
                  SoundPauseBetweenStimBlocksDur=0, LaserPrePauseDur=0, 
                  LaserPostPauseDur=0, LaserStimBlockNo=1, 
                  LaserPauseBetweenStimBlocksDur=0):
    """ Generate sound pulses in one channel and TTLs for sound and laser in 
    the other channel (Check ControlArduinoWithSoundBoard.ino code)."""
    
    with shelve.open(CalibrationFile) as Shelve:
        SBOutAmpF = Shelve['SBOutAmpF']
    
    print('Generating laser pulse...')
    LaserPulse = [round(LaserTTLVal/SBOutAmpF, 3)] * round(LaserPulseDur * Rate)
    LaserPulse[-1] = 0
    
    LaserPrePause = [0] * round(LaserPrePauseDur * Rate)
    LaserPostPause = [0] * round(LaserPostPauseDur * Rate)
    LaserUnit = LaserPrePause + LaserPulse + LaserPostPause
    LaserUnit = [LaserEl*TTLAmpF for LaserEl in LaserUnit]
    
    print('Generating Sound TTL...')
    SoundTTLPrePause = [0] * round(SoundPrePauseDur * Rate)
    SoundTTLPostPause = [0] * round(SoundPostPauseDur * Rate)
    SoundTTLPulse = [round(SoundTTLVal/SBOutAmpF, 3)] * round(SoundPulseDur * Rate)
    SoundTTLPulse[-1] = 0
    
    SoundTTLUnit = SoundTTLPrePause + SoundTTLPulse + SoundTTLPostPause
    SoundTTLUnit = [SoundTTLEl*TTLAmpF for SoundTTLEl in SoundTTLUnit]
    
    print('Summing sound TTL and laser units...')
    SoundTTLAndLaserUnit = [LaserUnit[_]+SoundTTLUnit[_] \
                            for _ in range(len(SoundTTLUnit))]
    
    print('Generating sound pulse...')
    SoundPrePause = [0] * round(Rate * SoundPrePauseDur)
    SoundPostPause = [0] * round(Rate * SoundPostPauseDur)
    SoundNoise = [random.random() for _ in range(0, round(Rate*SoundPulseDur))]
    SoundPulse = [SoundNoise[ElI]*2-1 for ElI,ElV in enumerate(SoundNoise)]
    
    # Preallocating memory
    SoundPulseFiltered = [0]*len(NoiseFrequency)
    SoundUnit = [0]*len(NoiseFrequency)
    SoundAndLaserList = [0]*len(NoiseFrequency)
    SoundAndLaser = [0]*len(NoiseFrequency)
    
    for Freq in range(len(NoiseFrequency)):
        Key= str(NoiseFrequency[Freq][0]) + '-' + str(NoiseFrequency[Freq][1])
        
        print('Filtering sound: ', NoiseFrequency[Freq], '...')
        passband = [NoiseFrequency[Freq][i]/(Rate/2) \
                    for i,j in enumerate(NoiseFrequency[Freq])]
        f2, f1 = signal.butter(4, passband, 'bandpass')
        SoundPulseFiltered[Freq] = signal.filtfilt(f2, f1, SoundPulse, \
                                                         padtype='odd', \
                                                         padlen=0)
        SoundPulseFiltered[Freq] = SoundPulseFiltered[Freq].tolist()
        SoundPulseFiltered[Freq][-1] = 0

        # Preallocating memory
        SoundUnit[Freq] = [0]*len(SoundAmpF[Key])
        SoundAndLaserList[Freq] = [0]*len(SoundAmpF[Key])
        SoundAndLaser[Freq] = [0]*len(SoundAmpF[Key])
    
        for AmpF in range(len(SoundAmpF[Key])):
            print('Applying amplification factor:', SoundAmpF[Key][AmpF], '...')
            SoundUnit[Freq][AmpF] = SoundPrePause + \
                                    SoundPulseFiltered[Freq] + \
                                    SoundPostPause
            Max = max(SoundUnit[Freq][AmpF])
            SoundUnit[Freq][AmpF] = [(SoundEl * (SBOutAmpF/Max))
                                     * SoundAmpF[Key][AmpF]
                                     for SoundEl in SoundUnit[Freq][AmpF]]
            
            # Preallocating memory
            SoundAndLaserList[Freq][AmpF] = [0]*(2*len(SoundUnit[Freq][AmpF]))
            SoundAndLaser[Freq][AmpF] = [0]*(2*len(SoundUnit[Freq][AmpF]))
            
            print('Interleaving channels...')
            for _ in range(len(SoundUnit[Freq][AmpF])):
                SoundAndLaserList[Freq][AmpF][_ *2] = SoundUnit[Freq][AmpF][_]
                SoundAndLaserList[Freq][AmpF][_ *2+1] = SoundTTLAndLaserUnit[_]
            
            SoundAndLaser[Freq][AmpF] = array.array('f')
            SoundAndLaser[Freq][AmpF].fromlist(SoundAndLaserList[Freq][AmpF])
            SoundAndLaser[Freq][AmpF] = bytes(SoundAndLaser[Freq][AmpF])
    
    print('Generating pause...')
    SoundPauseBetweenStimBlocksList = \
                         [0] * round(SoundPauseBetweenStimBlocksDur * Rate) * 2
    SoundAndLaserPauseBetweenStimBlocks = array.array('f')
    SoundAndLaserPauseBetweenStimBlocks.fromlist(SoundPauseBetweenStimBlocksList)
    SoundAndLaserPauseBetweenStimBlocks = bytes(SoundAndLaserPauseBetweenStimBlocks)
    
    print('Generating sound objects...')
    p = pyaudio.PyAudio()
    Stimulation = p.open(format=pyaudio.paFloat32,
                    channels=2,
                    rate=Rate,
                    output=True)
    
    class StartSoundAndLaser(threading.Thread):
        def run(self):
            for Freq in range(len(NoiseFrequency)):
                Key= str(NoiseFrequency[Freq][0]) + '-' \
                     + str(NoiseFrequency[Freq][1])
                
                for AmpF in range(len(SoundAmpF[Key])):
                    for OneBlock in range(SoundStimBlockNo):
                        for OnePulse in range(SoundPulseNo):
                            Stimulation.write(SoundAndLaser[Freq][AmpF])
                
                        Stimulation.write(SoundAndLaserPauseBetweenStimBlocks)
    
    
    print('Done generating sound and laser stimulus.')
    return(SoundAndLaser, SoundAndLaserPauseBetweenStimBlocks, StartSoundAndLaser)

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


def MicrOscilloscope(Rate, XLim, YLim, CalibrationFile, FramesPerBuf=512):
    """ Read data from sound board input and plot it until the windows is 
    closed. """
    
    with shelve.open(CalibrationFile) as Shelve:
        SBInAmpF = Shelve['SBInAmpF']
    
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


def MicrOscilloscopeRec(Rate, XLim, YLim, CalibrationFile, FramesPerBuf=512):
    """ Read data from sound board input, record it to a video and plot it 
    until the windows is closed (with a delay). """
    
    with shelve.open(CalibrationFile) as Shelve:
        SBInAmpF = Shelve['SBInAmpF']
    
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


def SoundCalOut(Rate):
    """ Generate 3s of 100Hz sine wave (1 to -1). """
    
    Freq = 10000; Time = 0.1
    
    print('Generating sound...')
    Pulse = [math.sin(2*math.pi*Freq*(_/Rate)) for _ in range(round(Rate*Time))]
    Pulse[-1] = 0
    
    print('Interleaving channels...')
    List = [0]*(2*len(Pulse))
    for _ in range(len(Pulse)):
        List[_ *2] = Pulse[_]
        List[_ *2+1] = 0
    
    Pulse = array.array('f')
    Pulse.fromlist(List)
    Pulse = bytes(Pulse)
    
    print('Generating sound objects...')
    p = pyaudio.PyAudio()
    Stimulation = p.open(format=pyaudio.paFloat32,
                    channels=2,
                    rate=Rate,
                    output=True)
    
    print('Playing...')
    for OnePulse in range(round(300/Time)):
        Stimulation.write(Pulse)


def SoundCalIn(Rate, SBOutAmpF):
    """ Generate 3s of 10KHz sine wave (1V to -1V) and read 1s of it. """
    Freq = 10000; Time = 0.1
    
    print('Generating laser pulse...')
    Pulse = [math.sin(2*math.pi*Freq*(_/Rate)) for _ in range(round(Rate*Time))]
    Pulse[-1] = 0
    
    print('Interleaving channels...')
    List = [0]*(2*len(Pulse))
    for _ in range(len(Pulse)):
        List[_ *2] = Pulse[_]
        List[_ *2+1] = 0
    
    Pulse = array.array('f')
    Pulse.fromlist(List)
    Pulse = bytes(Pulse)
    
    print('Generating sound objects...')
    p = pyaudio.PyAudio()
    Stimulation = p.open(format=pyaudio.paFloat32,
                    channels=2,
                    rate=Rate,
                    output=True)
    
    q = pyaudio.PyAudio()
    Reading = q.open(format=pyaudio.paFloat32,
                    channels=2,
                    rate=Rate,
                    input=True)
    
    class PlayThr(threading.Thread):
        def run(self):
            for OnePulse in range(round(30/Time)):
                Stimulation.write(Pulse)
    
    print('Playing...')
    PlayThr().start()
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