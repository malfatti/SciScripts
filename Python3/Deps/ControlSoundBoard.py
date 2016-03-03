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

This is a script that define functions allowing the use of a computer's sound
board as an analog I/O board.

"""

import array
import pyaudio
import random
import scipy.signal
import threading

def GenSound(Rate, SoundPulseDur, SoundPulseNo, SoundAmpF, NoiseFrequency, 
             TTLAmpF, SoundPrePauseDur=0, SoundPostPauseDur=0, 
             SoundStimBlockNo=1, SoundPauseBetweenStimBlocksDur=0):
    """ Generate sound pulses in one channel and TTLs in the other channel 
    (Check ControlArduinoWithSoundBoard.ino code)."""
    
    print('Generating Sound TTL...')
    SoundTTLPrePause = [0] * round(SoundPrePauseDur * Rate)
    SoundTTLPostPause = [0] * round(SoundPostPauseDur * Rate)
    if SoundPulseDur < 0.01:
        SoundTTLPulse = [0.6] * round(SoundPulseDur * Rate)
    else:
        Middle = [0]*round((SoundPulseDur-0.01) * Rate)
        Border = [0.6] * round(0.005 * Rate)
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
        print('Filtering sound: ', NoiseFrequency[Freq], '...')
        passband = [NoiseFrequency[Freq][i]/(Rate/2) \
                    for i,j in enumerate(NoiseFrequency[Freq])]
        f2, f1 = scipy.signal.butter(4, passband, 'bandpass')
        SoundPulseFiltered[Freq] = scipy.signal.filtfilt(f2, f1, SoundPulse, \
                                                         padtype='odd', \
                                                         padlen=0)
        SoundPulseFiltered[Freq] = SoundPulseFiltered[Freq].tolist()
        SoundPulseFiltered[Freq][-1] = 0
        
        # Preallocating memory
        SoundUnit[Freq] = [0]*len(SoundAmpF)
        SoundList[Freq] = [0]*len(SoundAmpF)
        Sound[Freq] = [0]*len(SoundAmpF)
        
        for AmpF in range(len(SoundAmpF)):
            print('Applying amplification factor:', SoundAmpF[AmpF], '...')
            SoundUnit[Freq][AmpF] = SoundPrePause + \
                                    SoundPulseFiltered[Freq] + \
                                    SoundPostPause
            SoundUnit[Freq][AmpF] = [SoundEl*SoundAmpF[AmpF] \
                                     for SoundEl in SoundUnit[Freq][AmpF]]
            
            # Preallocating memory
            SoundList[Freq][AmpF] = [0]*(2*len(SoundUnit[Freq][AmpF]))
            Sound[Freq][AmpF] = [0]*(2*len(SoundUnit[Freq][AmpF]))  
            
            print('Interleaving channels...')
            for _ in range(len(SoundUnit[Freq][AmpF])):
                SoundList[Freq][AmpF][_ *2] = SoundUnit[Freq][AmpF][_]
                SoundList[Freq][AmpF][_ *2+1] = SoundTTLUnit[_]
            
            Sound[Freq][AmpF] = array.array('f')
            Sound[Freq][AmpF].fromlist(SoundList[Freq][AmpF])
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
                for AmpF in range(len(SoundAmpF)):
                    for OneBlock in range(SoundStimBlockNo):
                        for OnePulse in range(SoundPulseNo):
                            Stimulation.write(Sound[Freq][AmpF])
                
                        Stimulation.write(SoundPauseBetweenStimBlocks)
    
    
    print('Done generating sound stimulus.')
    return(Sound, SoundPauseBetweenStimBlocks, StartSound)


def GenLaser(Rate, LaserPrePauseDur, LaserPulseDur, LaserPostPauseDur, \
             LaserPulseNo, LaserStimBlockNo, LaserPauseBetweenStimBlocksDur, \
             TTLAmpF):
    """ Generate square pulses in one channel that works as TTL for laser 
    (Check ControlArduinoWithSoundBoard.ino code)."""
           
    print('Generating laser pulse...')
    LaserPulse = [0.2] * round(LaserPulseDur * Rate)
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


def GenSoundLaser(Rate, SoundPrePauseDur, SoundPulseDur, SoundPostPauseDur, \
                  SoundPulseNo, SoundStimBlockNo, \
                  SoundPauseBetweenStimBlocksDur, SoundAmpF, NoiseFrequency, \
                  LaserPrePauseDur, LaserPulseDur, LaserPostPauseDur, \
                  LaserPauseBetweenStimBlocksDur, TTLAmpF):
    """ Generate sound pulses in one channel and TTLs for sound and laser in 
    the other channel (Check ControlArduinoWithSoundBoard.ino code)."""
           
    print('Generating laser pulse...')
    LaserPulse = [0.2] * round(LaserPulseDur * Rate)
    LaserPulse[-1] = 0
    
    LaserPrePause = [0] * round(LaserPrePauseDur * Rate)
    LaserPostPause = [0] * round(LaserPostPauseDur * Rate)
    LaserUnit = LaserPrePause + LaserPulse + LaserPostPause
    LaserUnit = [LaserEl*TTLAmpF for LaserEl in LaserUnit]
    
    print('Generating Sound TTL...')
    SoundTTLPrePause = [0] * round(SoundPrePauseDur * Rate)
    SoundTTLPostPause = [0] * round(SoundPostPauseDur * Rate)
    SoundTTLPulse = [0.6] * round(SoundPulseDur * Rate)
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
        print('Filtering sound: ', NoiseFrequency[Freq], '...')
        passband = [NoiseFrequency[Freq][i]/(Rate/2) \
                    for i,j in enumerate(NoiseFrequency[Freq])]
        f2, f1 = scipy.signal.butter(4, passband, 'bandpass')
        SoundPulseFiltered[Freq] = scipy.signal.filtfilt(f2, f1, SoundPulse, \
                                                         padtype='odd', \
                                                         padlen=0)
        SoundPulseFiltered[Freq] = SoundPulseFiltered[Freq].tolist()
        SoundPulseFiltered[Freq][-1] = 0

        # Preallocating memory
        SoundUnit[Freq] = [0]*len(SoundAmpF)
        SoundAndLaserList[Freq] = [0]*len(SoundAmpF)
        SoundAndLaser[Freq] = [0]*len(SoundAmpF)
    
        for AmpF in range(len(SoundAmpF)):
            print('Applying amplification factor:', SoundAmpF[AmpF], '...')
            SoundUnit[Freq][AmpF] = SoundPrePause + \
                                    SoundPulseFiltered[Freq] + \
                                    SoundPostPause
            SoundUnit[Freq][AmpF] = [SoundEl*SoundAmpF[AmpF] \
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
                for AmpF in range(len(SoundAmpF)):
                    for OneBlock in range(SoundStimBlockNo):
                        for OnePulse in range(SoundPulseNo):
                            Stimulation.write(SoundAndLaser[Freq][AmpF])
                
                        Stimulation.write(SoundAndLaserPauseBetweenStimBlocks)
    
    
    print('Done generating sound and laser stimulus.')
    return(SoundAndLaser, SoundAndLaserPauseBetweenStimBlocks, StartSoundAndLaser)


def GenSoundRec(Rate, SoundPulseDur, SoundPulseNo, SoundAmpF, NoiseFrequency, 
                TTLAmpF, SoundPrePauseDur=0, SoundPostPauseDur=0, 
                SoundStimBlockNo=1, SoundPauseBetweenStimBlocksDur=0):
    """ This will generate the sound pulses and the sound output objects.
    Remember to create input objects and call them accordingly. """
    
    print('Generating Sound TTL...')
    SoundTTLPrePause = [0] * round(SoundPrePauseDur * Rate)
    SoundTTLPostPause = [0] * round(SoundPostPauseDur * Rate)
    SoundTTLPulse = [0.6] * round(SoundPulseDur * Rate)
    SoundTTLPulse[-1] = 0
    
    SoundTTLUnit = SoundTTLPrePause + SoundTTLPulse + SoundTTLPostPause
    SoundTTLUnit = [SoundTTLEl*TTLAmpF for SoundTTLEl in SoundTTLUnit]
    
    print('Generating sound pulse...')
    SoundPrePause = [0] * round(Rate * SoundPrePauseDur)
    SoundPostPause = [0] * round(Rate * SoundPostPauseDur)        
    
    SoundNoise = [random.random() 
                  for _ in range(round(Rate*SoundPulseDur))]
    SoundPulse = [SoundNoise[ElI]*2-1 for ElI,ElV in enumerate(SoundNoise)]
    
    # Preallocating memory
    SoundPulseFiltered = [0]*len(NoiseFrequency)                                    
    SoundUnit = [0]*len(NoiseFrequency)                                             
    SoundList = [0]*len(NoiseFrequency)                                             
    Sound = [0]*len(NoiseFrequency)                                                 
    SoundRec = [0]*len(NoiseFrequency) 
    
    for Freq in range(len(NoiseFrequency)):
        print('Filtering sound: ', NoiseFrequency[Freq], '...')
        passband = [NoiseFrequency[Freq][i]/(Rate/2) 
                    for i,j in enumerate(NoiseFrequency[Freq])]
        f2, f1 = scipy.signal.butter(4, passband, 'bandpass')
        SoundPulseFiltered[Freq] = scipy.signal.filtfilt(f2, f1, SoundPulse, 
                                                         padtype='odd', 
                                                         padlen=0)
        SoundPulseFiltered[Freq] = SoundPulseFiltered[Freq].tolist()
        SoundPulseFiltered[Freq][-1] = 0
        
        # Preallocating memory
        SoundUnit[Freq] = [0]*len(SoundAmpF)
        SoundList[Freq] = [0]*len(SoundAmpF)
        Sound[Freq] = [0]*len(SoundAmpF)
        SoundRec[Freq] = [[] for _ in range(len(SoundAmpF))]
        
        for AmpF in range(len(SoundAmpF)):
            print('Applying amplification factor:', SoundAmpF[AmpF], '...')
            SoundUnit[Freq][AmpF] = (SoundPrePause + 
                                    SoundPulseFiltered[Freq] + 
                                    SoundPostPause)
                                    
            SoundUnit[Freq][AmpF] = [SoundEl*SoundAmpF[AmpF] 
                                     for SoundEl in SoundUnit[Freq][AmpF]]
            
            # Preallocating memory
            SoundList[Freq][AmpF] = [0]*(2*len(SoundUnit[Freq][AmpF]))
            Sound[Freq][AmpF] = [0]*(2*len(SoundUnit[Freq][AmpF]))
            
            print('Interleaving channels...')
            for _ in range(len(SoundUnit[Freq][AmpF])):
                SoundList[Freq][AmpF][_ *2] = SoundUnit[Freq][AmpF][_]
                SoundList[Freq][AmpF][_ *2+1] = SoundTTLUnit[_]
            
            Sound[Freq][AmpF] = array.array('f')
            Sound[Freq][AmpF].fromlist(SoundList[Freq][AmpF])
            Sound[Freq][AmpF] = bytes(Sound[Freq][AmpF])
    
    print('Generating pause...')
    SoundPauseBetweenStimBlocksList = \
                         [0] * round(SoundPauseBetweenStimBlocksDur * Rate) * 2
    SoundPauseBetweenStimBlocks = array.array('f')
    SoundPauseBetweenStimBlocks.fromlist(SoundPauseBetweenStimBlocksList)
    
    print('Generating sound objects...')
    p = pyaudio.PyAudio()
    Stimulation = p.open(format=pyaudio.paFloat32,
                         channels=2,
                         rate=Rate,
                         #frames_per_buffer = len(Sound[0][0]),
                         input=False,
                         output=True)

    return(Sound, SoundPauseBetweenStimBlocks, SoundRec, Stimulation)


def ReduceStim(SObj):
    """This function will reduce a stimulus that have frequencies at several
       AmpFs to a stimulus that has one AmpF per frequency, following the same
       index for NoiseFrequency and SoundAmpF. For example:
       
           SoundAmpF = [1, 0.5]
           NoiseFrequency = [[8000, 10000],[12000, 14000]]
       
           # create sound obj...
       
           Sound = ReduceStim(SObj)
       
       Sound will have NoiseFrequency[0] at SoundAmpF[0], NoiseFrequency[1] at 
       SoundAmpF[1], etc. So, obviously, len(SoundAmpF) must be equal to 
       len(NoiseFrequency)."""
    
    for _ in range(len(SObj)):
        SObj[_] = SObj[_][_]
    
    return(SObj)


def MicrOscilloscope(Rate, XLim, YLim, FramesPerBuf=512):
    """ Read data from sound board input and plot it until the windows is 
    closed. """
    
    import array
    import matplotlib.animation as animation
    import matplotlib.pyplot as plt
    import pyaudio

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
        Plot.set_ydata(Data)
        return Plot,
    
    Anim = animation.FuncAnimation(Fig, PltUp, frames=FramesPerBuf, interval=16, blit=False)


def MicrOscilloscopeRec(Rate, XLim, YLim, FramesPerBuf=512):
    """ Read data from sound board input, record it to a video and plot it 
    until the windows is closed (with a delay). """
    
    import array
    import matplotlib.animation as animation
    import matplotlib.pyplot as plt
    import pyaudio

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
        Plot.set_ydata(Data)
        return Plot,
    
    Anim = animation.FuncAnimation(Fig, PltUp, frames=FramesPerBuf, interval=16, blit=False)
    Anim.save('MicrOscilloscope.mp4', writer=Writer)