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

This is a simple script to generate pulses and send to the soundboard and to
an Arduino board. Basically it generates sound pulses, sound square pulses, 
and laser square pulses. The square pulses will be sent to the left channel 
and the sound pulses will be sent to the right channel. In our current setup, 
we have an arduino uno set to read an analog pin and when the voltage is 
within a range, it turns on a digital port (Check 
ControlArduinoWithSoundBoard.ino code). This way, you have control of one 
sound channel, one laser and "separate" TTLs for both.
"""

#%% Set parameters of the experiment

import array
import pyaudio
import random
import scipy.signal
import threading

def GenSound(Rate, SoundPrePauseDur, SoundPulseDur, SoundPostPauseDur, \
             SoundPulseNo, SoundStimBlockNo, SoundPauseBetweenStimBlocksDur, \
             SoundAmpF, NoiseFrequency, TTLAmpF):
    
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


def GenSoundRec(Rate, SoundPrePauseDur, SoundPulseDur, SoundPostPauseDur, 
                SoundPulseNo, SoundStimBlockNo, SoundPauseBetweenStimBlocksDur, 
                SoundAmpF, NoiseFrequency, TTLAmpF):
    # This will generate the sound pulses and the sound output objects.
    # Remember to create input objects and call them accordingly.
    
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

