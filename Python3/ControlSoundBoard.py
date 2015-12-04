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
an Arduino board. Basically it generates a sound pulse, a sound square pulse, 
and a laser square pulse. The square pulses will be sent to the left channel 
and the sound pulse will be sent to the right channel. In our current setup, 
we have an arduino uno set to read an analog pin and when the voltage is 
within a range, it turns on a digital port (Check 
ControlArduinoWithSoundBoard.ino code). This way, you have control of one 
sound channel, one laser and TTLs for both.

The first cell will set the experiment parameters. After, it creates an audio 
object for sound playback and set some attributes of it.

All the following cells send the stimulus to the board, each one with its own 
settings. 

This code is inspired in a Matlab code, written by Dr. Richardson Leão, PhD., 
with the same purpose (to control a sound board and a laser).
"""

#%% Set parameters of the experiment

import array
import pyaudio
import random
import scipy.signal
import threading

import matplotlib.pyplot as plt

#==========#==========#==========#==========#
#==========#==========#==========#==========#
Rate = 192000

## Fill all durations in SECONDS!

## Sound
# Silence before pulse
SoundPrePauseDur = 0.004
# Pulse duration
SoundPulseDur = 0.003
# Silence after pulse
SoundPostPauseDur = 0.093
# Amount of pulses per block
SoundPulseNo = 200
# Number of blocks
SoundStimBlockNo = 3
# Duration of pause between blocks
SoundPauseBetweenStimBlocksDur = 5
# Amplification factor. Max = 1
SoundAmpF = 0.1
NoiseFrequency = [14000, 16000]

## Laser
# Silence before pulse
LaserPrePauseDur = 0
# Pulse duration
LaserPulseDur = 0.01
# Silence after pulse
LaserPostPauseDur = 0.09
# Amount of pulses per block
LaserPulseNo = 200
# Number of blocks
LaserStimBlockNo = 3
# Duration of pause between blocks
LaserPauseBetweenStimBlocksDur = 5
# Amplification factor. Max = 1
LaserAmpF = 1

#==========#==========#==========#==========#
#==========#==========#==========#==========#

## Generate sound stimulus
SoundNoise = [random.random() for _ in range(0, round(Rate*SoundPulseDur))]
SoundPulse = [SoundNoise[ElI]*2-1 for ElI,ElV in enumerate(SoundNoise)]

passband = [NoiseFrequency[i]/(Rate/2) for i,j in enumerate(NoiseFrequency)]
f2, f1 = scipy.signal.butter(4, passband, 'bandpass')
SoundPulseFiltered = scipy.signal.filtfilt(f2, f1, SoundPulse, padtype='odd', padlen=0)
SoundPulseFiltered = SoundPulseFiltered.tolist()
SoundPulseFiltered[-1] = 0

SoundPrePause = [0] * round(Rate * SoundPrePauseDur)
SoundPostPause = [0] * round(Rate * SoundPostPauseDur)
SoundUnit = SoundPrePause + SoundPulseFiltered + SoundPostPause
SoundUnit = [SoundEl*SoundAmpF for SoundEl in SoundUnit]

# Generate Sound TTL
SoundTTLPulse = [0.6] * round(SoundPulseDur * Rate)
SoundTTLPulse[-1] = 0

SoundTTLPrePause = [0] * round(SoundPrePauseDur * Rate)
SoundTTLPostPause = [0] * round(SoundPostPauseDur * Rate)
SoundTTLUnit = SoundTTLPrePause + SoundTTLPulse + SoundTTLPostPause
SoundTTLUnit = [SoundTTLEl*LaserAmpF for SoundTTLEl in SoundTTLUnit]

SoundList = []
for _ in range(len(SoundUnit)):
    SoundList= SoundList + [SoundUnit[_]]
    SoundList= SoundList + [SoundTTLUnit[_]]


Sound = array.array('f')
Sound.fromlist(SoundList)

SoundPauseBetweenStimBlocksList = [0] * round(SoundPauseBetweenStimBlocksDur * Rate) * 2
SoundPauseBetweenStimBlocks = array.array('f')
SoundPauseBetweenStimBlocks.fromlist(SoundPauseBetweenStimBlocksList)


## Generate laser stimulus (No TTL needed, since it is a square already)
LaserPulse = [0.2] * round(LaserPulseDur * Rate)
LaserPulse[-1] = 0

LaserPrePause = [0] * round(LaserPrePauseDur * Rate)
LaserPostPause = [0] * round(LaserPostPauseDur * Rate)
LaserUnit = LaserPrePause + LaserPulse + LaserPostPause
LaserUnit = [LaserEl*LaserAmpF for LaserEl in LaserUnit]

LaserList = []
for _ in range(len(LaserUnit)):
    LaserList= LaserList + [0]
    LaserList= LaserList + [LaserUnit[_]]

Laser = array.array('f')
Laser.fromlist(LaserList)

LaserPauseBetweenStimBlocksList = [0] * round(LaserPauseBetweenStimBlocksDur * Rate)  * 2
LaserPauseBetweenStimBlocks = array.array('f')
LaserPauseBetweenStimBlocks.fromlist(LaserPauseBetweenStimBlocksList)


## Generate both stimulus together
SoundTTLAndLaserUnit = [LaserUnit[_]+SoundTTLUnit[_] for _ in range(len(SoundTTLUnit))]

SoundAndLaserList = []
for _ in range(len(SoundUnit)):
    SoundAndLaserList= SoundAndLaserList + [SoundUnit[_]]
    SoundAndLaserList= SoundAndLaserList + [SoundTTLAndLaserUnit[_]]

SoundAndLaser = array.array('f')
SoundAndLaser.fromlist(SoundAndLaserList)

SoundAndLaserPauseBetweenStimBlocksList = SoundPauseBetweenStimBlocksList + LaserPauseBetweenStimBlocksList
SoundAndLaserPauseBetweenStimBlocks = array.array('f')
SoundAndLaserPauseBetweenStimBlocks.fromlist(SoundAndLaserPauseBetweenStimBlocksList)


## Generate sound objects and define threads
p = pyaudio.PyAudio()
Stimulation = p.open(format=pyaudio.paFloat32,
                channels=2,
                rate=Rate,
                output=True)

class StartSound(threading.Thread):
    def run(self):
        for OneBlock in range(1, SoundStimBlockNo):
            for OnePulse in range(1, SoundPulseNo):
                Stimulation.write(bytes(Sound))
            
            Stimulation.write(bytes(SoundPauseBetweenStimBlocks))
            
            for OnePulse in range(1, SoundPulseNo):
                Stimulation.write(bytes(Sound))
            
        


class StartLaser(threading.Thread):
    def run(self):
        for OneBlock in range(1, LaserStimBlockNo):
            for OnePulse in range(1, LaserPulseNo):
                Stimulation.write(bytes(Laser))
            
            Stimulation.write(bytes(LaserPauseBetweenStimBlocks))
            
        for OnePulse in range(1, LaserPulseNo):
            Stimulation.write(bytes(Laser))
        
    


class StartSoundAndLaser(threading.Thread):
    def run(self):
        for OneBlock in range(1, SoundStimBlockNo):
            for OnePulse in range(1, SoundPulseNo):
                Stimulation.write(bytes(SoundAndLaser))
            
            Stimulation.write(bytes(SoundAndLaserPauseBetweenStimBlocks))
            
        for OnePulse in range(1, SoundPulseNo):
            Stimulation.write(bytes(SoundAndLaser))
        
    



#%% Run sound
StartSound().start()

    
#%% Run laser
StartLaser().start()


#%% Run sound and laser
StartSoundAndLaser().start()


#%%

import pyaudio
import array
import random

import matplotlib.pyplot as plt

Rate = 192000
SoundPulseDur = 4
SoundAmpF = 10
a = [random.uniform(-1,1) for _ in range(round(Rate*SoundPulseDur))]
b = [SoundEl*SoundAmpF for SoundEl in a]
#b = [0.0]*round(Rate*SoundPulseDur)

TestB = array.array('f')
TestB.fromlist(b)

p = pyaudio.PyAudio()
Stimulation = p.open(format=pyaudio.paFloat32,
                channels=2,
                rate=Rate,
                output=True)



#%%
Stimulation.write(bytes(TestB))

Stimulation.close()
p.terminate()

  if (SoundAndLaserTTLInV > 190 && SoundAndLaserTTLInV < 300) {
    digitalWrite(LaserOut, HIGH);
  } else {
    digitalWrite(LaserOut, LOW);
  }

  if (SoundAndLaserTTLInV > 320 && SoundAndLaserTTLInV < 800) {
    digitalWrite(SoundTTLOut, HIGH);
  } else {
    digitalWrite(SoundTTLOut, LOW);
  }

  if (SoundAndLaserTTLInV > 810) {
    digitalWrite(SoundTTLOut, HIGH);
    digitalWrite(LaserOut, HIGH);
  } else {
    digitalWrite(SoundTTLOut, LOW);
    digitalWrite(LaserOut, LOW);
  }
  

val = map(val, 0, 1023, 0, 63);  
#++++++++++++++++++++++++++++++++++++
  if (SoundAndLaserTTLInV < 50) {
    digitalWrite(LaserOut, LOW);
    digitalWrite(SoundTTLOut, LOW);
  }

  if (SoundAndLaserTTLInV > 60 && SoundAndLaserTTLInV < 150) {
    digitalWrite(LaserOut, HIGH);
  } else {
    digitalWrite(LaserOut, LOW);
  }

  if (SoundAndLaserTTLInV > 200 && SoundAndLaserTTLInV < 900) {
    digitalWrite(SoundTTLOut, HIGH);
  } else {
    digitalWrite(SoundTTLOut, LOW);
  }