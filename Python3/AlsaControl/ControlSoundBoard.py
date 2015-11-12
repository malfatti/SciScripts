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
an Arduino board.

The first cell will set the experiment parameters. After, it creates an ALSA 
pcm object for sound playback and set some attributes of it.

All the following cells send the stimulus to the board, each one with its own 
settings. 

This code is profoundly inspired in a Matlab code, written by Dr. Richardson 
Leão, PhD., with the same purpose (to control a sound board and a laser). 
"""


#%% Set parameters of the experiment

import alsaaudio
import array
#import audioop
#import time
#import numpy
#import os
import random
import scipy.signal
#import subprocess
#import serial
import threading

import matplotlib.pyplot as plt

#from decimal import Decimal

#==========#==========#==========#==========#
#==========#==========#==========#==========#
Rate = 128000

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
SoundStimBlockNo = 5
# Duration of pause between blocks
SoundPauseBetweenStimBlocksDur = 5
NoiseFrequency = [5000, 15000]

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
LaserStimBlockNo = 5
# Duration of pause between blocks
LaserPauseBetweenStimBlocksDur = 5

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

SoundList = []
for _ in range(len(SoundUnit)):
    SoundList= SoundList + [SoundUnit[_]]
    SoundList= SoundList + [0]

Sound = array.array('f')
Sound.fromlist(SoundList)

SoundPauseBetweenStimBlocksList = [0] * round(SoundPauseBetweenStimBlocksDur * Rate) * 2
SoundPauseBetweenStimBlocks = array.array('f')
SoundPauseBetweenStimBlocks.fromlist(SoundPauseBetweenStimBlocksList)


## Generate laser stimulus
LaserPulse = [1] * round(LaserPulseDur * Rate)
LaserPulse[-1] = 0

LaserPrePause = [0] * round(LaserPrePauseDur * Rate)
LaserPostPause = [0] * round(LaserPostPauseDur * Rate)
LaserUnit = LaserPrePause + LaserPulse + LaserPostPause

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
SoundAndLaserList = []
for _ in range(len(SoundUnit)):
    SoundAndLaserList= SoundAndLaserList + [SoundUnit[_]]
    SoundAndLaserList= SoundAndLaserList + [LaserUnit[_]]

SoundAndLaser = array.array('f')
SoundAndLaser.fromlist(SoundAndLaserList)

SoundAndLaserPauseBetweenStimBlocksList = SoundPauseBetweenStimBlocksList + LaserPauseBetweenStimBlocksList
SoundAndLaserPauseBetweenStimBlocks = array.array('f')
SoundAndLaserPauseBetweenStimBlocks.fromlist(SoundAndLaserPauseBetweenStimBlocksList)


## Generate sound and laser objects and define threads
SoundStim = alsaaudio.PCM(alsaaudio.PCM_PLAYBACK)
SoundStim.setchannels(2)
SoundStim.setrate(Rate)
SoundStim.setformat(alsaaudio.PCM_FORMAT_S32_LE)
SoundStim.setperiodsize(len(SoundList))

class StartSound(threading.Thread):
    def run(self):
        for OneBlock in range(1, SoundStimBlockNo):
            for OnePulse in range(1, SoundPulseNo):
                SoundStim.write(Sound)
            
            SoundStim.write(SoundPauseBetweenStimBlocks)
            
        for OnePulse in range(1, SoundPulseNo):
            SoundStim.write(Sound)
        
    


class StartLaser(threading.Thread):
    def run(self):
        for OneBlock in range(1, LaserStimBlockNo):
            for OnePulse in range(1, LaserPulseNo):
                SoundStim.write(Laser)
            
            SoundStim.write(LaserPauseBetweenStimBlocks)
            
        for OnePulse in range(1, LaserPulseNo):
            SoundStim.write(Laser)
        
    


class StartSoundAndLaser(threading.Thread):
    def run(self):
        for OneBlock in range(1, SoundStimBlockNo):
            for OnePulse in range(1, SoundPulseNo):
                SoundStim.write(SoundAndLaser)
            
            SoundStim.write(SoundAndLaserPauseBetweenStimBlocks)
            
        for OnePulse in range(1, SoundPulseNo):
            SoundStim.write(SoundAndLaser)
        
    



#%% Run sound
StartSound().start()


#%% Run laser
StartLaser().start()


#%% Run sound and laser
StartSoundAndLaser().start()

