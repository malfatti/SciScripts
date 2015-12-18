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

This is a script to generate a filtered sound pulse and send it to a sound 
board.

The first cell will set the experiment parameters and create an ALSA pcm 
object for sound playback. The following cell send the stimulus to the sound 
board.

This code is profoundly inspired in a Matlab code, written by Dr. Richardson 
Leão, PhD.
"""
#%% Set up everything

import alsaaudio
import array
#import numpy
import random
import scipy.signal

#==========#==========#==========#==========#
#==========#==========#==========#==========#
Rate = 128000
SoundPulseDur = 60 * 70 # in SECONDS!
NoiseFrequency = [10000, 14000]
#==========#==========#==========#==========#
#==========#==========#==========#==========#

# Generate sound stimulus

SoundNoise = [random.random() for _ in range(0, round(Rate))]
SoundPulse = [SoundNoise[ElI]*2-1 for ElI,ElV in enumerate(SoundNoise)]

passband = [NoiseFrequency[i]/(Rate/2) for i,j in enumerate(NoiseFrequency)]
f2, f1 = scipy.signal.butter(4, passband, 'bandpass')
SoundPulseFiltered = scipy.signal.filtfilt(f2, f1, SoundPulse, padtype='odd', padlen=0)
SoundPulseFiltered = SoundPulseFiltered.tolist()
SoundPulseFiltered[-1] = 0

Sound = array.array('f')                                                        
Sound.fromlist(SoundPulseFiltered)


## Generate sound objects

SoundStim = alsaaudio.PCM(alsaaudio.PCM_PLAYBACK)
SoundStim.setchannels(1)
SoundStim.setrate(Rate)
SoundStim.setformat(alsaaudio.PCM_FORMAT_S32_LE)
SoundStim.setperiodsize(128)

# Test:
#SoundStim.write(Sound))
#import matplotlib.pyplot as plt
#plt.plot(Sound)

#%% Run
for _ in range(SoundPulseDur):
    SoundStim.write(Sound)

