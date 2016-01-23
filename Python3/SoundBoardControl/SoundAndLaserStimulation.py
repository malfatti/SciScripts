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

This is a script to generate pulses and send to the soundboard, and then to a 
sound amplifier and an Arduino board. Basically it generates sound pulses, 
sound square pulses (TTLs), and laser square pulses. The square pulses will be 
sent to the left channel and the sound pulses will be sent to the right 
channel. In our current setup, we have an arduino uno set to read an analog 
pin and when the voltage is within a range, it turns on a digital port (Check 
ControlArduinoWithSoundBoard.ino code). This way, you can control one sound 
channel, one laser and "separate" TTLs for both.

The first cell will set the experiment parameters. 
The second, third and fourth cells create audio objects for sound playback and 
set some attributes of it. Then, threads are created so .write() commands can 
run on its own thread.

All the following cells send the stimulus to the sound board, each one with its 
own settings. 
"""
#%% Set Parameters

import ControlSoundBoard

#==========#==========#==========#==========#
Rate = 128000

# TTLs Amplification factor. DO NOT CHANGE unless you know what you're doing.
TTLAmpF = 1

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
SoundStimBlockNo = 1
# Duration of pause between blocks
SoundPauseBetweenStimBlocksDur = 5
# Amplification factor (consider using values <= 1). If using one value, keep it in a list.
SoundAmpF = [1.0, 0.9, 0.8, 0.7, 0.6, 0.5, 0.4, 0.3, 0.2, 0.1]
# Noise frequency. If using one freq., keep the list in a list, [[like this]].
NoiseFrequency = [[8000, 10000], [10000, 12000], [12000, 14000], [14000, 16000]]

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
LaserStimBlockNo = 1
# Duration of pause between blocks
LaserPauseBetweenStimBlocksDur = 5


#%% Prepare sound stimulation
Sound, SoundPauseBetweenStimBlocks, StartSound = \
    ControlSoundBoard.GenSound(Rate, SoundPulseDur, SoundPulseNo, SoundAmpF, 
                               NoiseFrequency, TTLAmpF, SoundPrePauseDur, 
                               SoundPostPauseDur, SoundStimBlockNo, 
                               SoundPauseBetweenStimBlocksDur)


#%% Prepare laser stimulation
Laser, LaserPauseBetweenStimBlocks, StartLaser = \
    ControlSoundBoard.GenLaser(Rate, LaserPrePauseDur, LaserPulseDur, 
                               LaserPostPauseDur, LaserPulseNo, 
                               LaserStimBlockNo, 
                               LaserPauseBetweenStimBlocksDur, TTLAmpF)


#%% Prepare sound and laser simultaneous stimulation
SoundAndLaser, SoundAndLaserPauseBetweenStimBlocks, StartSoundAndLaser = \
    ControlSoundBoard.GenSoundLaser(Rate, SoundPrePauseDur, SoundPulseDur, 
                                    SoundPostPauseDur, SoundPulseNo, 
                                    SoundStimBlockNo, 
                                    SoundPauseBetweenStimBlocksDur, SoundAmpF, 
                                    NoiseFrequency, LaserPrePauseDur, 
                                    LaserPulseDur, LaserPostPauseDur, 
                                    LaserPauseBetweenStimBlocksDur, TTLAmpF)


#%% Run sound
StartSound().start()


#%% Run laser
StartLaser().start()


#%% Run sound and laser
StartSoundAndLaser().start()

