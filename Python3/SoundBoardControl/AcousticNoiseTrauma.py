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

Rate = 128000
SoundDur = 60 * 70 # in SECONDS!
NoiseFrequency = [[10000, 14000]]
Intensity = [90]
TTLAmpF = 0

CalibrationFile = '/home/cerebro/Malfatti/Data/Test/' + \
                  '20160315202456-SoundMeasurement/SoundIntensity'
#CalibrationFile = '/home/malfatti/NotSynced/SoftwareTest/' + \
#                  'SoundMeasurements/20160125114052-SoundMeasurement/' + \
#                  'SoundIntensity'
#==========#==========#==========#==========#

import ControlSoundBoard
import shelve

with shelve.open(CalibrationFile) as Shelve:
    SoundIntensity = Shelve['SoundIntensity']
    SBInAmpF = Shelve['SBInAmpF']

SoundPulseDur = 0.5
SoundPulseNo = round(SoundDur/SoundPulseDur)


SoundAmpF = [[float(min(SoundIntensity[Hz].keys(), 
                                 key=lambda i: abs(SoundIntensity[Hz][i] - 
                                                   Intensity)))] 
                       for Hz in list(SoundIntensity)]


# Generate sound stimulus
Sound, SoundPauseBetweenStimBlocks, StartSound = ControlSoundBoard.GenSound(
                                                    Rate, SoundPulseDur, 
                                                    SoundPulseNo, SoundAmpF, 
                                                    NoiseFrequency, TTLAmpF, 
                                                    CalibrationFile)

#%% Run!!
""" To stop, close the console :) """
StartSound().start()