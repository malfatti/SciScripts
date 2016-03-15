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

Experiment: stimulation of the brainstem using light and sound, recording with
a silicon probe (16 channels) + 2 tungsten wires + reference screw.
"""
#%% Calibration

CalibrationFile = '/home/cerebro/Malfatti/Data/Test/' + \
                  '20160315153450-SoundMeasurement/SoundIntensity'
#CalibrationFile = '/home/malfatti/NotSynced/SoftwareTest/' + \
#                  'SoundMeasurements/20160125114052-SoundMeasurement/' + \
#                  'SoundIntensity'


#%% GPIAS

## Set experiment details

GPIASTimeBeforeTTL = 50    # in ms
GPIASTimeAfterTTL = 150    # in ms
FilterLow = 3       # High-pass frequency for bandpass filter
FilterHigh = 300     # Low-pass frequency
FilterOrder = 3       # butter order

#==========#==========#==========#==========#

import glob
import ControlSoundBoard

FileList = glob.glob('*.db'); FileList.sort()

ControlSoundBoard.GPIAS(FileList, CalibrationFile, GPIASTimeBeforeTTL, 
                        GPIASTimeAfterTTL, FilterLow, FilterHigh, FilterOrder)

ControlSoundBoard.PlotGPIAS(FileList)