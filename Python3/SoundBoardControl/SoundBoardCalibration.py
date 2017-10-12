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

This script can be used to calibrate the sound board in and out amplification 
factor.

For output:
    Write sine waves of amp 1 to sound board output. We use this to get the 
    sound board amp factor, that is, how much does the sound board being used 
    is increasing or decreasing the signal written to it.
    
For input:
    Read signal from sound board input. You have to apply a known amplitude 
    signal (a sine wave, for example) to the sound board input, so you can 
    check if the voltage you applied is the voltage being read by the sound 
    board.

It is very important to set the volume of the soundboard to unit level 
(usually 0dB, which is 100% of the intensity) so you know that no kind of 
frequency filter is being applied.
"""
#%% Set calibration
Rate = 192000; Freq = 10000; WaveDur = 10
SoundSystem = 'Jack-IntelOut-MackieIn-MackieOut-IntelIn'

from IO import SoundCard
import os, h5py
import numpy as np

Folder = os.environ['DATAPATH']+'/Tests/SoundMeasurements'
FileName = Folder + '/' + 'SoundMeasurements.hdf5'

#%% Output

SoundCard.SoundCalOut(Rate, Freq, WaveDur, Ch=1)
# SBOutAmpF is the generated signal divided by the measured signal
SBOutAmpF = 1/10

#%% Input
Repetitions = 4

SBInAmpF = np.zeros(Repetitions, dtype=np.float32)
for aa in range(Repetitions):
    Rec = SoundCard.SoundCalIn(Rate, Freq, WaveDur, SBOutAmpF, Ch=1)
    SBInAmpF[aa] = (max(Rec)-(min(Rec)))/2
    print(SBInAmpF[aa])

# Mean
SBInAmpF = sum(SBInAmpF)/len(SBInAmpF)
# SBInAmpF is the real amplitude divided by the measured amplitude
SBInAmpF = 1/SBInAmpF

print('SBInAmpF = ', str(SBInAmpF))

#%% NoiseRMS
import sounddevice as SD
NoiseRMS = SD.rec(Rate*10, Rate, 1, blocking=True)
NoiseRMS = ((np.mean(NoiseRMS**2))**0.5) * SBInAmpF


#%% Save
with h5py.File(FileName, 'w') as F:
    if SoundSystem not in F.keys(): F.create_group(SoundSystem)
    
    F[SoundSystem]['SBOutAmpF'] = SBOutAmpF
    F[SoundSystem]['SBInAmpF'] = SBInAmpF
    F[SoundSystem]['NoiseRMS'] = NoiseRMS

with h5py.File(FileName, 'r') as F:
    a = F[SoundSystem]['SBOutAmpF'][()]
    b = F[SoundSystem]['SBInAmpF'][()]
    c = F[SoundSystem]['NoiseRMS'][()]

