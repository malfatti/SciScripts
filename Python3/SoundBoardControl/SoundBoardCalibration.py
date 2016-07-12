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
    Write square pulses of amp 1 to sound board output. We use this to get the 
    sound board amp factor, that is, how much does the sound board being used 
    increase or decrease the signal written to it.
    
For input:
    Read signal from sound board input. You have to apply a known amplitude 
    signal (a sine wave, for example) sound board input, so you can check if 
    the voltage you applied is the voltage being read by the sound board. 
    For this, you can open 2 consoles, run the output cell in one and the 
    input cell in the other.

It is very important to set the volume of the soundboard to 0dB (which is 100%) 
so you know that no kind of frequency filter is being applied.
"""
#%% Set calibration
Rate = 128000; Freq = 1000; WaveDur = 10
SoundBoard = 'Intel_oAnalog-iAnalog'

import ControlSoundBoard
from datetime import datetime
import h5py

#%% Output

ControlSoundBoard.SoundCalOut(Rate, Freq, WaveDur)

#%% Input
Repetitions = 20
SBOutAmpF = 2

Data = [[] for _ in range(Repetitions)]
SBInAmpF = [[] for _ in range(Repetitions)]
for aa in range(Repetitions):
    Data[aa] = ControlSoundBoard.SoundCalIn(Rate, Freq, WaveDur, SBOutAmpF)
    SBInAmpF[aa] = (max(Data[aa])+(min(Data[aa])*-1))/2
    print(SBInAmpF[aa])

SBInAmpF = sum(SBInAmpF)/len(SBInAmpF)

print('SBInAmpF = ', str(SBInAmpF))

#%% Save
Date = datetime.now()
FileName = Date.strftime("%Y%m%d%H%M%S") + '-SBAmpFs.hdf5'
#FileName = '20160418173048-SBAmpFs.hdf5'
with h5py.File(FileName) as h5:
    h5.create_group(SoundBoard)
    h5[SoundBoard].attrs['SBOutAmpF'] = SBOutAmpF
    h5[SoundBoard].attrs['SBInAmpF'] = SBInAmpF

with h5py.File(FileName) as h5:
    a = h5[SoundBoard].attrs['SBOutAmpF']
    b = h5[SoundBoard].attrs['SBInAmpF']

"""
Malfatti = '/home/malfatti/Documents/PhD/Tests/20160712135926-SBAmpFs.hdf5'
"""