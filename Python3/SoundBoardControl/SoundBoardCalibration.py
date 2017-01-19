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
Rate = 192000; Freq = 1000; WaveDur = 10
Connections = 'IntelOut-MackieIn-MackieChOut-Speaker'

import ControlSoundBoard
from datetime import datetime
import h5py
import numpy as np

#%% Output

ControlSoundBoard.SoundCalOut(Rate, Freq, WaveDur)
# SBOutAmpF is the generated signal divided by the measured signal
SBOutAmpF = 1/1.8

#%% Input
Repetitions = 4

SBInAmpF = np.zeros(Repetitions, dtype=np.float32)
for aa in range(Repetitions):
    Rec = ControlSoundBoard.SoundCalIn(Rate, Freq, WaveDur, SBOutAmpF)
    SBInAmpF[aa] = (max(Rec)+(min(Rec)*-1))/2
    print(SBInAmpF[aa])

# Mean
SBInAmpF = sum(SBInAmpF)/len(SBInAmpF)
# SBInAmpF is the real amplitude divided by the measured amplitude
SBInAmpF = 1/SBInAmpF

print('SBInAmpF = ', str(SBInAmpF))

#%% Save
Date = datetime.now()
FileName = Date.strftime("%Y%m%d%H%M%S") + '-SBAmpFs.hdf5'
#FileName = '20161013123915-SBAmpFs.hdf5'
with h5py.File(FileName) as h5:
    h5.create_group(SoundBoard)
    h5[SoundBoard]['SBOutAmpF'] = SBOutAmpF
    h5[SoundBoard]['SBInAmpF'] = SBInAmpF

with h5py.File(FileName) as h5:
    a = h5[SoundBoard]['SBOutAmpF'][()]
    b = h5[SoundBoard]['SBInAmpF'][()]

"""
Malfatti = '/home/malfatti/Documents/PhD/Tests/20161013174053-SBAmpFs.hdf5'
Cerebro = '/home/cerebro/Malfatti/Test/20160418173048-SBAmpFs.hdf5'
"""