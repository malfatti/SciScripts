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
    Read signal from sound board input. You have to apply a known voltage to 
    sound board input, so you can check if the voltage you applied is the 
    voltage being read by the sound board. 
    We do this by putting two resistors (R1 = 9R2) in the positive and 
    negative of any 5V supply (Arduino, open usb cable...), so in the junction 
    R1-R2 the output is 0.5V.
    In the end, our circuit is this:
        
            |----900ohm ----|
        5V -|               |
                            |---------SoundBoardIn+
        Gnd-|               |
            |----100ohm ----|       |-SoundBoardIn-
                                    |
                                   Gnd

It is very important to set the volume of the soundboard to 0dB (which is 100%) 
so you know that no kind of frequency filter is being applied.
"""
#%% Set calibration
Rate = 128000

import ControlSoundBoard
from scipy import signal

#%% Output

ControlSoundBoard.SoundCalOut(Rate)

#%% Input

Data = ControlSoundBoard.SoundCalIn(Rate)
F, PxxSp = signal.welch(Data, Rate, nperseg=1024, scaling='spectrum')

RMS = (sum(PxxSp) * (F[1] - F[0]))**0.5

print(RMS)