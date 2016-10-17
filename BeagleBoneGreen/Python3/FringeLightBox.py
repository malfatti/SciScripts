#!/usr/bin/env python3
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

This is a script to simulate the light box from Fringe. The EEG signal will be 
filtered and the highest the power in that frequency, more leds will be turned 
off.
"""

import numpy as np
import sounddevice as SD

Rate=48000
Freq=1000

SD.default.samplerate=Rate
SD.default.channels=1

A = np.random.uniform(-1, 1, size=Rate)
B = [np.sin(2 * np.pi * Freq * (_/Rate))  for _ in range(Rate)]

#SD.play(A, blocking=True)

C = SD.playrec(B, blocking=True)