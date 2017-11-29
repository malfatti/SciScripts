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

This is a script to generate and control sound stimulation for gap-prepulse 
inhibition of acoustic startle reflex (GPIAS).
"""

#%% Settings
from Exps import GPIAS

Parameters = dict(
    AnimalName = 'RecoveryControl_04',
    StimType = ['Sound'],
    SoundCh = 3,
    TTLCh = 2,
    PiezoCh = [1],
    
    # Number of trials per freq. tested (1 trial = 1 stim w/ gap + 1 stim w/o gap)
    NoOfTrials = 9,
    
    # Fill all durations in SECONDS!
    SoundBGDur = 2.3,
    SoundGapDur = 0.04,
    SoundBGPrePulseDur = 0.1,
    SoundLoudPulseDur = 0.05,
    SoundBGAfterPulseDur = 0.51,
    SoundBetweenStimDur = [10, 20],
    NoiseFrequency = [[8000, 10000]],
    # NoiseFrequency = [[8000, 10000], [9000, 11000], [10000, 12000], 
    #                   [12000, 14000], [14000, 16000], [16000, 18000], [8000, 18000]],
    
    # Background and pulse intensities in dB. Supports float :)
    BGIntensity = [65],
    PulseIntensity = [105],
    
    
    ## Hardware parameters
    SoundSystem = 'Jack-Speaker-Mic',
    Setup = 'Test',
    # SoundSystem = 'Jack-IntelOut-MackieIn-MackieOut-IntelIn',
    # Setup = 'GPIAS',
)

GPIAS.Run(**Parameters)

