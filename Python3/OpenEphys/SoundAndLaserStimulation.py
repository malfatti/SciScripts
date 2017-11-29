#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: T. Malfatti
@year: 2015
@licence: GNU GPLv3 <https://raw.githubusercontent.com/malfatti/SciScripts/master/LICENSE>
@homepage: https://github.com/Malfatti/SciScripts
"""
#%% Settings
from Exps import SoundAndLaserStimulation

## === Experiment parameters === ##
Parameters = dict(
    AnimalName =        'Prevention_A5',
    StimType =          ['Sound'],
    # StimType =          ['Sound', 'Laser'],
    
    
    ## === Sound === ##
    Intensities =       [80, 70, 60, 50, 40],
    NoiseFrequency =    [[8000, 10000], 
                         [9000, 11000], 
                         [10000, 12000], 
                         [12000, 14000], 
                         [14000, 16000]],
    
    # Fill all durations in SECONDS!
    SoundPauseBeforePulseDur =      0.004,
    SoundPulseDur =                 0.003,
    SoundPauseAfterPulseDur =       0.093,
    SoundPulseNo =                  529,
    SoundPauseBetweenIntensities =  10,
    
    
    ## === Laser === ##
    LaserPauseBeforePulseDur =          0,
    LaserPulseDur =                     0.01,
    LaserPauseAfterPulseDur =           0.09,
    LaserPulseNo =                      529,
    LaserStimBlockNo =                  5,
    LaserPauseBetweenStimBlocksDur =    10,
    
    
    ## === Hardware === ##
    ABRCh =         [1,2],
    SoundCh =       26,
    TTLCh =         27,
    AnalogTTLs =    True,
    
    System =            'Jack-IntelOut-MackieIn-MackieOut-IntelIn',
    Setup =             'GPIAS',
)


SoundAndLaserStimulation.Run(**Parameters)
