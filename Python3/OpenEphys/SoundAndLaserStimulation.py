#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: T. Malfatti <malfatti@disroot.org>
@year: 2015
@license: GNU GPLv3 <https://raw.githubusercontent.com/malfatti/SciScripts/master/LICENSE>
@homepage: https://github.com/Malfatti/SciScripts
"""
#%% Sound and Laser stimulation
from Exps import SoundAndLaserStimulation

## === Experiment parameters === ##
# 3, 2, 4, 5, 1
# array([[ 8000, 10000],
#        [ 9000, 11000],
#        [14000, 16000],
#        [12000, 14000],
#        [10000, 12000]])
Parameters = dict(
    AnimalName          = 'RecoveryControl_02',
    StimType            = ['Sound'],
    # StimType            = ['Sound', 'Laser', 'SoundLaser'],
    
    
    ## === Sound === ##
    Intensities         = [80, 70, 60, 50, 40],
    NoiseFrequency      = [[8000, 10000], 
                           [9000, 11000], 
                           [10000, 12000], 
                           [12000, 14000], 
                           [14000, 16000]],
    
    # Fill all durations in SECONDS!
    SoundPauseBeforePulseDur    = 0.004,
    SoundPulseDur               = 0.003,
    SoundPauseAfterPulseDur     = 0.093,
    SoundPulseNo                = 529,
    PauseBetweenIntensities     = 10,
    
    
    ## === Laser === ##
    # 'Sq' for square pulses, 'Sin' for sin wave
    LaserType                       = 'Sq',
    
    # if LaserType == 'Sq'
    LaserPauseBeforePulseDur        = 0,
    LaserPulseDur                   = 0.01,
    LaserPauseAfterPulseDur         = 0.09,
    LaserPulseNo                    = 529,
    
    # if LaserType == 'Sin'
    LaserDur                        = 0.1*529,
    LaserFreq                       = 10,      # in Hz
    
    LaserStimBlockNo                = 5,
    LaserPauseBetweenStimBlocksDur  = 10,
    
    
    ## === Hardware === ##
    ABRCh           = [1],
    SoundCh         = 26,
    TTLCh           = 27,
    AnalogTTLs      = True,
    
    System          = 'Jack-IntelOut-MackieIn-MackieOut-IntelIn',
    Setup           = 'GPIAS',
)


Stimulation, InfoFile = SoundAndLaserStimulation.Prepare(**Parameters)

#%%
SoundAndLaserStimulation.Play(Stimulation, InfoFile, ['Sound'], DV='4330')
