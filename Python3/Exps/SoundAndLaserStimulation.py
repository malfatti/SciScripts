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
# 'B4', 'A4', 'B1', 
# 'A1', 'D1', 'A3', 
# 'A5', 'B2', 'B3', 
# 'A2', 'C3', 'D2', 
# 'C1', 'D4', 'C4'
# array([[10000, 12000],
#        [ 8000, 10000],
#        [ 9000, 11000],
#        [12000, 14000],
#        [14000, 16000]])
Parameters = dict(
    AnimalName          = 'A2',
    StimType            = ['Sound'],
    # StimType            = ['Sound', 'Laser', 'SoundLaser'],
    
    
    ## === Sound === ##
    Intensities         = [80, 75, 70, 65, 60, 55, 50, 45, 40, 35],
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
    
    
    ## === Probe === ##
    Remapped    = False,
    Probe       = None,     # None | 'A16'
    Adaptor     = 'A16OM16',     # None | 'CustomAdaptor' | 'RHAHeadstage' | 'A16OM16'
    
    
    ## === Hardware === ##
    ABRCh           = [1],
    SoundCh         = 2,
    TTLCh           = 3,
    AnalogTTLs      = True, 
    
    System          = 'Jack-IntelOut-Marantz-IntelIn',
    Setup           = 'UnitRec',
)


Stimulation, InfoFile = SoundAndLaserStimulation.Prepare(**Parameters)

#%%
SoundAndLaserStimulation.Play(Stimulation, InfoFile, ['Sound'], DV='4330')
