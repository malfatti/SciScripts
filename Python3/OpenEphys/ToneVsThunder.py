#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: T. Malfatti <malfatti@disroot.org>
@date: 2018-01-24
@license: GNU GPLv3 <https://raw.githubusercontent.com/malfatti/SciScripts/master/LICENSE>
@homepage: https://github.com/Malfatti/SciScripts
"""
#%% Sound and Laser stimulation
from Exps import ToneVsThunder

## === Experiment parameters === ##
Parameters = dict(
    AnimalName      = 'RecoveryControl_02',
    StimType        = ['Tone'],
    # StimType        = ['Tone', 'Thunder'],
    
    ## === Tone === ##
    ToneIntensities         = [80, 70, 60, 50, 40],
    ToneFrequency      = [[9000], [10000], [11000], [13000], [15000]],
    
    # Fill all durations in SECONDS!
    TonePauseBeforePulseDur    = 0.004,
    TonePulseDur               = 0.003,
    TonePauseAfterPulseDur     = 0.093,
    TonePulseNo                = 529,
    TonePauseBetweenIntensities     = 10,
    
    
    ## === Thunder === ##
    ThunderIntensities         = [80, 70, 60, 50, 40],
    ThunderFrequency      = [[4000, 22000]],
    
    # Fill all durations in SECONDS!
    ThunderPauseBeforePulseDur    = 0.004,
    ThunderPulseDur               = 0.050,
    ThunderPauseAfterPulseDur     = 0.046,
    ThunderPulseNo                = 529,
    ThunderPauseBetweenIntensities     = 10,
    
    
    ## === Hardware === ##
    SoundCh         = -1,
    TTLCh           = 0,
    AnalogTTLs      = True,
    
    System          = 'Jack-IntelOut-MackieIn-MackieOut-IntelIn',
    Setup           = 'GPIAS',
)


Stimulation, InfoFile = ToneVsThunder.Prepare(**Parameters)

#%%
ToneVsThunder.Play(Stimulation, InfoFile, ['Sound'], DV='4330')