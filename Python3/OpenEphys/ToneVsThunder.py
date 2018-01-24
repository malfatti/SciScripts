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
    StimType        = ['Tone', 'Thunder'],
    
    ## === Tone === ##
    ToneIntensities         = [90, 60],
    ToneFrequency           = [20000],
    
    # Fill all durations in SECONDS!
    TonePauseBeforePulseDur    = 0,
    TonePulseDur               = 0.05,
    TonePauseAfterPulseDur     = 4.95,
    TonePulseNo                = 30,
    TonePauseBetweenIntensities     = 10,
    
    
    ## === Thunder === ##
    ThunderIntensities         = [90, 60],
    ThunderFrequency      = [[8000, 18000]],
    
    # Fill all durations in SECONDS!
    ThunderPauseBeforePulseDur    = 0,
    ThunderPulseDur               = 0.05,
    ThunderPauseAfterPulseDur     = 4.95,
    ThunderPulseNo                = 30,
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
ToneVsThunder.Play(Stimulation, InfoFile, ['Thunder'], DV='4330', Ramp=False)
