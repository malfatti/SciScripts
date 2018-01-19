# -*- coding: utf-8 -*-
"""
@author: T. Malfatti <malfatti@disroot.org>
@year: 2015
@license: GNU GPLv3 <https://raw.githubusercontent.com/malfatti/SciScripts/master/LICENSE>
@homepage: https://github.com/Malfatti/SciScripts
"""

#%% GPIAS
from Exps import GPIAS

Parameters = dict(
    # array([4, 1, 3, 2, 5])
    AnimalName  = 'RecoveryControl_03',
    StimType    = ['Sound', 'CNO'],
    
    # Number of trials per freq. tested 
    # 1 trial is composed by 1 stim w/ gap and 1 stim w/o gap
    NoOfTrials  = 9,
    
    
    ## === Sound  === ##
    # Fill all durations in SECONDS!
    SoundBGDur              = 2.3,
    SoundGapDur             = 0.04,
    SoundBGPrePulseDur      = 0.1,
    SoundLoudPulseDur       = 0.05,
    SoundBGAfterPulseDur    = 0.51,
    SoundBetweenStimDur     = [10, 20],
    NoiseFrequency          = [[8000, 10000], 
                               [9000, 11000], 
                               [10000, 12000], 
                               [12000, 14000], 
                               [14000, 16000]],
    
    # Background and pulse intensities in dB
    BGIntensity     = [65],
    PulseIntensity  = [105],
    
    
    ## === Hardware  === ##
    SoundCh         = 3,
    TTLCh           = 2,
    PiezoCh         = [1],
    AnalogTTLs      = True,
    
    SoundSystem     = 'Jack-IntelOut-MackieIn-MackieOut-IntelIn',
    Setup           = 'GPIAS',
)

GPIAS.Run(**Parameters)

