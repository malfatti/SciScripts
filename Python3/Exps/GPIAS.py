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
    AnimalName  = 'A2',
    StimType    = ['Sound'],
    
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
    
    # NoiseFrequency          = [[9000, 11000],
    #                             [12000, 14000]],
    NoiseFrequency          = [[8000, 10000], 
                                [9000, 11000], 
                                [10000, 12000], 
                                [12000, 14000], 
                                [14000, 16000],
                                [8000, 18000]],
    
    # Background and pulse intensities in dB
    BGIntensity     = [60],
    PulseIntensity  = [105],
    
    
    ## === Hardware  === ##
    SoundCh         = 2,
    TTLCh           = 1,
    PiezoCh         = [3,4,5],
    AnalogTTLs      = True,
    
    SoundSystem     = 'Jack-IntelOut-Marantz-IntelIn',
    Setup           = 'GPIAS',
)

GPIAS.Run(**Parameters)
# array([3, 1, 4, 5, 2])
