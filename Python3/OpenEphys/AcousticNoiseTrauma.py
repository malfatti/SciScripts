#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: T. Malfatti <malfatti@disroot.org>
@date: 2017-10-04
@license: GNU GPLv3 <https://raw.githubusercontent.com/malfatti/SciScripts/master/LICENSE>
@homepage: https://github.com/Malfatti/SciScripts
"""

#%% Acoustic trauma
from Exps import AcousticNoiseTrauma

## === Experiment parameters === ##
Parameters = dict(
    AnimalName      = 'Prevention_A3_4_5',
    StimType        = ['Sound'],
    
    Intensities     = [95],
    NoiseFrequency  = [[9000, 11000]],
    SoundPulseDur   = 60,                 # in MINUTES!
    
    ## === Hardware === ##
    System  = 'Jack-IntelOut-MackieIn-MackieOut-IntelIn',
    Setup   = 'GPIAS',
)

AcousticNoiseTrauma.Run(**Parameters)
