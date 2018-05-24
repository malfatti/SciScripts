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
# 'A1', 'C4', 'A4', 
# 'C3', 'C1', 'D2', 
# 'D4', 'B1', 'B3', 
# 'A2', 'B4','D1', 
# 'B2', 'A5', 'A3'

Parameters = dict(
    AnimalName      = 'A1_C4_A4',
    StimType        = ['Sound'],
    
    Intensities     = [90],
    NoiseFrequency  = [[9000, 11000]],
    SoundPulseDur   = 60,                 # in MINUTES!
    
    ## === Hardware === ##
    System  = 'Jack-IntelOut-Marantz-IntelIn',
    Setup   = 'GPIAS',
)

AcousticNoiseTrauma.Run(**Parameters)
