#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Oct  4 11:43:58 2017

@author: malfatti
"""

#%% Acoustic trauma
from Exps import AcousticNoiseTrauma

## Experiment parameters
Parameters = dict(
    AnimalName =        'Prevention_A3_4_5',
    StimType =          ['Sound'],
    
    Intensities =       [95],
    NoiseFrequency =    [[8000, 12000]],
    SoundPulseDur =     90,                 # in MINUTES!
    
    System =            'Jack-IntelOut-MackieIn-MackieOut-IntelIn',
    Setup =             'GPIAS',
)

AcousticNoiseTrauma.Run(**Parameters)

