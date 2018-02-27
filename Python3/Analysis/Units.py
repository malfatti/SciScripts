#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: T. Malfatti <malfatti@disroot.org>
@date: 2018-02-26
@license: GNU GPLv3 <https://raw.githubusercontent.com/malfatti/SciScripts/master/LICENSE>
@homepage: https://github.com/Malfatti/SciScripts
"""
from DataAnalysis.Units import Klusta

Parameters = dict(
    Group = 'Recovery',
    
    TimeBeforeTTL = 50, 
    TimeAfterTTL = 50, 
    BinSize = 2,
    
    TTLCh = 17,
    AnalogTTLs = True,
    ProbeChSpacing = 50,
    SpksToPlot = 500,
    
    Show = False,
    Save = True,
    Ext = ['svg'],
)

Klusta.Units(**Parameters)


Group = 'Recovery'
TimeBeforeTTL = 50
TimeAfterTTL = 50
BinSize = 2

TTLCh = 17
AnalogTTLs = True
ProbeChSpacing = 50
SpksToPlot = 500

Show = False
Save = True
Ext = ['svg']

