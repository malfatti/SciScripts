# -*- coding: utf-8 -*-
"""
    Copyright (C) 2015  T. Malfatti
    
    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.
    
    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.
    
    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.

"""
#%% ABRs

## Set experiment details

ABRCh = [10]         # [RightChannel, LeftChannel], if order matters
ABRTimeBeforeTTL = 3    # in ms
ABRTimeAfterTTL = 9    # in ms
ABRTTLCh = 17            # TTL ch for ABR
FilterFreq = [300, 3000]         # frequency for filter
FilterOrder = 5         # butter order
AnalogTTLs = True
StimType = 'Sound'      # Stimulation type: 'Sound', 'Laser' or ['Sound', 'Laser']
Board = 'OE'

#==========#==========#==========#==========#

import glob
import KwikAnalysis

FileName = glob.glob('*.hdf5'); FileName = FileName[0]

KwikAnalysis.ABR(FileName, ABRCh, ABRTimeBeforeTTL, ABRTimeAfterTTL, ABRTTLCh, 
                 FilterFreq, FilterOrder, StimType, AnalogTTLs, Board)

KwikAnalysis.PlotABR2(FileName)


#%% GPIASs

## Set experiment details

PiezoCh = 1
GPIASTTLCh = 2
GPIASTimeBeforeTTL = 20    # in ms
GPIASTimeAfterTTL = 100    # in ms
FilterLow = 70       # High-pass frequency for bandpass filter
FilterHigh = 400     # Low-pass frequency
FilterOrder = 3       # butter order
RecFolder = 1

import glob
import KwikAnalysis

FileName = glob.glob('*.hdf5'); FileName.sort()
FileName = FileName[RecFolder-1]

KwikAnalysis.GPIASAnalogTTLs(RecFolder, FileName, GPIASTimeBeforeTTL, 
                             GPIASTimeAfterTTL, FilterLow, FilterHigh, 
                             FilterOrder, GPIASTTLCh, PiezoCh)

KwikAnalysis.PlotGPIAS2(FileName)

#%% TTLsLatencyTest

SoundCh = 1
SoundSqCh = 5
SoundTTLCh = 1

TimeBeforeTTL = 10   # in ms
TimeAfterTTL = 10    # in ms

FileName = glob.glob('*.hdf5'); FileName = FileName[0]