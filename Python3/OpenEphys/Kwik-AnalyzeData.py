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

ABRCh = [1, 2]         # [RightChannel, LeftChannel], if order matters
ABRTimeBeforeTTL = 3    # in ms
ABRTimeAfterTTL = 9    # in ms
ABRTTLCh = 2            # TTL ch for ABR
FilterLow = 300         # High-pass frequency for bandpass filter
FilterHigh = 3000       # Low-pass frequency
FilterOrder = 5         # butter order
StimType = 'Sound'

#==========#==========#==========#==========#

import glob
import KwikAnalysis

FileName = glob.glob('*.hdf5'); FileName = FileName[0]

KwikAnalysis.ABR(FileName, ABRCh, ABRTimeBeforeTTL, ABRTimeAfterTTL, ABRTTLCh, 
                 FilterLow, FilterHigh, FilterOrder, StimType)

KwikAnalysis.PlotABR(FileName)


#%% GPIASs

## Set experiment details

PiezoCh = 1
GPIASTTLCh = 1
GPIASTimeBeforeTTL = 50    # in ms
GPIASTimeAfterTTL = 150    # in ms
FilterLow = 3       # High-pass frequency for bandpass filter
FilterHigh = 300     # Low-pass frequency
FilterOrder = 3       # butter order

#==========#==========#==========#==========#

import glob
import KwikAnalysis

FileList = glob.glob('*.db'); FileList.sort()

KwikAnalysis.GPIAS(GPIASTimeBeforeTTL, GPIASTimeAfterTTL, FilterLow, 
                   FilterHigh, FilterOrder, GPIASTTLCh, PiezoCh)
KwikAnalysis.PlotGPIAS(FileList)

#%% TTLsLatencyTest

SoundCh = 1
SoundSqCh = 2
SoundTTLCh = 1

TimeBeforeTTL = 5   # in ms
TimeAfterTTL = 8    # in ms

FileName = glob.glob('*.hdf5'); FileName = FileName[0]