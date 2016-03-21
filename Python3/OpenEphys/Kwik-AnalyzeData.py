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

Experiment: stimulation of the brainstem using light and sound, recording with
a silicon probe (16 channels) + 2 tungsten wires + reference screw.
"""
#%% ABRs

## Set experiment details

ABRCh = [1, 2]         # [RightChannel, LeftChannel], if order matters
ABRTimeBeforeTTL = 0    # in ms
ABRTimeAfterTTL = 12    # in ms
ABRTTLCh = 1            # TTL ch for ABR
FilterLow = 300         # High-pass frequency for bandpass filter
FilterHigh = 3000       # Low-pass frequency
FilterOrder = 4         # butter order

#==========#==========#==========#==========#

import glob
import KwikAnalysis

FileName = glob.glob('*.db'); FileName = FileName[0][:-3]

KwikAnalysis.ABR(FileName, ABRCh, ABRTimeBeforeTTL, ABRTimeAfterTTL, ABRTTLCh, 
                 FilterLow, FilterHigh, FilterOrder)

FileName = glob.glob('*ABRs.db'); FileName = FileName[0][:-3]
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