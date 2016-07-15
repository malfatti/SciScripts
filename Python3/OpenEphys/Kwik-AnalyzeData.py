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
ABRTTLCh = 17            # TTL ch for ABR
ABRTimeBeforeTTL = 3    # in ms
ABRTimeAfterTTL = 9    # in ms
FilterFreq = [300, 3000]         # frequency for filter
FilterOrder = 5         # butter order
AnalogTTLs = True
StimType = 'Sound'      # Stimulation type: 'Sound', 'Laser' or ['Sound', 'Laser']
Board = 'OE'

#==========#==========#==========#==========#

from glob import glob
import KwikAnalysis

FileName = glob('*.hdf5'); FileName = FileName[0]

KwikAnalysis.ABR(FileName, ABRCh, ABRTTLCh, ABRTimeBeforeTTL, ABRTimeAfterTTL, 
                 FilterFreq, FilterOrder, StimType, AnalogTTLs, Board)

KwikAnalysis.PlotABR(FileName)


#%% GPIASs

## Set experiment details

GPIASCh = [1]
GPIASTTLCh = 2
GPIASTimeBeforeTTL = 20    # in ms
GPIASTimeAfterTTL = 100    # in ms
FilterFreq = [70, 400]     # frequency for filter
FilterOrder = 3       # butter order
AnalogTTLs = True
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

#%% Units
StimTTLCh = 17
PSTHTimeBeforeTTL = 0
PSTHTimeAfterTTL = 300
StimType = ['Sound_NaCl', 'Sound_CNO']
AnalogTTLs=True
Board='OE'
OverrideRec = 0

from glob import glob
import KwikAnalysis

FileName = glob('*.hdf5'); FileName = FileName[0]

KwikAnalysis.ClusterizeAll(FileName, StimType=['Sound'], AnalogTTLs=False, 
                           Board='OE', OverrideRec=[])
KwikAnalysis.UnitAnalysis(FileName, StimTTLCh, PSTHTimeBeforeTTL, 
                          PSTHTimeAfterTTL, StimType, AnalogTTLs, Board, 
                          OverrideRec)

#%% Clustering
StimTTLCh = 17
PSTHTimeBeforeTTL = 0
PSTHTimeAfterTTL = 300
StimType0 = ['Sound_NaCl']
StimType1 = ['Sound_CNO']
AnalogTTLs=True
Board='OE'
OverrideRec = 0

import KwikAnalysis
from glob import glob
from multiprocessing import Process

FileName = glob('*.hdf5'); FileName = FileName[0]

Clus_NaCl = Process(target=KwikAnalysis.ClusterizeAll, args=(FileName, StimType0, AnalogTTLs, Board, OverrideRec))
Clus_CNO = Process(target=KwikAnalysis.ClusterizeAll, args=(FileName, StimType1, AnalogTTLs, Board, OverrideRec))
Clus_NaCl.start(); Clus_CNO.start()
Clus_NaCl.join(); Clus_CNO.join()