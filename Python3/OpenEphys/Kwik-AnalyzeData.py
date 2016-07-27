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
StimType = ['Sound_NaCl', 'Sound_CNO']
#StimType = ['Sound']
Board = 'OE'

#==========#==========#==========#==========#

import KwikAnalysis
from glob import glob

FileName = glob('*.hdf5')[0]
KwikAnalysis.ABRAnalysis(FileName, ABRCh, ABRTTLCh, ABRTimeBeforeTTL, ABRTimeAfterTTL, 
                 FilterFreq, FilterOrder, StimType, AnalogTTLs, Board)

AnalysisFile = glob('../*.hdf5')[0]
KwikAnalysis.ABRPlot(AnalysisFile, FileName, Visible=True)
KwikAnalysis.ABRPlot3D(AnalysisFile, FileName, Visible=True)


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

#==========#==========#==========#==========#

import KwikAnalysis
from glob import glob

FileName = glob('*.hdf5'); FileName.sort(); FileName = FileName[RecFolder-1]
KwikAnalysis.GPIASAnalysis(RecFolder, FileName, GPIASCh, GPIASTTLCh, 
                             GPIASTimeBeforeTTL, GPIASTimeAfterTTL, FilterFreq, 
                             FilterOrder, AnalogTTLs)

AnalysisFile = glob('../*.hdf5')[0]
KwikAnalysis.GPIASPlot(AnalysisFile, FileName, Visible=True)


#%% TTLsLatencyTest

SoundCh = 1
SoundSqCh = 5
SoundTTLCh = 1

TimeBeforeTTL = 10   # in ms
TimeAfterTTL = 10    # in ms

FileName = glob.glob('*.hdf5')[0]


#%% Clustering
StimTTLCh = 17
StimType0 = ['Sound_NaCl']
StimType1 = ['Sound_CNO']
AnalogTTLs=True
Board='OE'
Override = {}
#Override = {'Rec'}

#==========#==========#==========#==========#

import KwikAnalysis
from glob import glob
from multiprocessing import Process

FileName = glob('*.hdf5')[0]

Clus_NaCl = Process(target=KwikAnalysis.ClusterizeSpks, args=(FileName, StimType0, AnalogTTLs, Board, Override))
Clus_CNO = Process(target=KwikAnalysis.ClusterizeSpks, args=(FileName, StimType1, AnalogTTLs, Board, Override))
Clus_NaCl.start(); Clus_CNO.start()
print('NaClPid =', str(Clus_NaCl.pid)); print('CNOPid =', str(Clus_CNO.pid))
Clus_NaCl.join(); Clus_CNO.join()


#%% Units
StimTTLCh = 17
PSTHTimeBeforeTTL = 0
PSTHTimeAfterTTL = 300
StimType = ['Sound_NaCl', 'Sound_CNO']
AnalogTTLs=True
Board='OE'
#Override = {'Rec':'0'}
Override0 = {'Rec':'0', 'Stim':StimType[0]}
Override1 = {'Rec':'0', 'Stim':StimType[1]}

#==========#==========#==========#==========#

import KwikAnalysis
from glob import glob
from multiprocessing import Process

FileName = glob('*.hdf5')[0]

UnitsPSTH_NaCl = Process(target=KwikAnalysis.UnitsPSTH, args=(
                                                           FileName, StimTTLCh, 
                                                           PSTHTimeBeforeTTL, 
                                                           PSTHTimeAfterTTL, 
                                                           StimType, AnalogTTLs, 
                                                           Board, Override0))
UnitsPSTH_CNO = Process(target=KwikAnalysis.UnitsPSTH, args=(
                                                           FileName, StimTTLCh, 
                                                           PSTHTimeBeforeTTL, 
                                                           PSTHTimeAfterTTL, 
                                                           StimType, AnalogTTLs, 
                                                           Board, Override1))
UnitsSpks_NaCl = Process(target=KwikAnalysis.UnitsSpks, args=(
                                                           FileName, StimType, 
                                                           Board, Override0))
UnitsSpks_CNO = Process(target=KwikAnalysis.UnitsSpks, args=(
                                                           FileName, StimType, 
                                                           Board, Override1))

UnitsPSTH_NaCl.start(); print('NaClPid =', str(UnitsPSTH_NaCl.pid))
UnitsPSTH_NaCl.join()

UnitsPSTH_CNO.start(); print('CNOPid =', str(UnitsPSTH_CNO.pid))
UnitsPSTH_CNO.join()


UnitsSpks_NaCl.start(); print('NaClPid =', str(UnitsSpks_NaCl.pid))
UnitsSpks_NaCl.join()

UnitsSpks_CNO.start(); print('CNOPid =', str(UnitsSpks_CNO.pid))
UnitsSpks_CNO.join()


AnalysisFile = glob('../*.hdf5')[0]
UnitsPlot_NaCl = Process(target=KwikAnalysis.UnitsSpksPSTH_ToSVG, args=(
                                                           AnalysisFile, 
                                                           FileName, Override0))
UnitsPlot_CNO = Process(target=KwikAnalysis.UnitsSpksPSTH_ToSVG, args=(
                                                           AnalysisFile, 
                                                           FileName, Override1))


UnitsPlot_NaCl.start(); UnitsPlot_CNO.start()
print('NaClPid =', str(UnitsPlot_NaCl.pid))
print('CNOPid =', str(UnitsPlot_CNO.pid))
UnitsPlot_NaCl.join(); UnitsPlot_CNO.join()

