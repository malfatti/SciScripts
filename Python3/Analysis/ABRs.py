#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Jul  8 16:40:55 2017

@author: malfatti
"""
#%% Single exp
from DataAnalysis import ABRs
from DataAnalysis.Plot import ABRs as ABRPlot
from glob import glob

Animal = 'CaMKIIahM4Dn06'; Exp = 'CaMKIIahM4Dn06-20160514-ABRs'
#AnalysisFile = Animal+'/'+Animal+'-Analysis.hdf5'
AnalysisFile = 'Test.hdf5'
AnalysisPath = '/'+Exp
InfoFile = Animal+'/'+Exp+'/20160514165646-CaMKIIahM4Dn06-SoundStim.hdf5'
Folders = glob(Animal+'/'+Exp+'/KwikFiles/*'); Folders.sort()
ABRCh=[1]; ABRTTLCh=21

ABRs.Analysis(Exp, Folders, InfoFile, AnalysisFile, ABRCh, ABRTTLCh)
ABRPlot.Traces(AnalysisPath, AnalysisFile, InfoFile, Save=False, Visible=True)
