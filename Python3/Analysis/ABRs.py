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

Animal = 'Prevention'; Exp = '20170623-Prevention_A2-ABRs'
AnalysisFile = Animal+'/'+Animal+'-Analysis.hdf5'
#AnalysisFile = 'Test.hdf5'
AnalysisPath = '/'+Exp
InfoFile = glob(Animal+'/'+Exp+'/*.hdf5')[0]
# Folders = glob(Animal+'/'+Exp+'/KwikFiles/*'); Folders.sort()
Folders = glob(Animal+'/'+Exp+'/2017-*'); Folders.sort()

ABRCh=[1]; ABRTTLCh=0

ABRs.Analysis(Exp, Folders, InfoFile, AnalysisFile, ABRCh, ABRTTLCh)
ABRPlot.Traces(AnalysisPath, AnalysisFile, InfoFile, Animal+'/Figs', Save=True, Visible=True)
