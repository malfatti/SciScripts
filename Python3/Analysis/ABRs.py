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

Group = 'Prevention'; Exp = '20170718-Prevention_A5-ABRs'
AnalysisFile = Group+'/'+Group+'-Analysis.hdf5'
#AnalysisFile = 'Test.hdf5'
AnalysisPath = '/'+Exp
InfoFile = glob(Group+'/'+Exp+'/*.dict')[0]
# Folders = glob(Animal+'/'+Exp+'/KwikFiles/*'); Folders.sort()
Folders = glob(Group+'/'+Exp+'/2017-*'); Folders.sort()

ABRCh=[1]; ABRTTLCh=0

ABRs.Analysis(Exp, Folders, InfoFile, AnalysisFile, ABRCh, ABRTTLCh, StimType=['NaCl', 'CNO'])
ABRPlot.Traces(AnalysisPath, AnalysisFile, InfoFile, Group+'/Figs', Save=True, Visible=True)
