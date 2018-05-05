#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jun 12 15:25:24 2017

@author: malfatti
"""
import subprocess

from glob import glob
from IO import Intan, OpenEphys


def DataLoader(Folder, Unit='uV', ChannelMap=[], AnalogTTLs=True):
    FilesExt = [F[-3:] for F in glob(Folder+'/*.*')]
    
    if 'kwd' in FilesExt: Data, Rate = OpenEphys.KwikLoad(Folder, Unit, ChannelMap)
    elif 'dat' in FilesExt: Data, Rate = OpenEphys.DatLoad(Folder, Unit, ChannelMap)
    elif 'ous' in FilesExt: Data, Rate = OpenEphys.OELoad(Folder, Unit, ChannelMap)
    elif 'int' in FilesExt: Data, Rate = Intan.Load(Folder, ChannelMap)
    else: print('Data format not supported.'); return(None)
    
    if not AnalogTTLs:
        if 'kwd' in FilesExt: 
            Kwds = glob(Folder+'/*.events')
            if len(Kwds) > 1: print('Multiple sessions not supported yet.'); return(None)
            
            EventsDict = 'ToBeContinued'
        else:
            EventsDict = EventsLoad(Folder)
        
        return(Data, Rate, EventsDict)
    else:
        return(Data, Rate)


def RunProcess(Cmd, LogFile=''):
    if LogFile == '': print('Logging disabled, outputting to STDOUT.')
    else: print('Check progress in file', LogFile)
    
    try:
        if LogFile == '': Log = subprocess.PIPE
        else:  Log = open(LogFile, 'w')
        
        P = subprocess.Popen(Cmd,
                             stdout=Log,
                             stderr=subprocess.STDOUT)
        
        print('Process id:', P.pid )
        P.communicate()[0]; ReturnCode = P.returncode
        if LogFile != '': Log.close()
    
    except Exception as e:
        ReturnCode = 1; print(e)
    
    return(ReturnCode)


