#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jun 12 15:25:24 2017

@author: malfatti
"""
import subprocess

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


