#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Feb  6 12:23:43 2018

@author: malfatti
"""
import os
import numpy as np

from datetime import datetime
from serial import Serial
from serial.tools.list_ports import comports
from time import time

def CreateObj(BaudRate):
    Port = comports()
    Arduino = Serial(Port[-1][0], BaudRate)
    
    return(Arduino)


def GetData(FramesPerBuf, ArduinoObj):
    Data = np.zeros((FramesPerBuf, 2), dtype='float32')
    for F in range(FramesPerBuf):
        Data[F,0] = ArduinoObj.readline()
        Data[F,1] = round(time(), 3)
    
    return(Data)


def Run(FramesPerBuf, FileName='', Plot=False):
    """
    To read the written file in matlab:
        
        fileID = fopen(filename);
        A = fread(fileID, shape, 'float')

    The shape will be in the filename.
    """
    #if Plot: PlotThread()
    
    ArduinoObj = CreateObj(115200)
    Date = datetime.now().strftime("%Y%m%d%H%M%S")
    DataLen = 0
            Shape = 'x'.join([str(_) for _ in Data.shape])
    
    try:
        while True:
            Data = GetData(FramesPerBuf, ArduinoObj)
            with open(Date+'.dat', 'ab') as File: File.write(Data.tobytes())
            DataLen += Data.shape[0]
    
    except KeyboardInterrupt:
        pass

    os.rename(Date+'.dat', FileName+'_'+DataLen+'x2.dat')
