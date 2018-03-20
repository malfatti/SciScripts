#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Feb  6 12:23:43 2018

@author: malfatti
"""
import os

from IO import Arduino
from datetime import datetime


def Run(FramesPerBuf, FileName='', Plot=False):
    """
    Grab serial data and continuously write to a .dat fileThe shape will be in the filename.
    """
    #if Plot: PlotThread()
    
    ArduinoObj = Arduino.CreateObj(115200)
    Date = datetime.now().strftime("%Y%m%d%H%M%S")
    DataLen = 0
    
    try:
        while True:
            Data = Arduino.GetSerialData(FramesPerBuf, ArduinoObj)
            with open(Date+'.dat', 'ab') as File: File.write(Data.tobytes())
            DataLen += Data.shape[0]
    
    except KeyboardInterrupt:
        pass

    os.rename(Date+'.dat', FileName+'_'+str(DataLen)+'x2.dat')

if __name__ == "__main__":
    Run(256, 'Test')
