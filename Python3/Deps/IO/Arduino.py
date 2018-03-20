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
     

This is a script to define functions allowing Arduino/Python integration.
"""

import numpy as np
import os
import time

from datetime import datetime
from serial import Serial
from serial.tools.list_ports import comports


## Level 0
def CreateObj(BaudRate):
    Port = comports()
    Arduino = Serial(Port[-1][0], BaudRate)
    
    return(Arduino)


def GetSerialData(FramesPerBuf, ArduinoObj):
    Data = np.zeros((FramesPerBuf, 2), dtype='float32')

    for F in range(FramesPerBuf):
        Line = ArduinoObj.readline()
        while Line in [b'\r\n', b'\n']:
            Line = ArduinoObj.readline()

        Data[F,0] = float(Line)
        Data[F,1] = time.clock()
        time.sleep(0.001)

    return(Data)


## Level 1
def CheckPiezoAndTTL(BaudRate=115200, XLim=(0, 128), YLim=(-5, 1028), 
                     FramesPerBuf=128):
    
    import matplotlib.animation as animation
    import matplotlib.pyplot as plt
    
    Arduino = CreateObj(BaudRate)
    
    Fig = plt.figure()
    Ax = plt.axes(xlim=XLim, ylim=YLim)
    
    Plots = [[], []]
    Plots[0] = Ax.plot([float('nan')]*FramesPerBuf, lw=1)[0]
    Plots[1] = Ax.plot([float('nan')]*FramesPerBuf, lw=1)[0]
    
    def AnimInit():
        for Plot in Plots:
            Plot.set_ydata([])
        return Plots
    
    def PltUp(n):
        Data = [[0]*FramesPerBuf, [0]*FramesPerBuf]
        for Frame in range(FramesPerBuf):
            Temp = Arduino.readline().decode(); Temp = Temp.split()
            if len(Temp) is not 2:
                Temp = [0, 0]
            Data[0][Frame] = Temp[0]; Data[1][Frame] = Temp[1]
        
        for Index, Plot in enumerate(Plots):
            Plot.set_ydata(Data[Index])
        
        return tuple(Plots)
    
    Anim = animation.FuncAnimation(Fig, PltUp, frames=FramesPerBuf, interval=10, blit=False)


def WriteSerialData(FramesPerBuf, FileName='', Plot=False):
    """
    Grab serial data and continuously write to a .dat file. The shape will be 
    in the filename.
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
    return(None)


def Oscilloscope(BaudRate=115200, XLim=(0, 128), YLim=(-5, 1028), 
                       FramesPerBuf=128):
    
    import matplotlib.animation as animation
    import matplotlib.pyplot as plt
    
    Arduino = CreateObj(BaudRate)
    
    Fig = plt.figure()
    Ax = plt.axes(xlim=XLim, ylim=YLim)
    Plot = Ax.plot([float('nan')]*FramesPerBuf, lw=1)[0]
    
    def AnimInit():
        Data = []
        Plot.set_ydata(Data)
        return Plot,
    
    def PltUp(n):
        Data = []
        for Frame in range(FramesPerBuf):
            Data.append(Arduino.readline())
        Plot.set_ydata(Data)
        return Plot,
    
    Anim = animation.FuncAnimation(Fig, PltUp, frames=FramesPerBuf, interval=10, blit=False)


