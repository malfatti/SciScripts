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

This is a script to generate sound stimulation for gap-prepulse inhibition of 
the acoustic startle reflex (GPIAS). It also records data from a sensor plugged 
in the sound board input.
"""
#%%
import matplotlib.animation as animation
import matplotlib.pyplot as plt
import serial
import serial.tools.list_ports

BaudRate = 19200
XLim = (0, BaudRate//50)
YLim = (-10, 50)
FramesPerBuf = BaudRate//50

Port = serial.tools.list_ports.comports()
Arduino = serial.Serial(Port[-1][0], 19200)

Fig = plt.figure()
Ax = plt.axes(xlim=XLim, ylim=YLim)
Plot, = Ax.plot([float('nan')]*FramesPerBuf, lw=1)

def AnimInit():
    Data = []
    Plot.set_ydata(Data)
    return Plot,

def PltUp(n):
    Data = []
    for Frame in range(FramesPerBuf):
        Data.append(Arduino.read())
    Plot.set_ydata(Data)
    return Plot,

Anim = animation.FuncAnimation(Fig, PltUp, frames=FramesPerBuf, interval=10, blit=False)