# -*- coding: utf-8 -*-
"""
Created on Mon Jan 25 14:30:21 2016

@author: malfatti
"""
#%%
import array
import matplotlib.animation as animation
import matplotlib.pyplot as plt
import serial
import serial.tools.list_ports

Rate = 128000
XLim = (0, 12800)
YLim = (-0.03, 0.03)
FramesPerBuf = 512

Port = serial.tools.list_ports.comports()
Arduino = serial.Serial(Port[-1][0], 19200)

Fig = plt.figure()
Ax = plt.axes(xlim=XLim, ylim=YLim)
Plot, = Ax.plot([], lw=1)

def AnimInit():
    Data = []
    Plot.set_ydata(Data)
    return Plot,

def PltUp(n):
    Data = []
    for Frame in range(FramesPerBuf):
        Data.append(array.array('f', Arduino.readlines()))
    Plot.set_ydata(Data)
    return Plot,

Anim = animation.FuncAnimation(Fig, PltUp, frames=FramesPerBuf, interval=5.12, blit=False)