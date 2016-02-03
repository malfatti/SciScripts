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

BaudRate = 19200
XLim = (0, BaudRate//50)
YLim = (-10, 100)
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

#%%

Port = serial.tools.list_ports.comports()
Arduino = serial.Serial(Port[-1][0], 19200)

# plot parameters
analogPlot = AnalogPlot(Port[0][-1], 100)

# set up animation
fig = plt.figure()
ax = plt.axes(xlim=(0, 100), ylim=(0, 1023))
a0, = ax.plot([], [])
a1, = ax.plot([], [])
anim = animation.FuncAnimation(fig, analogPlot.update, 
                               fargs=(a0, a1), 
                               interval=50)