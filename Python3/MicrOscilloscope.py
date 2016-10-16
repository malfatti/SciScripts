#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct 13 08:50:31 2016

@author: malfatti
"""
#%%
import Hdf5F
import sounddevice as SD
from queue import Queue, Empty
import numpy as np
from scipy import signal
import math

SBAmpFsFile = '/home/malfatti/Documents/PhD/Tests/20160712135926-SBAmpFs.hdf5'

SoundBoard = 'Intel_oAnalog-iAnalog'
Device = 'default'
Rate = 192000
DownSample = 10
Window = 200
Interval = 30
YLim = [0, 5.5**-14]
FreqBand = [20, 20000]
MicSens_dB = -47.46
MicSens_VPa = 10**(MicSens_dB/20)

def MicrOscilloscope(SoundBoard, Device, Rate, DownSample, Window, 
                     Interval, YLim, FreqBand):
    """ Read data from sound board input and plot it until the windows is 
        closed. """
    
    Params = {'backend': 'Qt5Agg'}
    from matplotlib import rcParams; rcParams.update(Params)
    from matplotlib.animation import FuncAnimation
    from matplotlib import pyplot as plt
    
    SBInAmpF = Hdf5F.SoundCalibration(SBAmpFsFile, SoundBoard, 'SBInAmpF')
    SoundQueue = Queue()
    DataLength = int(Window * Rate / (1000 * DownSample))
    DataPlot = np.zeros((DataLength, 1))
    
    def audio_callback(indata, outdata, frames, time, status):
        """This is called (from a separate thread) for each audio block."""
        if status:
            print(status, flush=True)
        # Fancy indexing with mapping creates a (necessary!) copy:
        SoundQueue.put(indata[::DownSample, 0])
    
    
    def PltUp(n):
        nonlocal DataPlot, SBInAmpF, FreqBand
        Block = True
        
        while True:
            try:
                Data = SoundQueue.get(block=Block)
                Data = Data * SBInAmpF
                
#                HWindow = signal.hanning(len(Data)//(Rate/1000))
                F, PxxSp = signal.welch(Data, Rate, nperseg=64, noverlap=0, 
                                        scaling='density')
                
#                Start = np.where(F > FreqBand[0])[0][0]-1
#                End = np.where(F > FreqBand[1])[0][0]-1
                BinSize = F[1] - F[0]
                RMS = sum(PxxSp * BinSize)**0.5
                dB = 20*(math.log(RMS/MicSens_VPa, 10)) + 94
                
                
            except Empty:
                break
#            Shift = len(Data)
#            DataPlot = np.roll(DataPlot, -Shift, axis=0)
#            DataPlot[-Shift:, 0] = Data
            Block = False
        
        for Col, Line in enumerate(Lines):
            print(dB, max(PxxSp))
            Line.set_xdata(F)
            Line.set_ydata(PxxSp)
        
        return(Lines)
    
    
    Fig, Ax = plt.subplots()
    Lines = Ax.plot(DataPlot)
    
#    if len(Channels) > 1:
#        Ax.legend(['Channel {}'.format(Channel) for Channel in Channels],
#                  loc='lower left', ncol=len(Channels))
#    
#    Ax.axis((0, len(DataPlot), -1, 1))
#    Ax.set_yticks([0])
    Ax.yaxis.grid(True)
#    Ax.tick_params(bottom='off', top='off', labelbottom='off',
#                   right='off', left='off', labelleft='off')
    Ax.set_xlim(FreqBand)
    Ax.set_ylim(YLim)
    Fig.tight_layout(pad=0)
    
    Stream = SD.Stream(device=Device, channels=1, blocksize=0, samplerate=Rate, callback=audio_callback, never_drop_input=True)
    Anim = FuncAnimation(Fig, PltUp, interval=Interval, blit=True)
    
    with Stream:
        plt.show()
        
    return(None)

#%%
MicrOscilloscope(SoundBoard, Device, Rate, DownSample, Window, Interval, YLim, FreqBand)

#%%
import ControlSoundBoard

SBAmpFsFile = '/home/malfatti/Documents/PhD/Tests/20161013174053-SBAmpFs.hdf5'

SoundBoard = 'Intel_oAnalog-iAnalog'
Rate = 48000
YLim = [0, 6**-14]
FramesPerBuf = 512
FreqBand = [20, 20000]
MicSens_dB = -47.46
MicSens_VPa = 10**(MicSens_dB/20)

ControlSoundBoard.MicrOscilloscope1(SoundBoard, Rate, YLim, FreqBand, MicSens_VPa, FramesPerBuf)