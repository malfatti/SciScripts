#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: T. Malfatti <malfatti@disroot.org>
@date: 2016-10-13
@license: GNU GPLv3 <https://raw.githubusercontent.com/malfatti/SciScripts/master/LICENSE>
@homepage: https://github.com/Malfatti/SciScripts
"""
#%%
import os
import numpy as np
import sounddevice as SD

from DataAnalysis.DataAnalysis import PSD, SignalIntensity
from IO import Hdf5
from queue import Queue, Empty


SBAmpFsFile = os.environ['DATAPATH']+'/Tests/SoundMeasurements/SoundMeasurements.hdf5'

def MicrOscilloscope(SoundBoard, Device, Rate, Window, Interval, YLim, 
                     FreqBand):
    """ Read data from sound board input and plot it until the windows is 
        closed. """
    
    Params = {'backend': 'TKAgg'}
    from matplotlib import rcParams; rcParams.update(Params)
    from matplotlib.animation import FuncAnimation
    from matplotlib import pyplot as plt
    
    SBInAmpF = Hdf5.SoundCalibration(SBAmpFsFile, SoundBoard, 'SBInAmpF')
    SoundQueue = Queue()
    DataLength = int(Window * Rate / (1000))
    DataPlot = np.zeros((DataLength, 1))
    
    def audio_callback(indata, outdata, frames, time, status):
        """This is called (from a separate thread) for each audio block."""
        if status:
            print(status, flush=True)
        # Fancy indexing with mapping creates a (necessary!) copy:
        SoundQueue.put(indata[:, 0])
    
    
    def PltUp(n):
        nonlocal DataPlot, SBInAmpF, FreqBand
        Block = True
        
        while True:
            try:
                Data = SoundQueue.get(block=Block)
                Data = Data * SBInAmpF *0.01
                
# #                HWindow = signal.hanning(len(Data)//(Rate/1000))
#                 F, PxxSp = signal.welch(Data, Rate, nperseg=64, noverlap=0, 
#                                         scaling='density')
                
# #                Start = np.where(F > FreqBand[0])[0][0]-1
# #                End = np.where(F > FreqBand[1])[0][0]-1
#                 BinSize = F[1] - F[0]
#                 RMS = sum(PxxSp * BinSize)**0.5
#                 dB = 20*(math.log(RMS/MicSens_VPa, 10)) + 94
                
                # F, PxxSp = PSD(Data, Rate)
                # Range = (F > FreqBand[0])*(F < FreqBand[1])
                # BinSize = F[1] - F[0]
                # RMS = (sum(PxxSp[Range]) * BinSize)**0.5
                # dB = 20*(np.log10((RMS/MicSens_VPa)/0.00002))
                
                SI = SignalIntensity(Data, Rate, FreqBand, MicSens_VPa, WindowSize=Window)
                F = SI['PSD'][0]
                PxxSp = SI['PSD'][1]
                dB = SI['dB']
                
                
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
    
    Stream = SD.Stream(samplerate=Rate,callback=audio_callback, never_drop_input=False)
    Anim = FuncAnimation(Fig, PltUp, interval=Interval, blit=True)
    
    with Stream:
        plt.show()
        
    return(None)

#%%
SoundBoard = 'Jack-IntelOut-IntelIn'
Device = 'default'
Rate = 192000
Window = 200
Interval = 30
Y = 5*10e-11
YLim = [-Y*0.1, Y]
FreqBand = [200, 4000]
MicSens_dB = -47.46
MicSens_VPa = 10**(MicSens_dB/20)

SD.default.device = Device
SD.default.samplerate = Rate
SD.default.channels = 1
SD.default.blocksize = 384

MicrOscilloscope(SoundBoard, Device, Rate, Window, Interval, YLim, FreqBand)

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