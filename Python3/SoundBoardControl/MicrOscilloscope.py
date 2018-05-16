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

from DataAnalysis.DataAnalysis import SignalIntensity
from IO import Hdf5
from queue import Queue, Empty


SBAmpFsFile = os.environ['DATAPATH']+'/Tests/SoundMeasurements/SoundMeasurements.hdf5'

def MicrOscilloscope(SoundBoard, Device, Rate, Window, Interval, YLim, 
                     FreqBand):
    """ Read data from sound board input and plot it until the windows is 
        closed. """
    
    Params = {'backend': 'Qt5Agg'}
    from matplotlib import rcParams; rcParams.update(Params)
    from matplotlib.animation import FuncAnimation
    from matplotlib import pyplot as plt
    
    # SD.default.device = Device
    SD.default.samplerate = Rate
    SD.default.channels = 1
    SD.default.blocksize = 192
    

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
                Data = Data * SBInAmpF
                
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
                
                SI, SIPSD = SignalIntensity(Data, Rate, FreqBand, MicSens_VPa)#, WindowSize=Window)
                F = SIPSD[0]
                PxxSp = SIPSD[1]
                dB = SI['dB']
                
                
            except Empty:
                break
            Shift = len(Data)
            DataPlot = np.roll(DataPlot, -Shift, axis=0)
            DataPlot[-Shift:, 0] = Data
            Block = False
        
        for Col, Line in enumerate(Lines):
            print(dB, max(PxxSp))
            Line.set_xdata(F)
            Line.set_ydata(PxxSp)
        
        return(Lines)
    
    
    Fig, Ax = plt.subplots()
    Lines = Ax.plot(DataPlot, lw=3)
    
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
    
    with Stream: plt.show()
        
    return(None)

#%%
SoundBoard = 'Jack-IntelOut-IntelIn'
Device = 'default'
Rate = 192000
Window = 256
Interval = 50
Y = 10e-9
YLim = [-Y*0.1, Y]
FreqBand = [200, 20000]
MicSens_dB = -47.46
MicSens_VPa = 10**(MicSens_dB/20)

MicrOscilloscope(SoundBoard, Device, Rate, Window, Interval, YLim, FreqBand)


#%%
import sounddevice as SD
from DataAnalysis.DataAnalysis import SignalIntensity

Device = 'system'
Rate = 192000
MicSens_dB = -47.46
MicSens_VPa = 10**(MicSens_dB/20)
FreqBand = [4000, 20000]

SD.default.device = Device
SD.default.samplerate = Rate
SD.default.channels = 1
SD.default.blocksize = 512

Rec = SD.rec(frames=3*Rate, blocking=True, channels=1)
SI, SIPSD = SignalIntensity(Rec[:,0], Rate, FreqBand, MicSens_VPa)#, WindowSize=Window)
