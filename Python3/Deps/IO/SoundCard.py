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

This is a script that defines functions allowing the use of a computer's sound
board as an analog I/O board.

"""
import numpy as np
import sounddevice as SD

from queue import Queue, Empty
from IO import Hdf5
from IO.SigGen import SineWave


def MicrOscilloscope(Rate, XLim, YLim, SoundBoard, FramesPerBuffer=512, Rec=False):
    """ Read data from sound board input, record it to a video and plot it 
    until the windows is closed (with a delay). """
    
    Params = {'backend': 'Qt5Agg'}
    from matplotlib import rcParams; rcParams.update(Params)
    from matplotlib.animation import FuncAnimation
    from matplotlib import pyplot as plt
    
#    SoundBoard = 'Intel_oAnalog-iAnalog'
#    Device = 'system'
#    Rate = 192000
    Channels = [0]
    DownSample = 10
    Window = 200
    Interval = 30
    SoundQueue = Queue()
    XLim=[0, 2000]
    YLim=[-0.05, 0.05]
    
    SD.default.device = 'system'
    SD.default.samplerate = Rate
#    SD.default.channels = 2
    
    def audio_callback(indata, outdata, frames, time, status):
        """This is called (from a separate thread) for each audio block."""
        if status: print(status, flush=True)
        SoundQueue.put(indata[::DownSample, Channels])
    
    
    def PltUp(n):
        global DataPlot, SBInAmpF
        Block = True
        
        while True:
            try:
                Data = SoundQueue.get(block=Block)
                Data = Data * SBInAmpF
            except Empty:
                break
            Shift = len(Data)
            DataPlot = np.roll(DataPlot, -Shift, axis=0)
            DataPlot[-Shift:, :] = Data
            Block = False
        
        for Col, Line in enumerate(Lines):
            Line.set_ydata(DataPlot[:, Col])
        
        return(Lines)
    
    
    SBInAmpF = Hdf5.SoundCalibration(SBAmpFsFile, SoundBoard, 'SBInAmpF')
    
    DataLength = int(Window * Rate / (1000 * DownSample))
    DataPlot = np.zeros((DataLength, len(Channels)))
    
    Fig, Ax = plt.subplots()
    Lines = Ax.plot(DataPlot)
    
    if len(Channels) > 1:
        Ax.legend(['Channel {}'.format(Channel) for Channel in Channels],
                  loc='lower left', ncol=len(Channels))
    
    Ax.axis((0, len(DataPlot), -1, 1))
    Ax.set_yticks([0])
    Ax.yaxis.grid(True)
    Ax.tick_params(bottom='off', top='off', labelbottom='off',
                   right='off', left='off', labelleft='off')
    Ax.set_xlim(XLim)
    Ax.set_ylim(YLim)
    Fig.tight_layout(pad=0)
    
    Stream = SD.Stream(channels=max(Channels)+1, blocksize=0, callback=audio_callback, never_drop_input=True)
    Anim = FuncAnimation(Fig, PltUp, interval=Interval, blit=False)
    
    with Stream:
        plt.show()
    
#    if Rec:
#        Writers = animation.writers['ffmpeg']
#        Writer = Writers(fps=15, metadata=dict(artist='Me'))
#        Anim.save('MicrOscilloscope.mp4', writer=Writer)
    
#    plt.show()
    
    return(None)


#def MicrOscilloscope1(SoundBoard, Rate, YLim, FreqBand, MicSens_VPa, FramesPerBuf=512, Rec=False):
#    import array
#    import Hdf5F
#    import math
#    import numpy as np
#    import pyaudio
#    from scipy import signal
#    
#    Params = {'backend': 'Qt5Agg'}
#    from matplotlib import rcParams; rcParams.update(Params)
#    import matplotlib.animation as animation
#    from matplotlib import pyplot as plt
#    
#    SBAmpFsFile = '/home/cerebro/Malfatti/Test/20170213142143-SBAmpFs.hdf5'
#    SBInAmpF = Hdf5F.SoundCalibration(SBAmpFsFile, SoundBoard,
#                                              'SBInAmpF')
#    
#    r = pyaudio.PyAudio()
#    
#    Plotting = r.open(format=pyaudio.paFloat32,
#                         channels=1,
#                         rate=Rate,
#                         input=True,
#                         output=False,
#                         input_device_index=17,
#                         frames_per_buffer=FramesPerBuf)
#                         #stream_callback=InCallBack)
#    
#    Fig = plt.figure()
#    Ax = plt.axes(xlim=FreqBand, ylim=YLim)
#    Plot, = Ax.plot([float('nan')]*(Rate//10), lw=1)
#    
#    def AnimInit():
#        Data = array.array('f', [])
#        Plot.set_ydata(Data)
#        return Plot,
#    
#    def PltUp(n):
##        Data = array.array('f', Plotting.read(Rate//10))
#        Data = array.array('f', Plotting.read(Rate//10, 
#                                              exception_on_overflow=False))
#        Data = [_ * SBInAmpF for _ in Data]
#        HWindow = signal.hanning(len(Data)//(Rate/1000))
#        F, PxxSp = signal.welch(Data, Rate, HWindow, nperseg=len(HWindow), noverlap=0, 
#                                scaling='density')
#        
#        Start = np.where(F > FreqBand[0])[0][0]-1
#        End = np.where(F > FreqBand[1])[0][0]-1
#        BinSize = F[1] - F[0]
#        RMS = sum(PxxSp[Start:End] * BinSize)**0.5
#        dB = 20*(math.log(RMS/MicSens_VPa, 10)) + 94
#        print(dB, max(PxxSp))
#        
#        Plot.set_xdata(F)
#        Plot.set_ydata(PxxSp)
#        return Plot,
#    
#    Anim = animation.FuncAnimation(Fig, PltUp, frames=FramesPerBuf, interval=16, 
#                                   blit=False)
#    
#    if Rec:
#        Writers = animation.writers['ffmpeg']
#        Writer = Writers(fps=15, metadata=dict(artist='Me'))
#        Anim.save('MicrOscilloscope.mp4', writer=Writer)
#    
#    plt.show()
#        
#    return(None)


def SoundCalOut(Rate, Freq, WaveDur, Ch=2):
    """ Generate a sine wave from 1 to -1 """
    Pulse = SineWave(Rate, Freq, 1, WaveDur)
    
    SD.default.device = 'system'
    SD.default.samplerate = Rate
    SD.default.channels = 1
    SD.default.blocksize = 384    
#    SD.default.latency = 'low'
    
#    Stream = SD.OutputStream()
#    
#    print('Playing... ', end='')
#    with Stream: 
#        for PulseNo in range(SoundPulseNo):
#            Stream.write(Pulse)
#   
    print('Playing...', end='')
    SD.play(Pulse, blocking=True, mapping=(Ch));
    print('Done.')
    
    return(None)


def SoundCalIn(Rate, Freq, WaveDur, SBOutAmpF, Ch=2):
    """ Generate sine wave (1V to -1V) and read 1s of it. """
    Pulse = SineWave(Rate, Freq, SBOutAmpF, WaveDur)
    
    SD.default.device = 'system'
    SD.default.samplerate = Rate
    SD.default.channels = 1
    SD.default.blocksize = 384
#    SD.default.latency = 'low'
    
    print('Measuring... ', end='')
    Rec = SD.playrec(Pulse, blocking=True, output_mapping=(Ch))
    print('Done.')
    return(Rec)

