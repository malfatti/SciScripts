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
#%%
import array
import Hdf5F
import KwikAnalysis
import math
import numpy as np
import pyaudio
from queue import Queue, Empty
import random
import sounddevice as SD
from scipy import signal
import threading

SoundTTLVal = 0.6; LaserTTLVal = 0.3
#SBAmpFsFile = '/home/cerebro/Malfatti/Test/20170213142143-SBAmpFs.hdf5'
SBAmpFsFile = '/home/malfatti/Documents/PhD/Tests/20170214093602-SBAmpFs.hdf5'

## Lower-level functions
def dBToAmpF(Intensities, CalibrationFile):
    print('Converting dB to AmpF...')
    SoundIntensity = Hdf5F.LoadSoundMeasurement(CalibrationFile, 
                                                    'SoundIntensity')
    
    SoundAmpF = {Hz: [float(min(SoundIntensity[Hz].keys(), 
                                key=lambda i: abs(SoundIntensity[Hz][i]-dB))) 
                      for dB in Intensities]
                 for Hz in list(SoundIntensity)}
    
    return(SoundAmpF)


def GenNoise(Rate, SoundPulseDur):
    print('Generating noise...')
    Noise = np.random.uniform(-1, 1, size=round(Rate*SoundPulseDur))
    Noise[-1] = 0
    Noise = np.array(Noise, dtype=np.float32)
    
    return(Noise)


def GenSineWave(Rate, Freq, AmpF, SoundPulseDur):
    print('Generating sine wave...')
    Pulse = [np.sin(2 * np.pi * Freq * (_/Rate)) * AmpF
             for _ in range(round(Rate*SoundPulseDur))]
    Pulse[-1] = 0
    Pulse = np.array(Pulse, dtype=np.float32)
    
    return(Pulse)


def BandpassFilterSound(SoundPulse, Rate, NoiseFrequency):
    # Preallocating memory
    SoundPulseFiltered = {}
    
    print('Filtering sound: ', end='')
    for Freq in range(len(NoiseFrequency)):
        FKey= str(NoiseFrequency[Freq][0]) + '-' + str(NoiseFrequency[Freq][1])
#        if len(NoiseFrequency[Freq]) == 1:
#            print('Filtering sound: ', NoiseFrequency[Freq], '...')
#            passband = [(NoiseFrequency[Freq][0]-1)/(Rate/2), 
#                        (NoiseFrequency[Freq][0]+1)/(Rate/2)]
#            f2, f1 = signal.butter(4, passband, 'bandpass')
#        else:        
        print(FKey, end='...')
        passband = [_/(Rate/2) for _ in NoiseFrequency[Freq]]
        
        f2, f1 = signal.butter(4, passband, 'bandpass')
        SoundPulseFiltered[FKey] = signal.filtfilt(f2, f1, SoundPulse, 
                                                   padtype='odd', padlen=0)
        SoundPulseFiltered[FKey][-1] = 0
    print(end='\n')
    return(SoundPulseFiltered)


def ApplySoundAmpF(SoundPulseFiltered, Rate, SoundAmpF, NoiseFrequency, 
                   SBOutAmpF, SoundPrePauseDur=0, SoundPostPauseDur=0):
    print('Applying amplification factors...')
    # Preallocating memory
    SoundPrePause = np.zeros(round(Rate * SoundPrePauseDur), dtype=np.float32)
    SoundPostPause = np.zeros(round(Rate * SoundPostPauseDur), dtype=np.float32)
    SoundUnit = {}
    
    for FKey in SoundPulseFiltered:
#    for Freq in range(len(NoiseFrequency)):
        SoundUnit[FKey] = {}
        
        for AmpF in range(len(SoundAmpF[FKey])):
            AKey = str(SoundAmpF[FKey][AmpF])
            SoundUnit[FKey][AKey] = np.concatenate([SoundPrePause, 
                                                    SoundPulseFiltered[FKey],
                                                    SoundPostPause])
            SoundUnit[FKey][AKey] = (SoundUnit[FKey][AKey]
                                     * SoundAmpF[FKey][AmpF]) * SBOutAmpF
    return(SoundUnit)


def GenTTL(Rate, PulseDur, TTLAmpF, TTLVal, SoundBoard, SBOutAmpF, PrePauseDur=0, 
           PostPauseDur=0):
    
    print('Generating Sound TTL...')
    TTLPulse = np.concatenate([
                  np.array([TTLVal] * round(Rate*PulseDur), dtype=np.float32),
                  np.array([TTLVal*-1] * round(Rate*PulseDur), dtype=np.float32)
                  ])
#    TTLPulse = [TTLVal] * round(Rate*PulseDur) + \
#               [TTLVal*-1] * round(Rate*PulseDur)
    TTLPulse[-1] = 0
    
    if PrePauseDur == 0:
        if PostPauseDur == 0:
            TTLUnit = TTLPulse
        else:
            TTLPostPause = np.zeros(round((PostPauseDur-PulseDur) * Rate), 
                            dtype=np.float32)
            TTLUnit = np.concatenate([TTLPulse, TTLPostPause])
    else:
        TTLPrePause = np.zeros(round(PrePauseDur * Rate), dtype=np.float32)
        if PostPauseDur == 0:
            TTLUnit = np.concatenate([TTLPrePause, TTLPulse])
        else:
            TTLPostPause = np.zeros(round((PostPauseDur-PulseDur) * Rate), 
                            dtype=np.float32)
            TTLUnit = np.concatenate([TTLPrePause, TTLPulse, TTLPostPause])
            
#    TTLPrePause = np.zeros(round(PrePauseDur * Rate), dtype=np.float32)
#    TTLPostPause = np.zeros(round((PostPauseDur-PulseDur) * Rate), 
#                            dtype=np.float32)
    
#    if SoundPulseDur < 0.01:
#        SoundTTLPulse = [round(SoundTTLVal/SBOutAmpF, 3)] * round(SoundPulseDur * Rate)
#    else:
#        Middle = [0]*round((SoundPulseDur-0.01) * Rate)
#        Border = [round(SoundTTLVal/SBOutAmpF, 3)] * round(0.005 * Rate)
#        SoundTTLPulse = Border + Middle + Border
    
#    TTLUnit = np.concatenate([TTLPrePause, TTLPulse, TTLPostPause])
    TTLUnit = (TTLUnit * TTLAmpF) * SBOutAmpF
    
    return(TTLUnit)


def InterleaveChannels(Right, Left, SoundAmpF, NoiseFrequency):
    print('Interleaving channels...')
    if Right == [0]:
        List = [0]*(2*len(Left))
        for _ in range(len(Left)):
            List[_ *2] = 0
            List[_ *2+1] = Left[_]
    
    elif Left == [0]:
        List = [0]*(2*len(Right))
        for _ in range(len(Right)):
            List[_ *2] = Right[_]
            List[_ *2+1] = 0
    
    else:
        List = [0]*len(NoiseFrequency)
        
        for Freq in range(len(NoiseFrequency)):
            Key= str(NoiseFrequency[Freq][0]) + '-' + \
                 str(NoiseFrequency[Freq][1])
            List[Freq] = [0]*len(SoundAmpF[Key])
            
            for AmpF in range(len(SoundAmpF[Key])):
                List[Freq][AmpF] = [0]*(2*len(Right[Freq][AmpF]))
                
                for _ in range(len(Right[Freq][AmpF])):
                    List[Freq][AmpF][_ *2] = Right[Freq][AmpF][_]
                    List[Freq][AmpF][_ *2+1] = Left[_]
    
    return(List)


def ListToByteArray(List, SoundAmpF, NoiseFrequency):
    print('Converting list to byte array...')
    if SoundAmpF == [0]:
        ByteArray = array.array('f', List)
        ByteArray = bytes(ByteArray)
    
    else:
        ByteArray = [0]*len(NoiseFrequency)
        
        for Freq in range(len(NoiseFrequency)):
            Key= str(NoiseFrequency[Freq][0]) + '-' + \
                 str(NoiseFrequency[Freq][1])
            ByteArray[Freq] = [0]*len(SoundAmpF[Key])
            
            for AmpF in range(len(SoundAmpF[Key])):
                ByteArray[Freq][AmpF] = array.array('f', List[Freq][AmpF])
                ByteArray[Freq][AmpF] = bytes(ByteArray[Freq][AmpF])
    
    return(ByteArray)


def GenPause(Rate, PauseBetweenStimBlocksDur):
    print('Generating pause...')
    PauseBetweenStimBlocksList = [0] * \
                                 round(PauseBetweenStimBlocksDur * Rate) * \
                                 2
    PauseBetweenStimBlocks = array.array('f', PauseBetweenStimBlocksList)
    PauseBetweenStimBlocks = bytes(PauseBetweenStimBlocks)
    
    return(PauseBetweenStimBlocks)


## Higher-level functions
def SoundStim(Rate, SoundPulseDur, SoundAmpF, NoiseFrequency, 
              TTLAmpF, SoundBoard, SoundPrePauseDur=0, SoundPostPauseDur=0, 
              SoundPauseBetweenStimBlocksDur=0, TTLs=True):
    """ Generate sound pulses in one channel and TTLs in the other channel 
    (Check ControlArduinoWithSoundBoard.ino code)."""
    
    SBOutAmpF = Hdf5F.SoundCalibration(SBAmpFsFile, SoundBoard,
                                               'SBOutAmpF')
    
    SoundPulse = GenNoise(Rate, SoundPulseDur)
    print('   ', end='')
    SoundPulseFiltered = BandpassFilterSound(SoundPulse, Rate, NoiseFrequency)
    print('   ', end='')
    SoundUnit = ApplySoundAmpF(SoundPulseFiltered, Rate, SoundAmpF, 
                               NoiseFrequency, SBOutAmpF, SoundPrePauseDur, 
                               SoundPostPauseDur)
    if TTLs:
        SoundTTLUnit = GenTTL(Rate, SoundPulseDur, TTLAmpF, SoundTTLVal, 
                              SoundBoard, SBOutAmpF, SoundPrePauseDur, 
                              SoundPostPauseDur)
    else:
        SoundTTLUnit = np.zeros((SoundPulse.size))
    
    Sound = {}
    for FKey in SoundUnit:
        Sound[FKey] = {}
        
        for AKey in SoundUnit[FKey]:
            Sound[FKey][AKey] = np.vstack((SoundTTLUnit, SoundUnit[FKey][AKey]))
    
    SoundPauseBetweenStimBlocks = np.zeros([2, 
                                            round(SoundPauseBetweenStimBlocksDur 
                                                  * Rate)])
    
    print('Done generating sound stimulus.')
    return(Sound, SoundPauseBetweenStimBlocks)


def LaserStim(Rate, LaserPulseDur, LaserPulseNo, TTLAmpF, SoundBoard, 
              Complexity='AllBlocks', LaserPrePauseDur=0, LaserPostPauseDur=0, 
              LaserStimBlockNo=1, LaserPauseBetweenStimBlocksDur=0):
    """ Generate square pulses in one channel that works as TTL for laser 
    (Check ControlArduinoWithSoundBoard.ino code)."""
    
    SBOutAmpF = Hdf5F.SoundCalibration(SBAmpFsFile, SoundBoard,
                                               'SBOutAmpF')
    
    LaserUnit = GenTTL(Rate, LaserPulseDur, TTLAmpF, LaserTTLVal, SoundBoard, 
                       SBOutAmpF, LaserPrePauseDur, LaserPostPauseDur)
    LaserList = InterleaveChannels([0], LaserUnit, [0], [0])
    Laser = ListToByteArray(LaserList, [0], [0])
    LaserPauseBetweenStimBlocks = GenPause(Rate, LaserPauseBetweenStimBlocksDur)    
    
    Stimulation = GenAudioObj(Rate, 'out')
    PlayLaser = RunSound(Laser, LaserPauseBetweenStimBlocks, LaserPulseNo, 
                         Stimulation, [0], [0], Complexity, LaserStimBlockNo, 
                         Multithread=False)
    
    print('Done generating laser stimulus.')
    return(Laser, PlayLaser)


def SoundLaserStim(Rate, SoundPulseDur, SoundPulseNo, SoundAmpF, 
                   NoiseFrequency, LaserPulseDur, LaserPulseNo, TTLAmpF, 
                   SoundBoard, Complexity='AllPulses', SoundPrePauseDur=0, 
                   SoundPostPauseDur=0, SoundStimBlockNo=1, 
                   SoundPauseBetweenStimBlocksDur=0, LaserPrePauseDur=0, 
                   LaserPostPauseDur=0, LaserStimBlockNo=1, 
                   LaserPauseBetweenStimBlocksDur=0):
    """ Generate sound pulses in one channel and TTLs for sound and laser in 
    the other channel (Check ControlArduinoWithSoundBoard.ino code)."""
    
    SBOutAmpF = Hdf5F.SoundCalibration(SBAmpFsFile, SoundBoard,
                                               'SBOutAmpF')
    
    LaserUnit = GenTTL(Rate, LaserPulseDur, TTLAmpF, LaserTTLVal, SoundBoard, 
                       SBOutAmpF, LaserPrePauseDur, LaserPostPauseDur)
    
    SoundTTLUnit = GenTTL(Rate, SoundPulseDur, TTLAmpF, SoundTTLVal, 
                          SoundBoard, SBOutAmpF, SoundPrePauseDur, 
                          SoundPostPauseDur)
    
    print('Summing sound TTL and laser units...')
    SoundTTLLaserUnit = [LaserUnit[_]+SoundTTLUnit[_] \
                         for _ in range(len(SoundTTLUnit))]
    
    SoundPulse = GenNoise(Rate, SoundPulseDur)
    print('   ', end='')
    SoundPulseFiltered = BandpassFilterSound(SoundPulse, Rate, NoiseFrequency)
    print('   ', end='')
    SoundUnit = ApplySoundAmpF(SoundPulseFiltered, Rate, SoundAmpF, 
                               NoiseFrequency, SBOutAmpF, SoundPrePauseDur, 
                               SoundPostPauseDur)
    
    SoundLaserList = InterleaveChannels(SoundUnit, SoundTTLLaserUnit, 
                                        SoundAmpF, NoiseFrequency)
    SoundLaser = ListToByteArray(SoundLaserList, SoundAmpF, NoiseFrequency)
    SoundLaserPauseBetweenStimBlocks = GenPause(Rate, 
                                                SoundPauseBetweenStimBlocksDur)
    
    Stimulation = GenAudioObj(Rate, 'out')
    PlaySoundLaser = RunSound(SoundLaser, SoundLaserPauseBetweenStimBlocks, 
                              SoundPulseNo, Stimulation, SoundAmpF, 
                              NoiseFrequency, Complexity, SoundStimBlockNo, 
                              Multithread=False)
    
    print('Done generating sound and laser stimulus.')
    return(SoundLaser, PlaySoundLaser)

def GPIAS(FileList, CalibrationFile, GPIASTimeBeforeTTL, GPIASTimeAfterTTL, FilterLow, 
          FilterHigh, FilterOrder):
    """ Analyze GPIAS recorded using sound board input. """
    
    for File in FileList:
        with shelve.open(CalibrationFile) as Shelve:
            SoundRec = Shelve['SoundRec']
            FakeTTLs = Shelve['FakeTTLs']
            DataInfo = Shelve['DataInfo']
        
        with shelve.open(CalibrationFile) as Shelve:
            SBInAmpF = Shelve['SBInAmpF']
        
        NoOfSamplesBefore = int(round((GPIASTimeBeforeTTL*DataInfo['Rate'])*10**-3))
        NoOfSamplesAfter = int(round((GPIASTimeAfterTTL*DataInfo['Rate'])*10**-3))
        NoOfSamples = NoOfSamplesBefore + NoOfSamplesAfter
        XValues = (range(-NoOfSamplesBefore, 
                         NoOfSamples-NoOfSamplesBefore)/DataInfo['Rate'])*10**3
        
        print('Preallocate memory...')
        GPIAS = [[0] for _ in range(len(DataInfo['NoiseFrequency']))]
        for Hz in range(len(DataInfo['NoiseFrequency'])):
            GPIAS[Hz] = [[0] for _ in range(DataInfo['NoOfTrials']*2)]
        
        AllTTLs = [[0] for _ in range(len(DataInfo['NoiseFrequency']))]
        for Hz in range(len(DataInfo['NoiseFrequency'])):
            AllTTLs[Hz] = [[0] for _ in range(DataInfo['NoOfTrials']*2)]
        
        for Hz in range(len(SoundRec)):        
            for Trial in range(DataInfo['NoOfTrials']*2):
                SoundStart = ((FakeTTLs[Hz][Trial][0] * 
                              len(SoundRec[Hz][Trial][0]))//4)-1
                SoundEnd = ((FakeTTLs[Hz][Trial][1] * 
                            len(SoundRec[Hz][Trial][0]))//4)-1
                
                Start = int(SoundStart-NoOfSamplesBefore)
                End = int(SoundEnd+NoOfSamplesAfter)
        
                GPIAS[Hz][Trial] = array.array('f', b''.join(SoundRec[Hz][Trial]))
                GPIAS[Hz][Trial] = GPIAS[Hz][Trial][Start:End]
                GPIAS[Hz][Trial] = [_/SBInAmpF for _ in GPIAS[Hz][Trial]]
                GPIAS[Hz][Trial] = abs(signal.hilbert(GPIAS[Hz][Trial]))
                
                passband = [FilterLow/(DataInfo['Rate']/2), 
                            FilterHigh/(DataInfo['Rate']/2)]
                f2, f1 = signal.butter(FilterOrder, passband, 'bandpass')
                GPIAS[Hz][Trial] = signal.filtfilt(f2, f1, GPIAS[Hz][Trial], 
                                                 padtype='odd', padlen=0)
                
                AllTTLs[Hz][Trial] = [SoundStart, SoundEnd]
            
            gData = GPIAS[Hz][:]
            NoGapAll = [gData[_] for _ in range(len(gData)) if _%2 == 0]
            GapAll = [gData[_] for _ in range(len(gData)) if _%2 != 0]
            NoGapSum = list(map(sum, zip(*NoGapAll)))
            GapSum = list(map(sum, zip(*GapAll)))
            
            gData = [0, 0]
            gData[0] = [_/DataInfo['NoOfTrials'] for _ in NoGapSum]
            gData[1] = [_/DataInfo['NoOfTrials'] for _ in GapSum]
            gData[0] = signal.savgol_filter(gData[0], 5, 2, mode='nearest')
            gData[1] = signal.savgol_filter(gData[1], 5, 2, mode='nearest')
            GPIAS[Hz] = gData[:]
            
            
            tData = AllTTLs[Hz][:]
            TTLNoGapAll = [tData[_] for _ in range(len(tData)) if _%2 == 0]
            TTLGapAll = [tData[_] for _ in range(len(tData)) if _%2 != 0]
            TTLNoGapSum = list(map(sum, zip(*TTLNoGapAll)))
            TTLGapSum = list(map(sum, zip(*TTLGapAll)))
            
            tData = [0, 0]
            tData[0] = [int(round(_/DataInfo['NoOfTrials'])) for _ in TTLNoGapSum]
            tData[1] = [int(round(_/DataInfo['NoOfTrials'])) for _ in TTLGapSum]
            AllTTLs[Hz] = tData[:]
            del(NoGapAll, GapAll, NoGapSum, GapSum, gData, tData)
        
        print('Saving data to ' + DataInfo['FileName'])
        with shelve.open(DataInfo['FileName']) as Shelve:
            Shelve['GPIAS'] = GPIAS
            Shelve['AllTTLs'] = AllTTLs
            Shelve['XValues'] = XValues
        print('Done.')
    
    return(None)


def PlotGPIAS(FileList):
    Params = KwikAnalysis.SetPlot(Backend='TkAgg', Params=True)
    from matplotlib import rcParams; rcParams.update(Params)
    from matplotlib import pyplot as plt
    
    for File in FileList:
        print('Loading data from ', File, ' ...')
        with shelve.open(File[:-3]) as Shelve:
            GPIAS = Shelve['GPIAS']
            AllTTLs = Shelve['AllTTLs']
            XValues = Shelve['XValues']
            DataInfo = Shelve['DataInfo']
        
        print('Plotting...')
        for Freq in range(len(DataInfo['NoiseFrequency'])):
            FigTitle = str(DataInfo['NoiseFrequency'][Freq]) + '\ Hz'
            Line0Label = 'No\ Gap'; Line1Label = 'Gap'
            SpanLabel = 'Sound\ Pulse'
            XLabel = 'time\ [ms]'; YLabel = 'voltage\ [mV]'
            
            plt.figure(Freq)
            plt.plot(XValues, GPIAS[Freq][0], 
                     color='r', label='$'+Line0Label+'$')
            plt.plot(XValues, GPIAS[Freq][1], 
                     color='b', label='$'+Line1Label+'$')
            plt.axvspan(XValues[AllTTLs[Freq][0][0]], XValues[AllTTLs[Freq][0][1]], 
                        color='k', alpha=0.5, lw=0, label='$'+SpanLabel+'$')
#            plt.axvspan(XValues[AllTTLs[Freq][1][0]], XValues[AllTTLs[Freq][1][1]], 
#                        color='b', alpha=0.5, lw=0, label='Sound pulse (Gap)')
            plt.suptitle('$'+FigTitle+'$')
            plt.ylabel('$'+YLabel+'$'); plt.xlabel('$'+XLabel+'$')
            plt.legend(loc='lower right')
            plt.locator_params(tight=True)
            plt.axes().spines['right'].set_visible(False)
            plt.axes().spines['top'].set_visible(False)
            plt.axes().yaxis.set_ticks_position('left')
            plt.axes().xaxis.set_ticks_position('bottom')
            plt.savefig('Figs/' + File[:-3] + '-' + 
                        str(DataInfo['NoiseFrequency'][Freq][0]) + '_' + 
                        str(DataInfo['NoiseFrequency'][Freq][1]) + '.svg', 
                        format='svg')
        print('Done.')
    
    return(None)


def MicrOscilloscope1(SoundBoard, Rate, YLim, FreqBand, MicSens_VPa, FramesPerBuf=512, Rec=False):
    Params = {'backend': 'Qt5Agg'}
    from matplotlib import rcParams; rcParams.update(Params)
    import matplotlib.animation as animation
    from matplotlib import pyplot as plt
    
    SBInAmpF = Hdf5F.SoundCalibration(SBAmpFsFile, SoundBoard,
                                              'SBInAmpF')
    
    r = pyaudio.PyAudio()
    
    Plotting = r.open(format=pyaudio.paFloat32,
                         channels=1,
                         rate=Rate,
                         input=True,
                         output=False,
                         input_device_index=18,
                         frames_per_buffer=FramesPerBuf)
                         #stream_callback=InCallBack)
    
    Fig = plt.figure()
    Ax = plt.axes(xlim=FreqBand, ylim=YLim)
    Plot, = Ax.plot([float('nan')]*(Rate//10), lw=1)
    
    def AnimInit():
        Data = array.array('f', [])
        Plot.set_ydata(Data)
        return Plot,
    
    def PltUp(n):
#        Data = array.array('f', Plotting.read(Rate//10))
        Data = array.array('f', Plotting.read(Rate//10, 
                                              exception_on_overflow=False))
        Data = [_ * SBInAmpF for _ in Data]
        HWindow = signal.hanning(len(Data)//(Rate/1000))
        F, PxxSp = signal.welch(Data, Rate, HWindow, nperseg=len(HWindow), noverlap=0, 
                                scaling='density')
        
        Start = np.where(F > FreqBand[0])[0][0]-1
        End = np.where(F > FreqBand[1])[0][0]-1
        BinSize = F[1] - F[0]
        RMS = sum(PxxSp[Start:End] * BinSize)**0.5
        dB = 20*(math.log(RMS/MicSens_VPa, 10)) + 94
        print(dB, max(PxxSp))
        
        Plot.set_xdata(F)
        Plot.set_ydata(PxxSp)
        return Plot,
    
    Anim = animation.FuncAnimation(Fig, PltUp, frames=FramesPerBuf, interval=16, 
                                   blit=False)
    
    if Rec:
        Writers = animation.writers['ffmpeg']
        Writer = Writers(fps=15, metadata=dict(artist='Me'))
        Anim.save('MicrOscilloscope.mp4', writer=Writer)
    
    plt.show()
        
    return(None)


def MicrOscilloscope(Rate, XLim, YLim, SoundBoard, FramesPerBuffer=512, Rec=False):
    """ Read data from sound board input, record it to a video and plot it 
    until the windows is closed (with a delay). """
    import Hdf5F
    import sounddevice as SD
    from queue import Queue, Empty
    import numpy as np
    SBAmpFsFile = '/home/malfatti/Documents/PhD/Tests/20160712135926-SBAmpFs.hdf5'

    SoundBoard = 'Intel_oAnalog-iAnalog'
    Device = 'default'
    Rate = 192000
    Channels = [0]
    DownSample = 10
    Window = 200
    Interval = 30
    SoundQueue = Queue()
    XLim=[0, 2000]
    YLim=[-0.05, 0.05]
    
#    SD.default.device = 'system'
#    SD.default.samplerate = Rate
#    SD.default.channels = 2
    
    def audio_callback(indata, outdata, frames, time, status):
        """This is called (from a separate thread) for each audio block."""
        if status:
            print(status, flush=True)
        # Fancy indexing with mapping creates a (necessary!) copy:
        SoundQueue.put(indata[::DownSample, Channels])
    
#    Fig = plt.figure()
#    Ax = plt.axes(xlim=XLim, ylim=YLim)
#    Plot, = Ax.plot([float('nan')]*(Rate//10), lw=1)
    
#    def AnimInit():
#        Data = np.array([np.nan])
#        Plot.set_ydata(Data)
#        return Plot,
    
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
    
#        Data = array.array('f', Plotting.read(Rate//10))
#        Data = SD.rec(Rate//10)
#        Data = Data * SBInAmpF
#        Plot.set_ydata(Data)
#        return Plot,
    
    Params = {'backend': 'Qt5Agg'}
    from matplotlib import rcParams; rcParams.update(Params)
    from matplotlib.animation import FuncAnimation
    from matplotlib import pyplot as plt
    
    SBInAmpF = Hdf5F.SoundCalibration(SBAmpFsFile, SoundBoard, 'SBInAmpF')
    
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
    
#    SD.default.
    
    Stream = SD.Stream(device=Device, channels=max(Channels)+1, blocksize=0, samplerate=Rate, callback=audio_callback, never_drop_input=True)
    Anim = FuncAnimation(Fig, PltUp, interval=Interval, blit=True)
    
    with Stream:
        plt.show()
    
#    Anim = FuncAnimation(Fig, PltUp, frames=FramesPerBuffer, interval=16, 
#                                   blit=False)
#    
#    if Rec:
#        Writers = animation.writers['ffmpeg']
#        Writer = Writers(fps=15, metadata=dict(artist='Me'))
#        Anim.save('MicrOscilloscope.mp4', writer=Writer)
    
#    plt.show()
    
    return(None)


def SoundCalOut(Rate, Freq, WaveDur):
    """ Generate a sine wave from 1 to -1 """
    Pulse = GenSineWave(Rate, Freq, 1, WaveDur)
    
    SD.default.device = 'system'
    SD.default.samplerate = Rate
    SD.default.channels = 1
    SD.default.blocksize = 0    
    SD.default.latency = 'low'
    
#    Stream = SD.OutputStream()
#    
#    print('Playing... ', end='')
#    with Stream: 
#        for PulseNo in range(SoundPulseNo):
#            Stream.write(Pulse)
#   
    print('Playing...', end='')
    SD.play(Pulse, blocking=True);
    print('Done.')
    
    return(None)


def SoundCalIn(Rate, Freq, WaveDur, SBOutAmpF):
    """ Generate sine wave (1V to -1V) and read 1s of it. """
    Pulse = GenSineWave(Rate, Freq, SBOutAmpF, WaveDur)
    
    SD.default.device = 'system'
    SD.default.samplerate = Rate
    SD.default.channels = 1
    SD.default.blocksize = 0    
    SD.default.latency = 'low'
    
    print('Measuring... ', end='')
    Rec = SD.playrec(Pulse, blocking=True)
    print('Done.')
    return(Rec)

