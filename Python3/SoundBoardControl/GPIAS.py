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

#%% Set parameters of the experiment

"""==========#==========#==========#=========="""
AnimalName = '1ODE'

Rate = 128000

## Fill all durations in SECONDS!

## Sound
# Background noise duration
SoundBackgroundDur = 2.3
# Gap duration
SoundGapDur = 0.04
# Duration of background noise before pulse
SoundBackgroundPrePulseDur = 0.1
# Pulse duration
SoundLoudPulseDur = 0.05
# Duration of background noise after pulse
SoundBackgroundAfterPulseDur = 0.51
# Duration of background noise between stimulation block
SoundBetweenStimDur = [10, 20]

# Background and pulse amplification factors for each frequency tested
SoundBackgroundAmpF = [0.03, 0.02, 0.015, 0.015]
#SoundPulseAmpF = [1, 0.9, 0.8, 0.7]    # Real pulse
SoundPulseAmpF = [0.05, 0.05, 0.05, 0.05]    # Small pulse, for habituation

# Freqs to test. If using one freq range, keep it in a list, [[like this]].
NoiseFrequency = [[8000, 10000], [10000, 12000], [12000, 14000], [14000, 16000]]

# Number of trials per freq. tested (1 trial = 1 stim w/ gap + 1 stim w/o gap)
NoOfTrials = 5

# TTLs Amplification factor. DO NOT CHANGE unless you know what you're doing.
TTLAmpF = 0
"""==========#==========#==========#=========="""

import array
import datetime
import ControlSoundBoard
import pyaudio
import random
import shelve

Date = datetime.datetime.now()
FileName = ''.join([Date.strftime("%Y%m%d%H%M%S"), '-', AnimalName])

## Prepare dict w/ experimental setup
DataInfo = dict((Name, eval(Name)) for Name in ['AnimalName', 'Rate', 
                                       'SoundBackgroundDur', 
                                       'SoundGapDur', 
                                       'SoundBackgroundPrePulseDur', 
                                       'SoundLoudPulseDur', 
                                       'SoundBackgroundAfterPulseDur', 
                                       'SoundBetweenStimDur', 
                                       'SoundBackgroundAmpF', 'SoundPulseAmpF', 
                                       'NoiseFrequency', 'NoOfTrials', 
                                       'TTLAmpF', 'FileName'])

print('Creating SoundBackground...')
SoundPulseDur = SoundBackgroundDur
SoundPulseNo = 1
SoundAmpF = SoundBackgroundAmpF
SoundBackground = ControlSoundBoard.GenSound(Rate, SoundPulseDur, SoundPulseNo, 
                                             SoundAmpF, NoiseFrequency, 
                                             TTLAmpF)[0]
del(SoundPulseDur, SoundPulseNo, SoundAmpF)
SoundBackground = ControlSoundBoard.ReduceStim(SoundBackground)


print('Creating SoundGap...')
SoundGap = [[], []]
SoundPulseDur = SoundGapDur
SoundPulseNo = 1
SoundAmpF = SoundBackgroundAmpF
SoundGap[0] = ControlSoundBoard.GenSound(Rate, SoundPulseDur,SoundPulseNo, 
                                        SoundAmpF, NoiseFrequency, 
                                        TTLAmpF=1)[0]
SoundGap[0] = ControlSoundBoard.ReduceStim(SoundGap[0])
SoundGap[1] = [0]*(round(Rate*SoundGapDur)*2)
SoundGap[1] = bytes(array.array('f',SoundGap[1]))
SoundGap[1] = [SoundGap[1]]*len(SoundGap[0])

del(SoundPulseDur, SoundPulseNo, SoundAmpF)


print('Creating SoundBackgroundPrePulse...')
SoundPulseDur = SoundBackgroundPrePulseDur
SoundPulseNo = 1
SoundAmpF = SoundBackgroundAmpF
SoundBackgroundPrePulse = ControlSoundBoard.GenSound(Rate, SoundPulseDur, 
                                                     SoundPulseNo, SoundAmpF, 
                                                     NoiseFrequency, TTLAmpF)[0]
del(SoundPulseDur, SoundPulseNo, SoundAmpF)
SoundBackgroundPrePulse = ControlSoundBoard.ReduceStim(SoundBackgroundPrePulse)


print('Creating SoundLoudPulse...')
SoundPulseDur = SoundLoudPulseDur
SoundPulseNo = 1
SoundAmpF = SoundPulseAmpF
SoundLoudPulse = ControlSoundBoard.GenSound(Rate, SoundPulseDur, SoundPulseNo, 
                                            SoundAmpF, NoiseFrequency, 
                                            TTLAmpF=1)[0]
del(SoundPulseDur, SoundPulseNo, SoundAmpF)
SoundLoudPulse = ControlSoundBoard.ReduceStim(SoundLoudPulse)


print('Creating SoundBackgroundAfterPulse...')
SoundPulseDur = SoundBackgroundAfterPulseDur
SoundPulseNo = 1
SoundAmpF = SoundBackgroundAmpF
SoundBackgroundAfterPulse = ControlSoundBoard.GenSound(Rate, SoundPulseDur, 
                                                     SoundPulseNo, SoundAmpF, 
                                                     NoiseFrequency, TTLAmpF)[0]
del(SoundPulseDur, SoundPulseNo, SoundAmpF)
SoundBackgroundAfterPulse = ControlSoundBoard.ReduceStim(SoundBackgroundAfterPulse)


print('Creating SoundBetweenStimDur...')
SBSUnitDur = 0.5
SoundPulseNo = round(SoundBetweenStimDur[1]/SBSUnitDur)
SoundAmpF = SoundBackgroundAmpF
SoundBetweenStim = ControlSoundBoard.GenSound(Rate, SoundPulseDur, 
                                                 SoundPulseNo, SoundAmpF, 
                                                 NoiseFrequency, TTLAmpF)[0]
del(SoundPulseNo, SoundAmpF)
SoundBetweenStim = ControlSoundBoard.ReduceStim(SoundBetweenStim)


print('Generating sound objects...')
p = pyaudio.PyAudio()
Stimulation = p.open(format=pyaudio.paFloat32,
                     channels=2,
                     rate=Rate,
                     input=False,
                     output=True)

q = pyaudio.PyAudio()
InOn = False
RecStop = False
def InCallBack(in_data, frame_count, time_info, status):
    if InOn:
        global SoundRec
        SoundRec[RealFreq][RealTrial].append(in_data)
        
    if RecStop:
        InFlag = pyaudio.paComplete
    else:
        InFlag = pyaudio.paContinue
    return(None, InFlag)

Reading = q.open(format=pyaudio.paFloat32,
                     channels=1,
                     rate=Rate,
                     input=True,
                     output=False,
                     stream_callback=InCallBack)


#%% Check sensor's signal
XLim = (0, Rate//10)
YLim = (-0.003, 0.003)
FramesPerBuf = 512
ControlSoundBoard.MicrOscilloscope(Rate, XLim, YLim)


#%% Run!!
print('Preallocating memory and pseudo-randomizing the experiment...')

Freqs = [In for In, El in enumerate(NoiseFrequency)]*NoOfTrials
random.shuffle(Freqs)

SoundRec = [[] for _ in range(len(NoiseFrequency))]
for _ in range(len(SoundRec)):
    SoundRec[_] = [[] for _ in range(NoOfTrials*2)]

FakeTTLs = [[] for _ in range(len(NoiseFrequency))]
for Freq in range(len(FakeTTLs)):
    FakeTTLs[Freq] = [[] for _ in range(NoOfTrials*2)]
    
    for Trial in range(len(FakeTTLs[Freq])):
        FakeTTLs[Freq][Trial] = [[], []]

FreqSlot = [0]*(len(NoiseFrequency)*NoOfTrials*2)
for FE in range(len(Freqs)):
    FreqSlot[FE*2] = (Freqs[0:FE+1].count(Freqs[FE])-1)*2
    FreqSlot[FE*2+1] = (Freqs[0:FE+1].count(Freqs[FE])-1)*2+1

# Play!!
for Freq in range(len(Freqs)):
    Trials = [0, 1]
    random.shuffle(Trials)
    
    for Trial in Trials:
        RealFreq = Freqs[Freq]; RealTrial = FreqSlot[Freq*2+Trial]
        SBSDur = random.randrange(SoundBetweenStimDur[0], SoundBetweenStimDur[1])
        NoOfPulses = round(SBSDur/SBSUnitDur)
        print('Playing ', str(NoiseFrequency[RealFreq]), ' trial ', str(Trial)) 
        
        for Pulse in range(NoOfPulses):
            Stimulation.write(SoundBetweenStim[RealFreq])
        
        Reading.start_stream()
        RecStop = False; InOn = True
        Stimulation.write(SoundBackground[RealFreq])
        Stimulation.write(SoundGap[Trial][RealFreq])
        Stimulation.write(SoundBackgroundPrePulse[RealFreq])
        FakeTTLs[RealFreq][RealTrial][0] = len(SoundRec[RealFreq][RealTrial])
        Stimulation.write(SoundLoudPulse[RealFreq])
        FakeTTLs[RealFreq][RealTrial][1] = len(SoundRec[RealFreq][RealTrial])
        Stimulation.write(SoundBackgroundAfterPulse[RealFreq])
        InOn = False; RecStop = True
        Reading.stop_stream()
print('Done.')

with shelve.open(FileName) as Shelve:
    Shelve['SoundRec'] = SoundRec
    Shelve['FakeTTLs'] = FakeTTLs
    Shelve['DataInfo'] = DataInfo

print('Data saved.')

#%% Analysis
#

#
#FileName = '20160214120627-1ODE'
#
#with shelve.open(FileName) as Shelve:
#    SoundRec = Shelve['SoundRec']
#    FakeTTLs = Shelve['FakeTTLs']
#    DataInfo = Shelve['DataInfo']
#
#for Key, Value in DataInfo.items():
#    exec(str(Key) + '=' + 'Value')
#del(Key, Value)


import array
import matplotlib.pyplot as plt
import numpy as np
import shelve
from scipy import signal


RecordingData = [0]*len(SoundRec)
SoundTTLs = [0]*len(SoundRec)
for Freq in range(len(NoiseFrequency)):
    RecordingData[Freq] = [0]*len(SoundRec[Freq])
    SoundTTLs[Freq] = [0]*len(SoundRec[Freq])
    
    for Trial in range(NoOfTrials*2):
        SoundTTLs[Freq][Trial] = [float('nan')]*((len(SoundRec[Freq][Trial])*
                                                  len(SoundRec[Freq][Trial][0]))
                                                  //4)
        SStart = ((FakeTTLs[Freq][Trial][0]*len(SoundRec[Freq][Trial][0]))//4)-1
        SEnd = ((FakeTTLs[Freq][Trial][1]*len(SoundRec[Freq][Trial][0]))//4)-1

        RecordingData[Freq][Trial] = array.array('f', 
                                                b''.join(SoundRec[Freq][Trial]))
        SoundTTLs[Freq][Trial][SStart:SEnd] = [max(RecordingData[Freq][Trial])
                                               *1.5]*len(range(SStart, SEnd))
        
        RecordingData[Freq][Trial] = RecordingData[Freq][Trial][
                                        SStart-(0.05*Rate):SStart+(0.190*Rate)]
        SoundTTLs[Freq][Trial] = SoundTTLs[Freq][Trial][SStart-(0.05*Rate):
                                                        SStart+(0.190*Rate)]
        
        RecordingData[Freq][Trial] = abs(signal.hilbert(
                                                RecordingData[Freq][Trial]))
        passband = [3/(Rate/2), 300/(Rate/2)]
        f2, f1 = signal.butter(1, passband, 'bandpass')
        RecordingData[Freq][Trial] = signal.filtfilt(f2, f1, 
                                         RecordingData[Freq][Trial], 
                                         padtype='odd', padlen=0)
    
    for Trial in range(NoOfTrials-1):
        RecordingData[Freq][0] = 
        RecordingData[Freq][1] =
