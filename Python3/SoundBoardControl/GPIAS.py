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

AnimalName = 'TestSetup01'
Rate = 128000

CalibrationFile = '/home/cerebro/Malfatti/Data/Test/' + \
                  '20160315153450-SoundMeasurement/SoundIntensity'
#CalibrationFile = '/home/malfatti/NotSynced/SoftwareTest/' + \
#                  'SoundMeasurements/20160125114052-SoundMeasurement/' + \
#                  'SoundIntensity'

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

# Background and pulse intensities in dB. Supports float :)
BackgroundIntensity = 60
PulseIntensity = 105

# Noise frequency. If using one freq., keep the list in a list, [[like this]].
# USE ONLY FREQUENCY BANDS THAT WERE CALIBRATED. To check the calibrated freqs, 
# just run the cell once and then list(SoundIntensity).
NoiseFrequency = [[8000, 10000], [10000, 12000], 
                  [12000, 14000], [14000, 16000]]

# Number of trials per freq. tested (1 trial = 1 stim w/ gap + 1 stim w/o gap)
NoOfTrials = 5

# TTLs Amplification factor. DO NOT CHANGE unless you know what you're doing.
TTLAmpF = 0
#==========#==========#==========#==========#

import array
import datetime
import ControlSoundBoard
import matplotlib.pyplot as plt
import pyaudio
import random
import shelve
from scipy import signal

with shelve.open(CalibrationFile) as Shelve:
    SoundIntensity = Shelve['SoundIntensity']
    SBInAmpF = Shelve['SBInAmpF']

Date = datetime.datetime.now()
FileName = ''.join([Date.strftime("%Y%m%d%H%M%S"), AnimalName, '-GPIAS'])

SoundBackgroundAmpF = [[float(min(SoundIntensity[Hz].keys(), 
                                 key=lambda i: abs(SoundIntensity[Hz][i] - 
                                                   BackgroundIntensity)))] 
                       for Hz in list(SoundIntensity)]

SoundPulseAmpF = [[float(min(SoundIntensity[Hz].keys(), 
                                 key=lambda i: abs(SoundIntensity[Hz][i] - 
                                                   BackgroundIntensity)))] 
                       for Hz in list(SoundIntensity)]

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
                                             TTLAmpF, CalibrationFile)[0]
del(SoundPulseDur, SoundPulseNo, SoundAmpF)

print('Creating SoundGap...')
SoundGap = [[], []]
SoundPulseDur = SoundGapDur
SoundPulseNo = 1
SoundAmpF = SoundBackgroundAmpF
SoundGap[0] = ControlSoundBoard.GenSound(Rate, SoundPulseDur,SoundPulseNo, 
                                        SoundAmpF, NoiseFrequency, 
                                        TTLAmpF, CalibrationFile)[0]
SoundGap[1] = [0, 0.6]*(round(Rate*SoundGapDur))
SoundGap[1][-1] = 0
SoundGap[1] = bytes(array.array('f',SoundGap[1]))
SoundGap[1] = [SoundGap[1]]*len(SoundGap[0])

del(SoundPulseDur, SoundPulseNo, SoundAmpF)


print('Creating SoundBackgroundPrePulse...')
SoundPulseDur = SoundBackgroundPrePulseDur
SoundPulseNo = 1
SoundAmpF = SoundBackgroundAmpF
SoundBackgroundPrePulse = ControlSoundBoard.GenSound(Rate, SoundPulseDur, 
                                                     SoundPulseNo, SoundAmpF, 
                                                     NoiseFrequency, TTLAmpF, 
                                                     CalibrationFile)[0]
del(SoundPulseDur, SoundPulseNo, SoundAmpF)

print('Creating SoundLoudPulse...')
SoundPulseDur = SoundLoudPulseDur
SoundPulseNo = 1
SoundAmpF = SoundPulseAmpF
SoundLoudPulse = ControlSoundBoard.GenSound(Rate, SoundPulseDur, SoundPulseNo, 
                                            SoundAmpF, NoiseFrequency, 
                                            TTLAmpF=1, CalibrationFile)[0]
del(SoundPulseDur, SoundPulseNo, SoundAmpF)

print('Creating SoundBackgroundAfterPulse...')
SoundPulseDur = SoundBackgroundAfterPulseDur
SoundPulseNo = 1
SoundAmpF = SoundBackgroundAmpF
SoundBackgroundAfterPulse = ControlSoundBoard.GenSound(Rate, SoundPulseDur, 
                                                     SoundPulseNo, SoundAmpF, 
                                                     NoiseFrequency, TTLAmpF, 
                                                     CalibrationFile)[0]
del(SoundPulseDur, SoundPulseNo, SoundAmpF)

print('Creating SoundBetweenStimDur...')
SBSUnitDur = 0.5; SoundPulseDur = SBSUnitDur
SoundPulseNo = round(SoundBetweenStimDur[1]/SBSUnitDur)
SoundAmpF = SoundBackgroundAmpF
SoundBetweenStim = ControlSoundBoard.GenSound(Rate, SoundPulseDur, 
                                                 SoundPulseNo, SoundAmpF, 
                                                 NoiseFrequency, TTLAmpF, 
                                                 CalibrationFile)[0]
del(SoundPulseDur, SoundPulseNo, SoundAmpF)

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
ControlSoundBoard.MicrOscilloscope(Rate, XLim, YLim, CalibrationFile)


#%% Run!!
print('Preallocating memory and pseudo-randomizing the experiment...')

Freqs = [In for In, El in enumerate(NoiseFrequency)]*NoOfTrials
random.shuffle(Freqs)

SoundRec = [[] for _ in range(len(NoiseFrequency))]
for _ in range(len(SoundRec)):
    SoundRec[_] = [[] for _ in range(NoOfTrials*2)]

FakeTTLs = [[] for _ in range(len(NoiseFrequency))]
for Hz in range(len(FakeTTLs)):
    FakeTTLs[Hz] = [[] for _ in range(NoOfTrials*2)]
    
    for Trial in range(len(FakeTTLs[Hz])):
        FakeTTLs[Hz][Trial] = [[], []]

FreqSlot = [0]*(len(NoiseFrequency)*NoOfTrials*2)
for FE in range(len(Freqs)):
    FreqSlot[FE*2] = (Freqs[0:FE+1].count(Freqs[FE])-1)*2
    FreqSlot[FE*2+1] = (Freqs[0:FE+1].count(Freqs[FE])-1)*2+1

# Play!!
for Hz in range(len(Freqs)):
    Trials = [0, 1]
    random.shuffle(Trials)
    
    for Trial in Trials:
        RealFreq = Freqs[Hz]; RealTrial = FreqSlot[Hz*2+Trial]
        SBSDur = random.randrange(SoundBetweenStimDur[0], SoundBetweenStimDur[1])
        NoOfPulses = round(SBSDur/SBSUnitDur)
        print('Playing ', str(NoiseFrequency[RealFreq]), ' trial ', str(Trial)) 
        
        for Pulse in range(NoOfPulses):
            Stimulation.write(SoundBetweenStim[RealFreq][0])
        
        Reading.start_stream()
        RecStop = False; InOn = True
        Stimulation.write(SoundBackground[RealFreq][0])
        Stimulation.write(SoundGap[Trial][RealFreq][0])
        Stimulation.write(SoundBackgroundPrePulse[RealFreq][0])
        FakeTTLs[RealFreq][RealTrial][0] = len(SoundRec[RealFreq][RealTrial])
        Stimulation.write(SoundLoudPulse[RealFreq][0])
        FakeTTLs[RealFreq][RealTrial][1] = len(SoundRec[RealFreq][RealTrial])
        Stimulation.write(SoundBackgroundAfterPulse[RealFreq][0])
        InOn = False; RecStop = True
        Reading.stop_stream()
print('Done.')

with shelve.open(FileName) as Shelve:
    Shelve['SoundRec'] = SoundRec
    Shelve['FakeTTLs'] = FakeTTLs
    Shelve['DataInfo'] = DataInfo

print('Data saved.')