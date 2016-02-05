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
the acoustic startle reflex (GPIAS) and record data from a sensor plugged in 
an Arduino board.
"""

#%% Set parameters of the experiment

"""==========#==========#==========#=========="""
Rate = 128000
BaudRate = 115200

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
SoundPulseAmpF = [1, 0.9, 0.8, 0.7]

# Freqs to test. If using one freq range, keep it in a list, [[like this]].
NoiseFrequency = [[8000, 10000]]#, [10000, 12000], [12000, 14000], [14000, 16000]]

# Number of trials per freq. tested (1 trial = 1 stim w/ gap + 1 stim w/o gap)
NoOfTrials = 1

# TTLs Amplification factor. DO NOT CHANGE unless you know what you're doing.
TTLAmpF = 0
"""==========#==========#==========#=========="""

import array
import datetime
import ControlArduino
import ControlSoundBoard
import matplotlib.pyplot as plt
import os
import pickle
import pyaudio
import random
from threading import Thread

Date = datetime.datetime.now()
Folder = ''.join([Date.strftime("%Y%m%d%H%M%S"), '-SoundMeasurement'])
os.makedirs(Folder)

## Prepare dict w/ experimental setup
DataInfo = [Rate, BaudRate, SoundBackgroundDur, SoundGapDur, 
            SoundBackgroundPrePulseDur, SoundLoudPulseDur, 
            SoundBackgroundAfterPulseDur, SoundBetweenStimDur, 
            SoundBackgroundAmpF, SoundPulseAmpF, NoiseFrequency, 
            NoiseFrequency, TTLAmpF, Folder]


print('Creating SoundBackground...')
#SoundPulseDur = 0.05
#SoundPulseNo = round(SoundBackgroundDur/0.05)
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
                                                     NoiseFrequency, 
                                                     TTLAmpF)[0]
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
                                                     NoiseFrequency, 
                                                     TTLAmpF)[0]
del(SoundPulseDur, SoundPulseNo, SoundAmpF)
SoundBackgroundAfterPulse = ControlSoundBoard.ReduceStim(SoundBackgroundAfterPulse)


print('Creating SoundBetweenStimDur...')
SoundPulseDur = 0.05
SoundPulseNo = round(SoundBetweenStimDur[1]/0.05)
SoundAmpF = SoundBackgroundAmpF
SoundBetweenStim = ControlSoundBoard.GenSound(Rate, SoundPulseDur, 
                                                 SoundPulseNo, SoundAmpF, 
                                                 NoiseFrequency, 
                                                 TTLAmpF)[0]
del(SoundPulseDur, SoundPulseNo, SoundAmpF)
SoundBetweenStim = ControlSoundBoard.ReduceStim(SoundBetweenStim)


print('Generating sound and arduino objects...')
p = pyaudio.PyAudio()
Stimulation = p.open(format=pyaudio.paFloat32,
                     channels=2,
                     rate=Rate,
                     input=False,
                     output=True)

Arduino = ControlArduino.CreateObj(BaudRate)

#class ClArduino(threading.Thread):
#    def __init__(self):
#        threading.Thread.__init__(self)
#        self.iterations = 0
#        self.daemon = True  # OK for main to exit even if instance is still running
#        self.paused = True  # start out paused
#        self.state = threading.Condition()
#        
#    def run(self):
#        while True:
#            PiezoRec[Freq][Trial].append(Arduino.read())
#    
#    def pause(self):
#        with self.state:
#            self.paused = True  # make self block and wait
#    
#    def resume(self):
#        with self.state:
#            self.paused = False
#            self.state.notify()  # unblock self if waiting
#    
#    def stop(self):
#        with self.state:
#            self.stopped = True

class ClArduino(Thread):
    def run(self):
        while True:
            Data = Arduino.readline(); Data = Data.split()
            if len(Data) is not 2:
                Data = [b'0', b'0']
            PiezoRec[RealFreq][RealTrial].append(Data[0])
            TTLs[RealFreq][RealTrial].append(Data[1])
    
    def pause(self):
        self.paused = True
    
    def resume(self):
        self.paused = False
    
    def stop(self):
        self.stopped = True


#%% Run!!
print('Preallocating memory and pseudo-randomizing the experiment...')

Freqs = [In for In, El in enumerate(NoiseFrequency)]*NoOfTrials
random.shuffle(Freqs)

PiezoRec = [[] for _ in range(len(NoiseFrequency))]
for _ in range(len(PiezoRec)):
    PiezoRec[_] = [[] for _ in range(NoOfTrials*2)]

TTLs = [[] for _ in range(len(NoiseFrequency))]
for _ in range(len(TTLs)):
    TTLs[_] = [[] for _ in range(NoOfTrials*2)]

FreqSlot = [0]*(len(NoiseFrequency)*NoOfTrials*2)
for FE in range(len(Freqs)):
    FreqSlot[FE*2] = (Freqs[0:FE+1].count(Freqs[FE])-1)*2
    FreqSlot[FE*2+1] = (Freqs[0:FE+1].count(Freqs[FE])-1)*2+1

# Play!!
ClArduino().start(); ClArduino().pause()
for Freq in range(len(Freqs)):
    Trials = [0, 1]
    random.shuffle(Trials)
    
    for Trial in Trials:
        RealFreq = Freqs[Freq]; RealTrial = FreqSlot[Freq*2+Trial]
        SBSDur = random.randrange(SoundBetweenStimDur[0], SoundBetweenStimDur[1])
        NoOfPulses = round(SBSDur/0.05)
        print('Playing ', str(NoiseFrequency[RealFreq]), ' trial ', str(Trial)) 
        
        for Pulse in range(NoOfPulses):
            Stimulation.write(SoundBetweenStim[RealFreq])
        
        Arduino.flush()
        ClArduino().resume()
        Stimulation.write(SoundBackground[RealFreq])
        Stimulation.write(SoundGap[Trial][RealFreq])
        Stimulation.write(SoundBackgroundPrePulse[RealFreq])
        Stimulation.write(SoundLoudPulse[RealFreq])
        Stimulation.write(SoundBackgroundAfterPulse[RealFreq])
        ClArduino().pause()
ClArduino().stop()
print('Done.')


File = open(Folder+'/'+'PiezoRec.pckl', 'wb')
pickle.dump(PiezoRec, File)
File.close()
del(File)

File = open(Folder+'/'+'TTLs.pckl', 'wb')
pickle.dump(TTLs, File)
File.close()
del(File)

File = open(Folder+'/'+'DataInfo.pckl', 'wb')
pickle.dump(DataInfo, File)
File.close()
del(File)
print('Data saved.')


#%% Analysis
RecordingData = [0]*len(PiezoRec)

for Freq in range(len(PiezoRec)):   
    RecordingData[Freq] = [0]*len(PiezoRec[Freq])
    
    for Trial in range(len(PiezoRec[Freq])):
        RecordingData[Freq][Trial] = [ord(PiezoRec[Freq][Trial][_]) for _ in range(len(PiezoRec[Freq][Trial]))]
    