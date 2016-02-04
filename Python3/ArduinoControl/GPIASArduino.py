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
BaudRate = 19200

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
NoiseFrequency = [[8000, 10000], [10000, 12000]]#, [12000, 14000], [14000, 16000]]

# Number of trials per freq. tested (1 trial = 1 stim w/ gap + 1 stim w/o gap)
NoOfTrials = 2

# TTLs Amplification factor. DO NOT CHANGE unless you know what you're doing.
TTLAmpF = 1
"""==========#==========#==========#=========="""

import array
import ControlArduino
import ControlSoundBoard
import pyaudio
import random
from threading import Thread

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
                                        SoundAmpF, NoiseFrequency, TTLAmpF)[0]
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
                                            TTLAmpF)[0]
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
SoundPulseDur = 0.05
SoundPulseNo = round(SoundBetweenStimDur[1]/0.05)
SoundAmpF = SoundBackgroundAmpF
SoundBetweenStim = ControlSoundBoard.GenSound(Rate, SoundPulseDur, 
                                                 SoundPulseNo, SoundAmpF, 
                                                 NoiseFrequency, TTLAmpF)[0]
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
            PiezoRec[ReadFreq][RealTrial].append(ord(Arduino.read()))
    
    def pause(self):
        self.paused = True
    
    def resume(self):
        self.paused = False
    
    def stop(self):
        self.stopped = True


#%% Check sensor's signal
XLim = (0, BaudRate//50)
YLim = (-10, 50)
FramesPerBuf = BaudRate//50
ControlArduino.Arduinoscilloscope(BaudRate, XLim, YLim, FramesPerBuf)


#%% Run!!
print('Preallocating memory and pseudo-randomizing the experiment...')

Freqs = [In for In, El in enumerate(NoiseFrequency)]*NoOfTrials
random.shuffle(Freqs)

PiezoRec = [[] for _ in range(len(NoiseFrequency))]
for _ in range(len(PiezoRec)):
    PiezoRec[_] = [[] for _ in range(NoOfTrials*2)]

FakeTTLs = [[] for _ in range(len(NoiseFrequency))]
for _ in range(len(FakeTTLs)):
    FakeTTLs[_] = [[] for _ in range(NoOfTrials*2)]

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
        NoOfPulses = round(SBSDur/0.05)
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

ClArduino().start(); ClArduino().pause()
for Freq in range(len(NoiseFrequency)):
    for Trial in range(NoOfTrials):
        ClArduino().resume()
        Stimulation.write(SoundBackground[Freq])
        ClArduino().pause()
ClArduino().stop()

#%% Analysis
RecordingData = [0]*len(SoundRec)

for Freq in range(len(SoundRec)):   
    RecordingData[Freq] = [0]*len(SoundRec[Freq])
    
    for Trial in range(len(SoundRec[Freq])):
        RecordingData[Freq][Trial] = array.array('f', 
                                                b''.join(SoundRec[Freq][Trial]))
    
##        passband = [50/(Rate/2), 300/(Rate/2)]
#        f2, f1 = scipy.signal.butter(4, 300/(Rate/2), 'lowpass')
#        RecordingData[Freq][Trial] = scipy.signal.filtfilt(f2, f1, RecordingData[Freq][Trial], 
#                                                           padtype='odd', 
#                                                           padlen=0)
#        del(f1,f2)

SoundTTLs = [0]*len(SoundRec); SoundTTLs[0] = [0]*len(SoundRec[0])
for Freq in range(len(NoiseFrequency)):
    for Trial in range(NoOfTrials*2):
        SoundTTLs[Freq][Trial] = [float('nan')]*((len(SoundRec[Freq][Trial])*
                                                  len(SoundRec[Freq][Trial][0]))//4)
        SStart = ((FakeTTLs[Freq][Trial][0]*len(SoundRec[Freq][Trial][0]))//4)-1
        SEnd = ((FakeTTLs[Freq][Trial][1]*len(SoundRec[Freq][Trial][0]))//4)-1
        SoundTTLs[Freq][Trial][SStart:SEnd] = [1]*len(range(SStart, SEnd))
