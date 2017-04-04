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
the acoustic startle reflex (GPIAS).
"""

#%% Set parameters of the experiment

AnimalName = 'Test01'

CalibrationFile = '/home/cerebro/Malfatti/Test/' + 'SoundMeasurements/' + \
                  'SoundMeasurements.hdf5'
#CalibrationFile = '/home/malfatti/Documents/PhD/Tests/' + \
#                  '20160419093139-SoundMeasurement/' + \
#                  '20160419093139-SoundMeasurement.hdf5'

# Sound board used
Setup = 'GPIAS'
System = 'Jack-IntelOut-MackieIn-MackieOut-IntelIn'

Rate = 192000
BaudRate = 115200
SoundCh = 6; SoundTTLCh = 3; PiezoCh = 8

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
BackgroundIntensity = [60]
PulseIntensity = [105]

# Noise frequency. If using one freq., keep the list in a list, [[like this]].
# USE ONLY FREQUENCY BANDS THAT WERE CALIBRATED. To check the calibrated freqs, 
# just run the cell once and then list(SoundIntensity).
NoiseFrequency = [[8000, 10000], [12000, 14000]]
#NoiseFrequency = [[8000, 10000], [9000, 11000], [10000, 12000], 
#                  [12000, 14000], [14000, 16000], [8000, 16000]]

# Number of trials per freq. tested (1 trial = 1 stim w/ gap + 1 stim w/o gap)
NoOfTrials = 9

# TTLs Amplification factor. DO NOT CHANGE unless you know what you're doing.
TTLAmpF = 0.07
#==========#==========#==========#==========#

import datetime
import ControlArduino
import ControlSoundBoard
import h5py
import Hdf5F
import numpy as np
import random
import sounddevice as SD

SoundIntensity = Hdf5F.LoadSoundMeasurement(CalibrationFile, Setup+'/'+System, 
                                            'SoundIntensity')

Date = datetime.datetime.now()
FileName = ''.join([Date.strftime("%Y%m%d%H%M%S"), '-', AnimalName, 
                    '-GPIAS.hdf5'])

SoundBackgroundAmpF = {Hz: [float(min(SoundIntensity[Hz].keys(), 
                            key=lambda i: abs(SoundIntensity[Hz][i]-dB))) 
                            for dB in BackgroundIntensity] 
                       for Hz in list(SoundIntensity)}

SoundPulseAmpF = {Hz: [float(min(SoundIntensity[Hz].keys(), 
                       key=lambda i: abs(SoundIntensity[Hz][i]-dB))) 
                       for dB in PulseIntensity] 
                  for Hz in list(SoundIntensity)}

## Prepare dict w/ experimental setup
DataInfo = dict((Name, eval(Name)) for Name in ['AnimalName', 'Rate', 'BaudRate',
                                       'Date', 'SoundBackgroundDur', 
                                       'SoundGapDur', 
                                       'SoundBackgroundPrePulseDur', 
                                       'SoundLoudPulseDur', 
                                       'SoundBackgroundAfterPulseDur', 
                                       'SoundBetweenStimDur',  
                                       'NoiseFrequency', 'NoOfTrials', 
                                       'TTLAmpF', 'FileName'])

print('Creating SoundBackground...')
SoundPulseDur = SoundBackgroundDur
SoundAmpF = SoundBackgroundAmpF
SoundBackground = ControlSoundBoard.SoundStim(Rate, SoundPulseDur, SoundAmpF, 
                                              NoiseFrequency, TTLAmpF, 
                                              System, TTLs=False)[0]
del(SoundPulseDur, SoundAmpF)

print('Creating SoundGap...')
SoundGap = {}
SoundPulseDur = SoundGapDur
SoundAmpF = SoundBackgroundAmpF
SoundGap['NoGap'] = ControlSoundBoard.SoundStim(Rate, SoundPulseDur, SoundAmpF, 
                                          NoiseFrequency, TTLAmpF, 
                                          System, TTLs=False)[0]

SoundGap['Gap'] = {FKey: {AKey: np.zeros(SoundGap['NoGap'][FKey][AKey].shape) 
                          for AKey in SoundGap['NoGap'][FKey]} 
                   for FKey in SoundGap['NoGap']}

del(SoundPulseDur, SoundAmpF)


print('Creating SoundBackgroundPrePulse...')
SoundPulseDur = SoundBackgroundPrePulseDur
SoundAmpF = SoundBackgroundAmpF
SoundBackgroundPrePulse = ControlSoundBoard.SoundStim(Rate, SoundPulseDur, 
                                                      SoundAmpF, NoiseFrequency, 
                                                      TTLAmpF, 
                                                      System, TTLs=False)[0]
del(SoundPulseDur, SoundAmpF)

print('Creating SoundLoudPulse...')
SoundPulseDur = SoundLoudPulseDur
SoundAmpF = SoundPulseAmpF

SoundLoudPulse = ControlSoundBoard.SoundStim(Rate, SoundPulseDur, SoundAmpF, 
                                             NoiseFrequency, TTLAmpF, System)[0]
del(SoundPulseDur, SoundAmpF)

print('Creating SoundBackgroundAfterPulse...')
SoundPulseDur = SoundBackgroundAfterPulseDur
SoundAmpF = SoundBackgroundAmpF
SoundBackgroundAfterPulse = ControlSoundBoard.SoundStim(Rate, SoundPulseDur, 
                                                        SoundAmpF, 
                                                        NoiseFrequency, TTLAmpF, 
                                                        System, TTLs=False)[0]
del(SoundPulseDur, SoundAmpF)

print('Creating SoundBetweenStim...')
SoundPulseDur = 1; SoundAmpF = SoundBackgroundAmpF
SoundBetweenStim = ControlSoundBoard.SoundStim(Rate, SoundPulseDur, SoundAmpF, 
                                               NoiseFrequency, TTLAmpF, 
                                               System)[0]
del(SoundPulseDur, SoundAmpF)

# Set Arduino and audio objects
Arduino = ControlArduino.CreateObj(BaudRate)
SD.default.device = 'system'
SD.default.samplerate = Rate
SD.default.blocksize = 384
SD.default.channels = 2
Map = (1, 2)

#%% Run!!
print('Preallocating memory and pseudo-randomizing the experiment...')
FreqsStr = ['-'.join([str(a) for a in b]) for b in NoiseFrequency]
TrialsStr = ['NoGap', 'Gap']
Freqs = [In for In, El in enumerate(NoiseFrequency)]*NoOfTrials
random.shuffle(Freqs); 

FreqSlot = [[0] for _ in range(len(Freqs)*2)]
for FE in range(len(Freqs)):
    FreqSlot[FE*2] = (Freqs[0:FE+1].count(Freqs[FE])-1)*2
    FreqSlot[FE*2+1] = (Freqs[0:FE+1].count(Freqs[FE])-1)*2+1

FreqOrder = [[0]]

# Play 3 trials only startle
for Pre in range(3):
    RealFreq = FreqsStr[-1]
    ABGKey = str(SoundBackgroundAmpF[RealFreq][0])
    APulseKey = str(SoundPulseAmpF[RealFreq][0])
    
    SBSDur = random.randrange(SoundBetweenStimDur[0], SoundBetweenStimDur[1])
    for Pulse in range(SBSDur):
        SD.play(SoundBetweenStim[RealFreq][ABGKey], blocking=True, mapping=Map)
    
    Arduino.write(b'P')
    SD.play(SoundBackground[RealFreq][ABGKey], blocking=True, mapping=Map)
    SD.play(SoundGap['NoGap'][RealFreq][ABGKey], blocking=True, mapping=Map)
    SD.play(SoundBackgroundPrePulse[RealFreq][ABGKey], blocking=True, mapping=Map)
    SD.play(SoundLoudPulse[RealFreq][APulseKey], blocking=True, mapping=Map)
    SD.play(SoundBackgroundAfterPulse[RealFreq][ABGKey], blocking=True, mapping=Map)
    Arduino.write(b'P')

# Play the test trials
for Hz in range(len(Freqs)):
    Trials = [0, 1]
    random.shuffle(Trials)
    print(str(Hz+1), 'of', str(len(Freqs)))
    
    for Trial in Trials:
        RealFreq = FreqsStr[Freqs[Hz]]; RealTrial = FreqSlot[Hz*2+Trial]
        FreqOrder[len(FreqOrder)-1] = [Freqs[Hz], RealTrial]
        FreqOrder.append([0])
        ABGKey = str(SoundBackgroundAmpF[RealFreq][0])
        APulseKey = str(SoundPulseAmpF[RealFreq][0])
        
        print('Playing ', RealFreq, ' trial ', TrialsStr[Trial])
        SBSDur = random.randrange(SoundBetweenStimDur[0], SoundBetweenStimDur[1])
        for Pulse in range(SBSDur):
            SD.play(SoundBetweenStim[RealFreq][ABGKey], blocking=True, mapping=Map)
        
        Arduino.write(b'P')
        SD.play(SoundBackground[RealFreq][ABGKey], blocking=True, mapping=Map)
        SD.play(SoundGap[TrialsStr[Trial]][RealFreq][ABGKey], blocking=True, mapping=Map)
        SD.play(SoundBackgroundPrePulse[RealFreq][ABGKey], blocking=True, mapping=Map)
        SD.play(SoundLoudPulse[RealFreq][APulseKey], blocking=True, mapping=Map)
        SD.play(SoundBackgroundAfterPulse[RealFreq][ABGKey], blocking=True, mapping=Map)
        Arduino.write(b'P')
        
FreqOrder.remove([0])

# Play 3 trials only startle
for Post in range(3):
    RealFreq = FreqsStr[-1]
    ABGKey = str(SoundBackgroundAmpF[RealFreq][0])
    APulseKey = str(SoundPulseAmpF[RealFreq][0])
    
    SBSDur = random.randrange(SoundBetweenStimDur[0], SoundBetweenStimDur[1])
    for Pulse in range(SBSDur):
        SD.play(SoundBetweenStim[RealFreq][ABGKey], blocking=True, mapping=Map)
    
    Arduino.write(b'P')
    SD.play(SoundBackground[RealFreq][ABGKey], blocking=True, mapping=Map)
    SD.play(SoundGap['NoGap'][RealFreq][ABGKey], blocking=True, mapping=Map)
    SD.play(SoundBackgroundPrePulse[RealFreq][ABGKey], blocking=True, mapping=Map)
    SD.play(SoundLoudPulse[RealFreq][APulseKey], blocking=True, mapping=Map)
    SD.play(SoundBackgroundAfterPulse[RealFreq][ABGKey], blocking=True, mapping=Map)
    Arduino.write(b'P')

print('Done.')

with h5py.File(FileName) as h5:
    h5.create_group('DataInfo')
    for Key, Value in DataInfo.items():
        h5['DataInfo'].attrs[Key] = Value
    
    h5['DataInfo'].create_group('SoundBackgroundAmpF')
    h5['DataInfo'].create_group('SoundPulseAmpF')
    
    for Key, Value in SoundBackgroundAmpF.items():
        h5['DataInfo']['SoundBackgroundAmpF'][Key] = Value
    
    for Key, Value in SoundPulseAmpF.items():
        h5['DataInfo']['SoundPulseAmpF'][Key] = Value
    
    h5['DataInfo']['Freqs'] = Freqs
    h5['DataInfo']['FreqOrder'] = FreqOrder
    h5['DataInfo']['FreqSlot'] = FreqSlot

print('Data saved.')