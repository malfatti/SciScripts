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

AnimalName = 'GPIAZon_NaCln02'

CalibrationFile = '/home/cerebro/Malfatti/Test/' + 'SoundMeasurements/' + \
                  'SoundMeasurements.hdf5'
#CalibrationFile = '/home/malfatti/Documents/PhD/Tests/' + \
#                  '20160419093139-SoundMeasurement/' + \
#                  '20160419093139-SoundMeasurement.hdf5'

# Sound setup and system used
Setup = 'GPIAS'
System = 'Jack-IntelOut-MackieIn-MackieOut-IntelIn'

Rate = 192000
BaudRate = 115200
SoundCh = 3; TTLCh = 1; PiezoCh = [8]

# Number of trials per freq. tested (1 trial = 1 stim w/ gap + 1 stim w/o gap)
NoOfTrials = 9

## Sound durations IN SECONDS
SoundBackgroundDur = 2.3
SoundGapDur = 0.04
SoundBackgroundPrePulseDur = 0.1
SoundLoudPulseDur = 0.05
SoundBackgroundAfterPulseDur = 0.51
SoundBetweenStimDur = [10, 20]

# Background and pulse intensities in dB. Supports float :)
BackgroundIntensity = [60]
PulseIntensity = [105]

# Noise frequency. If using one freq., keep the list in a list, [[like this]].
# USE ONLY FREQUENCY BANDS THAT WERE CALIBRATED. To check the calibrated freqs, 
# just run the cell once and then list(SoundIntensity).
#NoiseFrequency = [[8000, 16000]]
NoiseFrequency = [[8000, 10000], [10000, 12000], 
                  [12000, 14000], [14000, 16000], [8000, 16000]]

# TTLs Amplification factor. DO NOT CHANGE unless you know what you're doing.
TTLAmpF = 0.7
#==========#==========#==========#==========#

from datetime import datetime
import ControlArduino
import ControlSoundBoard
import Hdf5F
import numpy as np
import random
import sounddevice as SD

Date = datetime.now().strftime("%Y%m%d%H%M%S")
FileName = ''.join([Date, '-', AnimalName, '-GPIAS.hdf5'])

SoundBackgroundAmpF = ControlSoundBoard.dBToAmpF(BackgroundIntensity, 
                                                 CalibrationFile, 
                                                 Setup+'/'+System)

SoundPulseAmpF = ControlSoundBoard.dBToAmpF(PulseIntensity, 
                                            CalibrationFile, 
                                            Setup+'/'+System)


## Prepare dict w/ experimental setup
DataInfo = dict((Name, eval(Name)) for Name in ['AnimalName', 'Rate', 'BaudRate',
                                       'Date', 'SoundCh', 'TTLCh', 'PiezoCh', 
                                       'SoundBackgroundDur', 'SoundGapDur', 
                                       'SoundBackgroundPrePulseDur', 
                                       'SoundLoudPulseDur', 
                                       'SoundBackgroundAfterPulseDur', 
                                       'SoundBetweenStimDur',  
                                       'NoiseFrequency', 'NoOfTrials', 
                                       'TTLAmpF', 'FileName'])

SoundBackground, SoundGap, SoundBackgroundPrePulse, SoundLoudPulse, \
SoundBackgroundAfterPulse, SoundBetweenStim = ControlSoundBoard.GPIASStim(
        Rate, SoundBackgroundDur, SoundGapDur, SoundBackgroundPrePulseDur, 
        SoundLoudPulseDur, SoundBackgroundAfterPulseDur, SoundBetweenStimDur, 
        SoundBackgroundAmpF, SoundPulseAmpF, TTLAmpF, NoiseFrequency, 
        System)

# Set Arduino and audio objects
Arduino = ControlArduino.CreateObj(BaudRate)
SD.default.device = 'system'
SD.default.samplerate = Rate
SD.default.blocksize = 384
SD.default.channels = 2
Stim = SD.OutputStream(dtype='float32')


#%% Run!!
print('Preallocating memory and pseudo-randomizing the experiment...')
FreqsStr = ['-'.join([str(a) for a in b]) for b in NoiseFrequency]
TrialsStr = ['NoGap', 'Gap']
Freqs = [In for In, El in enumerate(NoiseFrequency)]*NoOfTrials
np.random.shuffle(Freqs); 

FreqSlot = [[0] for _ in range(len(Freqs)*2)]
for FE in range(len(Freqs)):
    FreqSlot[FE*2] = (Freqs[0:FE+1].count(Freqs[FE])-1)*2
    FreqSlot[FE*2+1] = (Freqs[0:FE+1].count(Freqs[FE])-1)*2+1

FreqOrder = [[0]]

print("Running...")
Stim.start()
# Play the Pre-trials
for Pre in range(3):
    RealFreq = FreqsStr[-1]
    FreqOrder[len(FreqOrder)-1] = [-1, -1]
    FreqOrder.append([0])
    ABGKey = str(SoundBackgroundAmpF[RealFreq][0])
    APulseKey = str(SoundPulseAmpF[RealFreq][0])
    SBSDur = random.randrange(SoundBetweenStimDur[0], SoundBetweenStimDur[1])
    
    print('Playing ', RealFreq, ' Pre-trial')
    Stim.write(SoundBetweenStim[RealFreq][ABGKey][:SBSDur*Rate, :])        
    Arduino.write(b'd')
    Stim.write(SoundBackground[RealFreq][ABGKey])
    Stim.write(SoundGap['NoGap'][RealFreq][ABGKey])
    Stim.write(SoundBackgroundPrePulse[RealFreq][ABGKey])
    Stim.write(SoundLoudPulse[RealFreq][APulseKey])
    Stim.write(SoundBackgroundAfterPulse[RealFreq][ABGKey])
    Arduino.write(b'w')

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
        SBSDur = random.randrange(SoundBetweenStimDur[0], SoundBetweenStimDur[1])
        
        print('Playing ', RealFreq, ' trial ', TrialsStr[Trial])
        Stim.write(SoundBetweenStim[RealFreq][ABGKey][:SBSDur*Rate, :])        
        Arduino.write(b'd')
        Stim.write(SoundBackground[RealFreq][ABGKey])
        Stim.write(SoundGap[TrialsStr[Trial]][RealFreq][ABGKey])
        Stim.write(SoundBackgroundPrePulse[RealFreq][ABGKey])
        Stim.write(SoundLoudPulse[RealFreq][APulseKey])
        Stim.write(SoundBackgroundAfterPulse[RealFreq][ABGKey])
        Arduino.write(b'w')

# Play the Pre-trials
for Pre in range(3):
    RealFreq = FreqsStr[-1]
    FreqOrder[len(FreqOrder)-1] = [-1, -1]
    FreqOrder.append([0])
    ABGKey = str(SoundBackgroundAmpF[RealFreq][0])
    APulseKey = str(SoundPulseAmpF[RealFreq][0])
    SBSDur = random.randrange(SoundBetweenStimDur[0], SoundBetweenStimDur[1])
    
    print('Playing ', RealFreq, ' Post-trial')
    Stim.write(SoundBetweenStim[RealFreq][ABGKey][:SBSDur*Rate, :])        
    Arduino.write(b'd')
    Stim.write(SoundBackground[RealFreq][ABGKey])
    Stim.write(SoundGap['NoGap'][RealFreq][ABGKey])
    Stim.write(SoundBackgroundPrePulse[RealFreq][ABGKey])
    Stim.write(SoundLoudPulse[RealFreq][APulseKey])
    Stim.write(SoundBackgroundAfterPulse[RealFreq][ABGKey])
    Arduino.write(b'w')

Stim.stop()
FreqOrder.remove([0])

DataInfo['FreqOrder'] = FreqOrder
DataInfo['FreqSlot'] = FreqSlot
DataInfo['Freqs'] = Freqs

Hdf5F.WriteGPIASInfo(DataInfo, SoundBackgroundAmpF, SoundPulseAmpF, Freqs, 
                     FreqOrder, FreqSlot, FileName)