# -*- coding: utf-8 -*-
"""
@author: T. Malfatti <malfatti@disroot.org>
@year: 2015
@license: GNU GPLv3 <https://raw.githubusercontent.com/malfatti/SciScripts/master/LICENSE>
@homepage: https://github.com/Malfatti/SciScripts

This script can be used to calibrate the sound board in and out amplification 
factor.

For output:
    Write sine waves of amp 1 to sound board output. We use this to get the 
    sound board amp factor, that is, how much does the sound board being used 
    is increasing or decreasing the signal written to it.
    
For input:
    Read signal from sound board input. You have to apply a known amplitude 
    signal (a sine wave, for example) to the sound board input, so you can 
    check if the voltage you applied is the voltage being read by the sound 
    board.

It is very important to set the volume of the soundboard to unit level 
(usually 0dB, which is 100% of the intensity) so you know that no kind of 
frequency filter is being applied.
"""
#%% Set calibration
Rate = 192000; Freq = 10000; WaveDur = 10
SoundSystem = 'Jack-IntelOut-IntelIn'

from IO import SoundCard
from IO.SigGen import SineWave

import os, h5py
import numpy as np

Folder = os.environ['DATAPATH']+'/Tests/SoundMeasurements'
FileName = Folder + '/' + 'SoundMeasurements.hdf5'

#%% Output
Pulse = SineWave(Rate, Freq, 1, WaveDur)
SoundCard.SoundCalOut(Pulse, Ch=2)
# SBOutAmpF is the generated signal divided by the measured signal
SBOutAmpF = 1/2

#%% Input
Repetitions = 4
Pulse = SineWave(Rate, Freq, SBOutAmpF, WaveDur)

SBInAmpF = np.zeros(Repetitions, dtype=np.float32)
for aa in range(Repetitions):
    Rec = SoundCard.SoundCalIn(Pulse, Ch=2)
    #SBInAmpF[aa] = (max(Rec)-(min(Rec)))/2
    SBInAmpF[aa] = (max(Rec[2000:])-(min(Rec[2000:])))/2
    print(SBInAmpF[aa])

# SBInAmpF is the real amplitude divided by the measured amplitude
SBInAmpF = 1/SBInAmpF.mean()

print('SBInAmpF = ', str(SBInAmpF))

#%% NoiseRMS
import sounddevice as SD
SD.default.device = 'system'
SD.default.samplerate = Rate
SD.default.channels = 1
SD.default.blocksize = 384
MicNoise = SD.rec(Rate*10, Rate, 1, blocking=True)


#%% Save
with h5py.File(FileName, 'w') as F:
    if SoundSystem not in F.keys(): F.create_group(SoundSystem)
    
    F[SoundSystem]['SBOutAmpF'] = SBOutAmpF
    F[SoundSystem]['SBInAmpF'] = SBInAmpF
    F[SoundSystem]['MicNoise'] = MicNoise

with h5py.File(FileName, 'r') as F:
    a = F[SoundSystem]['SBOutAmpF'][()]
    b = F[SoundSystem]['SBInAmpF'][()]
    c = F[SoundSystem]['MicNoise'][()]

