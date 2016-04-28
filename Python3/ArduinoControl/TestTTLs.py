# -*- coding: utf-8 -*-
"""
Created on Fri Mar 18 16:43:58 2016

@author: cerebro
"""

CalibrationFile = '/home/cerebro/Malfatti/Data/Test/' + \
                  '20160419093139-SoundMeasurement/' + \
                  '20160419093139-SoundMeasurement.hdf5'

SoundBoard = 'USBPre2_oAux-iAux'

import ControlArduino
import ControlSoundBoard
import pyaudio
import time

Sound = ControlSoundBoard.GenSound(128000, 0.003, 1, {'8000-10000': [0.004]}, 
                                   [[8000, 10000]], 1, CalibrationFile, 
                                   SoundBoard, SoundPostPauseDur=0.097)[0]

Laser = ControlSoundBoard.GenLaser(128000, 0.01, 1, 1, CalibrationFile, 
                                   SoundBoard, LaserPostPauseDur=0.09)[0]
                                   
SoundLaser = ControlSoundBoard.GenSoundLaser(128000, 0.003, 1, 
                                             {'8000-10000': [0.004]}, 
                                             [[8000, 10000]], 0.01, 1, 1, 
                                             CalibrationFile, SoundBoard,
                                             SoundPrePauseDur=0.004, 
                                             SoundPostPauseDur=0.093,
                                             LaserPostPauseDur=0.09,)[0]

p = pyaudio.PyAudio()
Stimulation = p.open(format=pyaudio.paFloat32,
                     channels=2,
                     rate=128000,
                     input=False,
                     output=True)

Arduino = ControlArduino.CreateObj(38400)

#%% Test OpenEphys
for _ in range(1000):
    Arduino.write(b'A')
    time.sleep(0.08)
    Arduino.write(b'B')
    time.sleep(0.08)
    Arduino.write(b'C')
    time.sleep(0.08)
    Arduino.write(b'D')
    time.sleep(0.08)
    Arduino.write(b'E')
    time.sleep(0.08)
    Arduino.write(b'F')
    time.sleep(0.08)
    Arduino.write(b'G')
    time.sleep(0.08)
    Arduino.write(b'P')
    time.sleep(3)

for _ in range(100):
    Stimulation.write(Sound)
time.sleep(3)

for _ in range(100):
    Stimulation.write(Laser)
time.sleep(3)

for _ in range(100):
    Stimulation.write(SoundLaser)
time.sleep(3)


#%% Test Intan RHA 120-145 40-70 180<
print('Testing serial TTLs...')
for _ in range(10):
    Arduino.write(b'A')
    time.sleep(0.08)
    Arduino.write(b'B')
    time.sleep(0.08)
    Arduino.write(b'C')
    time.sleep(0.08)
    Arduino.write(b'D')
    time.sleep(0.08)
    Arduino.write(b'E')
    time.sleep(0.08)
    Arduino.write(b'P')
    time.sleep(3)

print('Testing sound TTLs...')
for _ in range(100):
    Stimulation.write(Sound[0][0])
time.sleep(3)

for _ in range(100):
    Stimulation.write(Laser)
time.sleep(3)

for _ in range(100):
    Stimulation.write(SoundLaser[0][0])
time.sleep(3)
