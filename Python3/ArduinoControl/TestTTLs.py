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
"""

#CalibrationFile = '/home/cerebro/Malfatti/Data/Test/' + \
#                  '20160419093139-SoundMeasurement/' + \
#                  '20160419093139-SoundMeasurement.hdf5'

#SoundBoard = 'USBPre2_oAux-iAux'

from IO import Arduino
#import pyaudio
import time

#Sound = ControlSoundBoard.GenSound(128000, 0.003, 1, {'8000-10000': [0.04]}, 
#                                   [[8000, 10000]], 1, CalibrationFile, 
#                                   SoundBoard, SoundPostPauseDur=0.097)[0]
#
#Laser = ControlSoundBoard.GenLaser(128000, 0.01, 1, 1, CalibrationFile, 
#                                   SoundBoard, LaserPostPauseDur=0.09)[0]
#                                   
#SoundLaser = ControlSoundBoard.GenSoundLaser(128000, 0.003, 1, 
#                                             {'8000-10000': [0.04]}, 
#                                             [[8000, 10000]], 0.01, 1, 1, 
#                                             CalibrationFile, SoundBoard,
#                                             SoundPrePauseDur=0.004, 
#                                             SoundPostPauseDur=0.093,
#                                             LaserPostPauseDur=0.09,)[0]

#p = pyaudio.PyAudio()
#Stimulation = p.open(format=pyaudio.paFloat32,
#                     channels=2,
#                     rate=128000,
#                     input=False,
#                     output=True)

ArduinoObj = Arduino.CreateObj(115200)

#%% Test TTL connections
while True:
    ArduinoObj.write(b'A')
    time.sleep(0.08)
    ArduinoObj.write(b'B')
    time.sleep(0.08)
    ArduinoObj.write(b'C')
    time.sleep(0.08)
    ArduinoObj.write(b'D')
    time.sleep(0.08)
    ArduinoObj.write(b'E')
    time.sleep(0.08)
    ArduinoObj.write(b'F')
    time.sleep(0.08)
    ArduinoObj.write(b'G')
    time.sleep(0.08)
    ArduinoObj.write(b'P')
    time.sleep(3)

#%% Test Sound-Arduino TTLs
for _ in range(100):
    Stimulation.write(Sound[0][0])
time.sleep(3)

for _ in range(100):
    Stimulation.write(Laser)
time.sleep(3)

for _ in range(100):
    Stimulation.write(SoundLaser[0][0])
time.sleep(3)


#%% Test Intan RHA 120-145 40-70 180<
print('Testing serial TTLs...')
for _ in range(10):
    ArduinoObj.write(b'A')
    time.sleep(0.08)
    ArduinoObj.write(b'B')
    time.sleep(0.08)
    ArduinoObj.write(b'C')
    time.sleep(0.08)
    ArduinoObj.write(b'D')
    time.sleep(0.08)
    ArduinoObj.write(b'E')
    time.sleep(0.08)
    ArduinoObj.write(b'P')
    time.sleep(3)

#%% Test Sound-Arduino TTLs
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

#%% Test OpenEphys RecControl
i = 0
# ArduinoObj.write(b'd')
while True:
    i += 1; print(i)
    ArduinoObj.write(b'D')
    time.sleep(2)
    ArduinoObj.write(b'D')
    time.sleep(1)
# ArduinoObj.write(b'w')

# RecordEngine timeout
# JUCE Assertion failure in juce_MessageManager.cpp:194
# Write buffer overrun, resizing to90112

# JackAudioIODevice::errorCallback zombified - calling shutdown handler
# JackAudioIODevice::shutdown
# To try: nice -n -18 open-ephys