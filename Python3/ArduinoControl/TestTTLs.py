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

from IO import Arduino
import time

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
    ArduinoObj.write(b'd')
    time.sleep(3)
    ArduinoObj.write(b'w')
    time.sleep(2)
# ArduinoObj.write(b'w')

