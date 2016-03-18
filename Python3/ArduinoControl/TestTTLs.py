# -*- coding: utf-8 -*-
"""
Created on Fri Mar 18 16:43:58 2016

@author: cerebro
"""

import ControlArduino
import time

Arduino = ControlArduino.CreateObj(38400)

#%% Test OpenEphys
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
    Arduino.write(b'F')
    time.sleep(0.08)
    Arduino.write(b'G')
    time.sleep(0.08)
    Arduino.write(b'P')
    time.sleep(3)
