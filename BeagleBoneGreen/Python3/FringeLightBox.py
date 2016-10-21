#!/usr/bin/env python3
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

This is a script to simulate the light box from Fringe. The EEG signal will be 
filtered and the highest the power in that frequency, more leds will be turned 
off.
"""

import Adafruit_BBIO.GPIO as GPIO
import numpy as np
from queue import Queue, Empty
import sounddevice as SD
from scipy import signal

Rate = 48000
FreqBand = [8, 15]
FilterOrder = 4
Window = 256
Interval = 8
ButtonPower = 'P8_46'
ButtonPin = 'P8_40'
ButtonLed = 'P8_18'
LedPins = ['P9_'+str(_) for _ in range(21,32, 2)] + \
          ['P8_'+str(_) for _ in range(26,33, 2)] + ['P8_17', 'P8_27']

print('Checking pins...')
if ButtonPin in LedPins: 
    print('ButtonPin cannot be one of the LedPins.')
    raise SystemExit(1)

if ButtonLed in LedPins: 
    print('ButtonLed cannot be one of the LedPins.')
    raise SystemExit(1)

if ButtonPin == ButtonLed: 
    print('ButtonPin cannot be ButtonLed.')
    raise SystemExit(1)

print('Defining functions...')
def FilterSignal(Signal, Rate, Frequency, FilterOrder=4, Type='bandpass'):
    if Type not in ['bandpass', 'lowpass', 'highpass']:
        print("Choose 'bandpass', 'lowpass' or 'highpass'.")
        
    elif len(Frequency) not in [1, 2]:
        print('Frequency must have 2 elements for bandpass; or 1 element for \
        lowpass or highpass.')
        
    else:
        passband = [_/(Rate/2) for _ in Frequency]
        f2, f1 = signal.butter(FilterOrder, passband, Type)
        Signal = signal.filtfilt(f2, f1, Signal, padtype='odd', padlen=0)
        
        return(Signal)

def audio_callback(indata, outdata, frames, time, status):
    """This is called (from a separate thread) for each audio block."""
    if status:
        print(status, flush=True)
    # Fancy indexing with mapping creates a (necessary!) copy:
    SoundQueue.put(indata[:, 0])


def ProcessData(Rate, FreqBand):
#    Block = True
#    Data = np.array([])
    
#    while True:
#        try:
#            Data = SoundQueue.get(block=Block)
#            Data = np.append(Data, SoundQueue.get(block=Block))
#            
#        except Empty:
#            break
#        
#        Shift = len(Data)
#        SoundData = np.roll(SoundData, -Shift, axis=0)
#        SoundData[-Shift:, 0] = Data
#        Block = False
    
#    HWindow = signal.hanning(len(Data)//(Rate/1000))
#    F, PxxSp = signal.welch(Data, Rate, HWindow, nperseg=len(HWindow), 
#                            noverlap=0, scaling='density')
#    
#    Start = np.where(F > FreqBand[0])[0][0]-1
#    End = np.where(F > FreqBand[1])[0][0]-1
#    BinSize = F[1] - F[0]
#    RMS = sum(PxxSp[Start:End] * BinSize)**0.5
#    Max = max(PxxSp[Start:End])
#
    Data = SoundQueue.get(block=True)
    Data = FilterSignal(Data, Rate, FreqBand, FilterOrder, 'bandpass')
    RMS = np.mean(abs(Data))
    
    return(RMS)


print('Setting audio...')
SD.default.samplerate = Rate
SD.default.channels = 1
SD.default.blocksize = 0

SoundQueue = Queue()
Stream = SD.Stream(callback=audio_callback, never_drop_input=True)

print('Setting pins...')
for Pin in LedPins: GPIO.setup(Pin, GPIO.OUT)
GPIO.setup(ButtonPin, GPIO.IN)
GPIO.setup(ButtonLed, GPIO.OUT)
GPIO.setup(ButtonPower, GPIO.OUT)


with Stream:
    while True:
        print('Loop started!')
        print('Setting pins...')
        for Pin in LedPins: GPIO.output(Pin, GPIO.LOW)
        GPIO.output(ButtonPower, GPIO.HIGH)
        GPIO.output(ButtonLed, GPIO.HIGH)
        
        print('Press the microswitch to start.')
        GPIO.wait_for_edge(ButtonPin, GPIO.RISING)
        
        print('Calibrating RMS...')
        RefRMS = np.zeros([100])
        for _ in range(100):
            print(_)
            RefRMS[_] = ProcessData(Rate, FreqBand)
        RefRMS = np.mean(RefRMS)
        Slot = int(RefRMS/(len(LedPins)))
        print(RefRMS, Slot)
        
        print('Starting game...')
        GPIO.output(ButtonLed, GPIO.LOW)
        GPIO.output(ButtonPower, GPIO.LOW)
        for Pin in LedPins: GPIO.output(Pin, GPIO.HIGH)
        np.random.shuffle(LedPins)
        
        while True:
            print('Started Acq. loop!')
            Break = False
            RMS = ProcessData(Rate, FreqBand)
            print(RMS)
            
            for Limit in range(Slot, RefRMS, Slot):
                LowLimit = Limit - Slot;
                
                if RMS <= Limit and RMS > LowLimit:
                    print('Turn off', Limit/Slot, 'leds')
                    for Pin in LedPins[:int(Limit/Slot)]:
                        GPIO.output(Pin, GPIO.LOW)
                    
                    for Pin in LedPins[int(Limit/Slot):]:
                        GPIO.output(Pin, GPIO.HIGH)
                
                if RMS > RefRMS-Slot:
                    print('Turn off all leds and restart')
                    Break = True
                    break
            
            if Break: break
        
#            for (int Led = 0; Led < Limit/Slot; Led++) {
#              digitalWrite(LedPins[Led], LOW)
#            }
#            if (Limit/Slot == LastPin-FirstPin) { AllOff = true; }
#          }
#    
#          if (DataMean > (AmpMean + LimitP) && DataMean < (AmpMean + Limit) {
#            for (int Led = 0; Led < Limit/Slot; Led++) {
#                digitalWrite(LedPins[Led], LOW)
#            }
#            if (Limit/Slot == LastPin-FirstPin) { AllOff = true; }
#          }
#        }
#        if RMS
