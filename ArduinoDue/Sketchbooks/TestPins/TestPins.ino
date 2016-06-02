/*
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

  All pulses (uppercase chars) are 10ms long.
  A = pulse on pin 2
  B = pulse on pin 4
  C = pulse on pin 7
  D = pulse on pin 8
  E = pulse on pin 10
  F = pulse on pin 11
  G = pulse on pin 12
  P = pulse on pin 13
*/

const int sRate = 40; // in KHz
const int sPulseDur = 3; // in ms
const int sPulsePause = 97; // in ms
const int PulseNo = 200;

const long Pi = 3.1415; // in ms
const int sWaveDur = 100; // in ms
const int sWavePause = 900; // in ms
const int sWaveFreq = 1000; // in Hz
const int WaveNo = 1;

const int dPinNo = 8;
const int Pins[dPinNo] = {38, 41, 42, 45, 46, 49, 50, 51};

const int aPinNo = 12;
const int inPins[aPinNo] = {A0, A1, A2, A3, A4, A5, \
                            A6, A7, A8, A9, A10, A11
                           };

const int oAnalog[2] = {DAC0, DAC1};

void setup() {
  Serial.begin(38400);

  analogReadResolution(12);
  analogWriteResolution(12);
  randomSeed(analogRead(A11));
  analogWrite(DAC0, 0);
  analogWrite(DAC1, 0);

  for (int Pin = 0; Pin < aPinNo; Pin++) {
    pinMode(inPins[Pin], INPUT);
  }

  for (int Pin = 0; Pin < dPinNo; Pin++) {
    pinMode(Pins[Pin], OUTPUT);
    digitalWrite(Pins[Pin], LOW);
  }

  char ch = 0;
  int inPinV = 0;
}

void loop() {
  char ch = 0;
  int inPinV = 0;

  while (ch == 0) {
    ch = Serial.read();
  }

  if (ch == 'D') {
    for (int Pin = 0; Pin < dPinNo; Pin++) {
      digitalWrite(Pins[Pin], HIGH);
      delay(20);
      digitalWrite(Pins[Pin], LOW);
      delay(80);
    }
  }

  if (ch == 'a') {
    for (int Pin = 0; Pin < aPinNo; Pin++) {
      Serial.println(analogRead(inPins[Pin]));
    }
  }

  if (ch == 'A') {
    int Len = sRate * sPulseDur;
    int Delay = 1000 / sRate;

    long WN[Len];
    for (int El = 0; El < Len; El++) {
      WN[El] = random(530);
    }

    for ( int Pulse = 0; Pulse < PulseNo; Pulse++) {
      for (int El = 0; El < Len; El++) {
        dacc_set_channel_selection(DACC_INTERFACE, 0);
        dacc_write_conversion_data(DACC_INTERFACE, WN[El]);
        dacc_set_channel_selection(DACC_INTERFACE, 1);
        dacc_write_conversion_data(DACC_INTERFACE, WN[El]);
        //analogWrite(DAC0, WN[El]);
        //analogWrite(DAC1, WN[El]);
        delayMicroseconds(Delay);
      }
      analogWrite(DAC0, 0);
      analogWrite(DAC1, 0);
      delay(97);
    }
  }

  if (ch == 'S') {
    int Len = sRate * sWaveDur;
    int Delay = 1000 / sRate;

    long SW[Len];
    for (int El = 0; El < Len; El++) {
      SW[El] = sin(2*Pi*sWaveFreq*(El/sRate));
    }

    for ( int Pulse = 0; Pulse < PulseNo; Pulse++) {
      for (int El = 0; El < Len; El++) {
//        dacc_set_channel_selection(DACC_INTERFACE, 0);
//        dacc_write_conversion_data(DACC_INTERFACE, SW[El]);
          Serial.println(SW[El]);
        //analogWrite(DAC0, WN[El]);
        //analogWrite(DAC1, WN[El]);
        delayMicroseconds(Delay);
      }
      analogWrite(DAC0, 0);
      analogWrite(DAC1, 0);
      delay(97);
    }
  }
}
