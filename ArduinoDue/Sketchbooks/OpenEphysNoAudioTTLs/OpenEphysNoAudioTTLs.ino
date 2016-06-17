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
  A = pulse on pin 53
  B = pulse on pin 50
  C = pulse on pin 49
  D = pulse on pin 46
  E = pulse on pin 45
  F = pulse on pin 42
  G = pulse on pin 41
  P = pulse on pin 38
*/

const int inPin = A0;
const int PinNo = 8;
//const int Pins[PinNo] = {53, 50, 49, 46, 45, 42, 41, 38};
const int Pins[PinNo] = {38, 41, 42, 45, 46, 49, 50, 53};
const int Delay = 25;

void setup() {
  Serial.begin(38400);

//  // Set free running mode on ADC7 (pin A0)
//  ADC->ADC_MR |= 0x80;
//  ADC->ADC_CR = 2;
//  ADC->ADC_CHER = 0x80;

  pinMode(inPin, INPUT);
  for (int Pin = 0; Pin < PinNo; Pin++) {
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
//    while ((ADC->ADC_ISR & 0x80) == 0); // wait for conversion
//    inPinV = ADC->ADC_CDR[7];
//    
//    if (inPinV > 150) {
//      ch = -1;
//    }
  }

  if (ch == 'A') {
    digitalWrite(Pins[0], HIGH);
    delay(Delay);
    digitalWrite(Pins[0], LOW);
  }

  if (ch == 'a') {
    digitalWrite(Pins[0], HIGH);
    while (ch != 'z') {
      ch = Serial.read();
    }
    digitalWrite(Pins[0], LOW);
  }


  if (ch == 'B') {
    digitalWrite(Pins[1], HIGH);
    delay(Delay);
    digitalWrite(Pins[1], LOW);
  }

  if (ch == 'b') {
    digitalWrite(Pins[1], HIGH);
    while (ch != 'y') {
      ch = Serial.read();
    }
    digitalWrite(Pins[1], LOW);
  }

  if (ch == 'C') {
    digitalWrite(Pins[2], HIGH);
    delay(Delay);
    digitalWrite(Pins[2], LOW);
  }

  if (ch == 'c') {
    digitalWrite(Pins[2], HIGH);
    while (ch != 'x') {
      ch = Serial.read();
    }
    digitalWrite(Pins[2], LOW);
  }

  if (ch == 'D') {
    digitalWrite(Pins[3], HIGH);
    delay(Delay);
    digitalWrite(Pins[3], LOW);
  }

  if (ch == 'd') {
    digitalWrite(Pins[3], HIGH);
    while (ch != 'w') {
      ch = Serial.read();
    }
    digitalWrite(Pins[3], LOW);
  }

  if (ch == 'E') {
    digitalWrite(Pins[4], HIGH);
    delay(Delay);
    digitalWrite(Pins[4], LOW);
  }

  if (ch == 'e') {
    digitalWrite(Pins[4], HIGH);
    while (ch != 'v') {
      ch = Serial.read();
    }
    digitalWrite(Pins[4], LOW);
  }

  if (ch == 'F') {
    digitalWrite(Pins[5], HIGH);
    delay(Delay);
    digitalWrite(Pins[5], LOW);
  }

  if (ch == 'f') {
    digitalWrite(Pins[5], HIGH);
    while (ch != 'u') {
      ch = Serial.read();
    }
    digitalWrite(Pins[5], LOW);
  }

  if (ch == 'G') {
    digitalWrite(Pins[6], HIGH);
    delay(Delay);
    digitalWrite(Pins[6], LOW);
  }

  if (ch == 'g') {
    digitalWrite(Pins[6], HIGH);
    while (ch != 't') {
      ch = Serial.read();
    }
    digitalWrite(Pins[6], LOW);
  }

  if (ch == 'P') {
    digitalWrite(Pins[7], HIGH);
    delay(Delay);
    digitalWrite(Pins[7], LOW);
  }

//  if (inPinV < 200) {
//    digitalWrite(Pins[0], LOW);
//    digitalWrite(Pins[1], LOW);
//  }
//
//  if (inPinV >= 215 && inPinV < 275) {
//    digitalWrite(Pins[0], LOW);
//    digitalWrite(Pins[1], HIGH);
//  }
//
//  if (inPinV >= 550 && inPinV < 570) {
//    digitalWrite(Pins[0], HIGH);
//    digitalWrite(Pins[1], LOW);
//  }
//
//  if (inPinV >= 800) {
//    digitalWrite(Pins[0], HIGH);
//    digitalWrite(Pins[1], HIGH);
//  }
//
}
