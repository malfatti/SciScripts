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

const int inPin = 2;
const int PinNo = 8;
const int Pins[PinNo] = {2, 4, 7, 8, 10, 11, 12, 13};
const int Delay = 40;

void setup() {
  Serial.begin(38400);
  analogReference(INTERNAL);

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
    inPinV = analogRead(inPin);
    if (inPinV > 10) {
      ch = -1;
    }
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
  
  inPinV = map(inPinV, 0, 1023, 0, 255);

  if (inPinV >= 0 && inPinV < 30) {
    digitalWrite(Pins[0], LOW);
    digitalWrite(Pins[1], LOW);
  }

  if (inPinV >= 50 && inPinV < 80) {
    digitalWrite(Pins[0], LOW);
    digitalWrite(Pins[1], HIGH);
  }

  if (inPinV >= 125 && inPinV < 135) {
    digitalWrite(Pins[0], HIGH);
    digitalWrite(Pins[1], LOW);
  }

  if (inPinV >= 200) {
    digitalWrite(Pins[0], HIGH);
    digitalWrite(Pins[1], HIGH);
  }

}
