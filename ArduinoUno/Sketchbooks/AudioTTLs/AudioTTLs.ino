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

// defines for setting and clearing register bits
#ifndef cbi
#define cbi(sfr, bit) (_SFR_BYTE(sfr) &= ~_BV(bit))
#endif
#ifndef sbi
#define sbi(sfr, bit) (_SFR_BYTE(sfr) |= _BV(bit))
#endif

const int inSoundPin = 2;
const int outSoundPin = 2;
const int outLaserPin = 4;

void setup() {
  // set prescale to 16
  sbi(ADCSRA, ADPS2) ;
  cbi(ADCSRA, ADPS1) ;
  cbi(ADCSRA, ADPS0) ;

  Serial.begin(38400);
  analogReference(INTERNAL);

  pinMode(inSoundPin, INPUT);
  pinMode(outSoundPin, OUTPUT);
  pinMode(outLaserPin, OUTPUT);
  digitalWrite(outSoundPin, LOW);
  digitalWrite(outLaserPin, LOW);

  int inPinV = 0;
}

void loop() {
  cli();
  while (1) {
    int inPinV = analogRead(inSoundPin);

    if (inPinV >= 0 && inPinV < 30) {
      PORTD &= ~_BV (outSoundPin);
      PORTD &= ~_BV (outLaserPin);
      //    bitClear (PORTD, Pins[0]);
      //    bitClear (PORTD, Pins[1]);
      //    PORTD = B00000000;
      //    digitalWrite(Pins[0], LOW);
      //    digitalWrite(Pins[1], LOW);
    }

    if (inPinV >= 210 && inPinV < 180) {
      PORTD &= ~_BV (outSoundPin);
      PORTD |= _BV (outLaserPin);
      //    bitClear (PORTD, Pins[0]);
      //    bitSet (PORTD, Pins[1]);
      //    PORTD = B00010000;
      //    digitalWrite(Pins[0], LOW);
      //    digitalWrite(Pins[1], HIGH);
    }

    if (inPinV >= 420 && inPinV < 440) {
      //      PORTD |= _BV (outSoundPin);
      //      PORTD &= ~_BV (outLaserPin);
      //    bitSet (PORTD, Pins[0]);
      //    bitClear (PORTD, Pins[1]);
      PORTD = B00000100;
      //    digitalWrite(Pins[0], HIGH);
      //    digitalWrite(Pins[1], LOW);
    }

    if (inPinV >= 600) {
      PORTD |= _BV (outSoundPin);
      PORTD |= _BV (outLaserPin);
      //    bitSet (PORTD, Pins[0]);
      //    bitSet (PORTD, Pins[1]);
      //    PORTD = B00010100;
      //    digitalWrite(Pins[0], HIGH);
      //    digitalWrite(Pins[1], HIGH);
    }
  }
}
