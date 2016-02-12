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

Arduino will read an analog in and, depending on the voltage, it will turn on a
digital port, or both. This provides 5V TTLs even if you have analog small voltage.
*/
const int SoundAndLaserTTLIn = 3;
const int LaserOut =  13;
const int SoundTTLOut = 12;

void setup() {
  Serial.begin(19200);
  analogReference(INTERNAL);

  pinMode(SoundAndLaserTTLIn, INPUT);
  pinMode(LaserOut, OUTPUT);
  pinMode(SoundTTLOut, OUTPUT);

  digitalWrite(LaserOut, LOW);
  digitalWrite(SoundTTLOut, LOW);
}

void loop() {

  int SoundAndLaserTTLInV = analogRead(SoundAndLaserTTLIn);
  SoundAndLaserTTLInV = map(SoundAndLaserTTLInV, 0, 1023, 0, 255);

  if (SoundAndLaserTTLInV < 12) {
    digitalWrite(LaserOut, LOW);
    digitalWrite(SoundTTLOut, LOW);
  }

  if (SoundAndLaserTTLInV >= 14 && SoundAndLaserTTLInV < 30) {
    digitalWrite(LaserOut, HIGH);
  }

  if (SoundAndLaserTTLInV > 160 && SoundAndLaserTTLInV < 185) {
    digitalWrite(SoundTTLOut, HIGH);
  }

  if (SoundAndLaserTTLInV >= 253) {
    digitalWrite(LaserOut, HIGH);
    digitalWrite(SoundTTLOut, HIGH);
  }


}
