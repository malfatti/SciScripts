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

Sumarizing the commands:
 H = Helton's stimulus
 t = Test

 L = Full blue laser stimulus (5 blocks)
 l = Single blue laser stimulus block
 S = Full sound TTL (5 blocks)
 s = Single sound TTL block
 B = Full blue laser stim + sound TTL (5 blocks)
 b = Single blue laser stim + sound TTL block

 G = Full green laser stim + sound TTL (5 blocks)
 g = Single green laser stim + sound TTL block

 P = 1s light pulse
 O = On
 o = Off
 */

const int LaserPin =  13;
const int LaserTTLPin =  12;
const int SoundTTLPin =  8;

void setup() {
  Serial.begin(19200);
  pinMode(LaserPin, OUTPUT);
  pinMode(LaserTTLPin, OUTPUT);
  pinMode(SoundTTLPin, OUTPUT);

  digitalWrite(LaserPin, LOW);
  digitalWrite(LaserTTLPin, LOW);
  digitalWrite(SoundTTLPin, LOW);

  char ch = 0;
}

void loop() {

  char ch;

  while (ch == 0) {
    ch = Serial.read();
  }

  if (ch == 't') {
      digitalWrite(LaserPin, HIGH);
      digitalWrite(LaserTTLPin, HIGH);
      delay(5);
      digitalWrite(LaserPin, LOW);
      digitalWrite(LaserTTLPin, LOW);
      delay(95);
  }

  if (ch == 'H') {
    for (int i = 0; i < 1200; i++) {
      digitalWrite(LaserPin, HIGH);
      digitalWrite(LaserTTLPin, HIGH);
      delay(15);
      digitalWrite(LaserPin, LOW);
      digitalWrite(LaserTTLPin, LOW);
      delay(235);
    }
  }

  if (ch == 'L') {
    for (int i = 0; i < 4; i++) {
      for (int i = 0; i < 200; i++) {
        digitalWrite(LaserPin, HIGH);
        digitalWrite(LaserTTLPin, HIGH);
        delay(10);
        digitalWrite(LaserPin, LOW);
        digitalWrite(LaserTTLPin, LOW);
        delay(90);
      }
      delay(5000);
    }
    for (int i = 0; i < 200; i++) {
      digitalWrite(LaserPin, HIGH);
      digitalWrite(LaserTTLPin, HIGH);
      delay(10);
      digitalWrite(LaserPin, LOW);
      digitalWrite(LaserTTLPin, LOW);
      delay(90);
    }
  }

  if (ch == 'l') {
    for (int i = 0; i < 200; i++) {
      digitalWrite(LaserPin, HIGH);
      digitalWrite(LaserTTLPin, HIGH);
      delay(10);
      digitalWrite(LaserPin, LOW);
      digitalWrite(LaserTTLPin, LOW);
      delay(90);
    }
  }

  if (ch == 'S') {
    for (int i = 0; i < 4; i++) {
      for (int i = 0; i < 200; i++) {
        delay(4);
        digitalWrite(SoundTTLPin, HIGH);
        delay(3);
        digitalWrite(SoundTTLPin, LOW);
        delay(93);
        delayMicroseconds(80);
      }
      delay(5005);
      delayMicroseconds(500);
    }
    for (int i = 0; i < 200; i++) {
      delay(4);
      digitalWrite(SoundTTLPin, HIGH);
      delay(3);
      digitalWrite(SoundTTLPin, LOW);
      delay(93);
      delayMicroseconds(80);
    }
  }

  if (ch == 's') {
    for (int i = 0; i < 200; i++) {
      delay(4);
      digitalWrite(SoundTTLPin, HIGH);
      delay(3);
      digitalWrite(SoundTTLPin, LOW);
      delay(93);
      delayMicroseconds(80);
    }
  }

  if (ch == 'B') {
    for (int i = 0; i < 4; i++) {
      for (int i = 0; i < 200; i++) {
        digitalWrite(LaserPin, HIGH);
        digitalWrite(LaserTTLPin, HIGH);
        delay(4);
        digitalWrite(SoundTTLPin, HIGH);
        delay(3);
        digitalWrite(SoundTTLPin, LOW);
        delay(3);
        digitalWrite(LaserPin, LOW);
        digitalWrite(LaserTTLPin, LOW);
        delay(90);
        delayMicroseconds(55);
      }
      delay(5005);
      delayMicroseconds(500);
    }
    for (int i = 0; i < 200; i++) {
      digitalWrite(LaserPin, HIGH);
      digitalWrite(LaserTTLPin, HIGH);
      delay(4);
      digitalWrite(SoundTTLPin, HIGH);
      delay(3);
      digitalWrite(SoundTTLPin, LOW);
      delay(3);
      digitalWrite(LaserPin, LOW);
      digitalWrite(LaserTTLPin, LOW);
    }
  }

  if (ch == 'b') {
    for (int i = 0; i < 200; i++) {
      digitalWrite(LaserPin, HIGH);
      digitalWrite(LaserTTLPin, HIGH);
      delay(4);
      digitalWrite(SoundTTLPin, HIGH);
      delay(3);
      digitalWrite(SoundTTLPin, LOW);
      delay(3);
      digitalWrite(LaserPin, LOW);
      digitalWrite(LaserTTLPin, LOW);
      delay(90);
    }
  }

  if (ch == 'G') {
    for (int i = 0; i < 4; i++) {
      digitalWrite(LaserPin, HIGH);
      digitalWrite(LaserTTLPin, HIGH);
      for (int i = 0; i < 200; i++) {
        delay(4);
        digitalWrite(SoundTTLPin, HIGH);
        delay(3);
        digitalWrite(SoundTTLPin, LOW);
        delay(93);
      }
      digitalWrite(LaserPin, LOW);
      digitalWrite(LaserTTLPin, LOW);
      delay(5000);
    }
    digitalWrite(LaserPin, HIGH);
    digitalWrite(LaserTTLPin, HIGH);
    for (int i = 0; i < 200; i++) {
      delay(4);
      digitalWrite(SoundTTLPin, HIGH);
      delay(3);
      digitalWrite(SoundTTLPin, LOW);
      delay(93);
    }
    digitalWrite(LaserPin, LOW);
    digitalWrite(LaserTTLPin, LOW);
  }

  if (ch == 'g') {
    digitalWrite(LaserPin, HIGH);
    digitalWrite(LaserTTLPin, HIGH);
    for (int i = 0; i < 200; i++) {
      delay(4);
      digitalWrite(SoundTTLPin, HIGH);
      delay(3);
      digitalWrite(SoundTTLPin, LOW);
      delay(93);
    }
    digitalWrite(LaserPin, LOW);
    digitalWrite(LaserTTLPin, LOW);
  }

  if (ch == 'P') {
    digitalWrite(LaserPin, HIGH);
    digitalWrite(LaserTTLPin, HIGH);
    delay(1000);
    digitalWrite(LaserPin, LOW);
    digitalWrite(LaserTTLPin, LOW);
  }

  if (ch == 'O') {
    digitalWrite(LaserPin, HIGH);
    digitalWrite(LaserTTLPin, HIGH);
  }

  if (ch == 'o') {
    digitalWrite(LaserPin, LOW);
    digitalWrite(LaserTTLPin, LOW);
  }

}
