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

All pulses (uppercase chars) are 5ms long. 
 A = pulse on pin 2
 B = pulse on pin 4
 C = pulse on pin 7
 D = pulse on pin 8
 E = pulse on pin 12
 P = pulse on pin 13
*/

const int Pins[6] = {2, 4, 7, 8, 12, 13};
const int PinNo = 6;
const int Delay = 5;

void setup() {
  Serial.begin(38400);

  for (int Pin = 0; Pin < PinNo; Pin++) {
    pinMode(Pins[Pin], OUTPUT);
    digitalWrite(Pins[Pin], LOW);
  }
  
  char ch = 0;
}

void loop() {

  char ch;

  while (ch == 0) {
    ch = Serial.read();
  }

  if (ch == 'A') {
    digitalWrite(Pins[0], HIGH);
    delay(Delay);
    digitalWrite(Pins[0], LOW);
  }
  
  if (ch == 'a') {
    digitalWrite(Pins[0], HIGH);
  }

  if (ch == 'z') {
    digitalWrite(Pins[0], LOW);
  }


  if (ch == 'B') {
    digitalWrite(Pins[1], HIGH);
    delay(Delay);
    digitalWrite(Pins[1], LOW);
  }

  if (ch == 'C') {
    digitalWrite(Pins[2], HIGH);
    delay(Delay);
    digitalWrite(Pins[2], LOW);
  }

  if (ch == 'D') {
    digitalWrite(Pins[3], HIGH);
    delay(Delay);
    digitalWrite(Pins[3], LOW);
  }

  if (ch == 'E') {
    digitalWrite(Pins[4], HIGH);
    delay(Delay);
    digitalWrite(Pins[4], LOW);
  }
  
  if (ch == 'P') {
    digitalWrite(Pins[5], HIGH);
    delay(Delay);
    digitalWrite(Pins[5], LOW);
  }
  
}
