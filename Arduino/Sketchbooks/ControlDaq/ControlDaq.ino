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

In this program, Arduino will read serial input and, depending on the value, 
it will turn on and off (5ms delay) a digital pin. We use this to sinchronize 
DAQs in different computers.

Sumarizing the commands:
A = turn D13 on
*/

const int Pin13 =  13;

void setup() {
  Serial.begin(19200);
  pinMode(Pin13, OUTPUT);

  digitalWrite(Pin13, LOW);

  char ch = 0;
}

void loop() {

  char ch;

  while (ch == 0) {
    ch = Serial.read();
  }

  if (ch == 'A') {
      digitalWrite(Pin13, HIGH);
      delay(2)
      digitalWrite(Pin13, LOW);
  }

}
