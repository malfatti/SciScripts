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

*/

const int TTLs[6] = {13, 12, 8, 7, 4, 2};

void setup() {
  Serial.begin(19200);
  for (int i = 0; i < 6; i++) {
    pinMode(TTLs[i], OUTPUT);
    digitalWrite(TTLs[i], LOW);
  }
}

void loop() {
  
  digitalWrite(TTLs[0], HIGH);
  digitalWrite(TTLs[1], HIGH);
  digitalWrite(TTLs[2], HIGH);
  digitalWrite(TTLs[3], HIGH);
  digitalWrite(TTLs[4], HIGH);
  digitalWrite(TTLs[5], HIGH);
  delay(10);
  digitalWrite(TTLs[0], LOW);
  delay(40);
  digitalWrite(TTLs[1], LOW);
  delay(20);
  digitalWrite(TTLs[2], LOW);
  delay(30);
  digitalWrite(TTLs[3], LOW);
  delay(50);
  digitalWrite(TTLs[4], LOW);
  delay(50);
  digitalWrite(TTLs[5], LOW);
  delay(50);

}
