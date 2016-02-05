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

Arduino will read serial inputs and print them in serial,separated by tabs
and lines.
*/

const int Piezo = A0;
const int PiezoTTL = A3;

void setup() {
  Serial.begin(115200);
  analogReference(INTERNAL);
}


void loop() {
  Serial.print(analogRead(Piezo));
  Serial.print("\t");
  Serial.print(analogRead(PiezoTTL));
  Serial.println();
}
