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

The code to increase sampling rate to 38.5KHz (by setting prescale to 16) 
is from jmknapp. Source:
http://forum.arduino.cc/index.php?topic=6549.0
*/

// defines for setting and clearing register bits
#ifndef cbi
#define cbi(sfr, bit) (_SFR_BYTE(sfr) &= ~_BV(bit))
#endif
#ifndef sbi
#define sbi(sfr, bit) (_SFR_BYTE(sfr) |= _BV(bit))
#endif

// code:
const int Piezo = A0;
const int TTLs = A3;

void setup() {
  // set prescale to 16
  sbi(ADCSRA,ADPS2) ;
  cbi(ADCSRA,ADPS1) ;
  cbi(ADCSRA,ADPS0) ;

  Serial.begin(115200);
  analogReference(INTERNAL);
}


void loop() {
  int PiezoV = analogRead(Piezo);
  //delay(1);
  int TTLsV = analogRead(TTLs);
  
  Serial.print(PiezoV);
  Serial.print("\t");
  Serial.print(TTLsV);
  Serial.println();
}
