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

#include <Filters.h>

// Set led pins
const int FirstPin = 22;
const int LastPin = 52;
int LedPins[LastPin - FirstPin];

// Initialize variables
int Data;

// Set filter
const float HighPassFreq = 8.0;
const float LowPassFreq = 15.0;
FilterOnePole HighPassFilter( HIGHPASS, HighPassFreq );
FilterOnePole LowPassFilter( LOWPASS, LowPassFreq );


void setup()
{
  Serial.begin(38400);
  
  // Set free running mode on ADC7 (pin A0)
  ADC->ADC_MR |= 0x80;
  ADC->ADC_CR = 2;
  ADC->ADC_CHER = 0x80;


  // Set all led pins LOW
  for (int Pin = FirstPin; Pin < LastPin + 1; Pin++) {
    LedPins[Pin - FirstPin] = Pin;
    pinMode(LedPins[Pin - FirstPin], OUTPUT);
    digitalWrite(LedPins[Pin - FirstPin], LOW);
  }
}

void loop()
{
  while ((ADC->ADC_ISR & 0x80) == 0); // wait for conversion
  Data = ADC->ADC_CDR[7];
  FilteredData = 
  //  val = map(val, 0, 1023, 0, 255);

  Serial.println(Data);
}
