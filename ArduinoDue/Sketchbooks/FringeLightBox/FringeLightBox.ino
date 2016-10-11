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

// Set filter
const float HighPassFreq = 8.0;
const float LowPassFreq = 15.0;
FilterOnePole HighPassFilter( HIGHPASS, HighPassFreq );
FilterOnePole LowPassFilter( LOWPASS, LowPassFreq );

// Set Data
int Chunk = 1024;
float AmpMax = 1000;
float AmpMin = -1000;

float AmpMean = (AmpMax - AmpMin) / 2
                float Slot = ((AmpMax - AmpMin) / (LastPin - FirstPin)) / 2

int GenRandomInt(int Min, int Max, int Amount) {
  int Numbers[Amount];

  for (No = 0; No < Amount; No++) {
    Numbers[Amount] = random(Min, Max);
  }

  return Numbers;
}

int Shuffle(int Array) {
  // from Darth Hunterix
  // (http://stackoverflow.com/questions/32413209/shuffle-an-array-in-arduino-software#32417244)

  const size_t n = sizeof(Array) / sizeof(Array[0]);

  for (size_t i = 0; i < n - 1; i++) {
    size_t j = random(0, n - i);

    int t = Array[i];
    Array[i] = Array[j];
    Array[j] = t;
  }

  return Array;
}

void setup() {
  Serial.begin(38400);

  // Set free running mode on ADC7 (pin A0)
  ADC->ADC_MR |= 0x80;
  ADC->ADC_CR = 2;
  ADC->ADC_CHER = 0x80;


  // Set all led pins HIGH
  for (int Pin = FirstPin; Pin < LastPin + 1; Pin++) {
    LedPins[Pin - FirstPin] = Pin;
    pinMode(LedPins[Pin - FirstPin], OUTPUT);
    digitalWrite(LedPins[Pin - FirstPin], HIGH);
  }
}

void loop() {
  //Shuffle pins
  LedPins = Shuffle(LedPins);
  bool AllOff = false;
  
  while ( AllOff != true ) {
    // Set Data
    float Data[Chunk];
    float DataSum = 0.0;
    float DataMean = 0.0;
    float DataMin = 1023.0;
    float DataMax = -1023.0;

    // Acquire and filter data
    for (int Sample = 0; Sample < Chunk; Sample++) {
      while ((ADC->ADC_ISR & 0x80) == 0); // wait for conversion
      HighPassFilter.input(ADC->ADC_CDR[7]);
      LowPassFilter.input(HighPassFilter.output());

      Data[Sample] = LowPassFilter.output();

      if ( Data[Sample] > DataMax ) DataMax = Data[Sample];
      if ( Data[Sample] < minim ) DataMin = Data[Sample];
    }

    // Sum data samples
    for (Sample = 0; Sample < Chunk; Sample++) {
      DataSum += Data[Sample];
    }

    DataMean = DataSum / Chunk; // Mean
    DataPtP = DataMax - DataMin; // Peak-to-peak

    // Turn on leds depending on alpha intensity
    for (float Limit = Slot; Limit <= AmpMax; Limit = Limit + Slot) {
      float LimitP = Limit - slot;
      if (DataMean < (AmpMean - LimitP) && DataMean >= (AmpMean - Limit)) {
        for (int Led = 0; Led < Limit/Slot; Led++) {
          digitalWrite(LedPins[Led], LOW)
        }
        if (Limit/Slot == LastPin-FirstPin) { AllOff = true; }
      }

      if (DataMean > (AmpMean + LimitP) && DataMean < (AmpMean + Limit) {
        for (int Led = 0; Led < Limit/Slot; Led++) {
            digitalWrite(LedPins[Led], LOW)
        }
        if (Limit/Slot == LastPin-FirstPin) { AllOff = true; }
      }
    }
  }

  for (int Pin = FirstPin; Pin < LastPin + 1; Pin++) {
    digitalWrite(LedPins[Pin - FirstPin], HIGH);
  }
}
