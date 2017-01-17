LedNo = 8
Pins[8] = {6, 7, 8, 9, 10, 11, 12, 13};

void setup() {
  for(int Pin = 0; Pin < LedNo; Pin++) {
    digitalWrite(Pin, LOW)
  }

}

void loop() {
  for(int i = 0, x = 1; i < binLength; i++, x+=2) { 
      if(binNumber[i] == '0') state = LOW;
      if(binNumber[i] == '1') state = HIGH;
      digitalWrite(pins[i] + binLength - x, state);
    } 
}
