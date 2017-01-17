/*
// Blink 8 leds pseudorandomly

int Pins[8] = {6, 7, 8, 9, 10, 11, 12, 13};
int PinsNo = 8;

void setup() {
  for(int Pin = 0; Pin < PinsNo; Pin++) {
    pinMode(Pins[Pin], OUTPUT);
    digitalWrite(Pins[Pin], LOW); 
  }
  randomSeed(analogRead(A0));
}

void loop() {
  for(int Pin = 0; Pin < PinsNo; Pin++) {
    digitalWrite(Pins[Pin], random(2));
    delay(30);
  }
}
*/

// Generate random binary 8bit numbers and send to serial
void setup() {
  Serial.begin(38400);
  randomSeed(analogRead(A0));
}

void loop() {
  while (true) { 
    if(digitalRead(2) == 5) { break }
  }
  
  int Number = random(255);
  String Binary = String(Number, BIN);
  Serial.println(Binary);
}

