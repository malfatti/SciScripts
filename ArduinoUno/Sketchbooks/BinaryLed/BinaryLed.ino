int NI = 7;
int Led = 8;

void setup() {
  pinMode(NI, INPUT);
  pinMode(Led, OUTPUT);
  digitalWrite(Led, LOW);
  
}

void loop() {
  if(digitalRead(NI) == HIGH) {
    digitalWrite(Led, HIGH);
  }
  else {
    digitalWrite(Led, LOW);
  }
}
