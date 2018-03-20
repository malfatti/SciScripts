int val = 0;

void setup() {
  Serial.begin(115200);
//  Serial.println('a');
//  char a = 'b';
//  while (a != 'a') { a = Serial.read(); }
  
  pinMode(8, OUTPUT); digitalWrite(8, LOW);
  pinMode(11, OUTPUT);
}

void loop() {
//  while (Serial.available() == 0) {}
//  if (Serial.available() > 0) {
    int sensorvalue = analogRead(A0);
    float voltage = sensorvalue * (5.0 / 1023.0);
    sensorvalue = map(sensorvalue, 0, 1023, 0, 255);
    analogWrite(11, sensorvalue);
    Serial.println(voltage); Serial.println('\n');
    
//    val = Serial.read();
//    if (val == 'R') { Serial.println(voltage); }
//    if (val == 'P') { digitalWrite(8, HIGH); }
//    if (val == 'p') { digitalWrite(8, LOW); }
//  }
  delay(20);
}


