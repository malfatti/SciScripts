const int AI0 = A0;
int Data = 0;

void setup() {
  Serial.begin(19200);
  analogReference(INTERNAL);
}

void loop() {
  
  Data = analogRead(AI0);
  Serial.println(Data);
}
