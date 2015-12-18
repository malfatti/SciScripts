int analogPin = 3;
int val = 0;

void setup()
{
  Serial.begin(19200);
  analogReference(INTERNAL);
}

void loop()
{
  val = analogRead(analogPin);
  val = map(val, 0, 1023, 0, 255);

  Serial.println(val);
}
