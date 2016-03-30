int analogPin = 2;
int val = 0;

void setup()
{
  Serial.begin(38400);
  analogReference(INTERNAL);
}

void loop()
{
  val = analogRead(analogPin);
  val = map(val, 0, 1023, 0, 255);

  Serial.println(val);
}
