unsigned long val;

void setup()
{
  Serial.begin(38400);

  // Set free running mode on ADC7 (pin A0)
  ADC->ADC_MR |= 0x80;
  ADC->ADC_CR = 2;
  ADC->ADC_CHER = 0x80;
}

void loop()
{
  while ((ADC->ADC_ISR & 0x80) == 0); // wait for conversion
  val = ADC->ADC_CDR[7];
  //  val = map(val, 0, 1023, 0, 255);

  Serial.println(val);
}
