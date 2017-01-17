void setup() {
//  Reference: https://www.arduino.cc/en/Reference/PortManipulation
  DDRD = B11111111;
  PORTD = B00000000;

}

void loop() {
 PORTD = PORTD << 1;
}
