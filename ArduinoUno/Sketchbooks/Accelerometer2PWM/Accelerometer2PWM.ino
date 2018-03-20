#include <Wire.h>
#include <MMA8453_n0m1.h>

MMA8453_n0m1 accel;

void setup()
{
  accel.setI2CAddr(0x1C); //change your device address if necessary, default is 0x1C
  accel.dataMode(true, 2); //enable highRes 10bit, 2g range [2g,4g,8g]

  pinMode(9, OUTPUT);
  pinMode(10, OUTPUT);
  pinMode(11, OUTPUT);
}

void loop()
{
  accel.update();
  
  int X = map(accel.x(), 0, 1023, 0, 255);
  int Y = map(accel.y(), 0, 1023, 0, 255);
  int Z = map(accel.z(), 0, 1023, 0, 255);

  analogWrite(9, X);
  analogWrite(10, Y);
  analogWrite(11, Z);

  delay(20)

}
