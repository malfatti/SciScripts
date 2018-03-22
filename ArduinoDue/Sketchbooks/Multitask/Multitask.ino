#include <Wire.h>                                                               
#include <MMA8453_n0m1.h>                                                       
                                                                                
MMA8452Q accel;
int axes[3];
 
void setup() 
{ 
  pinMode(9, OUTPUT);
  pinMode(10, OUTPUT);
  pinMode(11, OUTPUT);
}

void loop()
{
  accel.axes(axes);
  
  int X = map(accel.x(), 0, 1023, 0, 255);
  int Y = map(accel.y(), 0, 1023, 0, 255);
  int Z = map(accel.z(), 0, 1023, 0, 255);

  analogWrite(9, X);
  analogWrite(10, Y);
  analogWrite(11, Z);
  
  delay(20);
}
