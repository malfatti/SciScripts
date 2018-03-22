#include <MMA8452Q.h>
#include <Wire.h> // so the Arduino IDE auto-detects the dependency

MMA8452Q accel;

int axes[3];

void setup() {
  Serial.begin(115200);
  Serial.print("working");
  if(accel.begin());
    while (1); // error
}

void loop() {
  // get and print raw axes values
  accel.axes(axes);

 Serial.print("x: ");
  Serial.print(axes[0]);
  Serial.print(", y: ");
  Serial.print(axes[1]);
  Serial.print(", z: ");
  Serial.println(axes[2]);

  delay(200);
}

