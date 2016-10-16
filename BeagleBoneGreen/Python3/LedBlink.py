#!/usr/bin/env python3
import AdafruitB_BIO.GPIO as GPIO
import time
GPIO.setup("USR3", GPIO.OUT)
while True:
    GPIO.output("USR3", GPIO.HIGH)
    time.sleep(0.5)
    GPIO.output("USR3", GPIO.LOW)
    time.sleep(0.5)