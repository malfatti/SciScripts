#!/bin/bash

Path=/sys/class/leds/beaglebone\:green\:usr3

for NUM in `seq 1 1 100`; do
    echo 255 > "$Path"/brightness
    sleep 1
    echo 0 > "$Path"/brightness
    sleep 1
done