#!/bin/bash

Set=$@

if [ \( "${Set,,}" == on \) ]; then
    echo "1" | sudo tee /sys/devices/system/cpu/cpufreq/boost > /dev/null

elif [ \( "${Set,,}" == off \) ]; then
    echo "0" | sudo tee /sys/devices/system/cpu/cpufreq/boost > /dev/null

elif [ \( "${Set,,}" == status \) ]; then
    cat /sys/devices/system/cpu/cpufreq/boost

else
    echo "Please use 'on', 'off' or 'status'."

fi
