#!/bin/bash

Governors="$(cat /sys/devices/system/cpu/cpu0/cpufreq/scaling_available_governors)"
Governor="$1"

if test "${Governors#*"$Governor"}" != "$Governors"; then
    for cpu in /sys/devices/system/cpu/cpu[0-9]*; do 
        echo "$Governor" > $cpu/cpufreq/scaling_governor; 
    done
elif [ "${Governor,,}" == current ]; then
    cat /sys/devices/system/cpu/cpu0/cpufreq/scaling_governor;
elif [ "${Governor,,}" == list ]; then
    echo "$Governors"
else
    echo "Usage: SetGovernor "$Governors""
fi
