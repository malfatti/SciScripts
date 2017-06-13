#!/bin/bash

Governors="$(cat /sys/devices/system/cpu/cpu0/cpufreq/scaling_available_governors)"
Governor="$1"

if test "${Governors#*"$Governor"}" != "$Governors"; then
    for cpu in /sys/devices/system/cpu/cpu[0-9]*; do 
       echo "$Governor" | sudo tee $cpu/cpufreq/scaling_governor > /dev/null;
    done
elif [ "${Governor,,}" == status ]; then
    cat /sys/devices/system/cpu/cpu0/cpufreq/scaling_governor;
elif [ "${Governor,,}" == list ]; then
    echo "$Governors"
else
    echo "Usage: "$0" ["$Governors"]"
fi
