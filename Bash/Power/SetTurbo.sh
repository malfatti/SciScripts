#!/bin/bash

Set=$@
TurboFile=/sys/devices/system/cpu/intel_pstate/no_turbo

if [ \( "${Set,,}" == on \) ]; then
    echo "0" | sudo tee "$TurboFile" > /dev/null

elif [ \( "${Set,,}" == off \) ]; then
    echo "1" | sudo tee "$TurboFile" > /dev/null

elif [ \( "${Set,,}" == status \) ]; then
    cat "$TurboFile"

else
    echo "Usage: "$0" [on | off | status]"

fi
