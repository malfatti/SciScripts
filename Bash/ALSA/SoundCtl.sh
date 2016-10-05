#!/bin/bash

if [ "${1,,}" == hdmi ]; then
    killall jackd &> /dev/null
    cp ~/.MyScripts/ALSA/HDMI ~/.asoundrc
    cp ~/.MyScripts/ALSA/HDMIVol.sh ~/.MyScripts/ALSA/VolControl.sh
elif [ "${1,,}" == analog ]; then
    killall jackd &> /dev/null
    cp ~/.MyScripts/ALSA/Analog ~/.asoundrc
    cp ~/.MyScripts/ALSA/PCHVol.sh ~/.MyScripts/ALSA/VolControl.sh
elif [ "${1,,}" == jack ]; then
    killall jackd &> /dev/null
    cp ~/.MyScripts/ALSA/Jack-ALSA ~/.asoundrc
    cp ~/.MyScripts/ALSA/PCHVol.sh ~/.MyScripts/ALSA/VolControl.sh
    bash ~/.MyScripts/ALSA/JackCtl.sh --rt r --rate 48000 --periodsize 1024 --periods 3 --priority 19
elif [ "${1,,}" == jackrt ]; then
    killall jackd &> /dev/null
    cp ~/.MyScripts/ALSA/Jack-ALSA ~/.asoundrc
    cp ~/.MyScripts/ALSA/PCHVol.sh ~/.MyScripts/ALSA/VolControl.sh
    bash ~/.MyScripts/ALSA/JackCtl.sh --rt R --rate 192000 --periodsize 384 --periods 3 --priority 89
elif [ "${1,,}" == spdif ]; then
    killall jackd &> /dev/null
    cp ~/.MyScripts/ALSA/SPDIf ~/.asoundrc
    cp ~/.MyScripts/ALSA/PCHVol.sh ~/.MyScripts/ALSA/VolControl.sh
else
    echo "Usage: SoundCtl [hdmi | analog | jack | jackrt | spdif]"
fi
