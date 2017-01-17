#!/bin/bash

if [ "${1,,}" == hdmi ]; then
    killall jackd &> /dev/null
    cp ~/Software/Git/Malfatti/SciScripts/Bash/ALSA/HDMI ~/.asoundrc
    cp ~/Software/Git/Malfatti/SciScripts/Bash/ALSA/HDMIVol.sh ~/Software/Git/Malfatti/SciScripts/Bash/ALSA/VolControl.sh
elif [ "${1,,}" == analog ]; then
    killall jackd &> /dev/null
    cp ~/Software/Git/Malfatti/SciScripts/Bash/ALSA/Analog ~/.asoundrc
    cp ~/Software/Git/Malfatti/SciScripts/Bash/ALSA/PCHVol.sh ~/Software/Git/Malfatti/SciScripts/Bash/ALSA/VolControl.sh
elif [ "${1,,}" == jack ]; then
    killall jackd &> /dev/null
    cp ~/Software/Git/Malfatti/SciScripts/Bash/ALSA/Jack-ALSA ~/.asoundrc
    cp ~/Software/Git/Malfatti/SciScripts/Bash/ALSA/PCHVol.sh ~/Software/Git/Malfatti/SciScripts/Bash/ALSA/VolControl.sh
    bash ~/Software/Git/Malfatti/SciScripts/Bash/ALSA/JackCtl.sh --rt r --rate 48000 --periodsize 1024 --periods 3 --priority 19
elif [ "${1,,}" == jackrt ]; then
    killall jackd &> /dev/null
    cp ~/Software/Git/Malfatti/SciScripts/Bash/ALSA/Jack-ALSA ~/.asoundrc
    cp ~/Software/Git/Malfatti/SciScripts/Bash/ALSA/PCHVol.sh ~/Software/Git/Malfatti/SciScripts/Bash/ALSA/VolControl.sh
    bash ~/Software/Git/Malfatti/SciScripts/Bash/ALSA/JackCtl.sh --rt R --rate 192000 --periodsize 128 --periods 9 --priority 89
elif [ "${1,,}" == spdif ]; then
    killall jackd &> /dev/null
    cp ~/Software/Git/Malfatti/SciScripts/Bash/ALSA/SPDIf ~/.asoundrc
    cp ~/Software/Git/Malfatti/SciScripts/Bash/ALSA/PCHVol.sh ~/Software/Git/Malfatti/SciScripts/Bash/ALSA/VolControl.sh
else
    echo "Usage: SoundCtl [hdmi | analog | jack | jackrt | spdif]"
fi
