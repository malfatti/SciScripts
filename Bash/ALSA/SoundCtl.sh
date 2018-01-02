#!/bin/bash

#JackPath="/home/malfatti/Software/Programs/Jack1/bin/"
JackPath="/usr/bin/"

if [ "${1,,}" == hdmi ]; then
    killall jackd &> /dev/null
    #cp ~/Software/Git/Malfatti/SciScripts/Bash/ALSA/HDMI ~/.asoundrc
    #cp ~/Software/Git/Malfatti/SciScripts/Bash/ALSA/HDMIVol.sh ~/Software/Git/Malfatti/SciScripts/Bash/ALSA/VolControl.sh
    cp ~/Software/Git/Malfatti/SciScripts/Bash/ALSA/Analog ~/.asoundrc
    cp ~/Software/Git/Malfatti/SciScripts/Bash/ALSA/HDMIVol.sh ~/Software/Git/Malfatti/SciScripts/Bash/ALSA/VolControl.sh
    bash ~/Software/Git/Malfatti/SciScripts/Bash/ALSA/JackCtl.sh --path "$JackPath" --rt r --card HDMI,3 --rate 48000 --periodsize 1024 --periods 3 --priority 19 --bridge false

elif [ "${1,,}" == analog ]; then
    killall jackd &> /dev/null
    cp ~/Software/Git/Malfatti/SciScripts/Bash/ALSA/Analog ~/.asoundrc
    cp ~/Software/Git/Malfatti/SciScripts/Bash/ALSA/PCHVol.sh ~/Software/Git/Malfatti/SciScripts/Bash/ALSA/VolControl.sh

elif [ "${1,,}" == jack ]; then
    killall jackd &> /dev/null
    cp ~/Software/Git/Malfatti/SciScripts/Bash/ALSA/Jack-ALSA ~/.asoundrc
    cp ~/Software/Git/Malfatti/SciScripts/Bash/ALSA/PCHVol.sh ~/Software/Git/Malfatti/SciScripts/Bash/ALSA/VolControl.sh
    bash ~/Software/Git/Malfatti/SciScripts/Bash/ALSA/JackCtl.sh --path "$JackPath" --rt r --card PCH --rate 48000 --periodsize 1024 --periods 3 --priority 19 --bridge true

elif [ "${1,,}" == jackrt ]; then
    killall jackd &> /dev/null
    cp ~/Software/Git/Malfatti/SciScripts/Bash/ALSA/Jack-ALSA ~/.asoundrc
    cp ~/Software/Git/Malfatti/SciScripts/Bash/ALSA/PCHVol.sh ~/Software/Git/Malfatti/SciScripts/Bash/ALSA/VolControl.sh
    bash ~/Software/Git/Malfatti/SciScripts/Bash/ALSA/JackCtl.sh --path "$JackPath" --rt R --card PCH --rate 192000 --periodsize 384 --periods 3 --priority 89 --bridge true

elif [ "${1,,}" == jackusbpre2 ]; then
    killall jackd &> /dev/null
    cp ~/Software/Git/Malfatti/SciScripts/Bash/ALSA/Jack-ALSA ~/.asoundrc
    cp ~/Software/Git/Malfatti/SciScripts/Bash/ALSA/PCHVol.sh ~/Software/Git/Malfatti/SciScripts/Bash/ALSA/VolControl.sh
    bash ~/Software/Git/Malfatti/SciScripts/Bash/ALSA/JackCtl.sh --path "$JackPath" --rt R --card USBPre2 --rate 192000 --periodsize 384 --periods 3 --priority 89 --bridge true

elif [ "${1,,}" == spdif ]; then
    killall jackd &> /dev/null
    cp ~/Software/Git/Malfatti/SciScripts/Bash/ALSA/SPDIf ~/.asoundrc
    cp ~/Software/Git/Malfatti/SciScripts/Bash/ALSA/PCHVol.sh ~/Software/Git/Malfatti/SciScripts/Bash/ALSA/VolControl.sh

else
    echo "Usage: SoundCtl [hdmi | analog | jack | jackrt | jackusbpre2 | spdif]"
fi
