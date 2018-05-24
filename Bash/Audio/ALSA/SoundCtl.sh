#!/bin/bash

#JackPath="/home/malfatti/Software/Programs/Jack1/bin/"
JackPath="/usr/bin/"

if [ "${1,,}" == hdmi ]; then
    killall jackd &> /dev/null
    #cp ~/Software/Git/Malfatti/SciScripts/Bash/Audio/ALSA/HDMI ~/.asoundrc
    #cp ~/Software/Git/Malfatti/SciScripts/Bash/Audio/ALSA/HDMIVol.sh ~/Software/Git/Malfatti/SciScripts/Bash/Audio/ALSA/VolControl.sh
    cp ~/Software/Git/Malfatti/SciScripts/Bash/Audio/ALSA/Analog ~/.asoundrc
    cp ~/Software/Git/Malfatti/SciScripts/Bash/Audio/ALSA/HDMIVol.sh ~/Software/Git/Malfatti/SciScripts/Bash/Audio/ALSA/VolControl.sh
    bash ~/Software/Git/Malfatti/SciScripts/Bash/Audio/ALSA/JackCtl.sh --path "$JackPath" --rt r --card HDMI,3 --rate 48000 --periodsize 1024 --periods 3 --priority 19 --bridge false

elif [ "${1,,}" == analog ]; then
    killall jackd &> /dev/null
    cp ~/Software/Git/Malfatti/SciScripts/Bash/Audio/ALSA/Analog ~/.asoundrc
    cp ~/Software/Git/Malfatti/SciScripts/Bash/Audio/ALSA/PCHVol.sh ~/Software/Git/Malfatti/SciScripts/Bash/Audio/ALSA/VolControl.sh

elif [ "${1,,}" == jack ]; then
    killall jackd &> /dev/null
    cp ~/Software/Git/Malfatti/SciScripts/Bash/Audio/ALSA/Jack-ALSA ~/.asoundrc
    cp ~/Software/Git/Malfatti/SciScripts/Bash/Audio/ALSA/PCHVol.sh ~/Software/Git/Malfatti/SciScripts/Bash/Audio/ALSA/VolControl.sh
    bash ~/Software/Git/Malfatti/SciScripts/Bash/Audio/ALSA/JackCtl.sh --path "$JackPath" --rt r --card PCH --rate 48000 --periodsize 1024 --periods 2 --priority 19 --bridge true

elif [ "${1,,}" == jackonly ]; then
    killall jackd &> /dev/null
    cp ~/Software/Git/Malfatti/SciScripts/Bash/Audio/ALSA/Jack-ALSA ~/.asoundrc
    cp ~/Software/Git/Malfatti/SciScripts/Bash/Audio/ALSA/PCHVol.sh ~/Software/Git/Malfatti/SciScripts/Bash/Audio/ALSA/VolControl.sh
    bash ~/Software/Git/Malfatti/SciScripts/Bash/Audio/ALSA/JackCtl.sh --path "$JackPath" --rt r --card PCH --rate 96000 --periodsize 1024 --periods 2 --priority 19 --bridge false

elif [ "${1,,}" == jackrt ]; then
    killall jackd &> /dev/null
    cp ~/Software/Git/Malfatti/SciScripts/Bash/Audio/ALSA/Jack-ALSA ~/.asoundrc
    cp ~/Software/Git/Malfatti/SciScripts/Bash/Audio/ALSA/PCHVol.sh ~/Software/Git/Malfatti/SciScripts/Bash/Audio/ALSA/VolControl.sh
    bash ~/Software/Git/Malfatti/SciScripts/Bash/Audio/ALSA/JackCtl.sh --path "$JackPath" --rt R --card PCH --rate 192000 --periodsize 384 --periods 3 --priority 89 --bridge true

elif [ "${1,,}" == jackonlyrt ]; then
    killall jackd &> /dev/null
    cp ~/Software/Git/Malfatti/SciScripts/Bash/Audio/ALSA/Jack-ALSA ~/.asoundrc
    cp ~/Software/Git/Malfatti/SciScripts/Bash/Audio/ALSA/PCHVol.sh ~/Software/Git/Malfatti/SciScripts/Bash/Audio/ALSA/VolControl.sh
    bash ~/Software/Git/Malfatti/SciScripts/Bash/Audio/ALSA/JackCtl.sh --path "$JackPath" --rt R --card PCH --rate 192000 --periodsize 512 --periods 4 --priority 89 --bridge false

elif [ "${1,,}" == jackonlyrtsu ]; then
    sudo killall jackd &> /dev/null
    sudo cp /home/malfatti/Software/Git/Malfatti/SciScripts/Bash/Audio/ALSA/Jack-ALSA /root/.asoundrc
    sudo bash /home/malfatti/Software/Git/Malfatti/SciScripts/Bash/Audio/ALSA/JackCtl.sh --path "$JackPath" --rt R --card PCH --rate 192000 --periodsize 512 --periods 4 --priority 89 --bridge false

elif [ "${1,,}" == jackusbpre2 ]; then
    killall jackd &> /dev/null
    cp ~/Software/Git/Malfatti/SciScripts/Bash/Audio/ALSA/Jack-ALSA ~/.asoundrc
    cp ~/Software/Git/Malfatti/SciScripts/Bash/Audio/ALSA/PCHVol.sh ~/Software/Git/Malfatti/SciScripts/Bash/Audio/ALSA/VolControl.sh
    bash ~/Software/Git/Malfatti/SciScripts/Bash/Audio/ALSA/JackCtl.sh --path "$JackPath" --rt R --card USBPre2 --rate 192000 --periodsize 384 --periods 3 --priority 89 --bridge true

elif [ "${1,,}" == spdif ]; then
    killall jackd &> /dev/null
    cp ~/Software/Git/Malfatti/SciScripts/Bash/Audio/ALSA/SPDIf ~/.asoundrc
    cp ~/Software/Git/Malfatti/SciScripts/Bash/Audio/ALSA/PCHVol.sh ~/Software/Git/Malfatti/SciScripts/Bash/Audio/ALSA/VolControl.sh

elif [ "${1,,}" == info ]; then
    if [ "${2,,}" == devsinuse ]; then
        fuser -fv /dev/snd/* /dev/dsp*
    elif [ "${2,,}" == jacklog ]; then
        less  ~/.log/JackSession.log
    elif [ "${2,,}" == rtthreads ]; then
        ps -eLo rtprio,pri,cgroup,class,pid,pcpu,%mem,user,comm --sort pri | less
    else
        less ~/.asoundrc
    fi
else
    echo "Usage: SoundCtl [hdmi | analog | jack | jackrt | jackusbpre2 | spdif | info]"
    echo ""
fi
