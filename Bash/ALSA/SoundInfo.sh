#!/bin/bash

if [ "${1,,}" == asoundrc ]; then
    less ~/.asoundrc
elif [ "${1,,}" == devsinuse ]; then
    fuser -fv /dev/snd/* /dev/dsp*
elif [ "${1,,}" == jacklog ]; then
    less  ~/.MyScripts/ALSA/JackSession.log
elif [ "${1,,}" == rtthreads ]; then
    ps -eLo rtprio,pri,cgroup,class,pid,pcpu,%mem,user,comm --sort pri | less
else
    echo "Usage: SoundInfo [asoundrc | devsinuse | jacklog | rtthreads]"
fi
