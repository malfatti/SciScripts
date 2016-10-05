#!/bin/bash

KernelName="$(uname -r)";
if test "${KernelName#*rt}" != "$KernelName"; then
    bash ~/Software/Git/Malfatti/SciScripts/Bash/ALSA/SoundCtl.sh jackrt
else
    bash ~/Software/Git/Malfatti/SciScripts/Bash/ALSA/SoundCtl.sh jack
fi
