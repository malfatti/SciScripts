#!/bin/bash

KernelName="$(uname -r)";
if test "${KernelName#*rt}" != "$KernelName"; then
    bash ~/.MyScripts/ALSA/SoundCtl.sh jackrt
else
    bash ~/.MyScripts/ALSA/SoundCtl.sh jack
fi
