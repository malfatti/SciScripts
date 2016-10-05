#!/bin/bash

if [ "${1,,}" == up ]; then
    amixer -c PCH sset Master,0 1+
elif [ "${1,,}" == down ]; then
    amixer -c PCH sset Master,0 1-
elif [ "${1,,}" == toggle ]; then
    amixer -c PCH sset Master,0 toggle
else
    echo "Usage: up, down or toggle"
fi
