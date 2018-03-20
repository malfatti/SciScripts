#!/bin/bash

Mode="$1"
In="$2"
Out="$3"
BS="$4"

if [ "${Mode,,}" == backup ]; then
    sudo dd if="$In" conv=sync,noerror bs="$BS" | gzip -c  > "$Out"
elif [ "${Mode,,}" == restore ]; then
    sudo gunzip -c "$In" | dd of="$Out" bs="$BS" status=progress
else
    echo "Usage:"
    echo "    DD <mode> <input> <output> <bs>"
fi
