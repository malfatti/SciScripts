#!/bin/bash

while getopts fs Option; do
    case "${Option}" in
        f)
            File=${OPTARG}
        ;;
        s)
            SHA=${OPTARG}
        ;;
    esac
done

SHACS=`grep -A 1 -i sha512 "$SHA" | awk '{ print $4 }' | sed -n 1p`
FileCS=`sha512sum "$File" | awk '{ print $1 }'

if [ SHACS == FileCS ]; then
    echo "File $File is valid";
else; 
    echo 'File is NOT ok!!';
    echo "File = $FILECS";
    echo "SHA  = $SHA";
fi
