#!/bin/bash

echo "Are you sure? [y/N]"; read Ans1

if [ \( "${Ans1,,}" == y \) -o \( "${Ans1,,}" == yes \) ]; then
    echo "ARE YOU COMPLETELY SURE? [y/N]"; read Ans

    if [ \( "${Ans,,}" == y \) -o \( "${Ans,,}" == yes \) ]; then
        echo "Removing files..."
        echo "$1" | xargs -n 1000 rm
        echo "Done."
    else
        echo "Aborted."
    fi
else
    echo "Aborted."
fi
