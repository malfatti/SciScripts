#!/bin/bash
Notebook=eDP1
HDMI=HDMI1
StateFile="$SCRIPTSPATH/Bash/Video/DisplayState"
State="$(cat $StateFile)"

if xrandr | grep "$HDMI connected"; then
    if [ $State == "OutOff" ]; then
        xrandr --output "$HDMI" --auto
        xrandr --output "$HDMI" --right-of "$Notebook"
        feh --bg-scale --randomize ~/Images/Others/Wallpapers/*
        echo "Extended" > $StateFile
    
    elif [ $State == "Extended" ]; then
        xrandr --output "$HDMI" --auto
        xrandr --output "$HDMI" --same-as "$Notebook"
        feh --bg-scale --randomize ~/Images/Others/Wallpapers/*
        echo "Mirror" > $StateFile
    
    #elif [ $State == "Mirror" ]; then
    #    xrandr --output "$HDMI" --auto
    #    xrandr --output "$Notebook" --off
    #    echo "InOff" > $StateFile
    
    else
        xrandr --output "$Notebook" --auto
        xrandr --output "$HDMI" --off
        feh --bg-scale --randomize ~/Images/Others/Wallpapers/*
        echo "OutOff" > $StateFile
    fi
else
    xrandr --output "$HDMI" --off
    feh --bg-scale --randomize ~/Images/Others/Wallpapers/*
    echo "OutOff" > $StateFile
fi
