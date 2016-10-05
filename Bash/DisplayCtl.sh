#!/bin/bash

## Change default display
Display="$1"

case "$Display" in
	"HdmiE") 
		xrandr --output HDMI1 --auto
		xrandr --output HDMI1 --right-of eDP1
		sleep 1
		feh --bg-scale --randomize ~/Images/Others/Wallpapers/*
		;;
	"HdmiM")
		xrandr --output HDMI1 --auto
                xrandr --output HDMI1 --same-as eDP1
                sleep 1
                feh --bg-scale --randomize ~/Images/Others/Wallpapers/*
                ;; 
	"HdmiOff")
		xrandr --output HDMI1 --off
                sleep 1                                                         
                feh --bg-scale --randomize ~/Images/Others/Wallpapers/*         
                ;; 
	*)
		printf "Available displays: \n  HdmiE, HdmiM, HdmiOff."
esac
