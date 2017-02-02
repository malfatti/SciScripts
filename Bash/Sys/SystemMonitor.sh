#!/bin/bash

# Path to battery /sys
BATPath=/sys/class/power_supply/BAT1/

# Low battery in %
BATLow=25

# Critical battery in % 
BATCritical=16

# Files to monitor size
Files=~/.log/*

# Maximum file size in bytes
MaxFileSize=5120000


BatteryCheck() {
#echo "BM loop starting..."
if [ -e $BATPath ]; then
	On=$(cat $BATPath/status)
    
	if [ $On == Discharging ]; then
		Current=$(cat $BATPath/capacity)
        
		if [ $Current -lt $BATCritical ]; then
            killall jackd; killall alsa_in
            sudo pm_suspend
			bash $SCRIPTSPATH/Bash/ALSA/SoundAutostart.sh
            
	    elif [[ $Current -gt $BATCritical && $Current -lt $BATLow ]]; then
#            echo "Battery below $BATLow"
            Brightness=$(xbacklight)

            if [ $Brightness -gt 30.0 ]; then
                xbacklight -set 30
            fi
            
			feh -xF $SCRIPTSPATH/Bash/Power/BatteryLow.jpg
			sleep 90
            
		elif [ $Current -ge $BATLow ]; then
#            echo "Battery above $BATLow"
			sleep 270
	    fi
    else
        sleep 40
	fi
fi
sleep 30
#echo "BM loop finished"
}

FileSizeCheck() {
#echo "FM loop starting..."
for File in $Files; do
    FileSize=$(stat -c "%s" "$File")
    if [ $FileSize -gt 5120000 ]; then
        echo ""$File" is bigger than 5MB. Truncating..."
        if [ "$File" == ~/.log/JackSession.log ]; then
            killall jackd; killall alsa_in; rm "$File"
            bash $SCRIPTSPATH/Bash/ALSA/SoundAutostart.sh
        else
            echo "" >  "$File"
        fi
    fi
done
sleep 60
#echo "FM loop finished"
}

while [ true ] ; do
    BatteryCheck &
    FileSizeCheck &
    wait
done
