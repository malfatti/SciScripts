#!/bin/bash

echo "Stopping Jack...";
killall jackd &> /dev/null; killall alsa_in &> /dev/null
#echo "Stopping Bumblebee..."; killall bumblebeed
echo "Stopping Alsa..."; sudo /etc/init.d/alsasound stop &> /dev/null
echo "Suspending..."; sudo pm-suspend

echo "Waking..."
echo "Locking screen..."; slock &
echo "Starting Alsa..."; sudo /etc/init.d/alsasound start &> /dev/null
#echo "Starting Bumblebee..."; sudo bumblebeed &>> ~/Software/Git/Malfatti/SciScripts/Bash/BumblebeeSession.log
echo "Starting Jack...";
bash ~/Software/Git/Malfatti/SciScripts/Bash/ALSA/SoundAutostart.sh &> /dev/null
echo "Awake."
