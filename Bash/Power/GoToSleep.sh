#!/bin/bash

echo "Stopping Jack..."; killall jackd; killall alsa_in
#echo "Stopping Bumblebee..."; killall bumblebeed
echo "Stopping Alsa..."; sudo /etc/init.d/alsasound stop
echo "Suspending..."; sudo pm-suspend

echo "Waking..."
echo "Locking screen..."; slock &
echo "Starting Alsa..."; sudo /etc/init.d/alsasound start
#echo "Starting Bumblebee..."; sudo bumblebeed &>> ~/.MyScripts/BumblebeeSession.log
echo "Starting Jack..."; bash ~/.MyScripts/ALSA/SoundAutostart.sh &> /dev/null
echo "Awake."
