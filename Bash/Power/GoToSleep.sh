#!/bin/bash

echo "Locking screen..."; slock &
echo "Stopping Jack...";
killall jackd &> /dev/null; killall alsa_in &> /dev/null
echo "Suspending..."; sudo pm-suspend

echo "Waking..."
echo "Starting Jack...";
bash ~/Software/Git/Malfatti/SciScripts/Bash/Audio/ALSA/SoundAutostart.sh &> /dev/null
echo "Awake."
