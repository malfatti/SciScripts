#!/bin/bash

while [[ $# -gt 1 ]]; do
    key="$1"

    case $key in
        --rt)
            RT="$2"
            shift
        ;;
        --priority)
            Prio="$2"
            shift
        ;;
        --rate)
            Rate="$2"
            shift
        ;;
        --periods)
            Periods="$2"
            shift
        ;;
        --periodsize)
            PeriodSize="$2"
            shift
        ;;
#        *)
#            echo "Usage: JackCtl --rt [ r | R ] --rate <rate> --periods <amount> --periodsize <period size>"
#        ;;
    esac
    shift
done

Catch=30000
ResamplingRate=48000
ResamplingQuality=1

echo "Jack options: -"$RT" -P"$Prio" -dalsa -dhw:PCH -r"$Rate" -p"$PeriodSize" -n"$Periods"" > ~/.MyScripts/ALSA/JackSession.log &

jackd -"$RT" -P"$Prio" -dalsa -dhw:PCH -r"$Rate" -p"$PeriodSize" -n"$Periods" &>> ~/.MyScripts/ALSA/JackSession.log &

sleep 2

echo "Building output bridge..." >> ~/.MyScripts/ALSA/JackSession.log &
echo "alsa_out options:  -j ALSAInput -dALSAOutput1 -f "$Catch" -q "$ResamplingQuality" -r "$ResamplingRate" -p "$PeriodSize" -n "$Periods"" &>> ~/.MyScripts/ALSA/JackSession.log &
/usr/bin/alsa_out -j ALSAInput -dALSAOutput1 -f "$Catch" -q "$ResamplingQuality" -r "$ResamplingRate" -p "$PeriodSize" -n "$Periods" &>> ~/.MyScripts/ALSA/JackSession.log &

echo "Building input bridge..." >> ~/.MyScripts/ALSA/JackSession.log &
echo "alsa_in options: -j ALSAOutput -dALSAInput1 -f "$Catch" -q "$ResamplingQuality" -r "$ResamplingRate" -p "$PeriodSize" -n "$Periods"" &>> ~/.MyScripts/ALSA/JackSession.log &
/usr/bin/alsa_in -j ALSAOutput -dALSAInput1 -f "$Catch" -q "$ResamplingQuality" -r "$ResamplingRate" -p "$PeriodSize" -n "$Periods" &>> ~/.MyScripts/ALSA/JackSession.log &

sleep 2

echo "Connecting bridge to output..." >> ~/.MyScripts/ALSA/JackSession.log &
jack_connect ALSAOutput:capture_1 system:playback_1
jack_connect ALSAOutput:capture_2 system:playback_2

echo "Connecting bridge to input..." >> ~/.MyScripts/ALSA/JackSession.log &
jack_connect system:capture_1 ALSAInput:playback_1
jack_connect system:capture_2 ALSAInput:playback_2
