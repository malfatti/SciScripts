#!/bin/bash

while [[ $# -gt 1 ]]; do
    key="$1"

    case $key in
        --path)
            Path="$2"
        ;;
        --rt)
            RT="$2"
            shift
        ;;
        --priority)
            Prio="$2"
            shift
        ;;
        --card)
            Card="$2"
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
        --bridge)
            if [ "${2,,}" == true ]; then
                Bridge="True"
            else
                Bridge="False"
            fi
            
            shift
        ;;
#        *)
#            echo "Usage: JackCtl --rt [ r | R ] --rate <rate> --periods <amount> --periodsize <period size>"
#        ;;
    esac
    shift
done

echo "Jack path: "$Path"" > ~/.log/JackSession.log &
echo "Jack options: jackd -"$RT" -P"$Prio" -dalsa -dhw:"$Card" -r"$Rate" -p"$PeriodSize" -n"$Periods"" >> ~/.log/JackSession.log &
echo "" >> ~/.log/JackSession.log &

"$Path"jackd -"$RT" -P"$Prio" -dalsa -dhw:"$Card" -r"$Rate" -p"$PeriodSize" -n"$Periods" &>> ~/.log/JackSession.log &

sleep 2

if [ $Bridge == True ]; then
    Catch=30000
    ResamplingRate=48000
    ResamplingQuality=1

    echo "Building output bridge..." >> ~/.log/JackSession.log &
    echo "alsa_out options:  -j ALSAInput -dALSAOutput1 -f "$Catch" -q "$ResamplingQuality" -r "$ResamplingRate" -p "$PeriodSize" -n "$Periods"" &>> ~/.log/JackSession.log &
    echo "" >> ~/.log/JackSession.log &

    "$Path"alsa_out -j ALSAInput -dALSAOutput1 -f "$Catch" -q "$ResamplingQuality" -r "$ResamplingRate" -p "$PeriodSize" -n "$Periods" &>> ~/.log/JackSession.log &
    #"$Path"alsa_out -j ALSAInput -dALSAOutput1 &>> ~/.log/JackSession.log &

    echo "Building input bridge..." >> ~/.log/JackSession.log &
    echo "alsa_in options: -j ALSAOutput -dALSAInput1 -f "$Catch" -q "$ResamplingQuality" -r "$ResamplingRate" -p "$PeriodSize" -n "$Periods"" &>> ~/.log/JackSession.log &
    echo "" >> ~/.log/JackSession.log &

    "$Path"alsa_in -j ALSAOutput -dALSAInput1 -f "$Catch" -q "$ResamplingQuality" -r "$ResamplingRate" -p "$PeriodSize" -n "$Periods" &>> ~/.log/JackSession.log &
    #"$Path"alsa_in -j ALSAOutput -dALSAInput1 &>> ~/.log/JackSession.log &

    sleep 2

    echo "Connecting bridge to output..." >> ~/.log/JackSession.log &
    "$Path"jack_connect ALSAOutput:capture_1 system:playback_1
    "$Path"jack_connect ALSAOutput:capture_2 system:playback_2

    echo "Connecting bridge to input..." >> ~/.log/JackSession.log &
    "$Path"jack_connect system:capture_1 ALSAInput:playback_1
    "$Path"jack_connect system:capture_2 ALSAInput:playback_2
fi

echo "" >> ~/.log/JackSession.log &
