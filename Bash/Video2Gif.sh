#!/bin/bash

while [[ $# -gt 1 ]]; do
    key="$1"

    case $key in
        --video)
            Vid="$2"
        ;;
        --start)
            Start="$2"
            shift
        ;;
        --dur)
            Dur="$2"
            shift
        ;;
        --fps)
            FPS="$2"
            shift
        ;;
        --wid)
            Wid="$2"
            shift
        ;;
    esac
    shift
done

VidName="${Vid%????}"

#echo "Separating frames..."
#ffmpeg -t $Dur -ss $Start -i $Vid "$VidName"-Temp-%04d.jpg
ffmpeg -t $Dur -ss $Start -i $Vid -vf scale=$Wid:-1 -r $FPS "$VidName".gif

#echo "Merging to gif..."
#convert -delay 1x"$FPS" -loop 0 "$VidName"-Temp*.jpg "$VidName".gif
#echo "Optimizing size..."
#convert -layers Optimize "$VidName".gif "$VidName"_small.gif
#echo "Cleaning..."
#rm "$VidName"-Temp*.jpg
#echo "Gif saved to "$VidName".gif."

