#!/bin/bash

Start=""; End=""

while [[ $# -gt 1 ]]; do
    key="$1"

    case $key in
        "--help")
            echo "Usage: "$0" [FileName] [DurationInSeconds]"
            exit 0
        ;;
        "--video")
            Video="$2"
            shift 2
        ;;
        "--start")
            Start="$2"
            shift 2
        ;;
        "--end")
            End="$2"
            shift 2
        ;;
        *)
            echo "Usage: "$0" [FileName] [DurationInSeconds]"
            exit 0
        ;;
    esac
done

NewFileName="${Video%????}"-Cut."${Video: -3}"

if [ "$End" == "" ]; then
    ffmpeg -i "$Video" -ss "$Start" -c copy "$NewFileName"
elif [ "$Start" == "" ]; then
    ffmpeg -i "$Video" -to "$End" -c copy "$NewFileName"
else
    ffmpeg -i "$Video" -ss "$Start" -to "$End" -c copy "$NewFileName"
fi

echo -en '\n'; echo -en '\n'; echo "Video saved as "$NewFileName""

