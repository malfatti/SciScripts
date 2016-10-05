#!/bin/bash

if [ "${1,,}" == --help ]; then
    echo "Usage: CutVideo [FileName] [DurationInSeconds]"
    exit 0
fi

Video="$1"
Start="$2"
NewFileName="${Video%????}"-Cut."${Video: -3}"

ffmpeg -i "$Video" -ss "$Start" -c copy "$NewFileName"
echo -en '\n'; echo -en '\n'; echo "Video saved as "$NewFileName""
