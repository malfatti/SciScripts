#!/bin/bash

Start="Empty"; End="Empty"

while getopts :v:s:e:h Option; do
    case "${Option}" in
        v)
            Video=${OPTARG}
        ;;
        s)
            Start=${OPTARG}
        ;;
        e)
            End=${OPTARG}
        ;;
        h) 
            echo -en '\n'
            echo "Usage: CutVideo -v <FileName> [-s <DurationInSeconds> | -e <DurationInSeconds>]"
            echo -en '\n'
            
            exit 0
        ;;
        :)
            echo "Option -"$OPTARG" needs a value. See 'CutVideo -h'."
            exit 0
        ;;
        \?)
            echo "Unknown option "$OPTARG". See 'CutVideo -h'."
            exit 0
        ;;
    esac
done

NewFileName="${Video%????}"-Cut."${Video: -3}"

if [[ "$Start" != "Empty" && "$End" != "Empty" ]]; then
    ffmpeg -i "$Video" -ss "$Start" -to "$End" -c copy "$NewFileName"
    echo -en '\n'; echo -en '\n'; echo "Video saved as "$NewFileName""
elif [[ "$Start" != "Empty" && "$End" == "Empty" ]]; then
    ffmpeg -i "$Video" -ss "$Start" -c copy "$NewFileName"
    echo -en '\n'; echo -en '\n'; echo "Video saved as "$NewFileName""
elif [[ "$Start" == "Empty" && "$End" != "Empty" ]]; then
    ffmpeg -i "$Video" -to "$End" -c copy "$NewFileName"
    echo -en '\n'; echo -en '\n'; echo "Video saved as "$NewFileName""
else
    echo "Choose a start time or an end time. See 'CutVideo -h'."
fi


