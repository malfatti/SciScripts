#!/bin/bash

Date=$(date +%Y%m%d-%H%M%S)

ffmpeg -f v4l2 \
       -framerate 25 \
       -video_size 1280x720 \
-i /dev/video1 Video-$Date.mp4
