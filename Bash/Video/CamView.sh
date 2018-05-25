#!/bin/bash

mplayer tv:// -tv driver=v4l2:width=1280:height=720:device=/dev/video1 -fps 30 -vf screenshot -vo x11
