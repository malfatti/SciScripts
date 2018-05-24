#!/bin/bash

## Power off storage device 

Device="$@"
echo "Syncing disks..."
sync

echo "Unmounting all partitions..."
for Part in "$Device"?; do udisksctl unmount --block-device $Part; done

echo "Powering off..."
udisksctl power-off --block-device "$Device"
