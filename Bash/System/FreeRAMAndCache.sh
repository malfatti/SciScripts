#!/bin/bash

echo "Clearing RAM cache..."
sudo su -c "echo 1 >'/proc/sys/vm/drop_caches'"
echo "Clearing Swap..."
sudo swapoff -a && sudo swapon -a
echo "Done."
