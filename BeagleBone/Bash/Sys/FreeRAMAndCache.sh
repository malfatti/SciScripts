#!/bin/bash

echo "Clearing RAM cache..."
echo 1 >'/proc/sys/vm/drop_caches'
echo "Done."
