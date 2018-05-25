#!/bin/bash

sed 's/[[:blank:]]*$//' "$1" > ."$1".tmp
mv ."$1".tmp "$1"
