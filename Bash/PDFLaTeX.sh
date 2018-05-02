#!/bin/bash

File="${1%.*}"

mkdir ."$File"
pdflatex -output-directory ."$File" "$1" #&>> ."$File"/"$File".log
mv ."$File"/"$File".pdf .

