#!/bin/bash

if [ "${1,,}" == --help ]; then
    echo "Usage: bash ./LaTex2Doc.sh [File name w/o .tex] [rtf | txt]"
fi

if [ -d 'DocVersion' ]; then
	echo "DocVersion Directory exists"
	else
		mkdir DocVersion
		echo "DocVersion Directory created"
fi

if [ $2 == 'txt' ]; then
	echo "Preparing $@.dvi..."
	latex -interaction=nonstopmode "$1" > /dev/null
	echo "Converting to text and removing formatting..."
	catdvi --debug=0 "$1".dvi > "$1"1.txt
	cat "$1"1.txt | echo -n `sed 's/^$/STARTPARA/'`|sed 's/STARTPARA/\n\n/g' > "$1"2.txt
	tr -d '\014' < "$1"2.txt > "$1".txt
	echo "Cleaning..."
	rm "$1"1.txt "$1"2.txt
	mv "$1".txt 'DocVersion/'
	echo "Done."
fi

if [ $2 == 'rtf' ]; then
	echo "Converting to rtf..."
	latex -interaction=nonstopmode "$1" > /dev/null
	bibtex "$1" > /dev/null
	latex2rtf -d0 "$1"
	mv "$1".rtf 'DocVersion/'
	echo "Done."
fi
